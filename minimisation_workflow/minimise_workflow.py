"""
A command line interface for running a minimisation workflow.

Workflow steps:
1. Load the ligands and create a perturbation map using the settings in the config file.
2. Following the Hybrid topology approach, create a hybrid topology for each edge in the perturbation map.
3. Minimise each end state using the pure and hybrid topologies with OpenMM.
4. Collect the positions of the minimised ligands and write them to an output SDF file.
5. Calculate the RMSD of the hybrid end states to the pure end states and calculate the relative energy difference using the pure topologies and save to CSV for analysis.
"""

import click
from pathlib import Path
import logging
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from gufe import SmallMoleculeComponent, LigandNetwork
from openfe.protocols.openmm_rfe._rfe_utils.relative import HybridTopologyFactory
from openfe.setup.ligand_network_planning import generate_lomap_network
from openfe.setup import KartografAtomMapper
from openfe.setup.atom_mapping.lomap_scorers import default_lomap_score
from openfe.protocols.openmm_utils import system_creation, settings_validation
from openfe.protocols.openmm_rfe import _rfe_utils
from functools import partial
import tqdm
import pandas as pd
from openff.toolkit import (
    RDKitToolkitWrapper, AmberToolsToolkitWrapper
)
from openff.toolkit.utils.toolkit_registry import (
    toolkit_registry_manager, ToolkitRegistry
)
from openmmforcefields.generators import SystemGenerator
from openmm import app
import openmm
from openmm import unit
from openff.units.openmm import to_openmm, ensure_quantity, from_openmm
from openff.units import unit as off_unit
from itertools import chain
from collections import defaultdict
from yammbs.analysis import get_internal_coordinate_rmsds
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("talk")

amber_rdkit = ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])


def internal_coord_plot(df, filename):
    """Make a boxplot of the internal coordinate RMSDs"""
    fig, ax = plt.subplots(1, 5, figsize=(12, 6))
    color_map = {
        "bond": "skyblue",
        "angle": "red",
        "dihedral": "orange",
        "improper": "lightgreen"
    }
    label_map = {
        "bond": "Bond (Ang)",
        "angle": "Angle (degrees)",
        "dihedral": "Dihedral (degrees)",
        "improper": "Improper (degrees)"
    }
    for i, interal_coord in enumerate(["bond", "angle", "dihedral", "improper"]):
        ax[i].boxplot(df[f"internal_rmsd_{interal_coord}"].values, boxprops={"facecolor": color_map[interal_coord]}, patch_artist=True, medianprops={"color": "k"})
        ax[i].set_xticks([1],[interal_coord.capitalize()])
        ax[i].set_ylabel(label_map[interal_coord])
    # add the all atom RMSD as well
    ax[-1].boxplot(df["rmsd_to_pure_ang"].values, boxprops={"facecolor": "purple"}, patch_artist=True, medianprops={"color": "k"})
    ax[-1].set_xticks([1], ["All atom RMSD"])
    ax[-1].set_ylabel("RMSD (Ang)")
    sns.despine()
    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close()

def plot_energy_difference(df, filename):
    """Make a histogram of the energy differences"""
    plt.figure(figsize=(8, 6))
    sns.histplot(df["energy_difference"], bins=30, kde=True, color="skyblue")
    plt.xlabel("Energy Difference (kcal/mol)")
    plt.ylabel("Count")
    sns.despine()
    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close()


def gen_charges(smc):
    """
    Generate AM1BCC partial charges for a SmallMoleculeComponent using
    the input conformer and antechamber as backend.

    Skip the charge generation if the component already has charges.
    """
    offmol = smc.to_openff()
    if offmol.partial_charges is None:
        print(f"INFO: generating partial charges for ligand {smc.name} -- this may be slow")
        with toolkit_registry_manager(amber_rdkit):
            offmol.assign_partial_charges(
                partial_charge_method="am1bcc",
                use_conformers=offmol.conformers
            )
    return SmallMoleculeComponent.from_openff(offmol)


def make_hybrid_factory(edge, system_generator):
    small_mols = [edge.componentA, edge.componentB]
    off_small_mols = {
        'stateA': [(edge.componentA, edge.componentA.to_openff())],
        'stateB': [(edge.componentB, edge.componentB.to_openff())],
        'both': [(m, m.to_openff()) for m in small_mols
                 if (m != edge.componentA and m != edge.componentB)]
    }
    # c. force the creation of parameters
    # This is necessary because we need to have the FF templates
    # registered ahead of solvating the system.
    for smc, mol in chain(off_small_mols['stateA'],
                          off_small_mols['stateB'],
                          off_small_mols['both']):
        system_generator.create_system(mol.to_topology().to_openmm(),
                                       molecules=[mol])

        # c. get OpenMM Modeller + a dictionary of resids for each component
        stateA_modeller, comp_resids = system_creation.get_omm_modeller(
            protein_comp=None,
            solvent_comp=None,
            small_mols=dict(chain(off_small_mols['stateA'],
                                  off_small_mols['both'])),
            omm_forcefield=None,
            solvent_settings=None,
        )

    stateA_topology = stateA_modeller.getTopology()
    stateA_positions = to_openmm(
        from_openmm(stateA_modeller.getPositions())
    )

    stateA_system = system_generator.create_system(
        stateA_modeller.topology,
        molecules=[m for _, m in chain(off_small_mols['stateA'],
                                       off_small_mols['both'])],
    )

    stateB_topology, stateB_alchem_resids = _rfe_utils.topologyhelpers.combined_topology(
        stateA_topology,
        # zeroth item (there's only one) then get the OFF representation
        off_small_mols['stateB'][0][1].to_topology().to_openmm(),
        exclude_resids=comp_resids[edge.componentA],
    )

    # b. get a list of small molecules for stateB
    stateB_system = system_generator.create_system(
        stateB_topology,
        molecules=[m for _, m in chain(off_small_mols['stateB'],
                                       off_small_mols['both'])],
    )

    #  c. Define correspondence mappings between the two systems
    ligand_mappings = _rfe_utils.topologyhelpers.get_system_mappings(
        edge.componentA_to_componentB,
        stateA_system, stateA_topology, comp_resids[edge.componentA],
        stateB_system, stateB_topology, stateB_alchem_resids,
        # These are non-optional settings for this method
        fix_constraints=True,
    )

    #  e. Finally get the positions
    stateB_positions = _rfe_utils.topologyhelpers.set_and_check_new_positions(
        ligand_mappings, stateA_topology, stateB_topology,
        old_positions=ensure_quantity(stateA_positions, 'openmm'),
        insert_positions=ensure_quantity(off_small_mols['stateB'][0][1].conformers[0], 'openmm'),
    )
    # 3. Create the hybrid topology
    # b. Get hybrid topology factory
    hybrid_factory = _rfe_utils.relative.HybridTopologyFactory(
        stateA_system, stateA_positions, stateA_topology,
        stateB_system, stateB_positions, stateB_topology,
        old_to_new_atom_map=ligand_mappings['old_to_new_atom_map'],
        old_to_new_core_atom_map=ligand_mappings['old_to_new_core_atom_map'],
        softcore_LJ_v2=True
    )
    return hybrid_factory


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


@click.command()
@click.option(
    "--output",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory to save the output files.",
)
@click.option(
    "--ligands",
    "-l",
    type=click.Path(exists=True, dir_okay=True, file_okay=True, path_type=Path),
    default=None,
    help="Directory containing input ligand files.",
)
@click.option(
    "--network",
    "-n",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    default=None,
    help="Optional JSON file containing a pre-defined perturbation map.",
)
def main(ligands: Path, output: Path, network: None | Path):
    """
    Command line interface for running a minimisation workflow.
    This can be configured using a YAML file.

    Workflow steps:
    1. Load the ligands and create a perturbation map.
    2. Following the Hybrid topology approach, create a hybrid topology for each edge in the perturbation map.
    3. Minimise each end state using the pure and hybrid topologies with OpenMM.
    4. Collect the positions of the minimised ligands and write them to an output SDF file.
    5. Calculate the RMSD of the hybrid end states to the pure end states and calculate the relative energy difference using the pure topologies and save to CSV for analysis.
    """
    if network is None:
        # load the ligands
        supplier = Chem.SDMolSupplier(ligands.as_posix(), removeHs=False)
        smcs = [SmallMoleculeComponent.from_rdkit(mol) for mol in supplier]
        logger.info(f"Loaded {len(smcs)} ligands from {ligands}")

        # Generate the partial charges
        logger.info("Generating partial charges for ligands")
        smcs = [gen_charges(smc) for smc in smcs]

        # create the perturbation map using the same settings as the industry benchmarks
        scorer = partial(default_lomap_score, charge_changes_score=0.1)
        mapper = KartografAtomMapper(map_hydrogens_on_hydrogens_only=True)
        ligand_network = generate_lomap_network(
            molecules=smcs, mappers=mapper, scorer=scorer
        )
        logger.info(f"Generated perturbation map with {len(ligand_network.edges)} edges")

        # create the output directory if it doesn't exist
        output.mkdir(parents=True, exist_ok=True)
        # save the perturbation map
        ligand_network.to_json(output / "ligand_network.json")
        logger.info(f"Saved perturbation map to {output / 'ligand_network.json'}")
    else:
        # load the perturbation map
        ligand_network = LigandNetwork.from_json(file=network)
        logger.info(f"Loaded perturbation map with {len(ligand_network.edges)} edges from {network}")
    
    # For each ligand in the perturbation map, create a pure end-state topology and minimise it
    pure_endstate_data = {}

    # create the system generator
    system_generator = SystemGenerator(
        small_molecule_forcefield="openff-2.2.0",
        nonperiodic_forcefield_kwargs={'nonbondedMethod': app.NoCutoff},
        cache=str(output / "ff_cache.json")
    )

    # Minimise each hybrid edge topologies
    hybrid_endstate_data = defaultdict(list)
    for edge in tqdm.tqdm(ligand_network.edges, desc="Minimising hybrid end-states"):
        htf = make_hybrid_factory(edge, system_generator)
        # c. Get the simulation object
        hybrid_system = htf.hybrid_system
        integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 1*unit.femtoseconds)
        simulation = app.Simulation(
            topology=htf.omm_hybrid_topology,
            system=hybrid_system,
            integrator=integrator
        )
        default_lambda = _rfe_utils.lambdaprotocol.LambdaProtocol()
        for i in [0, 1]:
            # set the positions to the hybrid positions
            simulation.context.setPositions(ensure_quantity(htf.hybrid_positions, "openmm"))
            # set all lambda values to the current end state
            for name, func in default_lambda.functions.items():
                val = func(i)
                simulation.context.setParameter(name, val)

            simulation.minimizeEnergy()
            state = simulation.context.getState(getPositions=True)
            positions = state.getPositions(asNumpy=True)
            # extract only the positions of the relevant component
            if i == 0:
                end_positions = htf.old_positions(positions).value_in_unit(unit.angstrom)
            else:
                end_positions = htf.new_positions(positions).value_in_unit(unit.angstrom)

            # set the name to correspond to the end state
            name = edge.componentA.name if i == 0 else edge.componentB.name
            hybrid_endstate_data[name].append(
                {
                    "edge": (edge.componentA.name, edge.componentB.name),
                    "positions": end_positions,
                }
            )
        # remove the simulation
        del simulation

    # create and minimise the pure end-states and save the positions and energy
    for smc in tqdm.tqdm(ligand_network.nodes, desc="Minimising pure end-states"):
        openff_mol = smc.to_openff()
        system = system_generator.create_system(openff_mol.to_topology().to_openmm(), molecules=[openff_mol])
        integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 1 * unit.femtoseconds)
        simulation = app.Simulation(
            topology=openff_mol.to_topology().to_openmm(),
            system=system,
            integrator=integrator
        )
        # make sure we wrap back to openmm units
        simulation.context.setPositions(ensure_quantity(openff_mol.conformers[0], "openmm"))
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        pure_energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        pure_endstate_data[smc.name] = {
            "positions": positions,
            "energy": pure_energy,
            "smc": smc
        }
        del state
        # create a ref rdkit mol from the pure endstate to use the RMSD calculations
        ref_mol = smc.to_rdkit()
        # set the positions
        conf = Chem.Conformer(ref_mol.GetNumAtoms())
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, (pos[0], pos[1], pos[2]))
        ref_mol.RemoveAllConformers()
        ref_mol.AddConformer(conf, assignId=True)

        # now loop over any hybrid end states that correspond to this ligand
        for hybrid_data in hybrid_endstate_data[smc.name]:
            hybrid_positions = hybrid_data["positions"]
            # create a prb rdkit mol from the hybrid endstate, make sure its a new mol
            prb_mol = Chem.Mol(smc.to_rdkit())
            conf = Chem.Conformer(prb_mol.GetNumAtoms())
            for i, pos in enumerate(hybrid_positions):
                conf.SetAtomPosition(i, (pos[0], pos[1], pos[2]))
            prb_mol.RemoveAllConformers()
            prb_mol.AddConformer(conf, assignId=True)
            # calculate the RMSD between the pure and hybrid positions using rdkit
            # we do no alignment first as the positions should be very similar
            rmsd = rdMolAlign.GetBestRMS(prb_mol, ref_mol)
            # save the RMSD to the hybrid data
            hybrid_data["rmsd_to_pure_ang"] = rmsd
            # calculate the energy at the hybrid positions using the pure topology
            simulation.context.setPositions(positions=hybrid_positions * unit.angstrom)
            state = simulation.context.getState(getEnergy=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            hybrid_data["energy_difference"] = energy - pure_energy
            del state

            # calculate the internal coordinate RMSD (bond, angle, dihedral) between the pure and hybrid positions do this highlight more differences
            offmol = smc.to_openff()
            internal_rmsd = get_internal_coordinate_rmsds(
                molecule=offmol,
                # use the pure end state positions as the reference
                reference=positions * off_unit.angstrom,
                # use the hybrid end state positions as the target
                target=hybrid_positions * off_unit.angstrom,

            )
            for key, value in internal_rmsd.items():
                hybrid_data[f"internal_rmsd_{key.lower()}"] = value



        # remove the simulation
        del simulation

    internal_coord_terms = ["Bond", "Angle", "Dihedral", "Improper"]
    # write out the results
    # 1.write out the pure end states to an SDF file
    writer = Chem.SDWriter((output / "pure_endstates.sdf").as_posix())
    for smc_name, data in pure_endstate_data.items():
        smc = data["smc"]
        rdkit_mol = smc.to_rdkit()
        conf = Chem.Conformer(rdkit_mol.GetNumAtoms())
        positions = data["positions"]
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, (pos[0], pos[1], pos[2]))
        rdkit_mol.RemoveAllConformers()
        rdkit_mol.AddConformer(conf, assignId=True)
        rdkit_mol.SetProp("Name", smc_name)
        rdkit_mol.SetProp("Pure_Energy_kcal_per_mol", str(data["energy"]))
        writer.write(rdkit_mol)
    writer.close()
    logger.info(f"Wrote pure end-states to {output / 'pure_endstates.sdf'}")
    # 2. write out the hybrid end states to an SDF file use one file per unique ligand
    for smc_name, hybrid_list in hybrid_endstate_data.items():
        writer = Chem.SDWriter((output / f"hybrid_endstates_{smc_name}.sdf").as_posix())
        for i, hybrid_data in enumerate(hybrid_list):
            rdkit_mol = pure_endstate_data[smc_name]["smc"].to_rdkit()
            conf = Chem.Conformer(rdkit_mol.GetNumAtoms())
            positions = hybrid_data["positions"]
            for j, pos in enumerate(positions):
                conf.SetAtomPosition(j, (pos[0], pos[1], pos[2]))
            rdkit_mol.RemoveAllConformers()
            rdkit_mol.AddConformer(conf, assignId=True)
            rdkit_mol.SetProp("Name", f"{smc_name}_from_edge_{hybrid_data['edge'][0]}_{hybrid_data['edge'][1]}")
            rdkit_mol.SetProp("RMSD_to_pure_ang", str(hybrid_data["rmsd_to_pure_ang"]))
            rdkit_mol.SetProp("Energy_Difference_kcal_per_mol", str(hybrid_data["energy_difference"]))
            for term in internal_coord_terms:
                rdkit_mol.SetProp(f"Internal_RMSD_{term}", str(hybrid_data.get(f"internal_rmsd_{term}", "NA")))
            writer.write(rdkit_mol)
        writer.close()
        logger.info(f"Wrote hybrid end-states to {output / f'hybrid_endstates_{smc_name}.sdf'}")
    # 3. write out a CSV file with the RMSD and energy differences for each edge
    rows = []
    for smc_name, hybrid_list in hybrid_endstate_data.items():
        for hybrid_data in hybrid_list:
            # format the data for the CSV
            edge = hybrid_data.pop("edge")
            hybrid_data["ligand_a"] = edge[0]
            hybrid_data["ligand_b"] = edge[1]
            hybrid_data["ligand"] = smc_name
            hybrid_data.pop("positions")
            rows.append(hybrid_data)
    df = pd.DataFrame(rows)
    df.to_csv(output / "hybrid_endstate_analysis.csv", index=False)
    logger.info(f"Wrote hybrid end-state analysis to {output / 'hybrid_endstate_analysis.csv'}")
    # plot the distributions of internal coord RMSD differences
    internal_coord_plot(df=df, filename=str(output / "internal_coordinate_rmsd.png"))
    logger.info(f"Wrote internal coordinate RMSD plot to {output / 'internal_coordinate_rmsd.png'}")
    plot_energy_difference(df=df, filename=str(output / "energy_difference_distribution.png"))
    logger.info(f"Wrote energy difference distribution plot to {output / 'energy_difference_distribution.png'}")



if __name__ == "__main__":
    main()









    
    


