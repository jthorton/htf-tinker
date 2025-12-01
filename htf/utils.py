from openff.units.openmm import to_openmm, ensure_quantity, from_openmm
from itertools import chain
from htf import DevelopmentHybridTopologyFactory
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit
import openmm
from gufe import LigandAtomMapping, ProteinComponent, SolventComponent
from openfe.protocols.openmm_rfe import _rfe_utils
from openfe.protocols.openmm_utils import system_creation
import ast
import copy
import math
import logging

logger = logging.getLogger(__name__)


def make_htf(mapping: LigandAtomMapping, settings, protein: ProteinComponent = None, solvent: SolventComponent = None) -> DevelopmentHybridTopologyFactory:
    """Code copied from the RBFE protocol to make an HTF."""

    system_generator = SystemGenerator(
        forcefields=settings.forcefield_settings.forcefields,
        small_molecule_forcefield=settings.forcefield_settings.small_molecule_forcefield,
        forcefield_kwargs={
            "constraints": app.HBonds,
            "rigidWater": True,
            "hydrogenMass": settings.forcefield_settings.hydrogen_mass * unit.amu,
            "removeCMMotion": settings.integrator_settings.remove_com
        },
        periodic_forcefield_kwargs={
            'nonbondedMethod': app.PME,
            'nonbondedCutoff': 0.9 * unit.nanometers,
        },
        barostat=openmm.MonteCarloBarostat(
            ensure_quantity(settings.thermo_settings.pressure, 'openmm'),
            ensure_quantity(settings.thermo_settings.temperature, 'openmm'),
            settings.integrator_settings.barostat_frequency.m,
        ),
        cache=None
    )
    small_mols = [mapping.componentA, mapping.componentB]
    # copy a lot of code from the RHT protocol
    off_small_mols = {
        'stateA': [(mapping.componentA, mapping.componentA.to_openff())],
        'stateB': [(mapping.componentB, mapping.componentB.to_openff())],
        'both': [(m, m.to_openff()) for m in small_mols
                 if (m != mapping.componentA and m != mapping.componentB)]
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
        # add the protein if passed
        protein_comp=protein,
        # add the solvent if passed
        solvent_comp=solvent,
        small_mols=dict(chain(off_small_mols['stateA'],
                              off_small_mols['both'])),
        omm_forcefield=system_generator.forcefield,
        solvent_settings=settings.solvation_settings,
    )
    # d. get topology & positions
    # Note: roundtrip positions to remove vec3 issues
    stateA_topology = stateA_modeller.getTopology()
    stateA_positions = to_openmm(
        from_openmm(stateA_modeller.getPositions())
    )

    # e. create the stateA System
    # Block out oechem backend in system_generator calls to avoid
    # any issues with smiles roundtripping between rdkit and oechem
    stateA_system = system_generator.create_system(
        stateA_modeller.topology,
        molecules=[m for _, m in chain(off_small_mols['stateA'],
                                       off_small_mols['both'])],
    )

    # 2. Get stateB system
    # a. get the topology
    stateB_topology, stateB_alchem_resids = _rfe_utils.topologyhelpers.combined_topology(
        stateA_topology,
        # zeroth item (there's only one) then get the OFF representation
        off_small_mols['stateB'][0][1].to_topology().to_openmm(),
        exclude_resids=comp_resids[mapping.componentA],
    )

    # b. get a list of small molecules for stateB
    # Block out oechem backend in system_generator calls to avoid
    stateB_system = system_generator.create_system(
        stateB_topology,
        molecules=[m for _, m in chain(off_small_mols['stateB'],
                                       off_small_mols['both'])],
    )

    #  c. Define correspondence mappings between the two systems
    ligand_mappings = _rfe_utils.topologyhelpers.get_system_mappings(
        mapping.componentA_to_componentB,
        stateA_system, stateA_topology, comp_resids[mapping.componentA],
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
    return DevelopmentHybridTopologyFactory(
        old_system=stateA_system,
        old_positions=stateA_positions,
        old_topology=stateA_topology,
        new_system=stateB_system,
        new_positions=stateB_positions,
        new_topology=stateB_topology,
        old_to_new_atom_map=ligand_mappings["old_to_new_atom_map"],
        old_to_new_core_atom_map=ligand_mappings["old_to_new_core_atom_map"],
        use_dispersion_correction=settings.alchemical_settings.use_dispersion_correction,
        softcore_alpha=settings.alchemical_settings.softcore_alpha,
        softcore_LJ_v2=True,
        softcore_LJ_v2_alpha=settings.alchemical_settings.softcore_alpha,
        interpolate_old_and_new_14s=settings.alchemical_settings.turn_off_core_unique_exceptions,
    )

def _parse_ghostly_output(ghostly_output_path: str) -> dict:
    """Parse the output file from Ghostly into a dictionary of applied corrections which can be processed further."""
    corrections = {
        "lambda_0":
            {"bridges":{}, "removed_angles": [], "removed_dihedrals": [], "stiffened_angles": [], "softened_angles": {}},
        "lambda_1":
            {"bridges":{}, "removed_angles": [], "removed_dihedrals": [], "stiffened_angles": [], "softened_angles": {}}
    }
    current_lambda = None
    bridges = False
    current_bridge = None
    bridge_id = None
    with open(ghostly_output_path, 'r') as f:
        for line in f.readlines():
            # split out the debugging prefix
            _, line = line.strip().split(" - ", 1)

            # control flags for bridges
            if "Ghost atom bridges at lambda = 0" in line:
                current_lambda = "lambda_0"
                bridges = True
                continue
            elif "Ghost atom bridges at lambda = 1" in line:
                current_lambda = "lambda_1"
                bridges = True
                continue
            elif "Bridge" in line and bridges:
                # Example output
                # Bridge 0: 1
                # start of a new bridge
                parts = line.strip().split()
                bridge_id = int(parts[1][0])
                current_bridge = int(parts[-1])
                corrections[current_lambda]["bridges"][bridge_id] = {"bridge_atom": current_bridge}
                continue
            elif "ghosts" in line and bridges:
                # Example output
                # ghosts: [0]
                parts = line.strip().split()
                ghosts = ast.literal_eval(parts[1])
                corrections[current_lambda]["bridges"][bridge_id]["ghosts"] = ghosts
                continue
            elif "physical" in line and bridges:
                # Example output
                # physical: [2,3,4,8]
                parts = line.strip().split()
                physicals = ast.literal_eval(parts[1])
                corrections[current_lambda]["bridges"][bridge_id]["physical"] = physicals
                continue
            elif "type" in line and bridges:
                # Example output
                # type: 4
                parts = line.strip().split()
                btype = int(parts[1])
                corrections[current_lambda]["bridges"][bridge_id]["type"] = btype
                continue
            # control flags for modifications
            elif "Applying modifications" in line:
                # Example output
                # Applying modifications to triple ghost junction at λ = 0:
                bridges = False
                if "λ = 0:" in line:
                    current_lambda = "lambda_0"
                elif "λ = 1:" in line:
                    current_lambda = "lambda_1"
                continue
            elif "Removing angle" in line and not bridges:
                # Example output
                # Removing angle: [2-1-8], 66.7563 [theta - 1.91878]^2
                parts = line.strip().split()
                angle = tuple(ast.literal_eval(parts[2][:-1].replace("-", ",")))
                corrections[current_lambda]["removed_angles"].append(angle)
                continue
            elif "Removing dihedral" in line and not bridges:
                # Example output
                # Removing dihedral: [5-2-1-8], 0.227995 cos(3 phi) + 0.227995
                parts = line.strip().split()
                dihedral = tuple(ast.literal_eval(parts[2][:-1].replace("-", ",")))
                corrections[current_lambda]["removed_dihedrals"].append(dihedral)
                continue
            elif "Stiffening angle" in line and not bridges:
                # angles are always stiffened to 100 kcal/mol/rad^2 and an equilibrium of 90 degrees
                # Example output
                # Stiffening angle: [4-1-8], 36.7766 [theta - 1.89114]^2 --> 100 [theta - 1.5708]^2
                parts = line.strip().split()
                angle = tuple(ast.literal_eval(parts[2][:-1].replace("-", ",")))
                corrections[current_lambda]["stiffened_angles"].append(angle)
                continue
            elif "Softening angle" in line and not bridges:
                # angles are later optimised by default in ghostly, so we just need to track which ones were softened
                # Example output
                # Softening angle: [17-19-30], 46.9 [theta - 1.8984]^2 --> 0.05 [theta - 1.8984]^2
                parts = line.strip().split()
                angle = tuple(ast.literal_eval(parts[2][:-1].replace("-", ",")))
                # add the default softening parameters k = 0.05 kcal/mol/rad^2 and theta stays the same
                corrections[current_lambda]["softened_angles"][angle] = {"k": 0.05}
                continue
            elif "Optimising angle" in line and not bridges:
                # if the softened angle is optimised, we need to update its parameters
                # Example output
                # Optimising angle: [17-19-30], 0.05 [theta - 1.8984]^2 --> 5 [theta - 1.66289]^2 (std err: 0.006 radian)
                parts = line.strip().split()
                angle = tuple(ast.literal_eval(parts[2][:-1].replace("-", ",")))
                new_params = line.strip().split("-->")[1].split()
                new_k = int(new_params[0])
                new_equ_theta = float(new_params[3].split("]")[0])
                corrections[current_lambda]["softened_angles"][angle] = {"k": new_k, "theta_eq": new_equ_theta}
                continue
            elif "junction" in line:
                continue
            else:
                raise RuntimeError(f"Could not parse line in Ghostly output: {line}")

    return corrections


def apply_ghostly_corrections(htf: DevelopmentHybridTopologyFactory, corrections: dict) -> DevelopmentHybridTopologyFactory:
    """Apply the ghostly corrections parsed from the output file to the HTF."""

    htf_corrected = copy.deepcopy(htf)

    new_hybrid_system = openmm.System()
    # add all the particles
    for i in range(htf_corrected.hybrid_system.getNumParticles()):
        new_hybrid_system.addParticle(htf_corrected.hybrid_system.getParticleMass(i))
    # add all constraints
    for i in range(htf_corrected.hybrid_system.getNumConstraints()):
        p1, p2, dist = htf_corrected.hybrid_system.getConstraintParameters(i)
        new_hybrid_system.addConstraint(p1, p2, dist)

    hybrid_forces = htf._hybrid_system_forces
    print(hybrid_forces)
    # copy all forces which do not need to be modified
    # We are only modifying angle and torsion forces with ghostly corrections
    # As the HTF stores all terms involving ghosts in the standard forces we can copy all others directly
    # The interpolated forces only contain terms for the core mapped atoms so we don't need to remove any
    forces_not_to_copy = ["standard_angle_force", "unique_atom_torsion_force"]
    for force_name, hybrid_force in hybrid_forces.items():
        if force_name not in forces_not_to_copy:
            new_force = copy.deepcopy(hybrid_force)
            new_hybrid_system.addForce(new_force)

    # now apply the ghostly corrections to the angle and torsion forces
    # first add a new standard angle and torsion force to the system
    new_harmonic_angle_force = openmm.HarmonicAngleForce()
    new_torsion_force = openmm.PeriodicTorsionForce()
    new_hybrid_system.addForce(new_harmonic_angle_force)
    new_hybrid_system.addForce(new_torsion_force)
    # get a quick lookup of the forces
    new_hybrid_forces = {force.getName(): force for force in new_hybrid_system.getForces()}
    # process angles
    custom_angle_force = new_hybrid_forces["CustomAngleForce"]
    old_hybrid_angle_force = hybrid_forces["standard_angle_force"]
    for i in range(old_hybrid_angle_force.getNumAngles()):
        p1, p2, p3, theta_eq, k = old_hybrid_angle_force.getAngleParameters(i)
        # check if we have one ghost atom for this angle
        angle = (p1, p2, p3)
        if len(htf_corrected._atom_classes["unique_old_atoms"].intersection(angle)) <= 2 or len(htf_corrected._atom_classes["unique_new_atoms"].intersection(angle)) <= 2:
            angle_reversed = (p3, p2, p1)
            # set up containers for the end state values
            lambda_0_k = k
            lambda_0_theta_eq = theta_eq
            lambda_1_k = k
            lambda_1_theta_eq = theta_eq
            # check for removed angles
            if angle in corrections["lambda_0"]["removed_angles"] or angle_reversed in corrections["lambda_0"]["removed_angles"]:
                lambda_0_k = 0.0 * unit.kilocalories_per_mole / (unit.radian ** 2)
            elif angle in corrections["lambda_1"]["removed_angles"] or angle_reversed in corrections["lambda_1"]["removed_angles"]:
                lambda_1_k = 0.0 * unit.kilocalories_per_mole / (unit.radian ** 2)
            # check for stiffened angles
            elif angle in corrections["lambda_0"]["stiffened_angles"] or angle_reversed in corrections["lambda_0"]["stiffened_angles"]:
                lambda_0_k = 100.0 * unit.kilocalories_per_mole / (unit.radian ** 2)  # default stiffening k value
                lambda_0_theta_eq = 0.5 * math.pi * unit.radian  # 90 degrees
            elif angle in corrections["lambda_1"]["stiffened_angles"] or angle_reversed in corrections["lambda_1"]["stiffened_angles"]:
                lambda_1_k = 100.0 * unit.kilocalories_per_mole / (unit.radian ** 2)  # default stiffening k value
                lambda_1_theta_eq = 0.5 * math.pi * unit.radian  # 90 degrees
            elif angle in corrections["lambda_0"]["softened_angles"] or angle_reversed in corrections["lambda_0"]["softened_angles"]:
                soften_params = corrections["lambda_0"]["softened_angles"].get(angle, corrections["lambda_0"]["softened_angles"].get(angle_reversed))
                lambda_0_k = soften_params["k"] * unit.kilocalories_per_mole / (unit.radian ** 2)
                if "theta_eq" in soften_params:
                    lambda_0_theta_eq = soften_params["theta_eq"] * unit.radian
            # add the term to the interpolated custom angle force
            print(f"Applying ghostly angle correction for angle {angle}: "
                  f"lambda_0 k = {lambda_0_k}, theta_eq = {lambda_0_theta_eq}; "
                  f"lambda_1 k = {lambda_1_k}, theta_eq = {lambda_1_theta_eq}")
            custom_angle_force.addAngle(p1, p2, p3,
                                        [lambda_0_theta_eq, lambda_0_k,
                                        lambda_1_theta_eq, lambda_1_k])

        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)
    print(openmm.XmlSerializer.serialize(custom_angle_force))
    print(openmm.XmlSerializer.serialize(new_harmonic_angle_force))

    # process torsions
    custom_torsion_force = new_hybrid_forces["CustomTorsionForce"]
    old_hybrid_torsion_force = hybrid_forces["unique_atom_torsion_force"]
    for i in range(old_hybrid_torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = old_hybrid_torsion_force.getTorsionParameters(i)
        # check if we have one ghost atom for this torsion
        torsion = (p1, p2, p3, p4)
        if len(htf_corrected._atom_classes["unique_old_atoms"].intersection(torsion)) <= 3 or len(htf_corrected._atom_classes["unique_new_atoms"].intersection(torsion)) <= 3:
            torsion_reversed = (p4, p3, p2, p1)
            # set up containers for the end state values
            lambda_0_k = k
            lambda_1_k = k
            # check for removed dihedrals
            if torsion in corrections["lambda_0"]["removed_dihedrals"] or torsion_reversed in corrections["lambda_0"]["removed_dihedrals"]:
                lambda_0_k = 0.0 * unit.kilocalories_per_mole
            elif torsion in corrections["lambda_1"]["removed_dihedrals"] or torsion_reversed in corrections["lambda_1"]["removed_dihedrals"]:
                lambda_1_k = 0.0 * unit.kilocalories_per_mole
            # add the term to the interpolated custom torsion force
            print(f"Applying ghostly torsion correction for torsion {torsion}: "
                  f"lambda_0 k = {lambda_0_k}; "
                  f"lambda_1 k = {lambda_1_k}")
            custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                            [periodicity, phase,
                                            lambda_0_k, periodicity,
                                             phase, lambda_1_k])
        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)
    print(openmm.XmlSerializer.serialize(custom_torsion_force))
    print(openmm.XmlSerializer.serialize(new_torsion_force))

    htf_corrected._hybrid_system = new_hybrid_system
    # set the hybrid system forces dict to the new one
    htf_corrected._hybrid_system_forces = {force.getName(): force for force in new_hybrid_system.getForces()}
    return htf_corrected
