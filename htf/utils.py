from openff.units.openmm import to_openmm, ensure_quantity, from_openmm
from itertools import chain
from htf import DevelopmentHybridTopologyFactory
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit
import openmm
from gufe import LigandAtomMapping, ProteinComponent, SolventComponent
from openfe.protocols.openmm_rfe import _rfe_utils
from openfe.protocols.openmm_utils import system_creation
import copy
import math
import logging
import json
import ast

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

def load_ghostly_corrections(ghostly_output_path: str) -> dict:
    """Parse the Ghostly modification output json file to extract the corrections to be applied to the HTF."""
    with open(ghostly_output_path, 'r') as f:
        corrections = json.load(f)
        # convert the string keys back to tuples
        for lambda_key in corrections.keys():
            for correction_type in corrections[lambda_key].keys():
                # these are stings for some reason, convert them back to tuples
                if correction_type in ["removed_angles", "removed_dihedrals"]:
                    corrections[lambda_key][correction_type] = [ast.literal_eval(tup_str) for tup_str in corrections[lambda_key][correction_type]]
                    # this is a list so we need to convert each to a tuple
                elif correction_type == "stiffened_angles":
                    corrections[lambda_key][correction_type] = [tuple(angle) for angle in corrections[lambda_key][correction_type]]
                elif correction_type == "softened_angles":
                    new_dict = {}
                    for tup_str, params in corrections[lambda_key][correction_type].items():
                        tup = ast.literal_eval(tup_str)
                        new_dict[tup] = params
                    corrections[lambda_key][correction_type] = new_dict
    return corrections


def apply_ghostly_corrections(htf: DevelopmentHybridTopologyFactory, corrections: dict) -> DevelopmentHybridTopologyFactory:
    """Apply the ghostly corrections parsed from the output file to the HTF."""

    htf_corrected = copy.deepcopy(htf)
    dummy_old_atoms = htf_corrected._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf_corrected._atom_classes["unique_new_atoms"]

    new_hybrid_system = openmm.System()
    # add all the particles
    for i in range(htf_corrected.hybrid_system.getNumParticles()):
        new_hybrid_system.addParticle(htf_corrected.hybrid_system.getParticleMass(i))
    # add all constraints
    for i in range(htf_corrected.hybrid_system.getNumConstraints()):
        p1, p2, dist = htf_corrected.hybrid_system.getConstraintParameters(i)
        new_hybrid_system.addConstraint(p1, p2, dist)

    hybrid_forces = htf._hybrid_system_forces
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
    new_hybrid_system.addForce(new_harmonic_angle_force)
    new_torsion_force = openmm.PeriodicTorsionForce()
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
        if 1<= len(dummy_old_atoms.intersection(angle)) < 3 or 1<= len(dummy_new_atoms.intersection(angle)) < 3:
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
                # check for softened angles
            elif (prob_angle:= angle) in corrections["lambda_0"]["softened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_0"]["softened_angles"]:
                soften_params = corrections["lambda_0"]["softened_angles"][prob_angle]
                lambda_0_k = soften_params["k"] * unit.kilocalories_per_mole / (unit.radian ** 2)
                lambda_0_theta_eq = soften_params["theta0"] * unit.radian
            elif (prob_angle:= angle) in corrections["lambda_1"]["softened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_1"]["softened_angles"]:
                soften_params = corrections["lambda_1"]["softened_angles"][prob_angle]
                lambda_1_k = soften_params["k"] * unit.kilocalories_per_mole / (unit.radian ** 2)
                lambda_1_theta_eq = soften_params["theta0"] * unit.radian

            # some angles involving dummy atoms need to be kept to ensure 3 redundant connections
            if lambda_0_k != lambda_1_k or lambda_0_theta_eq != lambda_1_theta_eq:
                # add the term to the interpolated custom angle force
                print(f"Applying ghostly angle correction for angle {angle}: "
                      f"lambda_0 k = {lambda_0_k}, theta_eq = {lambda_0_theta_eq}; "
                      f"lambda_1 k = {lambda_1_k}, theta_eq = {lambda_1_theta_eq}")
                custom_angle_force.addAngle(p1, p2, p3,
                                            [lambda_0_theta_eq, lambda_0_k,
                                            lambda_1_theta_eq, lambda_1_k])
            else:
                # both k and theta_eq values are the same, just add to the standard angle force
                new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

    # process torsions
    custom_torsion_force = new_hybrid_forces["CustomTorsionForce"]
    old_hybrid_torsion_force = hybrid_forces["unique_atom_torsion_force"]
    for i in range(old_hybrid_torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = old_hybrid_torsion_force.getTorsionParameters(i)
        # check if we have one ghost atom for this torsion
        torsion = (p1, p2, p3, p4)
        if 1<= len(dummy_old_atoms.intersection(torsion)) < 4 or 1<= len(dummy_new_atoms.intersection(torsion)) < 4:
            torsion_reversed = (p4, p3, p2, p1)
            # set up containers for the end state values
            lambda_0_k = k
            lambda_1_k = k
            # check for removed dihedrals
            if torsion in corrections["lambda_0"]["removed_dihedrals"] or torsion_reversed in corrections["lambda_0"]["removed_dihedrals"]:
                lambda_0_k = 0.0 * unit.kilocalories_per_mole
            elif torsion in corrections["lambda_1"]["removed_dihedrals"] or torsion_reversed in corrections["lambda_1"]["removed_dihedrals"]:
                lambda_1_k = 0.0 * unit.kilocalories_per_mole
            # some dihedrals involving ghost atoms need to be kept to ensure 3 redundant connections
            if lambda_0_k != lambda_1_k:
                # add the term to the interpolated custom torsion force
                print(f"Applying ghostly torsion correction for torsion {torsion}: "
                      f"lambda_0 k = {lambda_0_k}; "
                      f"lambda_1 k = {lambda_1_k}")
                custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                [periodicity, phase,
                                                lambda_0_k, periodicity,
                                                 phase, lambda_1_k])
            else:
                # both k values are the same, just add to the standard torsion force
                new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)
        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)

    htf_corrected._hybrid_system = new_hybrid_system
    # set the hybrid system forces dict to the new one
    htf_corrected._hybrid_system_forces = {force.getName(): force for force in new_hybrid_system.getForces()}
    return htf_corrected
