from openff.units.openmm import to_openmm, ensure_quantity, from_openmm
from itertools import chain, permutations
from htf import DevelopmentHybridTopologyFactory
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit
import openmm
from gufe import LigandAtomMapping, ProteinComponent, SolventComponent, SmallMoleculeComponent
from openfe.protocols.openmm_rfe import _rfe_utils
from openfe.protocols.openmm_utils import system_creation
import copy
import math
import logging
import json
import ast
from rdkit.Chem import Draw, AllChem
from collections import defaultdict

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

def _create_bond_lookup(htf: DevelopmentHybridTopologyFactory) -> dict:
    """
    Create a lookup dictionary of bonded atoms in the hybrid topology.

    Parameters
    ----------
    htf : DevelopmentHybridTopologyFactory
        The hybrid topology factory containing the hybrid topology.

    Returns
    -------
    dict
        A dictionary where keys are atom indices and values are sets of bonded atom indices.
    """
    # get all the bonds in the hybrid topology so we can check for an improper
    hybrid_bonds = list(htf.omm_hybrid_topology.bonds())
    # construct a lookup of bonded atoms
    bonded_atom_lookup = defaultdict(set)
    for bond in hybrid_bonds:
        a1 = bond[0].index
        a2 = bond[1].index
        bonded_atom_lookup[a1].add(a2)
        bonded_atom_lookup[a2].add(a1)
    return bonded_atom_lookup


def _find_dummy_junctions(htf: DevelopmentHybridTopologyFactory) -> dict:
    """Identify dummy-core atom junctions in the HTF and return a dictionary of their details."""
    junctions = {"lambda_0": {}, "lambda_1": {}}
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]

    # get all the bonds in the hybrid topology
    bonded_atom_lookup = _create_bond_lookup(htf=htf)

    def _collect_dummies(dummy_atoms: set[int], key: str):
        """Helper function to collect dummy-core junction atoms into the given key in junctions."""
        # track the number of junctions found
        junction_id = 0
        # track the dummy atoms already assigned to a junction
        assigned_dummies = set()
        for dummy_atom in dummy_atoms:
            # skip if already assigned as it is attached to the same core atom as another dummy
            if dummy_atom in assigned_dummies:
                continue
            bonded_physicals = [a for a in bonded_atom_lookup[dummy_atom] if a not in dummy_atoms]
            # if there are no bonded physical atoms this is part of a larger dummy group so skip
            if len(bonded_physicals) == 0:
                continue
            junction_atom = bonded_physicals[0]
            other_dummies = [a for a in bonded_atom_lookup[junction_atom] if a in dummy_atoms and a != dummy_atom]
            junctions[key][junction_id] = {
                "junction_atom": junction_atom,
                "dummies": [dummy_atom] + other_dummies,
                "physical": [a for a in bonded_atom_lookup[junction_atom] if a not in dummy_atoms],
            }
            assigned_dummies.update(junctions[key][junction_id]["dummies"])
            junction_id += 1

    _collect_dummies(dummy_new_atoms, "lambda_0")
    _collect_dummies(dummy_old_atoms, "lambda_1")
    return junctions


def _copy_hybrid_system(htf: DevelopmentHybridTopologyFactory) -> openmm.System:
    """
    Create a deep copy of the given HTF which is ready to be modified by scaling or ghostly corrections.

    Parameters
    ----------
    htf : DevelopmentHybridTopologyFactory
        The hybrid topology factory with the hybrid openmm system to copy.

    Returns
    -------
    openmm.System
        A minimal copy of the hybrid system ready for modification.

    Notes
    -----
    Things copied:
    - hybrid system particles
    - hybrid system constraints
    - barostat if present and box vectors
    """
    new_hybrid_system = openmm.System()
    # add all the particles
    for i in range(htf.hybrid_system.getNumParticles()):
        new_hybrid_system.addParticle(htf.hybrid_system.getParticleMass(i))
    # add all constraints
    for i in range(htf.hybrid_system.getNumConstraints()):
        p1, p2, dist = htf.hybrid_system.getConstraintParameters(i)
        new_hybrid_system.addConstraint(p1, p2, dist)
    # add the barostat and the box vectors
    for force in htf.hybrid_system.getForces():
        if isinstance(force, openmm.MonteCarloBarostat) or isinstance(force, openmm.MonteCarloMembraneBarostat):
            new_hybrid_system.addForce(copy.deepcopy(force))
            # also add the box vectors
            box_vectors = htf.hybrid_system.getDefaultPeriodicBoxVectors()
            new_hybrid_system.setDefaultPeriodicBoxVectors(*box_vectors)
            break
    return new_hybrid_system


def _scale_angles_and_torsions(htf: DevelopmentHybridTopologyFactory, scale_factor: float = 0.1, scale_angles: bool = True) -> DevelopmentHybridTopologyFactory:
    """
    Scale all angles and torsion force constants in the dummy-core junction by the given scale factor.

    Parameters
    ----------
    htf : DevelopmentHybridTopologyFactory
        The hybrid topology factory to modify.
    scale_factor : float, optional
        The factor by which to scale the angles and torsions (0 to 1), by default 0.1.
    scale_angles : bool, optional
        Whether to scale angles (True) or torsions (False), by default True.

    Note
    -----
    The HTF is edited inplace due to issues with deepcopying the HTF object.
    """
    assert 0 <= scale_factor <= 1, "Scale factor must be between 0 and 1."

    logger.info(f"Softening angles and torsions involving dummy atoms in the hybrid system by {(1.0 - scale_factor) * 100}%.")
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]

    softened_hybrid_system = _copy_hybrid_system(htf=htf)

    hybrid_forces = htf._hybrid_system_forces
    # copy all forces which do not need to be modified
    # We are only modifying angle and torsion forces which involve the bridge and dummy atoms
    # As the HTF stores all terms involving dummies in the standard forces we can copy all others directly
    # The interpolated forces only contain terms for the core mapped atoms so we don't need to remove any
    forces_not_to_copy = ["standard_angle_force", "unique_atom_torsion_force"]
    if not scale_angles:
        logger.info("Scaling of angles disabled. Only torsions will be softened.")
        forces_not_to_copy.remove("standard_angle_force")
    for force_name, hybrid_force in hybrid_forces.items():
        if force_name not in forces_not_to_copy:
            logger.info(f"Copying force {force_name} to new hybrid system without modification.")
            new_force = copy.deepcopy(hybrid_force)
            softened_hybrid_system.addForce(new_force)

    # now apply the softening to the angle and torsion forces
    # first add a new torsion force to the system
    softened_torsion_force = openmm.PeriodicTorsionForce()
    softened_hybrid_system.addForce(softened_torsion_force)

    # get a quick lookup of the forces
    new_hybrid_forces = {force.getName(): force for force in softened_hybrid_system.getForces()}

    # process angles
    if scale_angles:
        # if we scale angles add a new angle force to the system
        softened_harmonic_angle_force = openmm.HarmonicAngleForce()
        softened_hybrid_system.addForce(softened_harmonic_angle_force)
        logger.info("Processing dummy-core junction angles for softening.")
        logger.info("Adding softened angles to core_angle_force.")
        default_hybrid_angle_force = hybrid_forces["standard_angle_force"]
        softened_custom_angle_force = new_hybrid_forces["CustomAngleForce"]
        for i in range(default_hybrid_angle_force.getNumAngles()):
            p1, p2, p3, theta_eq, k = default_hybrid_angle_force.getAngleParameters(i)
            angle = (p1, p2, p3)
            # for the angle terms there must be at least one core atom and 1 or 2 dummy atoms
            # check lambda = 0 first
            if 1 <= len(dummy_new_atoms.intersection(angle)) < 3:
                # if we match a new unique atom the angle must be softened at lambda = 0
                # add the term to the interpolated custom angle force
                new_k = k * scale_factor
                logger.info(f"Softening angle {angle} at lambda=0: original k = {k}, new k = {new_k}")
                softened_custom_angle_force.addAngle(p1, p2, p3, [theta_eq, new_k, theta_eq, k])
            elif 1 <= len(dummy_old_atoms.intersection(angle)) < 3:
                # if we match an old unique atom the angle must be softened at lambda = 1
                # add the term to the interpolated custom angle force
                new_k = k * scale_factor
                logger.info(f"Softening angle {angle} at lambda=1: original k = {k}, new k = {new_k}")
                softened_custom_angle_force.addAngle(p1, p2, p3, [theta_eq, k, theta_eq, new_k])
            else:
                # the term does not involve any dummy atoms, so we can just copy it
                softened_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

    # process torsions
    logger.info("Processing dummy-core junction torsions for softening.")
    logger.info("Adding softened torsions to core_torsion_force.")
    default_hybrid_torsion_force = hybrid_forces["unique_atom_torsion_force"]
    softened_custom_torsion_force = new_hybrid_forces["CustomTorsionForce"]
    for i in range(default_hybrid_torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = default_hybrid_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        # for the torsion terms there must be at least one core atom and 1-3 dummy atoms
        # check lambda = 0 first
        if 1 <= len(dummy_new_atoms.intersection(torsion)) < 4:
            # if we match a new unique atom the torsion must be softened at lambda = 0
            # add the term to the interpolated custom torsion force
            new_k = k * scale_factor
            logger.info(f"Softening torsion {torsion} at lambda=0: original k = {k}, new k = {new_k}")
            softened_custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                     [periodicity, phase,
                                            new_k, periodicity,
                                             phase, k])
        elif 1 <= len(dummy_old_atoms.intersection(torsion)) < 4:
            # if we match an old unique atom the torsion must be softened at lambda = 1
            # add the term to the interpolated custom torsion force
            new_k = k * scale_factor
            logger.info(f"Softening torsion {torsion} at lambda=1: original k = {k}, new k = {new_k}")
            softened_custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                     [periodicity, phase,
                                            k, periodicity,
                                             phase, new_k])
        else:
            # the term does not involve any dummy atoms, so we can just copy it
            softened_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)

    htf._hybrid_system = softened_hybrid_system
    # set the hybrid system forces dict to the new one
    htf._hybrid_system_forces = {force.getName(): force for force in softened_hybrid_system.getForces()}
    return htf


def load_ghostly_corrections(ghostly_output_path: str) -> dict:
    """
    Parse the Ghostly modification output json file to extract the corrections to be applied to the HTF.

    Notes
    -----
    - The corrections are returned in a dictionary with the same structure as the Ghostly output but we use sets for the removed and stiffened angles/dihedrals for faster lookup.
    """
    with open(ghostly_output_path, 'r') as f:
        corrections = json.load(f)
        # convert the string keys back to tuples
        for lambda_key in corrections.keys():
            for correction_type in corrections[lambda_key].keys():
                # these are stings for some reason, convert them back to tuples
                if correction_type in ["removed_angles", "removed_dihedrals"]:
                    corrections[lambda_key][correction_type] = set([ast.literal_eval(tup_str) for tup_str in corrections[lambda_key][correction_type]])
                    # this is a list so we need to convert each to a tuple
                elif correction_type == "stiffened_angles":
                    corrections[lambda_key][correction_type] = set([tuple(angle) for angle in corrections[lambda_key][correction_type]])
                elif correction_type == "softened_angles":
                    new_dict = {}
                    for tup_str, params in corrections[lambda_key][correction_type].items():
                        tup = ast.literal_eval(tup_str)
                        new_dict[tup] = params
                    corrections[lambda_key][correction_type] = new_dict
    return corrections


def apply_ghostly_corrections(htf: DevelopmentHybridTopologyFactory, corrections: dict) -> DevelopmentHybridTopologyFactory:
    """
    Apply the ghostly corrections parsed from the output file to the HTF.

    Notes
    -----
    - The HTF is edited inplace due to issues with deepcopying the HTF object.
    - The method will track which corrections were applied and compare them to the supplied corrections.
    - The method will check that a correction is applied to all junctions involving dummy atoms identified using an internal method.

    Raises
    ------
    AssertionError
        If a parameter is changed by ghostly but we can not determine what type of correction it was.
    ValueError
        If a correction provided by ghostly is not applied to the HTF.
    """
    logger.info("Applying ghostly corrections to hybrid system.")
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]

    new_hybrid_system = _copy_hybrid_system(htf=htf)

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

    # track the applied corrections
    applied_corrections = {
        "lambda_0": {"removed_angles": set(), "stiffened_angles": set(), "softened_angles": set(), "removed_dihedrals": set()},
        "lambda_1": {"removed_angles": set(), "stiffened_angles": set(), "softened_angles": set(), "removed_dihedrals": set()}
       }

    # process angles
    custom_angle_force = new_hybrid_forces["CustomAngleForce"]
    old_hybrid_angle_force = hybrid_forces["standard_angle_force"]

    # set up the angle parameters for stiffening and zeroing
    ZERO_K = 0.0 * unit.kilocalories_per_mole / (unit.radian ** 2)
    STIFF_K = 100.0 * unit.kilocalories_per_mole / (unit.radian ** 2)
    STIFF_THETA = 0.5 * math.pi * unit.radian

    for angle_idx in range(old_hybrid_angle_force.getNumAngles()):
        p1, p2, p3, theta_eq, k = old_hybrid_angle_force.getAngleParameters(angle_idx)
        # check if we have one ghost atom for this angle
        angle = (p1, p2, p3)
        if 1 <= len(dummy_old_atoms.intersection(angle)) < 3 or 1<= len(dummy_new_atoms.intersection(angle)) < 3:
            angle_reversed = (p3, p2, p1)
            # set up containers for the end state values
            lambda_0_k = k
            lambda_0_theta_eq = theta_eq
            lambda_1_k = k
            lambda_1_theta_eq = theta_eq
            end_state, correction_type = None, None

            # check for removed angles
            if (prob_angle:= angle) in corrections["lambda_0"]["removed_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_0"]["removed_angles"]:
                lambda_0_k = ZERO_K
                end_state = 0
                correction_type = "removed_angles"
            elif (prob_angle:= angle) in corrections["lambda_1"]["removed_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_1"]["removed_angles"]:
                lambda_1_k = ZERO_K
                end_state = 1
                correction_type = "removed_angles"
            # check for stiffened angles
            elif (prob_angle:= angle) in corrections["lambda_0"]["stiffened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_0"]["stiffened_angles"]:
                lambda_0_k = STIFF_K  # default stiffening k value
                lambda_0_theta_eq = STIFF_THETA  # 90 degrees
                end_state = 0
                correction_type = "stiffened_angles"
            elif (prob_angle:= angle) in corrections["lambda_1"]["stiffened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_1"]["stiffened_angles"]:
                lambda_1_k = STIFF_K  # default stiffening k value
                lambda_1_theta_eq = STIFF_THETA  # 90 degrees
                end_state = 1
                correction_type = "stiffened_angles"
                # check for softened angles
            elif (prob_angle:= angle) in corrections["lambda_0"]["softened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_0"]["softened_angles"]:
                soften_params = corrections["lambda_0"]["softened_angles"][prob_angle]
                lambda_0_k = soften_params["k"] * unit.kilocalories_per_mole / (unit.radian ** 2)
                lambda_0_theta_eq = soften_params["theta0"] * unit.radian
                end_state = 0
                correction_type = "softened_angles"
            elif (prob_angle:= angle) in corrections["lambda_1"]["softened_angles"] or (prob_angle:= angle_reversed) in corrections["lambda_1"]["softened_angles"]:
                soften_params = corrections["lambda_1"]["softened_angles"][prob_angle]
                lambda_1_k = soften_params["k"] * unit.kilocalories_per_mole / (unit.radian ** 2)
                lambda_1_theta_eq = soften_params["theta0"] * unit.radian
                end_state = 1
                correction_type = "softened_angles"

            # some angles involving dummy atoms need to be kept to ensure 3 redundant connections
            if lambda_0_k != lambda_1_k or lambda_0_theta_eq != lambda_1_theta_eq:
                # add the term to the interpolated custom angle force
                print(f"Applying ghostly angle correction for angle {angle}: "
                      f"lambda_0 k = {lambda_0_k}, theta_eq = {lambda_0_theta_eq}; "
                      f"lambda_1 k = {lambda_1_k}, theta_eq = {lambda_1_theta_eq}")
                logger.info(f"Applying ghostly angle correction for angle {angle}: "
                      f"lambda_0 k = {lambda_0_k}, theta_eq = {lambda_0_theta_eq}; "
                      f"lambda_1 k = {lambda_1_k}, theta_eq = {lambda_1_theta_eq}")
                custom_angle_force.addAngle(p1, p2, p3,
                                            [lambda_0_theta_eq, lambda_0_k,
                                            lambda_1_theta_eq, lambda_1_k])
                # log this as a correction applied
                assert correction_type is not None, "Correction type should not be None if k or theta_eq differ!"
                applied_corrections[f"lambda_{end_state}"][correction_type].add(prob_angle)

            else:
                # both k and theta_eq values are the same, just add to the standard angle force
                new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

    # process torsions
    custom_torsion_force = new_hybrid_forces["CustomTorsionForce"]
    old_hybrid_torsion_force = hybrid_forces["unique_atom_torsion_force"]

    # set up the torsion parameters for zeroing
    TORSION_ZERO_K = 0.0 * unit.kilocalories_per_mole

    # get all the bonds in the hybrid topology so we can check for an improper
    bonded_atom_lookup = _create_bond_lookup(htf=htf)

    for torsion_idx in range(old_hybrid_torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = old_hybrid_torsion_force.getTorsionParameters(torsion_idx)
        # check if we have one ghost atom for this torsion
        torsion = (p1, p2, p3, p4)
        if 1<= len(dummy_old_atoms.intersection(torsion)) < 4 or 1<= len(dummy_new_atoms.intersection(torsion)) < 4:
            torsion_reversed = (p4, p3, p2, p1)
            # check if we have an improper torsion (central atoms bonded)
            if not (p1 in bonded_atom_lookup[p2] and p2 in bonded_atom_lookup[p3] and p3 in bonded_atom_lookup[p4]):
                # this is an improper with a dummy atom and should be skipped
                # generate all permutations of the other atoms to check for removal
                central_atom = p1
                other_atoms = [p2, p3, p4]
                torsion_variants = set()
                for perm in permutations(other_atoms):
                    torsion_variants.add((central_atom, perm[0], perm[1], perm[2]))
                    # add the reverse as well as ghostly may list either
                    torsion_variants.add((perm[2], perm[1], perm[0], central_atom))
            else:
                torsion_variants = {torsion, torsion_reversed}

            # set up containers for the end state values
            lambda_0_k = k
            lambda_1_k = k
            end_state = None

            # check for removed dihedrals
            if matched:= corrections["lambda_0"]["removed_dihedrals"].intersection(torsion_variants):
                lambda_0_k = TORSION_ZERO_K
                end_state = 0
            elif matched:= corrections["lambda_1"]["removed_dihedrals"].intersection(torsion_variants):
                lambda_1_k = TORSION_ZERO_K
                end_state = 1
            # some dihedrals involving ghost atoms need to be kept to ensure 3 redundant connections
            if lambda_0_k != lambda_1_k:
                # add the term to the interpolated custom torsion force
                print(f"Applying ghostly torsion correction for torsion {torsion}: "
                      f"lambda_0 k = {lambda_0_k}; "
                      f"lambda_1 k = {lambda_1_k}")
                logger.info(f"Applying ghostly torsion correction for torsion {torsion}: "
                      f"lambda_0 k = {lambda_0_k}; "
                      f"lambda_1 k = {lambda_1_k}")
                custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                [periodicity, phase,
                                                lambda_0_k, periodicity,
                                                 phase, lambda_1_k])
                # log this as a correction applied
                assert end_state is not None, "End state should not be None if k values differ!"
                applied_corrections[f"lambda_{end_state}"]["removed_dihedrals"].update(matched)
            else:
                # both k values are the same, just add to the standard torsion force
                new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)
        else:
            # the term does not involve any ghost atoms, so we can just copy it
            new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)


    # compare the supplied and applied corrections
    for lambda_key in corrections.keys():
        for correction_type in corrections[lambda_key].keys():
            supplied = set()
            applied = applied_corrections[lambda_key][correction_type]
            if correction_type in ["removed_angles", "stiffened_angles"]:
                supplied = corrections[lambda_key][correction_type]
            elif correction_type == "softened_angles":
                supplied = set([tuple(tup) for tup in corrections[lambda_key][correction_type].keys()])
            elif correction_type == "removed_dihedrals":
                supplied = corrections[lambda_key][correction_type]
            not_applied = supplied - applied
            if len(not_applied) > 0 and correction_type == "removed_dihedrals":
                # in some cases dihedrals are listed to be removed but are not present in the HTF these involve linear nitrile groups for example
                # check if these missed dihedrals are this type
                dummy_group = dummy_new_atoms if lambda_key == "lambda_1" else dummy_old_atoms
                for dihedral in list(not_applied):
                    # if this is a linear group torsion then atom 2 or 3 will be bonded to only 2 other atoms
                    a2_bonds = bonded_atom_lookup[dihedral[1]] - dummy_group
                    a3_bonds = bonded_atom_lookup[dihedral[2]] - dummy_group
                    if len(a2_bonds) <= 2 or len(a3_bonds) <= 2:
                        not_applied.remove(dihedral)
            if len(not_applied) > 0:
                raise ValueError(f"The following {correction_type} corrections for {lambda_key} were not applied: {not_applied}")


    htf._hybrid_system = new_hybrid_system
    # set the hybrid system forces dict to the new one
    htf._hybrid_system_forces = {force.getName(): force for force in new_hybrid_system.getForces()}
    return htf


def draw_ghostly_modifications(molecule: SmallMoleculeComponent, htf:DevelopmentHybridTopologyFactory, modifications: dict, filename: str, end_state: int):
    """
    Draw the ghostly modifications for an end state on the molecule and save to the given output path.

    Parameters
    ----------
    molecule : SmallMoleculeComponent
        The molecule to draw the modifications on, this should be one end state of the HTF.
    htf :DevelopmentHybridTopologyFactory
        The hybrid topology factory containing the hybrid topology we will use this to find the remaining terms.
    modifications : dict
        The ghostly modifications parsed from the ghostly output file use this to highlight stiffened and removed terms.
    filename : str
        The output file path to save the image to.
    end_state : int
        The end state to draw (0 or 1).

    Notes
    -----
    The output is an SVG image saved to the given filename.
    The atom highlighting colours are:
    - Light green: core atoms
    - Light grey: dummy atoms
    The valence term highlighting colours are:
    - Green: normal term crossing the junction
    - Red: stiffened term crossing the junction
    - Blue: softened term crossing the junction
    """
    # define some atom colours
    core_atom_colour = (0.7, 1.0, 0.7) # light green
    dummy_atom_colour = (0.5, 0.5, 0.5) # light grey
    normal_term_colour = (0.0, 1.0, 0.0) # green
    stiffened_term_colour = (1.0, 0.0, 0.0) # red
    softened_term_colour = (0.0, 0.0, 1.0) # blue

    # find all the dummy junctions in the HTF
    # we are plotting the molecule at end state 0/1 so we need to find the corrections at the opposite lambda to plot on that molecule
    junctions = _find_dummy_junctions(htf)[f"lambda_{1 - end_state}"]
    # generate a list of all valence forces in the hybrid system that we need to check
    forces_by_valence = {
        "constraints": [htf.hybrid_system],
        "bonds": [htf._hybrid_system_forces[force] for force in ["CustomBondForce", "HarmonicBondForce"]],
        "angles": [htf._hybrid_system_forces[force] for force in ["CustomAngleForce", "HarmonicAngleForce"]],
        "dihedrals": [htf._hybrid_system_forces[force] for force in ["CustomTorsionForce", "PeriodicTorsionForce"]],
    }

    # construct a lookup of bonded atoms
    bonded_atom_lookup = _create_bond_lookup(htf=htf)

    # get the core atoms to colour
    core_atoms = htf._atom_classes["core_atoms"]
    # get the endstate mapping to map hybrid indices to end state molecule indices
    end_state_mapping = htf._hybrid_to_new_map if end_state == 1 else htf._hybrid_to_old_map
    dummies = htf._atom_classes["unique_old_atoms"] if end_state == 0 else htf._atom_classes["unique_new_atoms"]
    # get the rdkit molecule for drawing and a 2D depiction
    rd_mol = molecule.to_rdkit()
    AllChem.Compute2DCoords(rd_mol)

    to_draw_by_junction = defaultdict(list)

    for junction_id, junction_info in junctions.items():
        junction_atom = junction_info["junction_atom"]
        junction_dummies = junction_info["dummies"]

        # find all valence terms which involve the junction atom and at least one dummy atom
        for valence_type, forces in forces_by_valence.items():
            for force in forces:
                if valence_type == "bonds":
                    for i in range(force.getNumBonds()):
                        terms = force.getBondParameters(i)
                        bond = (terms[0], terms[1])
                        if junction_atom in bond and any(dummy in bond for dummy in junction_dummies):
                            # found a bond crossing the junction
                            # now we need the atom index in the end state molecule
                            atoms_to_highlight = [end_state_mapping[atom] for atom in bond]

                            # find the bond index in the rdkit molecule for these atoms
                            bond_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight).GetIdx()

                            # add all the data for this junction type
                            to_draw_by_junction[tuple(sorted(atoms_to_highlight))].append(
                                {
                                    "legend": "Bond" + " " + repr([int(x) for x in atoms_to_highlight]),
                                    "Bonds": [bond_idx],
                                    "Bond Colours": {bond_idx: normal_term_colour},
                                    "Atoms": atoms_to_highlight,
                                }
                            )

                # H-atoms which are dummies will be in the constraints section not bonds
                elif valence_type == "constraints":
                    for i in range(force.getNumConstraints()):
                        p1, p2, _ = force.getConstraintParameters(i)
                        constraint = (p1, p2)
                        if junction_atom in constraint and any(dummy in constraint for dummy in junction_dummies):
                            # found a constraint crossing the junction
                            atoms_to_highlight = [end_state_mapping[atom] for atom in constraint]
                            bond_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight).GetIdx()
                            to_draw_by_junction[tuple(sorted(atoms_to_highlight))].append(
                                {
                                    "legend": "Constraint" + " " + repr([int(x) for x in atoms_to_highlight]),
                                    "Bonds": [bond_idx],
                                    "Bond Colours": {bond_idx: normal_term_colour},
                                    "Atoms": atoms_to_highlight,
                                }
                            )
                elif valence_type == "angles":
                    for i in range(force.getNumAngles()):
                        terms = force.getAngleParameters(i)
                        angle = tuple(terms[:3])
                        if junction_atom in angle and any(dummy in angle for dummy in junction_dummies):
                            # found an angle crossing the junction
                            # check if this angle is modified
                            if angle in modifications[f"lambda_{1 - end_state}"]["removed_angles"] or tuple(
                                    reversed(angle)) in modifications[f"lambda_{1 - end_state}"]["removed_angles"]:
                                continue  # skip drawing removed angles

                            atoms_to_highlight = [end_state_mapping[atom] for atom in angle]
                            # workout the junction bond atoms
                            junction_bond = [atom for atom in angle if atom == junction_atom or atom in junction_dummies]
                            # this should always be 2 atoms unless it is a dummy only angle which spans a core atom
                            if len(junction_bond) == 3:
                                junction_bond = junction_bond[:2]
                            assert len(junction_bond) == 2, "Junction bond in angle should be 2 atoms."
                            # convert to end state indices
                            junction_bond = tuple(sorted([end_state_mapping[atom] for atom in junction_bond]))
                            # find the bonds for this angle
                            bond1_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight[:2]).GetIdx()
                            bond2_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight[1:]).GetIdx()

                            # colour depending on the modification type
                            if angle in modifications[f"lambda_{1 - end_state}"]["stiffened_angles"] or tuple(reversed(angle)) in modifications[f"lambda_{1 - end_state}"]["stiffened_angles"]:
                                bond_colour = {bond_idx: stiffened_term_colour for bond_idx in [bond1_idx, bond2_idx]}
                                legend = "Stiffened Angle" + " " + repr([int(x) for x in atoms_to_highlight])
                            elif angle in modifications[f"lambda_{1 - end_state}"]["softened_angles"] or tuple(reversed(angle)) in modifications[f"lambda_{1 - end_state}"]["softened_angles"]:
                                bond_colour = {bond_idx: softened_term_colour for bond_idx in [bond1_idx, bond2_idx]}
                                legend = "Softened Angle" + " " + repr([int(x) for x in atoms_to_highlight])
                            else:
                                bond_colour = {bond_idx: normal_term_colour for bond_idx in [bond1_idx, bond2_idx]}
                                legend = "Angle" + " " + repr([int(x) for x in atoms_to_highlight])
                            to_draw_by_junction[junction_bond].append(
                                {
                                    "legend": legend,
                                    "Bonds": [bond1_idx, bond2_idx],
                                    "Bond Colours": bond_colour,
                                    "Atoms": atoms_to_highlight,
                                }
                            )
                else:
                    for i in range(force.getNumTorsions()):
                        terms = force.getTorsionParameters(i)
                        torsion = tuple(terms[:4])
                        if junction_atom in torsion and any(dummy in torsion for dummy in junction_dummies):
                            # check if we have an improper torsion
                            if not (torsion[0] in bonded_atom_lookup[torsion[1]] and torsion[1] in bonded_atom_lookup[torsion[2]] and torsion[2] in bonded_atom_lookup[torsion[3]]):
                                # this is an improper with a dummy atom and should be skipped
                                # generate all permutations of the other atoms to check for removal
                                central_atom = torsion[0]
                                other_atoms = [torsion[1], torsion[2], torsion[3]]
                                torsion_variants = []
                                for perm in permutations(other_atoms):
                                    torsion_variants.append((central_atom, perm[0], perm[1], perm[2]))
                            else:
                                torsion_variants = [torsion, tuple(reversed(torsion))]
                            # check if this torsion is modified
                            if any(tor in modifications[f"lambda_{1 - end_state}"]["removed_dihedrals"] for tor in torsion_variants):
                                continue  # skip drawing removed dihedrals
                            # found a torsion crossing the junction
                            atoms_to_highlight = tuple([end_state_mapping[atom] for atom in torsion])
                            # find the junction atoms for this torsion
                            junction_bond = [atom for atom in torsion if atom == junction_atom or atom in junction_dummies]
                            # this should always be 2 atoms
                            assert len(junction_bond) == 2, "Junction bond in torsion should be 2 atoms."
                            # convert to end state indices
                            junction_bond = tuple(sorted([end_state_mapping[atom] for atom in junction_bond]))
                            try:
                                bond1_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight[:2]).GetIdx()
                                bond2_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight[1:3]).GetIdx()
                                bond3_idx = rd_mol.GetBondBetweenAtoms(*atoms_to_highlight[2:]).GetIdx()
                                torsion_label = "Dihedral"
                            except AttributeError:
                                # this is probably an improper dihedral with the central atom listed first so change how we get the bonds
                                bond1_idx = rd_mol.GetBondBetweenAtoms(atoms_to_highlight[1], atoms_to_highlight[0]).GetIdx()
                                bond2_idx = rd_mol.GetBondBetweenAtoms(atoms_to_highlight[0], atoms_to_highlight[2]).GetIdx()
                                bond3_idx = rd_mol.GetBondBetweenAtoms(atoms_to_highlight[0], atoms_to_highlight[3]).GetIdx()
                                torsion_label = "Improper"
                                # if we have an improper we want to list the central atom first and sort the others
                                other_atoms = sorted(atoms_to_highlight[1:])
                                atoms_to_highlight = (atoms_to_highlight[0],) + tuple(other_atoms)

                            # due to the way torsions are stored in openmm with perodicities stored in different terms we may find the same torsion multiple times
                            # to avoid duplicates we check if we have already added this torsion
                            for highlight_data in to_draw_by_junction[junction_bond]:
                                if atoms_to_highlight == highlight_data["Atoms"]:
                                    continue

                            to_draw_by_junction[junction_bond].append(
                                {
                                    "legend": torsion_label + " " + repr([int(x) for x in atoms_to_highlight]),
                                    "Bonds": [bond1_idx, bond2_idx, bond3_idx],
                                    "Bond Colours": {bond_idx: normal_term_colour for bond_idx in [bond1_idx, bond2_idx, bond3_idx]},
                                    "Atoms": atoms_to_highlight,
                                }
                            )


    # set up containers for drawing
    molecules = []
    highlight_atoms = []
    highlight_bonds = []
    bond_colours = []
    atom_colours = []
    legends = []

    # set up the atom colours
    atom_to_color = {end_state_mapping[idx]: core_atom_colour for idx in core_atoms}  # set the junction atom colour
    for dummy in dummies:
        atom_to_color[end_state_mapping[dummy]] = dummy_atom_colour  # set the dummy atom colours on all dummies

    # add a copy of the molecule with the atoms involved in the term annotated with their map numbers
    for highlight_junction in to_draw_by_junction.values():
        for highlight_data in highlight_junction:
            highlight_atoms.append(highlight_data["Atoms"])
            highlight_bonds.append(highlight_data["Bonds"])
            bond_colours.append(highlight_data["Bond Colours"])
            atom_colours.append(atom_to_color)
            legends.append(highlight_data["legend"])
            mol_copy = copy.deepcopy(rd_mol)
            atoms = highlight_data["Atoms"]
            for atom_idx in atoms:
                atom = mol_copy.GetAtomWithIdx(atom_idx)
                atom.SetProp("molAtomMapNumber", str(atom_idx))
            molecules.append(mol_copy)



    res = Draw.MolsToGridImage(
        molecules,
        molsPerRow=4,
        highlightAtomLists=highlight_atoms,
        highlightBondLists=highlight_bonds,
        highlightBondColors=bond_colours,
        highlightAtomColors=atom_colours,
        subImgSize=(350, 350),
        legends=legends,
        useSVG=True,
    )
    with open(filename, "w") as fh:
        fh.write(res)
