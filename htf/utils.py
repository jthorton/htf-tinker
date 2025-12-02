from openff.units.openmm import to_openmm, ensure_quantity, from_openmm
from itertools import chain
from htf import DevelopmentHybridTopologyFactory
from openmmforcefields.generators import SystemGenerator
from openmm import app, unit
import openmm
from gufe import SmallMoleculeComponent, LigandAtomMapping, ProteinComponent, SolventComponent
from openfe.protocols.openmm_rfe import _rfe_utils
from openfe.protocols.openmm_utils import system_creation
import copy
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

def _find_dummy_junctions(htf: DevelopmentHybridTopologyFactory) -> dict:
    """Identify dummy-core atom junctions in the HTF and return a dictionary of their details."""
    junctions = {"lambda_0": {}, "lambda_1": {}}
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]
    junction_id = 0
    # get all the bonds in the hybrid topology
    hybrid_bonds = list(htf.omm_hybrid_topology.bonds())
    # construct a lookup of bonded atoms
    bonded_atom_lookup = {}
    for bond in hybrid_bonds:
        a1 = bond[0].index
        a2 = bond[1].index
        if a1 not in bonded_atom_lookup:
            bonded_atom_lookup[a1] = set()
        if a2 not in bonded_atom_lookup:
            bonded_atom_lookup[a2] = set()
        bonded_atom_lookup[a1].add(a2)
        bonded_atom_lookup[a2].add(a1)

    # find bridges for lambda_0, these are the unique new atoms which should not be present at lambda=0
    for dummy_atom in dummy_new_atoms:
        # find bonded physical atoms
        bonded_physicals = [a for a in bonded_atom_lookup[dummy_atom] if a not in dummy_new_atoms]
        if len(bonded_physicals) == 1:
            # this is a bridge atom
            junction_atom = bonded_physicals[0]
            # find other dummy atoms bonded to the same physical atom
            other_dummies = [a for a in bonded_atom_lookup[junction_atom] if a in dummy_new_atoms and a != dummy_atom]
            junctions["lambda_0"][junction_id] = {
                "junction_atom": junction_atom,
                "dummies": [dummy_atom] + other_dummies,
                "physical": [a for a in bonded_atom_lookup[junction_atom] if a not in dummy_new_atoms],
            }
            junction_id += 1
    # find bridges for lambda_1, these are the unique old atoms which should not be present at lambda=1
    junction_id = 0
    for dummy_atom in dummy_old_atoms:
        # find bonded physical atoms
        bonded_physicals = [a for a in bonded_atom_lookup[dummy_atom] if a not in dummy_old_atoms]
        if len(bonded_physicals) == 1:
            # this is a bridge atom
            junction_atom = bonded_physicals[0]
            # find other dummy atoms bonded to the same physical atom
            other_dummies = [a for a in bonded_atom_lookup[junction_atom] if a in dummy_old_atoms and a != dummy_atom]
            junctions["lambda_1"][junction_id] = {
                "junction_atom": junction_atom,
                "dummies": [dummy_atom] + other_dummies,
                "physical": [a for a in bonded_atom_lookup[junction_atom] if a not in dummy_old_atoms],
            }
            junction_id += 1
    return junctions


def _scale_angles_and_torsions(htf: DevelopmentHybridTopologyFactory, scale_factor: float = 0.1) -> DevelopmentHybridTopologyFactory:
    """
    Scale all angles and torsion force constants in the dummy-core junction by the given scale factor.

    Parameters
    ----------
    htf : DevelopmentHybridTopologyFactory
        The hybrid topology factory to modify.
    scale_factor : float, optional
        The factor by which to scale the angles and torsions (0 to 1), by default 0.1.

    Returns
    -------
    DevelopmentHybridTopologyFactory
        A new HTF with softened angles and torsions crossing the dummy core junctions in the hybrid system.
    """
    assert 0 <= scale_factor <= 1, "Scale factor must be between 0 and 1."

    logger.info(f"Softening angles and torsions involving dummy atoms in the hybrid system by {(1.0 - scale_factor) * 100}%.")
    htf_softened = copy.deepcopy(htf)
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]

    new_hybrid_system = openmm.System()
    # add all the particles
    logger.info("Copying particles and constraints to new hybrid system.")
    for i in range(htf.hybrid_system.getNumParticles()):
        new_hybrid_system.addParticle(htf.hybrid_system.getParticleMass(i))
    # add all constraints
    for i in range(htf.hybrid_system.getNumConstraints()):
        p1, p2, dist = htf.hybrid_system.getConstraintParameters(i)
        new_hybrid_system.addConstraint(p1, p2, dist)

    hybrid_forces = htf._hybrid_system_forces
    # copy all forces which do not need to be modified
    # We are only modifying angle and torsion forces which involve the bridge and dummy atoms
    # As the HTF stores all terms involving dummies in the standard forces we can copy all others directly
    # The interpolated forces only contain terms for the core mapped atoms so we don't need to remove any
    forces_not_to_copy = ["standard_angle_force", "unique_atom_torsion_force"]
    for force_name, hybrid_force in hybrid_forces.items():
        if force_name not in forces_not_to_copy:
            logger.info(f"Copying force {force_name} to new hybrid system without modification.")
            new_force = copy.deepcopy(hybrid_force)
            new_hybrid_system.addForce(new_force)

    # now apply the softening to the angle and torsion forces
    # first add a new standard angle and torsion force to the system
    new_harmonic_angle_force = openmm.HarmonicAngleForce()
    new_torsion_force = openmm.PeriodicTorsionForce()
    new_hybrid_system.addForce(new_harmonic_angle_force)
    new_hybrid_system.addForce(new_torsion_force)

    # get a quick lookup of the forces
    new_hybrid_forces = {force.getName(): force for force in new_hybrid_system.getForces()}

    # process angles
    logger.info("Processing dummy-core junction angles for softening.")
    old_hybrid_angle_force = hybrid_forces["standard_angle_force"]
    new_custom_angle_force = new_hybrid_forces["CustomAngleForce"]
    for i in range(old_hybrid_angle_force.getNumAngles()):
        p1, p2, p3, theta_eq, k = old_hybrid_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # for the angle terms there must be at least one core atom and 1 or 2 dummy atoms
        # check lambda = 0 first
        if 1 <= len(dummy_new_atoms.intersection(angle)) < 3:
            # if we match a new unique atom the angle must be softened at lambda = 0
            # add the term to the interpolated custom angle force
            new_k = k * scale_factor
            logger.info(f"Softening angle {angle} at lambda=0: original k = {k}, new k = {new_k}")
            new_custom_angle_force.addAngle(p1, p2, p3, [theta_eq, new_k, theta_eq, k])
        elif 1 <= len(dummy_old_atoms.intersection(angle)) < 3:
            # if we match an old unique atom the angle must be softened at lambda = 1
            # add the term to the interpolated custom angle force
            new_k = k * scale_factor
            logger.info(f"Softening angle {angle} at lambda=1: original k = {k}, new k = {new_k}")
            new_custom_angle_force.addAngle(p1, p2, p3, [theta_eq, k, theta_eq, new_k])
        else:
            # the term does not involve any dummy atoms, so we can just copy it
            new_harmonic_angle_force.addAngle(p1, p2, p3, theta_eq, k)

    # process torsions
    logger.info("Processing dummy-core junction torsions for softening.")
    old_hybrid_torsion_force = hybrid_forces["unique_atom_torsion_force"]
    new_custom_torsion_force = new_hybrid_forces["CustomTorsionForce"]
    for i in range(old_hybrid_torsion_force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = old_hybrid_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        # for the torsion terms there must be at least one core atom and 1-3 dummy atoms
        # check lambda = 0 first
        if 1 <= len(dummy_new_atoms.intersection(torsion)) < 4:
            # if we match a new unique atom the torsion must be softened at lambda = 0
            # add the term to the interpolated custom torsion force
            new_k = k * scale_factor
            logger.info(f"Softening torsion {torsion} at lambda=0: original k = {k}, new k = {new_k}")
            new_custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                [periodicity, phase,
                                            new_k, periodicity,
                                             phase, k])
        elif 1 <= len(dummy_old_atoms.intersection(torsion)) < 4:
            # if we match an old unique atom the torsion must be softened at lambda = 1
            # add the term to the interpolated custom torsion force
            new_k = k * scale_factor
            logger.info(f"Softening torsion {torsion} at lambda=1: original k = {k}, new k = {new_k}")
            new_custom_torsion_force.addTorsion(p1, p2, p3, p4,
                                                [periodicity, phase,
                                            k, periodicity,
                                             phase, new_k])
        else:
            # the term does not involve any dummy atoms, so we can just copy it
            new_torsion_force.addTorsion(p1, p2, p3, p4, periodicity, phase, k)

    htf_softened._hybrid_system = new_hybrid_system
    # set the hybrid system forces dict to the new one
    htf_softened._hybrid_system_forces = {force.getName(): force for force in new_hybrid_system.getForces()}
    return htf_softened