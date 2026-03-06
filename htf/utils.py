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
from collections import defaultdict
from rdkit import Chem
from typing import Iterable

logger = logging.getLogger(__name__)

def make_htf(mapping: LigandAtomMapping, settings, protein: ProteinComponent = None, solvent: SolventComponent = None, corrections = None) -> DevelopmentHybridTopologyFactory:
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
        valence_correction_terms=corrections
    )

def _find_dummy_junctions(htf: DevelopmentHybridTopologyFactory) -> dict:
    """Identify dummy-core atom junctions in the HTF and return a dictionary of their details."""
    junctions = {"lambda_0": {}, "lambda_1": {}}
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]

    # get all the bonds in the hybrid topology
    hybrid_bonds = list(htf.omm_hybrid_topology.bonds())
    # construct a lookup of bonded atoms
    bonded_atom_lookup = defaultdict(set)
    for bond in hybrid_bonds:
        a1 = bond[0].index
        a2 = bond[1].index
        bonded_atom_lookup[a1].add(a2)
        bonded_atom_lookup[a2].add(a1)


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

    softened_hybrid_system = openmm.System()
    # add all the particles
    logger.info("Copying particles and constraints to new hybrid system.")
    for i in range(htf.hybrid_system.getNumParticles()):
        softened_hybrid_system.addParticle(htf.hybrid_system.getParticleMass(i))
    # add all constraints
    for i in range(htf.hybrid_system.getNumConstraints()):
        p1, p2, dist = htf.hybrid_system.getConstraintParameters(i)
        softened_hybrid_system.addConstraint(p1, p2, dist)

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

    htf_softened._hybrid_system = softened_hybrid_system
    # set the hybrid system forces dict to the new one
    htf_softened._hybrid_system_forces = {force.getName(): force for force in softened_hybrid_system.getForces()}
    return htf_softened


def _derive_dummy_junction_corrections(mapping: LigandAtomMapping, force_field: str) -> dict[str, dict[str, set[tuple[int, ...]]]]:
    """
    For the given mapping and forcefield, derive dummy junction corrections based on the best practices of Fleck at al.

    Notes
    -----
    Currently only works for smirnoff force fields
    """
    from openff.toolkit import ForceField
    all_corrections = {
        "lambda_0": _get_correction_dict(),
        "lambda_1": _get_correction_dict()
    }
    # load the force field to label the molecules with parameters
    if not force_field.endswith(".offxml"):
        force_field = force_field + ".offxml"
    ff = ForceField(force_field)


    for state_key in ["lambda_0", "lambda_1"]:
        print(f"Finding corrections for {state_key}")
        # the corrections for the end state are derived from the ligand at the opposite end state
        if state_key == "lambda_0":
            smc = mapping.componentB
            core_atoms = set(mapping.componentA_to_componentB.values())
        else:
            smc = mapping.componentA
            core_atoms = set(mapping.componentA_to_componentB.keys())
        print(f"Ligand at this state is: {smc} with core atoms {core_atoms}")
        rdkit_mol = smc.to_rdkit()
        dummy_atoms = set([a.GetIdx() for a in rdkit_mol.GetAtoms() if a.GetIdx() not in core_atoms])
        ff_labels = ff.label_molecules(smc.to_openff().to_topology())[0]
        rotor_atoms = _find_free_rotors(rdkit_mol)
        print(f"Found rotor atoms: {rotor_atoms}")


        # make a lookup for the bonds
        bond_look_up = defaultdict(list)
        for bond in rdkit_mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            bond_look_up[a1].append(a2)
            bond_look_up[a2].append(a1)

        # now find the dummy core junctions
        assigned_dummies = set()
        dummy_junctions = []
        for dummy_atom in dummy_atoms:
            if dummy_atom in assigned_dummies:
                continue
            bonded_physicals = [a for a in bond_look_up[dummy_atom] if a not in dummy_atoms]
            # if there are no bonded physical atoms this is part of a larger dummy group so skip
            if len(bonded_physicals) == 0:
                continue
            junction_atom = bonded_physicals[0]
            other_dummies = [a for a in bond_look_up[junction_atom] if a in dummy_atoms and a != dummy_atom]
            dummy_junctions.append({
                "junction_atom": junction_atom,
                "dummies": set([dummy_atom] + other_dummies),
                "physical": set([a for a in bond_look_up[junction_atom] if a not in dummy_atoms]),
            })
            assigned_dummies.update([dummy_atom] + other_dummies)

        # now derive the corrections for each junction
        for junction in dummy_junctions:
            # dispatch the junction to the correction
            physicals = len(junction["physical"])
            match physicals:
                case 1:
                    print(f"Deriving corrections for terminal junction with physical atom {junction['physical']} and dummy atoms {junction['dummies']}")
                    corrections = _derive_terminal_corrections(junction, rdkit_mol, ff_labels, core_atoms, rotor_atoms)
                    print(f"Derived corrections for terminal junction: {corrections}")
                case 2:
                    print(f"Deriving corrections for dual junction with physical atoms {junction['physical']} and dummy atoms {junction['dummies']}")
                    corrections = _derive_dual_corrections(junction, rdkit_mol, ff_labels, core_atoms, rotor_atoms)
                    print(f"Derived corrections for dual junction: {corrections}")
                case 3:
                    print(f"Deriving corrections for triple junction with physical atoms {junction['physical']} and dummy atoms {junction['dummies']}")
                    corrections = _derive_triple_corrections(junction, rdkit_mol, ff_labels, core_atoms, rotor_atoms)
                    print(f"Derived corrections for triple junctions: {corrections}")
                case _:
                    raise NotImplementedError(f"Higher order junctions not currently supported junction is order {physicals}.")

            # validate the corrections
            if not _check_dummy_junction(junction=junction, corrections=corrections, force_field_labels=ff_labels, core_atoms=core_atoms):
                raise ValueError(f"Corrections failed validation for junction with physical atoms {junction['physical']} and dummy atoms {junction['dummies']}. Corrections were: {corrections}")
            # add these corrections to the overall corrections dict
            for correction_type, terms in corrections.items():
                all_corrections[state_key][correction_type].update(terms)

            # now remove all improper torsions involving a dummy and physical atom
            corrections = _find_improper_corrections(junction, ff_labels)
            all_corrections[state_key]["removed_impropers"].update(corrections["removed_impropers"])
    return all_corrections



def _find_free_rotors(rdkit_mold) -> set[int]:
    """
    Get the atom indices of any free rotor terminal atoms in the molecule.

    Note
    ----
    We define free rotors as any terminal atoms in the following groups:
        - CH3
        - NH2
        - OH
        - SH
    """
    from rdkit import Chem
    rotor_indices = set()
    # smarts wrote to match the terminal atom and the connected atom
    for rotor_smarts in ["[H]-[CX4H3]", "[H]-[NX3H2]", "[H]-[OX2H]", "[H]-[SX2H]"]:
        rotor_query = Chem.MolFromSmarts(rotor_smarts)
        matches = rdkit_mold.GetSubstructMatches(rotor_query)
        # get the terminal atom in each match and add to the set
        for match in matches:
            terminal_atoms = [a for a in match if rdkit_mold.GetAtomWithIdx(a).GetDegree() == 1]
            rotor_indices.update(terminal_atoms)
    return rotor_indices

def _get_heaviest_terminal_atom(dihedrals: Iterable[tuple[int, int, int, int]], excluded_atoms: set[int], rdkit_mol: Chem.Mol) -> int:
    """
    Get the heaviest terminal atom from a list of dihedrals, excluding any atoms in the excluded_atoms set.

    Notes
    -----
    - This is used to determine which dihedral to keep in many junctions,
        as we want to keep the one which terminates in the heaviest atom to minimize the impact of the constraint.
    """
    terminal_atoms_with_mass = []
    for dihedral in dihedrals:
        terminal_atom = [a for a in dihedral if a not in excluded_atoms][0]
        mass = rdkit_mol.GetAtomWithIdx(terminal_atom).GetMass()
        terminal_atoms_with_mass.append((terminal_atom, mass))
    # sort by mass then index, ensure deterministic if multiple terminal atoms have the same mass by sorting by index as a tiebreaker
    terminal_atoms_with_mass.sort(key=lambda x: (x[1], x[0]), reverse=True)
    heaviest_terminal_atom = terminal_atoms_with_mass[0][0]
    return heaviest_terminal_atom


def _derive_terminal_corrections(junction: dict[str, set[int]], rdkit_mol, force_field_labels, core_atoms, free_rotors) -> dict[str, set[frozenset[int]]]:
    """
    Derive the corrections for a terminal dummy junction as per the best practices of Fleck et al.

    Notes
    -----
    - This results in a single bond, angle and dihedral term being kept per terminal junction
    - Terms between terminal junction ghost atoms are all kept
    - All dihedrals will terminate in a single heavy atom to avoid single and dual anchor dihedral constraints
    - If the dihedral terminates in a free rotor we stiffen the dihedral and change the periodicity to 1 and the phase to 0 to
        try and uncouple the dummy from the rapid rotations of the free rotor.
    """
    # find all dihedrals which involve a dummy atom and the physical junction atom
    junction_atom = junction["junction_atom"]
    dummy_atoms = junction["dummies"]
    physical_atoms = junction["physical"]
    dummy_dihedrals = set()
    for dihedral in force_field_labels["ProperTorsions"].keys():
        if junction_atom in dihedral and dummy_atoms.intersection(dihedral) and physical_atoms.intersection(dihedral):
            dummy_dihedrals.add(dihedral)

    # now we need to determine which dihedral to keep,
    # select the one which terminates in the heaviest atom to minimize the impact of the constraint
    # get the set of atoms to ignore when finding the terminal atom
    excluded_atoms = set(dummy_atoms) | set(physical_atoms) | {junction_atom}
    heaviest_terminal_atom = _get_heaviest_terminal_atom(
        dihedrals=dummy_dihedrals, excluded_atoms=excluded_atoms, rdkit_mol=rdkit_mol
    )
    corrections = _get_correction_dict()
    # remove all dihedrals which do not terminate in the heaviest terminal atom
    for dihedral in dummy_dihedrals:
        if heaviest_terminal_atom not in dihedral:
            corrections["removed_dihedrals"].add(frozenset(dihedral))
        elif heaviest_terminal_atom in dihedral and heaviest_terminal_atom in free_rotors:
            # if the dihedral terminates in a free rotor we stiffen the dihedral to
            # try and uncouple the dummy from the rapid rotations of the free rotor.
            corrections["stiffened_dihedrals"].add(frozenset(dihedral))
    return corrections


def _derive_dual_corrections(junction: dict[str, set[int]], rdkit_mol, force_field_labels, core_atoms, rotor_atoms) -> dict[str, set[tuple[int, ...]]]:
    """
    Derive the corrections for a dual dummy junction as per the best practices of Fleck et al.

    Notes
    -----
    - In the case of a dual junction with two dummy groups we do not remove any valence terms between them as we use the
        single bond-angle-dihedral method to avoid the groups overlapping.
    - If the physical junction atom is planar we stiffen the dihedral we retain to prevent the junction flapping, we use a k value
        of 100 kcal/mol a periodicity of 1 and a phase of PI to se the equilibrium value to 0.
    """
    # find all dihedrals which involve a dummy atom and the physical junction atom
    junction_atom = junction["junction_atom"]
    dummy_atoms = junction["dummies"]
    physical_atoms = junction["physical"]
    # find torsions involving dummy atoms
    # track those which terminate in the dummy junction atom (bonded to the physical junction atom)
    dummy_junction_dihedrals = set()
    # also track dihedrals which terminate at one of the physical atoms and originate from within one of the dummy groups
    dummy_group_dihedrals = set()
    other_core_atoms = core_atoms - physical_atoms - {junction_atom}
    for dihedral in force_field_labels["ProperTorsions"].keys():
        # do checks to classify the dihedral
        junction_in_dihedral = junction_atom in dihedral
        dummy_in_dihedral = dummy_atoms.intersection(dihedral)
        physical_in_dihedral = physical_atoms.intersection(dihedral)
        other_core_in_dihedral = other_core_atoms.intersection(dihedral)

        # check if the dihedral terminates at the dummy junction atom and involves no other dummy atoms
        if junction_in_dihedral and dummy_in_dihedral and physical_in_dihedral and other_core_in_dihedral:
            dummy_junction_dihedrals.add(dihedral)
        # check if the dihedral originates from the dummy group but not the dummy junction atom
        elif physical_in_dihedral and dummy_in_dihedral and junction_in_dihedral and not other_core_in_dihedral:
            dummy_group_dihedrals.add(dihedral)
    corrections = _get_correction_dict()

    # find the dihedral to terminate in by mass
    excluded_atoms = set(dummy_atoms) | set(physical_atoms) | {junction_atom}
    heaviest_terminal_atom = _get_heaviest_terminal_atom(
        dihedrals=dummy_junction_dihedrals, excluded_atoms=excluded_atoms, rdkit_mol=rdkit_mol
    )

    # we need to check if the junction is planar and if we need to stiffen the dihedral
    junction_hybridisation = rdkit_mol.GetAtomWithIdx(junction_atom).GetHybridization()
    if junction_hybridisation == Chem.HybridizationType.SP2:
        junction_is_planar = True
    else:
        junction_is_planar = False

    # remove all dihedrals which do not terminate in the heaviest terminal atom
    for dihedral in dummy_junction_dihedrals:
        if heaviest_terminal_atom not in dihedral:
            corrections["removed_dihedrals"].add(frozenset(dihedral))
        # we need to stiffen the dihedral if we want to keep the dummy group in plan with the rest of the molecule
        # or if we have to anchor using a free rotor
        elif heaviest_terminal_atom in dihedral and (heaviest_terminal_atom in rotor_atoms or junction_is_planar):
            corrections["stiffened_dihedrals"].add(frozenset(dihedral))

    # find the physical atom the kept dihedrals pass through
    heavy_atom = rdkit_mol.GetAtomWithIdx(heaviest_terminal_atom)
    physical_atom_to_keep = [a.GetIdx() for a in heavy_atom.GetNeighbors() if a.GetIdx() in physical_atoms][0]

    # remove redundant angles not running through the kept dihedral
    for angle in force_field_labels["Angles"].keys():
        if dummy_atoms.intersection(angle) and junction_atom in angle and physical_atoms.intersection(angle) and physical_atom_to_keep not in angle:
            # if the angle involves the junction atom, a dummy atom and a physical atom it might be redundant
            # dual junctions have two possible branches
            corrections["removed_angles"].add(frozenset(angle))

    # remove single and dual anchor dihedral constraints if the dummy group is larger than a single atom
    if dummy_group_dihedrals:
        # find the heaviest terminal atom to anchor dihedrals from the dummy group to
        # we expect this to be one of the physical atom
        all_dummies = set([a.GetIdx() for a in rdkit_mol.GetAtoms() if a.GetIdx() not in core_atoms])
        excluded_atoms = all_dummies | {junction_atom}
        dummy_group_terminal_atom = _get_heaviest_terminal_atom(
            dihedrals=dummy_group_dihedrals, excluded_atoms=excluded_atoms, rdkit_mol=rdkit_mol
        )

        # remove any dihedrals which terminate in this core atom
        for dihedral in dummy_group_dihedrals:
            if dummy_group_terminal_atom not in dihedral:
                corrections["removed_dihedrals"].add(frozenset(dihedral))

    return corrections


def _derive_triple_corrections(junction: dict[str, set[int]], rdkit_mol, force_field_labels, core_atoms, rotor_atoms) -> dict[str, set[tuple[int, ...]]]:
    """
    Derive the corrections for a triple dummy junction as per the best practices of Fleck et al.

    Notes
    -----
    - Non-planar triple junctions can not be fully separated, so we just soften all angles and remove dihedrals terminating in the
        dummy junction atom originating from the physical system.
    - Non-planar triple junctions with large dummy groups will have dihedrals removed to avoid single and dual anchor dihedral constraints.
    """
    corrections = _get_correction_dict()

    junction_atom = junction["junction_atom"]
    physical_atoms = junction["physical"]
    dummy_atoms = junction["dummies"]

    # determine the nature of the junction by the hybridisation of the physical junction atom
    junction_hybridisation = rdkit_mol.GetAtomWithIdx(junction_atom).GetHybridization()

    if junction_hybridisation == Chem.HybridizationType.SP2:
        print("Junction is planar, applying planar junction corrections.")
        raise NotImplementedError("Planar triple junctions not supported")
    else:
        print("Junction is non-planar, applying non-planar junction corrections.")
        # first we need to find the three angles to soften
        for angle in force_field_labels["Angles"].keys():
            junction_in_angle = junction_atom in angle
            # we only want a single dummy atom in the angle
            dummy_in_angle = len(dummy_atoms.intersection(angle)) == 1
            physical_in_angle = len(physical_atoms.intersection(angle)) == 1
            if junction_in_angle and dummy_in_angle and physical_in_angle:
                corrections["softened_angles"].add(frozenset(angle))

        # make sure we found 3 angles to soften as expected
        if len(corrections["softened_angles"]) != 3:
            raise ValueError(f"Expected to find 3 angles to soften for a triple junction, but found {len(corrections['softened_angles'])}. "
                             f"Found the following softened angles: {corrections['softened_angles']}")

        # remove all dihedrals which terminate in the dummy junction atom and originate from the physical system
        # also catch any dihedrals which terminate in a physical atom and originate from within the dummy group
        other_core_atoms = core_atoms - physical_atoms - {junction_atom}
        dummy_group_dihedrals = set()
        for dihedral in force_field_labels["ProperTorsions"].keys():
            junction_in_dihedral = junction_atom in dihedral
            dummy_in_dihedral = len(dummy_atoms.intersection(dihedral)) == 1
            physical_in_dihedral = len(physical_atoms.intersection(dihedral)) == 1
            other_core_in_dihedral = other_core_atoms.intersection(dihedral)
            # if the dihedral originates from the physical and terminates in the dummy junction we need to remove
            if dummy_in_dihedral and physical_in_dihedral and junction_in_dihedral and other_core_in_dihedral:
                corrections["removed_dihedrals"].add(frozenset(dihedral))
            # if the dihedral originates from the dummy group and terminates in the physical atom we need to collect for
            # anchor corrections
            elif dummy_in_dihedral and junction_in_dihedral and physical_in_dihedral and not other_core_in_dihedral:
                dummy_group_dihedrals.add(dihedral)

        if dummy_group_dihedrals:
            # find the heaviest terminal atom to anchor dihedrals from the dummy group to
            # we expect this to be one of the physical atom
            all_dummies = set([a.GetIdx() for a in rdkit_mol.GetAtoms() if a.GetIdx() not in core_atoms])
            excluded_atoms = all_dummies | {junction_atom}
            dummy_group_terminal_atom = _get_heaviest_terminal_atom(
                dihedrals=dummy_group_dihedrals, excluded_atoms=excluded_atoms, rdkit_mol=rdkit_mol
            )

            # remove any dihedrals which terminate in this core atom
            for dihedral in dummy_group_dihedrals:
                if dummy_group_terminal_atom not in dihedral:
                    corrections["removed_dihedrals"].add(frozenset(dihedral))

    return corrections

def _find_improper_corrections(junction: dict[str, set[int]], force_field_labels) -> dict[str, set[tuple[int, ...]]]:
    """
    Find any improper torsions which need to be removed for a given junction.

    Notes
    -----
    - We remove any improper torsions which involve the junction atom and a dummy atom as recommended by Fleck at al
    """
    corrections = _get_correction_dict()
    junction_atom = junction["junction_atom"]
    dummy_atoms = junction["dummies"]
    for dihedral in force_field_labels["ImproperTorsions"].keys():
        # if any dummy and core atom are in the improper make sure to remove it
        if junction_atom in dihedral and dummy_atoms.intersection(dihedral):
            corrections["removed_impropers"].add(frozenset(dihedral))
    return corrections


def _check_dummy_junction(junction: dict[str, set[int]], corrections: dict[str, set[tuple[int, ...]]], force_field_labels, core_atoms) -> bool:
    """
    Check the dummy junction atom is correctly connected via only 3 or less non-redundant valence terms.
    """
    junction_atom = junction["junction_atom"]
    dummy_atoms = junction["dummies"]
    physical_atoms = junction["physical"]
    other_core_atoms = core_atoms - physical_atoms - {junction_atom}
    # check there is a single bond/constraint between the dummy and physical junction atoms
    for dummy_atom in dummy_atoms:
        valence_terms = {
            "bonds": set(),
            "angles": set(),
            "dihedrals": set(),
        }
        # check bond terms
        junction_bond = {dummy_atom, junction_atom}
        for bond in force_field_labels["Bonds"].keys():
            if set(bond) == junction_bond:
                valence_terms["bonds"].add(bond)
        if len(valence_terms["bonds"]) != 1:
            return False

        # check angle terms
        for angle in force_field_labels["Angles"].keys():
            angle = frozenset(angle)
            junction_in_angle = junction_atom in angle
            dummy_in_angle = dummy_atom in angle
            physical_in_angle = len(physical_atoms.intersection(angle)) == 1
            if junction_in_angle and dummy_in_angle and physical_in_angle:
                if angle not in corrections["removed_angles"]:
                    valence_terms["angles"].add(angle)

        # now check dihedral terms
        for dihedral in force_field_labels["ProperTorsions"].keys():
            dihedral = frozenset(dihedral)
            junction_in_dihedral = junction_atom in dihedral
            dummy_in_dihedral = dummy_atom in dihedral
            physical_in_dihedral = len(physical_atoms.intersection(dihedral)) == 1
            other_core_in_dihedral = len(other_core_atoms.intersection(dihedral)) == 1
            if junction_in_dihedral and dummy_in_dihedral and physical_in_dihedral and other_core_in_dihedral:
                if dihedral not in corrections["removed_dihedrals"]:
                    valence_terms["dihedrals"].add(dihedral)

        # now we need to check the valence terms are correct
        total_terms = len(valence_terms["bonds"]) + len(valence_terms["angles"]) + len(valence_terms["dihedrals"])
        if total_terms > 3 or total_terms == 0:
            # check if this is a triple junction with softened angles then this is expected
            if total_terms == 4 and all(angle in corrections["softened_angles"] for angle in valence_terms["angles"]):
                continue
            return False

    return True


def _get_correction_dict():
    return {
        "removed_angles": set(),
        "softened_angles": set(),
        "removed_dihedrals": set(),
        "stiffened_dihedrals": set(),
        "removed_impropers": set(),
    }

