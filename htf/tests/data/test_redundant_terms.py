import pytest
from collections import defaultdict


def _find_junction_terms(htf, junction_atom, end_state: int) -> dict:
    """
    Find all valence terms in the HTF that cross the given junction atom and are present at the given end state.

    Parameters
    ----------
    htf : HTF
        The HTF instance to check.
    junction_atom : int
        The index of the core junction atom to check around.
    end_state : int
        The end state to check (0 or 1).

    Returns
    -------
    dict
        A dictionary with keys "Constraints", "Bonds", "Angles", "Dihedrals", "Impropers" and values as sets of tuples of the atom indices involved in each term that crosses the junction and is present at the given end state.
    """
    end_state_dummy_atoms = htf._atom_classes["unique_old_atoms"] if end_state == 0 else htf._atom_classes["unique_new_atoms"]
    core_atoms = htf._atom_classes["core_atoms"]
    # get all the bonds in the hybrid topology
    hybrid_bonds = list(htf.omm_hybrid_topology.bonds())
    # construct a lookup of bonded atoms
    bonded_atom_lookup = defaultdict(set)
    for bond in hybrid_bonds:
        a1 = bond[0].index
        a2 = bond[1].index
        bonded_atom_lookup[a1].add(a2)
        bonded_atom_lookup[a2].add(a1)
    # for the junction atom check its bonded to dummies and get all the bonded dummies there should be one dummy in this case
    dummy_atoms = [a for a in bonded_atom_lookup[junction_atom] if a in end_state_dummy_atoms]

    # track the found connections
    found_connections = {"Constraints": set(), "Bonds": set(), "Angles": set(), "Dihedrals": set(), "Impropers": set()}

    # now check the bonds
    hybrid_forces = htf._hybrid_system_forces
    junction_bond = {junction_atom, dummy_atoms[0]}
    # check if this is a constraint
    hybrid_system = htf.hybrid_system
    for constraint_idx in range(hybrid_system.getNumConstraints()):
        a1, a2, _ = hybrid_system.getConstraintParameters(constraint_idx)
        if {a1, a2} == junction_bond:
            found_connections["Constraints"].add((a1, a2))

    for bond_force_name in ["core_bond_force", "standard_bond_force"]:
        bond_force = hybrid_forces[bond_force_name]
        for bond_idx in range(bond_force.getNumBonds()):
            terms = bond_force.getBondParameters(bond_idx)
            bond = {*terms[:2]}
            if bond == junction_bond:
                # no forces play with bonds so leave this
                found_connections["Bonds"].add(tuple([*terms[:2]]))

    # now check the angles
    for angle_force_name in ["core_angle_force", "standard_angle_force"]:
        angle_force = hybrid_forces[angle_force_name]
        for angle_idx in range(angle_force.getNumAngles()):
            terms = angle_force.getAngleParameters(angle_idx)
            angle = {*terms[:3]}
            if junction_atom in angle and any(a in end_state_dummy_atoms for a in angle):
                # is this term interpolated off, this will be in the custom angle force
                if angle_force_name == "core_angle_force":
                    # workout which endstate we should be checking
                    if end_state == 0:
                        k = terms[6]
                    else:
                        k = terms[4]
                    if k == 0.0:
                        # this term is interpolated off so we can skip it
                        continue
                found_connections["Angles"].add(tuple([*terms[:3]]))

    # now check the dihedrals
    for dihedral_force_name in ["custom_torsion_force", "unique_atom_torsion_force"]:
        dihedral_force = hybrid_forces[dihedral_force_name]
        for dihedral_idx in range(dihedral_force.getNumTorsions()):
            terms = dihedral_force.getTorsionParameters(dihedral_idx)
            a1, a2, a3, a4 = terms[:4]
            dihedral = {a1, a2, a3, a4}
            if junction_atom in dihedral and any(a in end_state_dummy_atoms for a in dihedral):
                # check if the term is interpolated off, this will be in the custom torsion force
                if dihedral_force_name == "custom_torsion_force":
                    # workout which endstate we should be checking
                    if end_state == 0:
                        k = terms[9]
                    else:
                        k = terms[6]
                    if k == 0.0:
                        # this term is interpolated off so we can skip it
                        continue
                # check if each atom is bonded to the next
                if a1 in bonded_atom_lookup[a2] and a2 in bonded_atom_lookup[a3] and a3 in bonded_atom_lookup[a4]:
                    found_connections["Dihedrals"].add(tuple([a1, a2, a3, a4]))
                else:
                    found_connections["Impropers"].add(tuple([a1, a2, a3, a4]))
    return found_connections


def test_dual_junction_single_branch(htf_toluene_pyridine):
    """
    Test for redundant terms following the dual junction with a single branch rule in the toluene-pyridine hybrid topology.

    Dual junction rules:
    - There should be a single bond crossing the junction and it should not be a constraint.
    - There should be 2 angles crossing the junction which terminate at the dummy junction atom (these should be stiffened though not checked).
    - There should be 3 dihedrals crossing the dummy core junction and they should all terminate at the same fore atom.
    - There should be no impropers crossing the junction which involve the dummy atoms.
    """
    htf = htf_toluene_pyridine["htf"]
    core_junction_atom = 1
    dummy_junction_atom = 0
    core_atoms = htf._atom_classes["core_atoms"]
    found_connections = _find_junction_terms(htf, junction_atom=core_junction_atom, end_state=0)
    print(found_connections)
    # make sure we have a single bond crossing the junction and that it is not a constraint
    assert not found_connections["Constraints"]
    assert len(found_connections["Bonds"]) == 1
    assert found_connections["Bonds"] == {(dummy_junction_atom, core_junction_atom)}
    # there should be two angles crossing the junction which terminate at the dummy junction atom
    dummy_terminal_angles = []
    for angle in found_connections["Angles"]:
        # terminal angles should terminate at the dummy atom and have 2 core atoms
        if dummy_junction_atom in angle and len(set(angle).intersection(core_atoms)) == 2:
            dummy_terminal_angles.append(angle)
    assert len(dummy_terminal_angles) == 2
    # there should be 3 other angles crossing the junction which terminate at the core junction atom
    core_terminal_angles = []
    for angle in found_connections["Angles"]:
        # terminal angles should terminate at the core junction atom and have 1 core atom
        if core_junction_atom in angle and len(set(angle).intersection(core_atoms)) == 1:
            core_terminal_angles.append(angle)
    assert len(core_terminal_angles) == 3
    # there should be 3 dihedrals crossing the junction, and they should all terminate at the same core atom
    core_dihedrals = []
    for dihedral in found_connections["Dihedrals"]:
        # terminal dihedrals should have 1 core atom and 3 dummy atoms
        if {*dihedral[1:3]} == {core_junction_atom, dummy_junction_atom}:
            core_dihedrals.append(dihedral)
    assert len(core_dihedrals) == 3
    core_dihedral_terminal_atoms = set()
    for dihedral in core_dihedrals:
        if dihedral[0] in core_atoms:
            core_dihedral_terminal_atoms.add(dihedral[0])
        else:
            core_dihedral_terminal_atoms.add(dihedral[3])
    assert len(core_dihedral_terminal_atoms) == 1

    # there should be no impropers involving the dummy atoms
    assert not found_connections["Impropers"]







