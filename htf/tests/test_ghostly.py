import openmm

from htf.utils import load_ghostly_corrections, apply_ghostly_corrections
from openff.units import unit as offunit
from openmm import unit
import math


def test_parse_ghostly_output_chloro(ghostly_output_chloroethane_to_ethane, chloroethane_ethane_ghostly_modifications):
    """
    Make sure we can correctly parse ghostly output files and convert the tuples for a simple chloroethane to ethane case.
    """
    corrections = load_ghostly_corrections(ghostly_output_chloroethane_to_ethane)
    expected_corrections = chloroethane_ethane_ghostly_modifications
    assert corrections == expected_corrections

def test_parse_ghostly_output_hsp90(hsp90_11_to_12_ghostly_modifications, ghostly_output_hsp90_11_to_12):
    """
    Make sure we can correctly parse ghostly output files for a more complex HSP90 case with removed and stiffened angles.
    """
    corrections = load_ghostly_corrections(ghostly_output_hsp90_11_to_12)
    expected_corrections = hsp90_11_to_12_ghostly_modifications
    assert corrections == expected_corrections


def test_apply_ghostly_corrections_mass_and_constraints(chloroethane_ethane_ghostly_modifications, htf_chloro_ethane):
    """
    Make sure that applying the ghostly corrections does not change the number of particles, masses, or constraints
    in the hybrid system for chloroethane to ethane.
    """

    corrections = chloroethane_ethane_ghostly_modifications
    htf = htf_chloro_ethane["htf"]
    corrected_htf = apply_ghostly_corrections(htf, corrections)

    # we need to make sure we have the same number of particles and masses
    original_hybrid_system = htf.hybrid_system
    corrected_hybrid_system = corrected_htf.hybrid_system
    assert original_hybrid_system.getNumParticles() == corrected_hybrid_system.getNumParticles()
    for i in range(original_hybrid_system.getNumParticles()):
        original_mass = original_hybrid_system.getParticleMass(i)
        corrected_mass = corrected_hybrid_system.getParticleMass(i)
        assert original_mass == corrected_mass
    # now check that the constraints are the same
    original_constraints = original_hybrid_system.getNumConstraints()
    corrected_constraints = corrected_hybrid_system.getNumConstraints()
    assert original_constraints == corrected_constraints
    for i in range(original_constraints):
        original_constraint = original_hybrid_system.getConstraintParameters(i)
        corrected_constraint = corrected_hybrid_system.getConstraintParameters(i)
        assert original_constraint == corrected_constraint


def test_apply_ghostly_corrections_angles_chloro(chloroethane_ethane_ghostly_modifications, htf_chloro_ethane):
    """
    Make sure that the ghostly corrections are correctly applied to the htf object angles for chloroethane to ethane.
    """
    corrections = chloroethane_ethane_ghostly_modifications
    htf = htf_chloro_ethane["htf"]
    corrected_htf = apply_ghostly_corrections(htf, corrections)
    chloro_labels = htf_chloro_ethane["chloro_labels"]
    ethane_labels = htf_chloro_ethane["ethane_labels"]
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]
    mapping = htf_chloro_ethane["mapping"]

    corrected_hybrid_system = corrected_htf.hybrid_system
    corrected_forces = {force.getName(): force for force in corrected_hybrid_system.getForces()}

    # based on the corrections check that the angle and torsion force contains the expected parameters
    # there should be 0 standard angle force terms (non-interpolated)
    # as we now soften the unique angles they should have been moved to the custom angle force
    standard_angle_force = corrected_forces["HarmonicAngleForce"]
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 0

    # there should then be 15 interpolated angle terms
    # included a mix of unique and fully mapped terms as we now soften unique angles in this junction type
    custom_angle_force = corrected_forces["CustomAngleForce"]
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"

    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 15
    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # check if this angle is expected to be in the correction terms
        if 1<= len(dummy_new_atoms.intersection(angle)) < 3:
            # this angle or the reverse should be in the lambda_0 softened angles
            assert (prob_angle:= angle) in corrections["lambda_0"]["softened_angles"] or (prob_angle:= angle[::-1]) in corrections["lambda_0"]["softened_angles"]
            # check that the parameters are as expected
            # we should be turning on the angle so it should be soft at lambda_0 and on at lambda_1
            assert params[0] == corrections["lambda_0"]["softened_angles"][prob_angle]["theta0"]  # theta0
            assert params[1] == (corrections["lambda_0"]["softened_angles"][prob_angle]["k"] * offunit.kilocalorie_per_mole).m_as(offunit.kilojoule_per_mole)  # k at lambda_0
            # at lambda_1 it should be fully on so check against the ethane parameters
            # map the angle to ethane indices
            e1 = mapping.componentA_to_componentB[p1] if p1 !=8 else 0
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentB_to_componentA[p3] if p3 !=8 else 0
            ethane_angle = ethane_labels["Angles"][(e1, e2, e3)]
            assert params[2] == ethane_angle.angle.m_as(offunit.radian)  # theta0
            assert params[3] == ethane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)

        elif 1<= len(dummy_old_atoms.intersection(angle)) < 3:
            # this angle or the reverse should be in the lambda_1 softened angles
            assert (prob_angle:= angle) in corrections["lambda_1"]["softened_angles"] or (prob_angle:= angle[::-1]) in corrections["lambda_1"]["softened_angles"]
            # check that the parameters are as expected
            # we should be turning off the angle so it should be on at lambda_0 and soft at lambda_1
            # at lambda_0 it should be fully on so check against the chloroethane parameters
            chloro_angle = chloro_labels["Angles"][angle]
            assert params[0] == chloro_angle.angle.m_as(offunit.radian)  # theta0
            assert params[1] == chloro_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)
            assert params[2] == corrections["lambda_1"]["softened_angles"][prob_angle]["theta0"]  # theta0
            assert params[3] == (corrections["lambda_1"]["softened_angles"][prob_angle]["k"] * offunit.kilocalorie_per_mole).m_as(offunit.kilojoule_per_mole)  # k at lambda_1
        else:
            # this is a fully mapped angle with no dummies so it should not be in the corrections and should be
            # using the normal interpolated parameters
            assert angle not in corrections["lambda_0"]["softened_angles"] and angle[::-1] not in corrections["lambda_0"]["softened_angles"]
            assert angle not in corrections["lambda_1"]["softened_angles"] and angle[::-1] not in corrections["lambda_1"]["softened_angles"]
            chloro_angle = chloro_labels["Angles"][angle]
            e1 = mapping.componentA_to_componentB[p1]
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentB_to_componentA[p3]
            ethane_angle = ethane_labels["Angles"][(e1, e2, e3)]
            # lambda_0 angle
            assert params[0] == chloro_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == chloro_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)
            # lambda_1 angle
            assert params[2] == ethane_angle.angle.m_as(offunit.radian)
            # lambda_1 k
            assert params[3] == ethane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)


def test_apply_ghostly_corrections_torsions_chloro(chloroethane_ethane_ghostly_modifications, htf_chloro_ethane):
    """
    Make sure that the ghostly corrections are correctly applied to the htf object torsions for chloroethane to ethane.
    """
    corrections = chloroethane_ethane_ghostly_modifications
    htf = htf_chloro_ethane["htf"]
    corrected_htf = apply_ghostly_corrections(htf, corrections)
    chloro_labels = htf_chloro_ethane["chloro_labels"]
    ethane_labels = htf_chloro_ethane["ethane_labels"]
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]
    mapping = htf_chloro_ethane["mapping"]

    corrected_hybrid_system = corrected_htf.hybrid_system
    corrected_forces = {force.getName(): force for force in corrected_hybrid_system.getForces()}

    # there should be 9 interpolated torsion terms only involving ghost atoms
    # there should be 3 terms for the chloroethane unique torsions with 2 peroidicities each
    # and 3 terms for the ethane unique torsions with 1 periodicity
    custom_torsion_force = corrected_forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"
    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 9

    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        # check if this torsion is expected to be in the correction terms
        if 1<= len(dummy_new_atoms.intersection(torsion)) < 4:
            # this torsion or the reverse should be in the lambda_0 removed torsions
            assert torsion in corrections["lambda_0"]["removed_dihedrals"] or torsion[::-1] in corrections["lambda_0"]["removed_dihedrals"]
            # check that the parameters are as expected
            # we should be turning on the torsion so it should be soft at lambda_0 and on at lambda_1
            # so we need to check against the ethane parameters
            e1 = mapping.componentA_to_componentB[p1] if p1 !=8 else 0
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentB_to_componentA[p3]
            e4 = mapping.componentB_to_componentA[p4] if p4 !=8 else 0
            ethane_torsion = ethane_labels["ProperTorsions"][(e1, e2, e3, e4)]
            assert params[0] in ethane_torsion.periodicity
            term_index = ethane_torsion.periodicity.index(params[0])
            assert params[1] == ethane_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[2] == 0.0  # k
            assert params[3] == ethane_torsion.periodicity[term_index]
            assert params[4] == ethane_torsion.phase[term_index].m_as(offunit.radian)
            assert params[5] == ethane_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)

        elif 1<= len(dummy_old_atoms.intersection(torsion)) < 4:
            # this torsion or the reverse should be in the lambda_1 removed torsions
            assert torsion in corrections["lambda_1"]["removed_dihedrals"] or torsion[::-1] in corrections["lambda_1"]["removed_dihedrals"]
            # check that the parameters are as expected
            # we should be turning off the torsion so it should be on at lambda_0 and soft at lambda_1
            # so we need to check against the chloroethane parameters
            chloro_torsion = chloro_labels["ProperTorsions"][torsion]
            assert params[0] in chloro_torsion.periodicity
            term_index = chloro_torsion.periodicity.index(params[0])
            assert params[1] == chloro_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[2] == chloro_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)  # k
            assert params[3] == chloro_torsion.periodicity[term_index]  # periodicity
            assert params[4] == chloro_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[5] == 0.0  # k

        else:
            assert False, f"All torsions in this test should involve at least one dummy atom, but got {torsion}"

    # check the standard torsion force has the correct number of terms
    standard_torsion_force = corrected_forces["PeriodicTorsionForce"]
    num_standard_torsions = standard_torsion_force.getNumTorsions()
    # there should be 6 fully mapped torsions which each have a single periodicity
    assert num_standard_torsions == 6
    # make sure the terms are correct
    for i in range(num_standard_torsions):
        p1, p2, p3, p4, periodicity, phase, k = standard_torsion_force.getTorsionParameters(i)
        chloro_torsion = chloro_labels["ProperTorsions"][(p1, p2, p3, p4)]
        e1 = mapping.componentA_to_componentB[p1]
        e2 = mapping.componentA_to_componentB[p2]
        e3 = mapping.componentA_to_componentB[p3]
        e4 = mapping.componentA_to_componentB[p4]
        ethane_torsion = ethane_labels["ProperTorsions"][(e1, e2, e3, e4)]
        # check against chloroethane parameters
        assert periodicity in chloro_torsion.periodicity
        term_index = chloro_torsion.periodicity.index(periodicity)
        assert phase == chloro_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == chloro_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole
        # check against ethane parameters which should be the same
        assert periodicity in ethane_torsion.periodicity
        term_index = ethane_torsion.periodicity.index(periodicity)
        assert phase == ethane_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == ethane_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole


def test_apply_ghostly_corrections_hsp90_angles(hsp90_11_to_12_ghostly_modifications, htf_hsp90_11_to_12):
    """
    Make sure that the ghostly corrections are correctly applied to the htf object angles for HSP90 11 to 12.
    """
    corrections = hsp90_11_to_12_ghostly_modifications
    htf = htf_hsp90_11_to_12["htf"]
    hsp90_11_labels = htf_hsp90_11_to_12["hsp90_11_labels"]
    hsp90_12_labels = htf_hsp90_11_to_12["hsp90_12_labels"]
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]
    corrected_htf = apply_ghostly_corrections(htf, corrections)

    corrected_hybrid_system = corrected_htf.hybrid_system
    corrected_forces = {force.getName(): force for force in corrected_hybrid_system.getForces()}

    # based on the corrections check that the angle force contains the expected parameters
    # there should be 6 standard angle force terms (non-interpolated) corresponding to the methyl
    # which involves only dummy atoms in the lambda 1 state and 2 dummies and 1 core atom
    standard_angle_force = corrected_forces["HarmonicAngleForce"]
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 6
    # check that all of the angles involve only new dummy atoms
    for i in range(num_angles):
        p1, p2, p3, theta0, k = standard_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        assert len(dummy_new_atoms.intersection(angle)) >= 2
        # map to new indices
        h1 = htf._hybrid_to_new_map[p1]
        h2 = htf._hybrid_to_new_map[p2]
        h3 = htf._hybrid_to_new_map[p3]
        hsp90_12_angle = hsp90_12_labels["Angles"][(h1, h2, h3)]
        # check the terms are correct
        assert theta0 == hsp90_12_angle.angle.m_as(offunit.radian) * unit.radian
        assert k == hsp90_12_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2) * unit.kilojoule_per_mole / unit.radian ** 2

    # there should then be 64 interpolated angle terms
    # included a mix of unique and fully mapped terms as we now soften unique angles in this junction type
    custom_angle_force = corrected_forces["CustomAngleForce"]
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"

    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 60

    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # check if this angle is expected to be in the correction terms
        if len(dummy_new_atoms.intersection(angle)) == 1:
            # this angle or the reverse should be in the lambda_0 angles removed or stiffened angles
            removed = None
            if angle in corrections["lambda_0"]["removed_angles"] or angle[::-1] in corrections["lambda_0"]["removed_angles"]:
                removed = True
            elif  angle in corrections["lambda_0"]["stiffened_angles"] or angle[::-1] in corrections["lambda_0"]["stiffened_angles"]:
                removed = False
            else:
                assert False, f"Angle {angle} with one dummy atom should be in either removed or stiffened angles in lambda_0 corrections"
            # now check the parameters but we have to map from the hybrid index to the end state index
            h1 = htf._hybrid_to_new_map[p1]
            h2 = htf._hybrid_to_new_map[p2]
            h3 = htf._hybrid_to_new_map[p3]
            hsp90_12_angle = hsp90_12_labels["Angles"][(h1, h2, h3)]
            if removed:
                lambda_0_k = 0 * offunit.kilojoule_per_mole / offunit.radian ** 2
                lambda_0_theta0 = hsp90_12_angle.angle
            else:
                # the angle is stiffened at lambda_0
                lambda_0_k = 100 * offunit.kilocalorie_per_mole / offunit.radian ** 2
                lambda_0_theta0 = 0.5 * math.pi * offunit.radian
            assert params[0] == lambda_0_theta0.m_as(offunit.radian)  # theta0
            assert params[1] == lambda_0_k.m_as(offunit.kilojoule_per_mole/ offunit.radian ** 2)  # k at lambda_0
            # at lambda_1 it should be fully on so check against the hsp90 12 parameters
            assert params[2] == hsp90_12_angle.angle.m_as(offunit.radian)
            assert params[3] == hsp90_12_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)

        # same checks for old atoms
        elif len(dummy_old_atoms.intersection(angle)) == 1:
            # this angle or the reverse should be in the lambda_1 angles removed or stiffened angles
            removed = None
            if angle in corrections["lambda_1"]["removed_angles"] or angle[::-1] in corrections["lambda_1"]["removed_angles"]:
                removed = True
            elif  angle in corrections["lambda_1"]["stiffened_angles"] or angle[::-1] in corrections["lambda_1"]["stiffened_angles"]:
                removed = False
            else:
                assert False, f"Angle {angle} with one dummy atom should be in either removed or stiffened angles in lambda_1 corrections"
            # now check the parameters but we have to map from the hybrid index to the end state index
            h1 = htf._hybrid_to_old_map[p1]
            h2 = htf._hybrid_to_old_map[p2]
            h3 = htf._hybrid_to_old_map[p3]
            hsp90_11_angle = hsp90_11_labels["Angles"][(h1, h2, h3)]
            if removed:
                lambda_1_k = 0 * offunit.kilojoule_per_mole / offunit.radian ** 2
                lambda_1_theta0 = hsp90_11_angle.angle
            else:
                # the angle is stiffened at lambda_1
                lambda_1_k = 100 * offunit.kilocalorie_per_mole / offunit.radian ** 2
                lambda_1_theta0 = 0.5 * math.pi * offunit.radian
            # at lambda_0 it should be fully on so check against the hsp90 11 parameters
            assert params[0] == hsp90_11_angle.angle.m_as(offunit.radian)
            assert params[1] == hsp90_11_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)
            assert params[2] == lambda_1_theta0.m_as(offunit.radian)  # theta0
            assert params[3] == lambda_1_k.m_as(offunit.kilojoule_per_mole/ offunit.radian ** 2)  # k at lambda_1
        else:
            # this is a fully mapped angle with no dummies so it should not be in the corrections and should be
            # using the normal interpolated parameters
            assert angle not in corrections["lambda_0"]["removed_angles"] and angle[::-1] not in corrections["lambda_0"]["removed_angles"]
            assert angle not in corrections["lambda_0"]["stiffened_angles"] and angle[::-1] not in corrections["lambda_0"]["stiffened_angles"]
            assert angle not in corrections["lambda_1"]["removed_angles"] and angle[::-1] not in corrections["lambda_1"]["removed_angles"]
            assert angle not in corrections["lambda_1"]["stiffened_angles"] and angle[::-1] not in corrections["lambda_1"]["stiffened_angles"]
            # now check the parameters but we have to map from the hybrid index to the end state index
            h1_old = htf._hybrid_to_old_map[p1]
            h2_old = htf._hybrid_to_old_map[p2]
            h3_old = htf._hybrid_to_old_map[p3]
            hsp90_11_angle = hsp90_11_labels["Angles"][(h1_old, h2_old, h3_old)]
            h1_new = htf._hybrid_to_new_map[p1]
            h2_new = htf._hybrid_to_new_map[p2]
            h3_new = htf._hybrid_to_new_map[p3]
            hsp90_12_angle = hsp90_12_labels["Angles"][(h1_new, h2_new, h3_new)]
            # lambda_0 angle
            assert params[0] == hsp90_11_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == hsp90_11_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)
            # lambda_1 angle
            assert params[2] == hsp90_12_angle.angle.m_as(offunit.radian)
            # lambda_1 k
            assert params[3] == hsp90_12_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian ** 2)


def test_apply_ghostly_corrections_hsp90_torsions(hsp90_11_to_12_ghostly_modifications, htf_hsp90_11_to_12):
    """
    Make sure that the ghostly corrections are correctly applied to the htf object torsions for HSP90 11 to 12.
    TODO fix this test, this currently keeps impropers due to a mistake in the mapping passed to ghostly
    """
    corrections = hsp90_11_to_12_ghostly_modifications
    htf = htf_hsp90_11_to_12["htf"]
    hsp90_11_labels = htf_hsp90_11_to_12["hsp90_11_labels"]
    hsp90_12_labels = htf_hsp90_11_to_12["hsp90_12_labels"]
    dummy_old_atoms = htf._atom_classes["unique_old_atoms"]
    dummy_new_atoms = htf._atom_classes["unique_new_atoms"]
    corrected_htf = apply_ghostly_corrections(htf, corrections)

    corrected_hybrid_system = corrected_htf.hybrid_system
    corrected_forces = {force.getName(): force for force in corrected_hybrid_system.getForces()}

    # Most of the torsions in this molecule should remain the same due to the small localised change
    standard_torsion_force = corrected_forces["PeriodicTorsionForce"]
    num_torsions = standard_torsion_force.getNumTorsions()
    assert num_torsions == 128
    # check that all of the standard torsions have at most 1 dummy atom
    for i in range(num_torsions):
        p1, p2, p3, p4, periodicity, phase, k = standard_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        if len(dummy_new_atoms.intersection(torsion)) == 0 and len(dummy_old_atoms.intersection(torsion)) == 0:
            # map to old and new indices
            h1_old = htf._hybrid_to_old_map[p1]
            h2_old = htf._hybrid_to_old_map[p2]
            h3_old = htf._hybrid_to_old_map[p3]
            h4_old = htf._hybrid_to_old_map[p4]
            torsion_old = (h1_old, h2_old, h3_old, h4_old)
            if torsion_old in hsp90_11_labels["ProperTorsions"]:
                improper_scale = 1.0
                state_a_torsion = hsp90_11_labels["ProperTorsions"][torsion_old]
                h1_new = htf._hybrid_to_new_map[p1]
                h2_new = htf._hybrid_to_new_map[p2]
                h3_new = htf._hybrid_to_new_map[p3]
                h4_new = htf._hybrid_to_new_map[p4]
                state_b_torsion = hsp90_12_labels["ProperTorsions"][(h1_new, h2_new, h3_new, h4_new)]
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                torsion_old = (h2_old, h1_old, h3_old, h4_old)
                state_a_torsion = hsp90_11_labels["ImproperTorsions"][torsion_old]
                h1_new = htf._hybrid_to_new_map[p1]
                h2_new = htf._hybrid_to_new_map[p2]
                h3_new = htf._hybrid_to_new_map[p3]
                h4_new = htf._hybrid_to_new_map[p4]
                torsion_new = (h2_new, h1_new, h3_new, h4_new)
                state_b_torsion = hsp90_12_labels["ImproperTorsions"][torsion_new]
                improper_scale = 1 / 3
        elif 1<= len(dummy_new_atoms.intersection(torsion)) < 4:
            # this torsion is held fixed and used as a nonredundant connection and so it should be set to the lambda 1 state
            h1 = htf._hybrid_to_new_map[p1]
            h2 = htf._hybrid_to_new_map[p2]
            h3 = htf._hybrid_to_new_map[p3]
            h4 = htf._hybrid_to_new_map[p4]
            if (h1, h2, h3, h4) in hsp90_12_labels["ProperTorsions"]:
                improper_scale = 1.0
                state_a_torsion = hsp90_12_labels["ProperTorsions"][(h1, h2, h3, h4)]
                state_b_torsion = state_a_torsion
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                improper_scale = 1 / 3
                state_a_torsion = hsp90_12_labels["ImproperTorsions"][(h2, h1, h3, h4)]
                state_b_torsion = state_a_torsion
        elif 1<= len(dummy_old_atoms.intersection(torsion)) < 4:
            # this torsion is held fixed and used as a nonredundant connection and so it should be set to the lambda 0 state
            h1 = htf._hybrid_to_old_map[p1]
            h2 = htf._hybrid_to_old_map[p2]
            h3 = htf._hybrid_to_old_map[p3]
            h4 = htf._hybrid_to_old_map[p4]
            if (h1, h2, h3, h4) in hsp90_11_labels["ProperTorsions"]:
                improper_scale = 1.0
                state_a_torsion = hsp90_11_labels["ProperTorsions"][(h1, h2, h3, h4)]
                state_b_torsion = state_a_torsion
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                improper_scale = 1 / 3
                state_a_torsion = hsp90_11_labels["ImproperTorsions"][(h2, h1, h3, h4)]
                state_b_torsion = state_a_torsion
        else:
            assert False, f"Torsion {torsion} should not have more than 1 dummy atom in the standard torsion force"

        # check the terms are correct against both states
        assert periodicity in state_a_torsion.periodicity
        term_index = state_a_torsion.periodicity.index(periodicity)
        assert phase == state_a_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == state_a_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole * improper_scale
        # make sure the parameters are the same at the other end state
        assert periodicity in state_b_torsion.periodicity
        term_index = state_b_torsion.periodicity.index(periodicity)
        assert phase == state_b_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == state_b_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole * improper_scale

    # we should then have a subset of torsions which are interpolated involving ghost atoms which are removed
    # and some which are mapped but change
    custom_torsion_force = corrected_forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"

    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 7

    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        # check if this torsion is expected to be in the correction terms
        if 1<= len(dummy_new_atoms.intersection(torsion)) < 4:
            # this torsion or the reverse should be in the lambda_0 removed dihedrals
            assert torsion in corrections["lambda_0"]["removed_dihedrals"] or torsion[::-1] in corrections["lambda_0"]["removed_dihedrals"]
            # check that the parameters are as expected
            # we should be turning on the torsion so it should be soft at lambda_0 and on at lambda_1
            # so we need to check against the hsp90 12 parameters
            h1 = htf._hybrid_to_new_map[p1]
            h2 = htf._hybrid_to_new_map[p2]
            h3 = htf._hybrid_to_new_map[p3]
            h4 = htf._hybrid_to_new_map[p4]
            if (h1, h2, h3, h4) in hsp90_12_labels["ProperTorsions"]:
                improper_scale = 1.0
                state_b_torsion = hsp90_12_labels["ProperTorsions"][(h1, h2, h3, h4)]
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                improper_scale = 1 / 3
                state_b_torsion = hsp90_12_labels["ImproperTorsions"][(h2, h1, h3, h4)]
            assert params[0] in state_b_torsion.periodicity
            term_index = state_b_torsion.periodicity.index(params[0])
            assert params[1] == state_b_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[2] == 0.0  # k
            assert params[3] == state_b_torsion.periodicity[term_index]
            assert params[4] == state_b_torsion.phase[term_index].m_as(offunit.radian)
            assert params[5] == state_b_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale

        elif 1<= len(dummy_old_atoms.intersection(torsion)) < 4:
            # this torsion or the reverse should be in the lambda_1 removed dihedrals
            assert torsion in corrections["lambda_1"]["removed_dihedrals"] or torsion[::-1] in corrections["lambda_1"]["removed_dihedrals"]
            # check that the parameters are as expected
            # we should be turning off the torsion so it should be on at lambda_0 and soft at lambda_1
            # so we need to check against the hsp90 11 parameters
            h1 = htf._hybrid_to_old_map[p1]
            h2 = htf._hybrid_to_old_map[p2]
            h3 = htf._hybrid_to_old_map[p3]
            h4 = htf._hybrid_to_old_map[p4]
            if (h1, h2, h3, h4) in hsp90_11_labels["ProperTorsions"]:
                improper_scale = 1.0
                state_a_torsion = hsp90_11_labels["ProperTorsions"][(h1, h2, h3, h4)]
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                improper_scale = 1 / 3
                state_a_torsion = hsp90_11_labels["ImproperTorsions"][(h2, h1, h3, h4)]
            assert params[0] in state_a_torsion.periodicity
            term_index = state_a_torsion.periodicity.index(params[0])
            assert params[1] == state_a_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[2] == state_a_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale  # k
            assert params[3] == state_a_torsion.periodicity[term_index]  # periodicity
            assert params[4] == state_a_torsion.phase[term_index].m_as(offunit.radian)  # phase
            assert params[5] == 0.0  # k
        else:
            assert False, f"All torsions in this test should involve at least one dummy atom, but got {torsion}"
