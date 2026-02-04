import re

from htf.utils import load_ghostly_corrections, apply_ghostly_corrections, draw_ghostly_modifications
from openff.units import unit as offunit
from openmm import unit
import pytest


def test_parse_ghostly_output_chloro(ghostly_output_chloroethane_to_ethane, chloroethane_ethane_ghostly_modifications):
    """
    Make sure we can correctly parse ghostly output files and convert the tuples for a simple chloroethane to ethane case.
    """
    corrections = load_ghostly_corrections(ghostly_output_chloroethane_to_ethane)
    expected_corrections = chloroethane_ethane_ghostly_modifications
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


def test_apply_ghostly_corrections_unused_correction(chloroethane_ethane_ghostly_modifications, htf_chloro_ethane):
    """Make sure an error is raised if a correction is provided for an angle that is not used as it's not a dummy in the htf"""
    corrections = chloroethane_ethane_ghostly_modifications
    corrections["lambda_0"]["removed_angles"].add((1, 2, 6))  # add an angle that does not involve a dummy
    htf = htf_chloro_ethane["htf"]
    with pytest.raises(ValueError, match=re.escape("The following removed_angles corrections for lambda_0 were not applied: {(1, 2, 6)}")):
        _ = apply_ghostly_corrections(htf, corrections)


def test_draw_ghostly_corrections_chloro(chloroethane_ethane_ghostly_modifications, htf_chloro_ethane, tmp_path):
    """
    Test that we can draw the ghostly corrections for chloroethane to ethane without error.
    """
    htf = htf_chloro_ethane["htf"]
    corrected_htf = apply_ghostly_corrections(htf, chloroethane_ethane_ghostly_modifications)
    draw_ghostly_modifications(
        molecule=htf_chloro_ethane["chloroethane"],
        htf=corrected_htf,
        modifications=chloroethane_ethane_ghostly_modifications,
        filename=(tmp_path / "chloroethane_ghostly.svg").as_posix(),
        end_state=0
    )
    draw_ghostly_modifications(
        molecule=htf_chloro_ethane["ethane"],
        htf=corrected_htf,
        modifications=chloroethane_ethane_ghostly_modifications,
        filename=(tmp_path / "ethane_ghostly.svg").as_posix(),
        end_state=1
    )


def test_draw_ghostly_corrections_toluene_to_pyridine(toluene_pyridine_ghostly_modifications, htf_toluene_pyridine, tmp_path):
    """
    Test that we can draw the ghostly corrections for toluene to pyridine without error.
    """
    htf = htf_toluene_pyridine["htf"]
    corrected_htf = apply_ghostly_corrections(htf, toluene_pyridine_ghostly_modifications)
    draw_ghostly_modifications(
        molecule=htf_toluene_pyridine["toluene"],
        htf=corrected_htf,
        modifications=toluene_pyridine_ghostly_modifications,
        filename=(tmp_path / "toluene_ghostly.svg").as_posix(),
        end_state=0
    )
    draw_ghostly_modifications(
        molecule=htf_toluene_pyridine["pyridine"],
        htf=corrected_htf,
        modifications=toluene_pyridine_ghostly_modifications,
        filename=(tmp_path / "pyridine_ghostly.svg").as_posix(),
        end_state=1
    )


def test_draw_ghostly_corrections_propane_to_dimethyl_ether(propane_dimethylether_ghostly_modifications, htf_propane_dimethyl_ether, tmp_path):
    """
    Test that we can draw the ghostly corrections for propane to dimethyl ether without error.
    """
    htf = htf_propane_dimethyl_ether["htf"]
    corrected_htf = apply_ghostly_corrections(htf, propane_dimethylether_ghostly_modifications)
    draw_ghostly_modifications(
        molecule=htf_propane_dimethyl_ether["propane"],
        htf=corrected_htf,
        modifications=propane_dimethylether_ghostly_modifications,
        filename=(tmp_path / "propane_ghostly.svg").as_posix(),
        end_state=0
    )
    draw_ghostly_modifications(
        molecule=htf_propane_dimethyl_ether["dimethyl_ether"],
        htf=corrected_htf,
        modifications=propane_dimethylether_ghostly_modifications,
        filename=(tmp_path / "dimethyl_ether_ghostly.svg").as_posix(),
        end_state=1
    )