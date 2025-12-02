"""Test the dihedral and angle scaling factor correction method in the HTF for a range of transformations"""
import copy

from htf.utils import _scale_angles_and_torsions
import openmm
from openmm import unit, app
from openfe.protocols.openmm_rfe import _rfe_utils
import pytest
from openff.units import unit as offunit


def test_chloro_ethane_no_scale_energy(htf_chloro_ethane):
    """Make sure there is no change to the system energy when the softening factor is 0.0"""
    htf = htf_chloro_ethane["htf"]

    softened_htf = _scale_angles_and_torsions(htf=htf, scale_factor=1.0)
    # we need to make sure we have the same number of particles in the systems
    original_hybrid_system = htf.hybrid_system
    softened_hybrid_system = softened_htf.hybrid_system
    assert original_hybrid_system.getNumParticles() == softened_hybrid_system.getNumParticles()
    for i in range(original_hybrid_system.getNumParticles()):
        original_mass = original_hybrid_system.getParticleMass(i)
        softened_mass = softened_hybrid_system.getParticleMass(i)
        assert original_mass == softened_mass
    # now check that the constraints are the same
    original_constraints = original_hybrid_system.getNumConstraints()
    softened_constraints = softened_hybrid_system.getNumConstraints()
    assert original_constraints == softened_constraints
    for i in range(original_constraints):
        original_constraint = original_hybrid_system.getConstraintParameters(i)
        softened_constraint = softened_hybrid_system.getConstraintParameters(i)
        assert original_constraint == softened_constraint

    # now check that the single point energies are the same for the hybrid systems at the end states
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = openmm.Platform.getPlatformByName("CPU")
    default_lambda = _rfe_utils.lambdaprotocol.LambdaProtocol()

    # set the nonbonded method to NoCutoff to avoid any cutoff issues
    for system in [original_hybrid_system, softened_hybrid_system]:
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
            if isinstance(force, openmm.CustomNonbondedForce):
                force.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)

    energies = []
    for end_state in [0.0, 1.0]:
        for system in [original_hybrid_system, softened_hybrid_system]:
            hybrid_simulation = app.Simulation(
                topology=htf.omm_hybrid_topology,
                system=system,
                integrator=copy.deepcopy(integrator),
                platform=platform
            )
            # set the lambda parameters
            for name, func in default_lambda.functions.items():
                val = func(end_state)
                hybrid_simulation.context.setParameter(name, val)

            hybrid_simulation.context.setPositions(htf.hybrid_positions)
            hybrid_state = hybrid_simulation.context.getState(getEnergy=True)
            hybrid_energy = hybrid_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            energies.append(hybrid_energy)
        # now compare the energies
        assert energies[0] == pytest.approx(energies[1])


def test_chloro_ethane_scale_factor_angles(htf_chloro_ethane):
    """Make sure that the scaling factor is applied correctly to the angle terms"""
    htf = htf_chloro_ethane["htf"]
    scale_factor = 0.1
    mapping = htf_chloro_ethane["mapping"]
    chloro_labels = htf_chloro_ethane["chloro_labels"]
    ethane_labels = htf_chloro_ethane["ethane_labels"]
    softened_htf = _scale_angles_and_torsions(htf=htf, scale_factor=scale_factor)
    softened_hybrid_system = softened_htf.hybrid_system
    forces = {force.getName(): force for force in softened_hybrid_system.getForces()}

    # there should be 0 standard angle force terms (non-interpolated)
    # as we now soften the unique angles they should have been moved to the custom angle force
    standard_angle_force = forces["HarmonicAngleForce"]
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 0

    # there should then be 15 interpolated angle terms
    # included a mix of unique and fully mapped terms
    custom_angle_force = forces["CustomAngleForce"]
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"

    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 15
    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        # p1, p2, p3 are the index in chloroethane get the expected parameters from the labels
        if p1 == 0 or p3 == 0:
            # this is the chlorine atom which goes to a dummy in ethane
            chloro_angle = chloro_labels["Angles"][(p1, p2, p3)]
            # lambda_0 angle
            assert params[0] == chloro_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == chloro_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            # lambda_1 angle stays the same
            assert params[2] == chloro_angle.angle.m_as(offunit.radian)
            # lambda_1 k is scaled by (1 - 0.9) = 0.1
            expected_k = chloro_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2) * scale_factor
            assert params[3] == expected_k
        elif p1 == 8 or p3 == 8:
            # this is the hydrogen atom which goes to a dummy in chloroethane
            # map the terms
            e1 = mapping.componentA_to_componentB[p1] if p1 != 8 else 0
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3] if p3 != 8 else 0
            ethane_angle = ethane_labels["Angles"][(e1, e2, e3)]
            # lambda_0 should have the same angle
            assert params[0] == ethane_angle.angle.m_as(offunit.radian)
            # lambda_0 k is scaled by (1 - 0.9) = 0.1
            expected_k = ethane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2) * scale_factor
            assert params[1] == expected_k
            # lambda_1 angle
            assert params[2] == ethane_angle.angle.m_as(offunit.radian)
            # lambda_1 k
            assert params[3] == ethane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
        else:
            # fully mapped angle
            chloro_angle = chloro_labels["Angles"][(p1, p2, p3)]
            e1 = mapping.componentA_to_componentB[p1]
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3]
            ethane_angle = ethane_labels["Angles"][(e1, e2, e3)]
            # lambda_0 angle
            assert params[0] == chloro_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == chloro_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            # lambda_1 angle
            assert params[2] == ethane_angle.angle.m_as(offunit.radian)
            # lambda_1 k
            assert params[3] == ethane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)


def test_chloro_ethane_scale_factor_torsions(htf_chloro_ethane):
    htf = htf_chloro_ethane["htf"]
    scale_factor = 0.1
    mapping = htf_chloro_ethane["mapping"]
    chloro_labels = htf_chloro_ethane["chloro_labels"]
    ethane_labels = htf_chloro_ethane["ethane_labels"]
    softened_htf = _scale_angles_and_torsions(htf=htf, scale_factor=scale_factor)
    softened_hybrid_system = softened_htf.hybrid_system
    forces = {force.getName(): force for force in softened_hybrid_system.getForces()}

    # there should be 9 interpolated torsion terms only involving ghost atoms
    # there should be 3 terms for the chloroethane unique torsions with 2 peroidicities each
    # and 3 terms for the ethane unique torsions with 1 periodicity
    custom_torsion_force = forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"
    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 9

    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        # p1, p2, p3, p4 are the index in chloroethane get the expected parameters from the labels
        if p1 == 0 or p4 == 0:
            # this is the chlorine atom which goes to a dummy in ethane
            chloro_torsion = chloro_labels["ProperTorsions"][(p1, p2, p3, p4)]
            # lambda_0 periodicity
            assert params[0] in chloro_torsion.periodicity
            term_index = chloro_torsion.periodicity.index(params[0])
            # lambda_0 phase
            assert params[1] == chloro_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_0 k
            assert params[2] == chloro_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
            # lambda_1 periodicity stays the same
            assert params[3] in chloro_torsion.periodicity
            # lambda_1 phase stays the same
            assert params[4] == chloro_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_1 k is scaled by (1 - 0.9) = 0.1
            expected_k = chloro_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * scale_factor
            assert params[5] == expected_k
        elif p1 == 8 or p4 == 8:
            # this is the hydrogen atom which goes to a dummy in chloroethane
            # map the terms
            e1 = mapping.componentA_to_componentB[p1] if p1 != 8 else 0
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3]
            e4 = mapping.componentA_to_componentB[p4] if p4 != 8 else 0
            ethane_torsion = ethane_labels["ProperTorsions"][(e1, e2, e3, e4)]
            # lambda_0 periodicity
            assert params[0] in ethane_torsion.periodicity
            term_index = ethane_torsion.periodicity.index(params[0])
            # lambda_0 phase
            assert params[1] == ethane_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_0 k is scaled by (1 - 0.9) = 0.1
            expected_k = ethane_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * scale_factor
            assert params[2] == expected_k
            # lambda_1 periodicity
            assert params[3] in ethane_torsion.periodicity
            # lambda_1 phase
            assert params[4] == ethane_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_1 k
            assert params[5] == ethane_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
        else:
            assert False, f"All torsions in this test should involve a transforming atom but found the following torsion: {(p1, p2, p3, p4)}"

    # check the standard torsion force has the correct number of terms
    standard_torsion_force = forces["PeriodicTorsionForce"]
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
        # check against ethane parameters
        assert periodicity in ethane_torsion.periodicity
        term_index = ethane_torsion.periodicity.index(periodicity)
        assert phase == ethane_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == ethane_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole

