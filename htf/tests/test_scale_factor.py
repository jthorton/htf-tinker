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


def test_toluene_to_pyridine_scale_factor_angles(htf_toluene_pyridine):
    """Make sure that the scaling factor is applied correctly to the angle terms in a more complex system involving angles
    with two dummy atoms.
    """
    htf = htf_toluene_pyridine["htf"]
    scale_factor = 0.1
    mapping = htf_toluene_pyridine["mapping"]
    toluene_labels = htf_toluene_pyridine["toluene_labels"]
    pyridine_labels = htf_toluene_pyridine["pyridine_labels"]
    softened_htf = _scale_angles_and_torsions(htf=htf, scale_factor=scale_factor)
    softened_hybrid_system = softened_htf.hybrid_system
    forces = {force.getName(): force for force in softened_hybrid_system.getForces()}

    # there should be 3 standard angle force terms (non-interpolated) corresponding to the 3 angles
    # formed by the unique old dummy atoms in the methyl of toluene
    standard_angle_force = forces["HarmonicAngleForce"]
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 3
    for i in range(num_angles):
        p1, p2, p3, angle_eq, k = standard_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # make sure all the atoms are in the dummy region of toluene
        assert set(angle).issubset(htf._atom_classes["unique_old_atoms"])
        # now compare the parameters to the toluene labels
        toluene_angle = toluene_labels["Angles"][angle]
        # angle_eq
        assert angle_eq == toluene_angle.angle.m_as(offunit.radian) * unit.radian
        # k
        assert k == toluene_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2) * unit.kilojoule_per_mole / unit.radian**2

    # there should be 21 interpolated angle terms covering fully mapped and scaled angles
    custom_angle_force = forces["CustomAngleForce"]
    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 21
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"
    # track how many are scaled we expect 5
    scaled_angles = 0
    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        # p1, p2, p3 are the index in toluene/pyridine get the expected parameters from the labels
        if p1 in htf._atom_classes["unique_old_atoms"] or p3 in htf._atom_classes["unique_old_atoms"]:
            # this angle involves at least one old dummy atom from toluene
            toluene_angle = toluene_labels["Angles"][(p1, p2, p3)]
            # lambda_0 angle
            assert params[0] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == toluene_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            # lambda_1 angle stays the same
            assert params[2] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_1 k is scaled by 0.1
            expected_k = toluene_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2) * scale_factor
            assert params[3] == expected_k
            scaled_angles += 1
        # there are no dummy atoms in pyridine so all others must be fully mapped
        else:
            # fully mapped angle
            toluene_angle = toluene_labels["Angles"][(p1, p2, p3)]
            e1 = mapping.componentA_to_componentB[p1]
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3]
            pyridine_angle = pyridine_labels["Angles"][(e1, e2, e3)]
            # lambda_0 angle should be the toluene angle
            assert params[0] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_0 k should be the toluene k
            assert params[1] == toluene_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            # lambda_1 angle should be the pyridine angle
            assert params[2] == pyridine_angle.angle.m_as(offunit.radian)
            # lambda_1 k should be the pyridine k
            assert params[3] == pyridine_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)

    assert scaled_angles == 5

def test_toluene_to_pyridine_scale_factor_torsions(htf_toluene_pyridine):
    """Make sure that the scaling factor is applied correctly to the torsion terms in a more complex system involving
    torsions with two dummy atoms and improper torsions.
    """
    htf = htf_toluene_pyridine["htf"]
    scale_factor = 0.1
    mapping = htf_toluene_pyridine["mapping"]
    toluene_labels = htf_toluene_pyridine["toluene_labels"]
    pyridine_labels = htf_toluene_pyridine["pyridine_labels"]
    softened_htf = _scale_angles_and_torsions(htf=htf, scale_factor=scale_factor)
    softened_hybrid_system = softened_htf.hybrid_system
    forces = {force.getName(): force for force in softened_hybrid_system.getForces()}

    # check the standard torsion force has the correct number of terms
    standard_torsion_force = forces["PeriodicTorsionForce"]
    num_standard_torsions = standard_torsion_force.getNumTorsions()
    # Note the 5 impropers should be conserved which should give us 15 improper potentials in total but
    # due to the degenerate order in which the impropers match only 2 (6 terms total) are stored in this force the rest
    # are in the interpolated force
    # So we have 22 torsions in total
    assert num_standard_torsions == 22
    for i in range(num_standard_torsions):
        p1, p2, p3, p4, periodicity, phase, k = standard_torsion_force.getTorsionParameters(i)
        torsion = (p1, p2, p3, p4)
        # make sure this is a fully mapped torsion
        assert len(htf._atom_classes["unique_old_atoms"].intersection(torsion)) == 0
        # now compare the parameters to the toluene and pyridine labels they should be the same at both end states
        # check if we have a proper or improper torsion
        if torsion in toluene_labels["ProperTorsions"]:
            toluene_torsion = toluene_labels["ProperTorsions"][torsion]
            e1 = mapping.componentA_to_componentB[p1]
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3]
            e4 = mapping.componentA_to_componentB[p4]
            pyridine_torsion = pyridine_labels["ProperTorsions"][(e1, e2, e3, e4)]
            # used to account for smirnoff improper see below
            improper_scale = 1.0
        else:
            # if this is an improper openff expects the central atom to be index 1
            # but openmm stores it as index 0 so change the order
            torsion = (p2, p1, p3, p4)
            # smirnoff improper are also applied 3 times to the system with different permutations of the connected
            # atoms so we need to account for that here
            improper_scale = 1 / 3
            toluene_torsion = toluene_labels["ImproperTorsions"][torsion]
            e1 = mapping.componentA_to_componentB[p1]
            e2 = mapping.componentA_to_componentB[p2]
            e3 = mapping.componentA_to_componentB[p3]
            e4 = mapping.componentA_to_componentB[p4]
            pyridine_torsion = pyridine_labels["ImproperTorsions"][(e2, e1, e3, e4)]
        # check against toluene parameters
        assert periodicity in toluene_torsion.periodicity
        term_index = toluene_torsion.periodicity.index(periodicity)
        assert phase == toluene_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole * improper_scale
        # check against pyridine parameters
        assert periodicity in pyridine_torsion.periodicity
        term_index = pyridine_torsion.periodicity.index(periodicity)
        assert phase == pyridine_torsion.phase[term_index].m_as(offunit.radian) * unit.radian
        assert k == pyridine_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * unit.kilojoule_per_mole * improper_scale

    # check the interpolated terms - this force has the scaled torsions and impropers
    # Note some impropers are incorrectly scaled here due to the degenerate ordering mentioned above
    # this should have no effect on the energy as the same value is used for both end states
    custom_torsion_force = forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"
    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 39

    # track the number of scaled torsions we expect 13
    scaled_torsions = 0
    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        # p1, p2, p3, p4 are the index in toluene/pyridine get the expected parameters from the labels
        if 1<= len({p1, p2, p3, p4}.intersection(htf._atom_classes["unique_old_atoms"])) < 4:
            # this is a torsion involving at least one old dummy atom from toluene and should be scaled
            if (p1, p2, p3, p4) in toluene_labels["ProperTorsions"]:
                toluene_torsion = toluene_labels["ProperTorsions"][(p1, p2, p3, p4)]
                improper_scale = 1.0
            else:
                # if this is an improper openff expects the central atom to be index 1
                # but openmm stores it as index 0 so change the order
                torsion = (p2, p1, p3, p4)
                # smirnoff improper are also applied 3 times to the system with different permutations of the connected
                # atoms so we need to account for that here
                improper_scale = 1 / 3
                toluene_torsion = toluene_labels["ImproperTorsions"][torsion]
            # lambda_0 periodicity
            assert params[0] in toluene_torsion.periodicity
            term_index = toluene_torsion.periodicity.index(params[0])
            # lambda_0 phase
            assert params[1] == toluene_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_0 k
            assert params[2] == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale
            # lambda_1 periodicity stays the same
            assert params[3] in toluene_torsion.periodicity
            # lambda_1 phase stays the same
            assert params[4] == toluene_torsion.phase[term_index].m_as(offunit.radian)
            # lambda_1 k is scaled by 0.1
            expected_k = toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale * scale_factor
            assert params[5] == expected_k
            scaled_torsions += 1

        else:
            # this is a fully mapped torsion which is changing between toluene and pyridine
            # however the HTF only allows scaling to or from zero potentials not between potentials so the check
            # is complicated
            if params[0] == 0.0 and params[1] == 0.0 and params[2] == 0.0:
                # this is a pyridine torsion which is zeroed at lambda_0
                # map the torsion
                e1 = mapping.componentA_to_componentB[p1]
                e2 = mapping.componentA_to_componentB[p2]
                e3 = mapping.componentA_to_componentB[p3]
                e4 = mapping.componentA_to_componentB[p4]
                if (e1, e2, e3, e4) in pyridine_labels["ProperTorsions"]:
                    pyridine_torsion = pyridine_labels["ProperTorsions"][(e1, e2, e3, e4)]
                    improper_scale = 1.0
                else:
                    # if this is an improper openff expects the central atom to be index 1
                    # but openmm stores it as index 0 so change the order
                    torsion = (e2, e1, e3, e4)
                    # smirnoff improper are also applied 3 times to the system with different permutations of the connected
                    # atoms so we need to account for that here
                    improper_scale = 1 / 3
                    pyridine_torsion = pyridine_labels["ImproperTorsions"][torsion]
                # lambda_1 periodicity
                assert params[3] in pyridine_torsion.periodicity
                term_index = pyridine_torsion.periodicity.index(params[3])
                # lambda_1 phase
                assert params[4] == pyridine_torsion.phase[term_index].m_as(offunit.radian)
                # lambda_1 k
                assert params[5] == pyridine_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale
            elif params[3] == 0.0 and params[4] == 0.0 and params[5] == 0.0:
                # this is a toluene torsion which is zeroed at lambda_1
                if (p1, p2, p3, p4) in toluene_labels["ProperTorsions"]:
                    toluene_torsion = toluene_labels["ProperTorsions"][(p1, p2, p3, p4)]
                    improper_scale = 1.0
                else:
                    # if this is an improper openff expects the central atom to be index 1
                    # but openmm stores it as index 0 so change the order
                    torsion = (p2, p1, p3, p4)
                    # smirnoff improper are also applied 3 times to the system with different permutations of the connected
                    # atoms so we need to account for that here
                    improper_scale = 1 / 3
                    toluene_torsion = toluene_labels["ImproperTorsions"][torsion]
                # lambda_0 periodicity
                assert params[0] in toluene_torsion.periodicity
                term_index = toluene_torsion.periodicity.index(params[0])
                # lambda_0 phase
                assert params[1] == toluene_torsion.phase[term_index].m_as(offunit.radian)
                # lambda_0 k
                assert params[2] == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole) * improper_scale

    assert scaled_torsions == 13
