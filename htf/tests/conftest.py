import openmm
from openmm import unit
import pytest
from importlib import resources
from gufe import SmallMoleculeComponent, LigandAtomMapping, ProteinComponent
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from htf.utils import make_htf
from htf import DevelopmentHybridTopologyFactory
from openff.toolkit import ForceField
from openff.units import unit as offunit


@pytest.fixture(scope="module")
def chloroethane():
    """Load chloroethane with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(f / "chloroethane.sdf")


@pytest.fixture(scope="module")
def fluoroethane():
    """Load fluoroethane with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(f / "fluoroethane.sdf")


@pytest.fixture(scope="module")
def ethane():
    """Load ethane with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(f / "ethane.sdf")


@pytest.fixture(scope="module")
def chlorobenzene():
    """Load chlorobenzene with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(
            f / "t4_lysozyme_data" / "chlorobenzene.sdf"
        )


@pytest.fixture(scope="module")
def fluorobenzene():
    """Load fluorobenzene with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(
            f / "t4_lysozyme_data" / "fluorobenzene.sdf"
        )


@pytest.fixture(scope="module")
def benzene():
    """Load benzene with partial charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(
            f / "t4_lysozyme_data" / "benzene.sdf"
        )


@pytest.fixture(scope="module")
def hsp90_11():
    """Load HSP90 ligand 11 with charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(f / "ghostly" / "hsp90_11.sdf")


@pytest.fixture(scope="module")
def hsp90_12():
    """Load HSP90 ligand 12 with charges from sdf file."""
    with resources.files("htf.tests.data") as f:
        return SmallMoleculeComponent.from_sdf_file(f / "ghostly" / "hsp90_12.sdf")


@pytest.fixture(scope="module")
def chloroethane_to_fluoroethane_mapping(chloroethane, fluoroethane):
    """Return a mapping from chloroethane to fluoroethane."""
    return LigandAtomMapping(
        componentA=chloroethane,
        componentB=fluoroethane,
        componentA_to_componentB={
            # perfect one-to-one mapping
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
        },
    )


@pytest.fixture(scope="module")
def chloroethane_to_ethane_mapping(chloroethane, ethane):
    """Return a mapping from chloroethane to ethane."""
    return LigandAtomMapping(
        componentA=chloroethane,
        componentB=ethane,
        componentA_to_componentB={
            # Cl-H not mapped, all others one-to-one
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
        },
    )


@pytest.fixture(scope="module")
def chlorobenzene_to_fluorobenzene_mapping(chlorobenzene, fluorobenzene):
    """Return a mapping from chlorobenzene to fluorobenzene."""
    return LigandAtomMapping(
        componentA=chlorobenzene,
        componentB=fluorobenzene,
        componentA_to_componentB={
            # perfect one-to-one mapping
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 10,
            11: 11,
        },
    )


@pytest.fixture(scope="module")
def chlorobenzene_to_benzene_mapping(chlorobenzene, benzene):
    """Return a mapping from chlorobenzene to benzene."""
    return LigandAtomMapping(
        componentA=chlorobenzene,
        componentB=benzene,
        componentA_to_componentB={
            # Cl-H not mapped, all others one-to-one
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 10,
            11: 11,
        },
    )


@pytest.fixture(scope="module")
def hsp90_11_to_12_mapping(hsp90_11, hsp90_12):
    """Return a mapping from HSP90 ligand 11 to ligand 12."""
    return LigandAtomMapping(
        componentA=hsp90_11,
        componentB=hsp90_12,
        componentA_to_componentB={
            20: 20,
            21: 21,
            22: 22,
            23: 23,
            24: 24,
            25: 25,
            26: 26,
            27: 27,
            28: 28,
            29: 29,
            30: 30,
            31: 35,
            0: 0,
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
            8: 8,
            9: 9,
            10: 10,
            11: 11,
            12: 12,
            13: 13,
            14: 14,
            15: 15,
            16: 16,
            17: 17,
            18: 18,
            19: 19,
        },
    )


@pytest.fixture(scope="module")
def t4_lysozyme_solvated():
    """Load the T4 lysozyme L99A structure and solvent from the pdb file."""
    with resources.files("htf.tests.data") as f:
        return ProteinComponent.from_pdb_file(
            (f / "t4_lysozyme_data" / "t4_lysozyme_solvated.pdb").as_posix()
        )


@pytest.fixture(scope="module")
def htf_chloro_fluoroethane(
    chloroethane, fluoroethane, chloroethane_to_fluoroethane_mapping
):
    """Generate the htf for chloroethane to fluoroethane."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure we interpolate the 1-4 exceptions involving dummy atoms if present
    settings.alchemical_settings.turn_off_core_unique_exceptions = True
    small_ff = settings.forcefield_settings.small_molecule_forcefield
    if ".offxml" not in small_ff:
        small_ff += ".offxml"
    ff = ForceField(small_ff)
    chloro_openff = chloroethane.to_openff()
    chloro_charges = chloro_openff.partial_charges.m_as(offunit.elementary_charge)
    chloro_labels = ff.label_molecules(chloro_openff.to_topology())[0]
    fluoro_openff = fluoroethane.to_openff()
    fluoro_charges = fluoro_openff.partial_charges.m_as(offunit.elementary_charge)
    fluoro_labels = ff.label_molecules(fluoro_openff.to_topology())[0]
    htf = make_htf(mapping=chloroethane_to_fluoroethane_mapping, settings=settings)
    hybrid_system = htf.hybrid_system
    forces = {force.getName(): force for force in hybrid_system.getForces()}

    return {
        "htf": htf,
        "hybrid_system": hybrid_system,
        "forces": forces,
        "chloro_labels": chloro_labels,
        "fluoro_labels": fluoro_labels,
        "mapping": chloroethane_to_fluoroethane_mapping,
        "chloroethane": chloroethane,
        "fluoroethane": fluoroethane,
        "chloro_charges": chloro_charges,
        "fluoro_charges": fluoro_charges,
        "electrostatic_scale": ff.get_parameter_handler("Electrostatics").scale14,
        "vdW_scale": ff.get_parameter_handler("vdW").scale14,
        "force_field": ff,
    }


@pytest.fixture(scope="module")
def htf_chloro_ethane(chloroethane, ethane, chloroethane_to_ethane_mapping):
    """Generate the htf for chloroethane to ethane with interpolate 1-4s on!"""
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure we interpolate the 1-4 exceptions involving dummy atoms
    settings.alchemical_settings.turn_off_core_unique_exceptions = True
    small_ff = settings.forcefield_settings.small_molecule_forcefield
    if ".offxml" not in small_ff:
        small_ff += ".offxml"
    ff = ForceField(small_ff)

    chloro_openff = chloroethane.to_openff()
    chloro_charges = chloro_openff.partial_charges.m_as(offunit.elementary_charge)
    chloro_labels = ff.label_molecules(chloro_openff.to_topology())[0]
    ethane_openff = ethane.to_openff()
    ethane_charges = ethane_openff.partial_charges.m_as(offunit.elementary_charge)
    ethane_labels = ff.label_molecules(ethane_openff.to_topology())[0]
    htf = make_htf(mapping=chloroethane_to_ethane_mapping, settings=settings)
    hybrid_system = htf.hybrid_system
    forces = {force.getName(): force for force in hybrid_system.getForces()}

    return {
        "htf": htf,
        "hybrid_system": hybrid_system,
        "forces": forces,
        "chloro_labels": chloro_labels,
        "ethane_labels": ethane_labels,
        "mapping": chloroethane_to_ethane_mapping,
        "chloroethane": chloroethane,
        "ethane": ethane,
        "chloro_charges": chloro_charges,
        "ethane_charges": ethane_charges,
        "electrostatic_scale": ff.get_parameter_handler("Electrostatics").scale14,
        "vdW_scale": ff.get_parameter_handler("vdW").scale14,
        "force_field": ff,
    }


def apply_box_vectors_and_fix_nn_force(
    hybrid_topology_factory: DevelopmentHybridTopologyFactory, force_field: ForceField
):
    """
    Edit the systems in the hybrid topology factory to have the correct box vectors and nonbonded force settings for the T4 lysozyme system.
    """
    hybrid_system = hybrid_topology_factory.hybrid_system
    # as we use a pre-solvated system, we need to correct the nonbonded methods and cutoffs and set the box vectors
    box_vectors = [
        openmm.vec3.Vec3(x=6.90789161545809, y=0.0, z=0.0) * unit.nanometer,
        openmm.vec3.Vec3(x=0.0, y=6.90789161545809, z=0.0) * unit.nanometer,
        openmm.vec3.Vec3(x=3.453945807729045, y=3.453945807729045, z=4.88461700499211)
        * unit.nanometer,
    ]
    hybrid_system.setDefaultPeriodicBoxVectors(*box_vectors)
    for force in hybrid_system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setNonbondedMethod(openmm.NonbondedForce.PME)
            force.setCutoffDistance(
                force_field.get_parameter_handler("Electrostatics").cutoff.m_as(
                    offunit.nanometer
                )
                * unit.nanometer
            )
            force.setUseDispersionCorrection(False)
            force.setUseSwitchingFunction(False)
        elif isinstance(force, openmm.CustomNonbondedForce):
            force.setCutoffDistance(
                force_field.get_parameter_handler("Electrostatics").cutoff.m_as(
                    offunit.nanometer
                )
                * unit.nanometer
            )
            force.setNonbondedMethod(force.CutoffPeriodic)
            force.setUseLongRangeCorrection(False)
            force.setUseSwitchingFunction(False)

    # make sure both end state systems have the same cutoff method and distance
    for end_state in [
        hybrid_topology_factory._old_system,
        hybrid_topology_factory._new_system,
    ]:
        end_state.setDefaultPeriodicBoxVectors(*box_vectors)
        for force in end_state.getForces():
            if isinstance(force, openmm.NonbondedForce):
                force.setNonbondedMethod(openmm.NonbondedForce.PME)
                force.setCutoffDistance(
                    force_field.get_parameter_handler("Electrostatics").cutoff.m_as(
                        offunit.nanometer
                    )
                    * unit.nanometer
                )
                force.setUseDispersionCorrection(False)
                force.setUseSwitchingFunction(False)


@pytest.fixture(scope="module")
def htf_chlorobenzene_fluorobenzene(
    chlorobenzene,
    fluorobenzene,
    chlorobenzene_to_fluorobenzene_mapping,
    t4_lysozyme_solvated,
):
    """Generate the htf for chlorobenzene to fluorobenzene."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure we interpolate the 1-4 exceptions involving dummy atoms if present
    settings.alchemical_settings.turn_off_core_unique_exceptions = True
    small_ff = settings.forcefield_settings.small_molecule_forcefield
    if ".offxml" not in small_ff:
        small_ff += ".offxml"
    ff = ForceField(small_ff)
    chloro_openff = chlorobenzene.to_openff()
    chloro_charges = chloro_openff.partial_charges.m_as(offunit.elementary_charge)
    chloro_labels = ff.label_molecules(chloro_openff.to_topology())[0]
    fluoro_openff = fluorobenzene.to_openff()
    fluoro_charges = fluoro_openff.partial_charges.m_as(offunit.elementary_charge)
    fluoro_labels = ff.label_molecules(fluoro_openff.to_topology())[0]
    htf = make_htf(
        mapping=chlorobenzene_to_fluorobenzene_mapping,
        settings=settings,
        protein=t4_lysozyme_solvated,
    )
    hybrid_system = htf.hybrid_system

    apply_box_vectors_and_fix_nn_force(hybrid_topology_factory=htf, force_field=ff)

    forces = {force.getName(): force for force in hybrid_system.getForces()}

    return {
        "htf": htf,
        "hybrid_system": hybrid_system,
        "forces": forces,
        "chloro_labels": chloro_labels,
        "fluoro_labels": fluoro_labels,
        "mapping": chlorobenzene_to_fluorobenzene_mapping,
        "chlorobenzene": chlorobenzene,
        "fluorobenzene": fluorobenzene,
        "chloro_charges": chloro_charges,
        "fluoro_charges": fluoro_charges,
        "electrostatic_scale": ff.get_parameter_handler("Electrostatics").scale14,
        "vdW_scale": ff.get_parameter_handler("vdW").scale14,
        "force_field": ff,
    }


@pytest.fixture(scope="module")
def htf_chlorobenzene_benzene(
    chlorobenzene, benzene, chlorobenzene_to_benzene_mapping, t4_lysozyme_solvated
):
    """Generate the htf for chlorobenzene to benzene with interpolate 1-4s on!"""
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure we interpolate the 1-4 exceptions involving dummy atoms
    settings.alchemical_settings.turn_off_core_unique_exceptions = True
    small_ff = settings.forcefield_settings.small_molecule_forcefield
    if ".offxml" not in small_ff:
        small_ff += ".offxml"
    ff = ForceField(small_ff)

    chloro_openff = chlorobenzene.to_openff()
    chloro_charges = chloro_openff.partial_charges.m_as(offunit.elementary_charge)
    chloro_labels = ff.label_molecules(chloro_openff.to_topology())[0]
    benzene_openff = benzene.to_openff()
    benzene_charges = benzene_openff.partial_charges.m_as(offunit.elementary_charge)
    benzene_labels = ff.label_molecules(benzene_openff.to_topology())[0]
    htf = make_htf(
        mapping=chlorobenzene_to_benzene_mapping,
        settings=settings,
        protein=t4_lysozyme_solvated,
    )
    hybrid_system = htf.hybrid_system

    apply_box_vectors_and_fix_nn_force(hybrid_topology_factory=htf, force_field=ff)

    forces = {force.getName(): force for force in hybrid_system.getForces()}

    return {
        "htf": htf,
        "hybrid_system": hybrid_system,
        "forces": forces,
        "chloro_labels": chloro_labels,
        "benzene_labels": benzene_labels,
        "mapping": chlorobenzene_to_benzene_mapping,
        "chlorobenzene": chlorobenzene,
        "benzene": benzene,
        "chloro_charges": chloro_charges,
        "benzene_charges": benzene_charges,
        "electrostatic_scale": ff.get_parameter_handler("Electrostatics").scale14,
        "vdW_scale": ff.get_parameter_handler("vdW").scale14,
        "force_field": ff,
    }


@pytest.fixture
def ghostly_output_chloroethane_to_ethane():
    with resources.files("htf.tests.data") as f:
        return (f / "ghostly" / "chloroethane_ethane.json").as_posix()


@pytest.fixture
def ghostly_output_hsp90_11_to_12():
    with resources.files("htf.tests.data") as f:
        return (f / "ghostly" / "hsp90_11_hsp90_12.json").as_posix()


@pytest.fixture()
def chloroethane_ethane_ghostly_modifications():
    return {
        "lambda_0": {
            "removed_angles": [],
            "removed_dihedrals": [
                (6, 2, 1, 8),
                (7, 2, 1, 8),
                (5, 2, 1, 8),
            ],
            "stiffened_angles": [],
            "softened_angles": {
                (4, 1, 8): {"k": 5.0, "theta0": 2.1329879928506985},
                (2, 1, 8): {"k": 5.0, "theta0": 1.9101780962868387},
                (3, 1, 8): {"k": 5.0, "theta0": 1.6792925223888766},
            },
        },
        "lambda_1": {
            "removed_angles": [],
            "removed_dihedrals": [(0, 1, 2, 6), (0, 1, 2, 7), (0, 1, 2, 5)],
            "stiffened_angles": [],
            "softened_angles": {
                (0, 1, 4): {"k": 5.0, "theta0": 2.230554843104804},
                (0, 1, 2): {"k": 5.0, "theta0": 1.8982497804952694},
                (0, 1, 3): {"k": 5.0, "theta0": 1.543215972886785},
            },
        },
    }


@pytest.fixture()
def hsp90_11_to_12_ghostly_modifications():
    return {
        "lambda_0": {
            "removed_angles": [(12, 19, 33)],
            "removed_dihedrals": [
                (12, 19, 33, 36),
                (12, 19, 33, 34),
                (12, 19, 33, 35),
                (9, 12, 19, 33),
                (14, 12, 19, 33),
            ],
            "stiffened_angles": [(18, 19, 33)],
            "softened_angles": {},
        },
        "lambda_1": {
            "removed_angles": [(12, 19, 32)],
            "removed_dihedrals": [(9, 12, 19, 32), (14, 12, 19, 32)],
            "stiffened_angles": [(18, 19, 32)],
            "softened_angles": {},
        },
    }


@pytest.fixture()
def htf_hsp90_11_to_12(
    hsp90_11,
    hsp90_12,
    hsp90_11_to_12_mapping,
):
    """Generate the htf for HSP90 ligand 11 to ligand 12."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure we interpolate the 1-4 exceptions involving dummy atoms if present
    settings.alchemical_settings.turn_off_core_unique_exceptions = True
    small_ff = settings.forcefield_settings.small_molecule_forcefield
    if ".offxml" not in small_ff:
        small_ff += ".offxml"
    ff = ForceField(small_ff)
    hsp90_11_openff = hsp90_11.to_openff()
    hsp90_11_charges = hsp90_11_openff.partial_charges.m_as(offunit.elementary_charge)
    hsp90_11_labels = ff.label_molecules(hsp90_11_openff.to_topology())[0]
    hsp90_12_openff = hsp90_12.to_openff()
    hsp90_12_charges = hsp90_12_openff.partial_charges.m_as(offunit.elementary_charge)
    hsp90_12_labels = ff.label_molecules(hsp90_12_openff.to_topology())[0]
    htf = make_htf(
        mapping=hsp90_11_to_12_mapping,
        settings=settings,
    )
    hybrid_system = htf.hybrid_system

    apply_box_vectors_and_fix_nn_force(hybrid_topology_factory=htf, force_field=ff)

    forces = {force.getName(): force for force in hybrid_system.getForces()}

    return {
        "htf": htf,
        "hybrid_system": hybrid_system,
        "forces": forces,
        "hsp90_11_labels": hsp90_11_labels,
        "hsp90_12_labels": hsp90_12_labels,
        "mapping": hsp90_11_to_12_mapping,
        "hsp90_11": hsp90_11,
        "hsp90_12": hsp90_12,
        "hsp90_11_charges": hsp90_11_charges,
        "hsp90_12_charges": hsp90_12_charges,
        "electrostatic_scale": ff.get_parameter_handler("Electrostatics").scale14,
        "vdW_scale": ff.get_parameter_handler("vdW").scale14,
        "force_field": ff,
    }
