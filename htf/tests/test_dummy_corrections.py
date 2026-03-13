import pathlib

from htf.utils import _derive_dummy_junction_corrections, make_htf, _draw_dummy_corrections
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol


def test_terminal_corrections(ejm_50_to_ejm_55_mapping):
    corrections = _derive_dummy_junction_corrections(ejm_50_to_ejm_55_mapping, "openff-2.0.0.offxml")
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure to use the same forcefield
    settings.forcefield_settings.small_molecule_forcefield = "openff-2.0.0.offxml"
    htf = make_htf(ejm_50_to_ejm_55_mapping, settings=settings, corrections=corrections)
    # print(corrections)
    # _draw_dummy_corrections(
    #     mapping=ejm_50_to_ejm_55_mapping,
    #     corrections=corrections,
    #     force_field="openff-2.0.0.offxml",
    #     output_dir=pathlib.Path("ejm_50_to_ejm_55_corrections")
    # )


def test_toluene_to_pyridine_corrections(toluene_to_pyridine_mapping):
    corrections = _derive_dummy_junction_corrections(toluene_to_pyridine_mapping, "openff-2.0.0.offxml")
    settings = RelativeHybridTopologyProtocol.default_settings()
    # make sure to use the same forcefield
    settings.forcefield_settings.small_molecule_forcefield = "openff-2.0.0.offxml"
    htf = make_htf(toluene_to_pyridine_mapping, settings=settings, corrections=corrections)
    # print(corrections)
    # _draw_dummy_corrections(
    #     mapping=toluene_to_pyridine_mapping,
    #     corrections=corrections,
    #     force_field="openff-2.0.0.offxml",
    #     output_dir=pathlib.Path("toluene_to_pyridine_corrections")
    # )


def test_propane_to_dimethyl_ether_corrections(propane_to_dimethyl_ether_mapping):
    corrections = _derive_dummy_junction_corrections(propane_to_dimethyl_ether_mapping, "openff-2.0.0.offxml")
    # print(corrections)


def test_triple_corrections_tyk2(ejm_31_to_ejm_42_mapping):
    corrections = _derive_dummy_junction_corrections(ejm_31_to_ejm_42_mapping, "openff-2.0.0.offxml")
    # print(corrections)