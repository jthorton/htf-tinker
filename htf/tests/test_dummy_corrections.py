import pathlib
from dataclasses import fields

import pytest

from htf.utils import (
    _derive_dummy_junction_corrections,
    make_htf,
    _draw_dummy_corrections,
    _find_dummy_junctions_rdkit,
    _find_free_rotors,
    _find_improper_corrections,
    JunctionData,
    _prune_multiple_path_dihedrals,
    CorrectionData,
    _get_heaviest_dihedral_anchor,
    _derive_uniform_valence_pruning,
)
from openff.toolkit import Molecule, ForceField
from openff.units import unit as offunit
from openmm import unit as ommunit
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe import SolventComponent


def test_single_terminal_corrections(ejm_50_to_ejm_42_mapping, tmp_path):
    """Test single dummy group terminal corrections for a single and triple branch type."""
    corrections = _derive_uniform_valence_pruning(
        ejm_50_to_ejm_42_mapping, "openff-2.0.0.offxml"
    )
    # check the single terminal dummy group correction first ejm_50
    state_1 = corrections["lambda_1"]
    # check the blank corrections first
    assert not state_1.removed_angles
    assert not state_1.removed_impropers
    assert not state_1.softened_angles
    assert not state_1.stiffened_angles
    assert not state_1.stiffened_dihedrals
    # there should be 2 redundant dihedrals removed
    assert state_1.removed_dihedrals == {
        # both dihedrals should terminate in the dummy atom 32
        frozenset({32, 19, 29, 31}),
        frozenset({32, 19, 30, 31}),
    }

    # here we have the same type of correction but with three dummy groups attached to the core
    state_0 = corrections["lambda_0"]
    # check the blank corrections first
    assert not state_0.removed_angles
    assert not state_0.softened_angles
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_impropers
    # there should 6 redundant dihedrals removed 2 per dummy group atom {32, 33, 34}
    assert state_0.removed_dihedrals == {
        frozenset({32, 19, 29, 31}),  # dummy group 1 redundant dihedrals
        frozenset({32, 19, 30, 31}),
        frozenset({33, 19, 29, 31}),  # dummy group 2 redundant dihedrals
        frozenset({33, 19, 30, 31}),
        frozenset({34, 19, 29, 31}),  # dummy group 3 redundant dihedrals
        frozenset({34, 19, 30, 31}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_50_to_ejm_42_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_50_to_ejm_42_corrections",
    )


def test_tyk2_dummy_group_interaction_corrections(ejm_50_to_ejm_55_mapping, tmp_path):
    """Test corrections which should leave interactions between two dummy groups not on the same core junction atom (ejm_50)."""
    corrections = _derive_uniform_valence_pruning(
        ejm_50_to_ejm_55_mapping, "openff-2.0.0.offxml"
    )
    # # check the simple non-corrected end state first ejm_55
    # # we need no corrections for this triple terminal junction due to the geometry of the junction and lack of dummy groups
    # state_0 = corrections["lambda_0"]
    # # make sure all corrections are blank
    # assert not state_0.softened_angles
    # assert not state_0.stiffened_angles
    # assert not state_0.stiffened_dihedrals
    # assert not state_0.removed_impropers
    # assert not state_0.removed_dihedrals
    # assert not state_0.removed_angles
    #
    # # check the more complicated interacting dummy group case ejm_50
    # # this has a terminal and dual anchor double branch dummy group which are connected
    # state_1 = corrections["lambda_1"]
    # # first check the blank corrections
    # assert not state_1.stiffened_dihedrals
    # assert not state_1.stiffened_angles
    # assert not state_1.softened_angles
    # assert not state_1.removed_impropers
    # # make sure 2 redundant angles are removed around the dual anchor junction involving dummy atoms {29, 30}
    # assert state_1.removed_angles == {
    #     frozenset({19, 29, 31}),  # dummy atom 1 on dual anchor junction
    #     frozenset({19, 30, 31}),  # dummy atom 2 on dual anchor junction
    # }
    # # make sure 2 redundant dihedrals are removed which terminate in the dummy atoms {29, 30}
    # assert state_1.removed_dihedrals == {
    #     frozenset({16, 17, 19, 29}),  # redundant dihedral on dummy atom 1
    #     frozenset({16, 17, 19, 30}),  # redundant dihedral on dummy atom 2
    # }
    print(corrections)
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_50_to_ejm_55_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=pathlib.Path("ejm_50_to_ejm_55_corrections_unified"),
    )


def test_tyk2_dual_junction_large_single_branch_corrections(
    ejm_31_to_ejm_49_mapping, tmp_path
):
    """Test the dual junction single branch corrections for a large dummy group which has an improper in the junction (ejm_49).
    Here we want to make sure the improper spanning the core and the dummy junction atom is removed while the improper
    spanning the dummy group and the core junction atom is kept.
    """
    corrections = _derive_uniform_valence_pruning(
        ejm_31_to_ejm_49_mapping, "openff-2.0.0.offxml"
    )
    # check the simple dual anchor single branch correction first ejm_31
    state_1 = corrections["lambda_1"]
    # check blank corrections first
    assert not state_1.softened_angles
    assert not state_1.stiffened_angles
    assert not state_1.stiffened_dihedrals
    # make sure the improper spanning the core atoms and the dummy junction atom is removed
    assert state_1.removed_impropers == {
        frozenset({16, 17, 18, 19})  # dummy junction improper dummy atom 19
    }
    # make sure the redundant angle is removed
    assert state_1.removed_angles == {frozenset({17, 18, 19})}
    # make sure the 4 redundant dihedrals are removed
    assert state_1.removed_dihedrals == {
        frozenset({16, 17, 19, 28}),  # dummy junction redundant dihedral
        frozenset({17, 18, 19, 29}),  # dummy group redundant dihedrals
        frozenset({17, 18, 19, 30}),
        frozenset({17, 18, 19, 31}),
    }

    # check the dual anchor single branch large group ejm_49
    state_0 = corrections["lambda_0"]
    # check the blank corrections first
    assert not state_0.softened_angles
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    # make sure the redundant angle is removed linking the dummy junction atom 31
    assert state_0.removed_angles == {frozenset({17, 18, 31})}
    # make sure only 1 improper is removed which spanns the core and dummy junction
    assert state_0.removed_impropers == {frozenset({16, 17, 18, 31})}
    # make sure the 3 redundant dihedrals are removed
    assert state_0.removed_dihedrals == {
        frozenset({16, 17, 27, 31}),  # dummy junction redundant dihedral
        frozenset({17, 18, 30, 31}),  # dummy group redundant dihedrals
        frozenset({32, 17, 18, 31}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_31_to_ejm_49_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_31_to_ejm_49_corrections",
    )


def test_find_toluene_to_pyridine_corrections(toluene_to_pyridine_mapping, tmp_path):
    """Simple test case with corrections at one end (toluene), dual anchor single branch correction type."""
    corrections = _derive_uniform_valence_pruning(
        toluene_to_pyridine_mapping, "openff-2.0.0.offxml"
    )
    # there should be no corrections in lambda_0
    state_0 = corrections["lambda_0"]
    for f in fields(state_0):
        assert not getattr(state_0, f.name)

    state_1 = corrections["lambda_1"]
    # check blank corrections first
    assert not state_1.softened_angles
    assert not state_1.stiffened_angles
    assert not state_1.stiffened_dihedrals
    # check removed angle
    assert state_1.removed_angles == {frozenset({0, 1, 2})}
    # check the removed improper on the junction
    assert state_1.removed_impropers == {frozenset({0, 1, 2, 6})}
    # check removed dihedrals
    assert state_1.removed_dihedrals == {
        frozenset({0, 1, 2, 3}),  # dummy anchor constraint
        frozenset({0, 1, 2, 10}),  # dummy anchor constraint
        frozenset({0, 1, 6, 14}),  # dummy anchor constraint
        frozenset({0, 1, 2, 8}),  # dummy group dual anchor constraint
        frozenset({0, 1, 2, 9}),  # dummy group dual anchor constraint
        frozenset({0, 1, 2, 7}),  # dummy group dual anchor constraint
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=toluene_to_pyridine_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "toluene_to_pyridine_corrections",
    )


def test_propane_to_dimethyl_ether_corrections(
    propane_to_dimethyl_ether_mapping, tmp_path
):
    """Free rotor test case with corrections at one end (dimethyl ether), dual anchor two branch correction type."""
    corrections = _derive_uniform_valence_pruning(
        propane_to_dimethyl_ether_mapping, "openff-2.0.0.offxml"
    )
    # there should be no corrections at lambda_0
    state_0 = corrections["lambda_0"]
    for f in fields(state_0):
        assert not getattr(state_0, f.name)

    state_1 = corrections["lambda_1"]
    # check the blank corrections first
    assert not state_1.softened_angles
    assert not state_1.removed_impropers
    assert not state_1.stiffened_angles
    # check we have two removed angles one for each dummy atom
    assert state_1.removed_angles == {
        frozenset({0, 1, 9}),  # extra dummy angle constraint
        frozenset({0, 1, 10}),  # extra other dummy angle constraint
    }
    # check we have 10 removed in total 5 for each group
    assert state_1.removed_dihedrals == {
        frozenset({0, 9, 5, 7}),  # dummy group 1 redundant dihedrals
        frozenset({0, 1, 2, 9}),
        frozenset({0, 1, 3, 9}),
        frozenset({0, 9, 5, 6}),
        frozenset({0, 1, 4, 9}),
        frozenset({0, 10, 5, 6}),  # dummy group 2 redundant dihedrals
        frozenset({0, 1, 2, 10}),
        frozenset({0, 1, 10, 4}),
        frozenset({0, 10, 5, 7}),
        frozenset({0, 1, 10, 3}),
    }
    # check we have two stiffened dihedrals as the anchor is a free rotor
    # they should also terminate in the same atom
    assert state_1.stiffened_dihedrals == {
        # note both dihedrals terminate in atom 8
        frozenset({8, 0, 5, 9}),  # dummy group 1
        frozenset({8, 0, 10, 5}),  # dummy group 2
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=propane_to_dimethyl_ether_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "propane_to_dimethyl_ether_corrections",
    )


def test_triple_corrections_tyk2(ejm_31_to_ejm_42_mapping, tmp_path):
    """Triple non-planar junction correction type at both ends, end stateB (lambda_0 corrections) is a dummy group with pruned group anchor dihedral constraints."""
    corrections = _derive_uniform_valence_pruning(
        ejm_31_to_ejm_42_mapping, "openff-2.0.0.offxml"
    )
    # check simple case first which is for ejm_31 and lambda_1
    state_1 = corrections["lambda_1"]
    # check the blank corrections first
    assert not state_1.stiffened_dihedrals
    assert not state_1.removed_impropers
    assert not state_1.stiffened_angles
    assert not state_1.removed_angles
    # check the softened angles
    assert state_1.softened_angles == {
        # should all involve the dummy atom 30
        frozenset({19, 30, 31}),
        frozenset({19, 29, 30}),
        frozenset({17, 19, 30}),
    }
    # check that the 2 possible dihedrals coupling the dummy are removed
    assert state_1.removed_dihedrals == {
        # both dihedrals terminate in 30
        frozenset({17, 18, 19, 30}),
        frozenset({16, 17, 19, 30}),
    }

    state_0 = corrections["lambda_0"]
    # check the blank corrections first
    assert not state_0.stiffened_dihedrals
    assert not state_0.stiffened_angles
    assert not state_0.removed_impropers
    assert not state_0.removed_angles
    # check the softened angle
    assert state_0.softened_angles == {
        # should all involve the dummy atom 31
        frozenset({17, 19, 31}),
        frozenset({19, 30, 31}),
        frozenset({19, 29, 31}),
    }
    # check the removed dihedrals
    assert state_0.removed_dihedrals == {
        frozenset({16, 17, 19, 31}),  # dummy junction redundant dihedrals
        frozenset({17, 18, 19, 31}),
        frozenset({32, 19, 30, 31}),  # dummy group redundant dihedrals anchor 1
        frozenset({33, 19, 30, 31}),
        frozenset({34, 19, 30, 31}),
        frozenset({32, 19, 29, 31}),  # dummy group redundant dihedrals anchor 2
        frozenset({33, 19, 29, 31}),
        frozenset({34, 19, 29, 31}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_31_to_ejm_42_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_31_to_ejm_42_corrections",
    )


def test_terminal_co_linear_tyk2_corrections(jmc_28_to_jmc_30_mapping, tmp_path):
    """Triple terminal junction to a terminal co-linear junction (jmc_30)"""
    corrections = _derive_uniform_valence_pruning(
        jmc_28_to_jmc_30_mapping, "openff-2.0.0.offxml"
    )
    # start with the simple single co-linear terminal group for jmc_30
    state_0 = corrections["lambda_0"]
    # check the blank corrections
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_impropers
    assert not state_0.stiffened_angles
    assert not state_0.removed_angles
    assert not state_0.softened_angles
    # make sure 2 redundant dihedrals are removed which involve 35
    assert state_0.removed_dihedrals == {
        frozenset({34, 35, 36, 21}),
        frozenset({35, 34, 19, 21}),
    }

    # check the triple terminal junction with dummy atoms {36, 37, 38}
    state_1 = corrections["lambda_1"]
    # check the blank corrections
    assert not state_1.stiffened_angles
    assert not state_1.removed_impropers
    assert not state_1.stiffened_dihedrals
    assert not state_1.softened_angles
    assert not state_1.removed_angles
    # make sure 6 redundant dihedrals are removed which involve {36, 37, 38} 2 per dummy atom
    assert state_1.removed_dihedrals == {
        frozenset({35, 19, 36, 21}),  # dummy atom 36 redundant dihedrals
        frozenset({34, 35, 36, 21}),
        frozenset({37, 35, 19, 21}),  # dummy atom 37 redundant dihedrals
        frozenset({37, 34, 35, 21}),
        frozenset({35, 19, 21, 38}),  # dummy atom 38 redundant dihedrals
        frozenset({34, 35, 21, 38}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=jmc_28_to_jmc_30_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "jmc_28_to_jmc_30_corrections",
    )


def test_tyk2_triple_junction_cyclopropane_corrections(
    jmc_30_to_ejm_46_mapping, tmp_path
):
    """Make sure we correctly prune dummy junction anchors in 3 membered rings with multiple possible ways around the ring."""
    corrections = _derive_uniform_valence_pruning(
        jmc_30_to_ejm_46_mapping, "openff-2.0.0.offxml"
    )
    # check ejm_30 and make sure that all dihedrals in the cyclopropane group are removed
    state_1 = corrections["lambda_1"]
    # check the blank terms first
    assert not state_1.stiffened_dihedrals
    assert not state_1.removed_impropers
    assert not state_1.stiffened_angles
    assert not state_1.removed_angles
    # check the 3 softened angles which should involve dummy atom 34
    assert state_1.softened_angles == {
        frozenset({34, 19, 21}),
        frozenset({34, 36, 21}),
        frozenset({34, 20, 21}),
    }
    assert state_1.removed_dihedrals == {
        # dummy anchor dihedrals must all be removed terminating in 34
        frozenset({34, 19, 20, 21}),
        frozenset({34, 19, 21, 31}),
        frozenset({32, 34, 20, 21}),
        frozenset({33, 34, 20, 21}),
        frozenset({17, 34, 19, 21}),
        frozenset({35, 34, 19, 21}),  # remove the dummy group rotor constraints
        frozenset({34, 35, 36, 21}),
    }
    # check the ejm_46 corrections which a slightly simpler
    state_0 = corrections["lambda_0"]
    # check the blank corrections
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_impropers
    assert not state_0.stiffened_angles
    assert not state_0.removed_angles
    # check the 3 softened angles which should involve dummy atom 35
    assert state_0.softened_angles == {
        frozenset({35, 20, 21}),
        frozenset({34, 35, 21}),
        frozenset({35, 19, 21}),
    }
    assert state_0.removed_dihedrals == {
        # remove all possible torsions terminating in 35
        # there are 6 possible but 2 involve the same atoms in opposite directions around the ring so we
        # have just 5 entries
        frozenset({35, 19, 21, 31}),
        frozenset({35, 19, 20, 21}),
        frozenset({33, 35, 20, 21}),
        frozenset({32, 35, 20, 21}),
        frozenset({17, 19, 21, 35}),
    }
    _draw_dummy_corrections(
        mapping=jmc_30_to_ejm_46_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "jmc_30_to_ejm_46_corrections",
    )


def test_triple_junction_interactions_hif2a(hif2a_155_to_231_mapping):
    """Test removing triple junction dihedrals between dummy groups."""
    correction = _derive_dummy_junction_corrections(
        hif2a_155_to_231_mapping, "openff-2.0.0.offxml"
    )
    print(correction)
    _draw_dummy_corrections(
        mapping=hif2a_155_to_231_mapping,
        corrections=correction,
        force_field="openff-2.0.0.offxml",
        output_dir=pathlib.Path("hif2a_155_to_231_corrections"),
    )

def test_faah_26_to_28(faah_26_to_28_mapping):
    corrections = _derive_uniform_valence_pruning(faah_26_to_28_mapping, "openff-2.0.0.offxml")
    print(corrections)
    _draw_dummy_corrections(
        mapping=faah_26_to_28_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=pathlib.Path("faah_26_to_28_corrections_uniform"),
    )

def test_shp2_triple_junction_dummy_group_interaction_corrections(
    shp2_099_1_ex7_to_ex9, tmp_path
):
    """Test corrections for a triple non-planar junction which also has interactions with a terminal junction near a ring (example-9)."""
    corrections = _derive_uniform_valence_pruning(
        shp2_099_1_ex7_to_ex9, "openff-2.0.0.offxml"
    )
    print(corrections)
    # start with the example-9 case
    state_0 = corrections["lambda_0"]
    # first check the blank corrections
    assert not state_0.removed_impropers
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_angles
    # check softened angles should all involve the dummy atom 41
    assert state_0.softened_angles == {
        frozenset({41, 3, 4}),
        frozenset({41, 3, 21}),
        frozenset({41, 2, 3}),
    }
    # now check the removed dihedrals
    assert state_0.removed_dihedrals == {
        # all dihedrals terminating in 41 should be removed
        frozenset({41, 3, 4, 27}),
        frozenset({41, 3, 4, 26}),
        frozenset({41, 3, 4, 5}),
        frozenset({41, 3, 2, 24}),
        frozenset({41, 3, 2, 25}),
        frozenset({41, 3, 2, 1}),
        # remove dual rotor constraints in the terminal group
        frozenset({38, 21, 3, 2}),
        frozenset({36, 21, 3, 2}),
        frozenset({37, 21, 3, 2}),
    }
    # now check the other state which is a triple junction and a terminal group
    state_1 = corrections["lambda_1"]
    # check the blank corrections first
    assert not state_1.removed_impropers
    assert not state_1.stiffened_angles
    assert not state_1.stiffened_dihedrals
    assert not state_1.removed_angles
    assert state_1.softened_angles == {
        # all angles should terminate in the dummy atom 22
        frozenset({22, 3, 4}),
        frozenset({22, 3, 2}),
        frozenset({22, 3, 21}),
    }
    assert state_1.removed_dihedrals == {
        # all dihedrals terminating in 22 should be removed not coupled to the other dummy group
        frozenset({22, 3, 2, 25}),
        frozenset({22, 3, 2, 26}),
        frozenset({22, 3, 2, 1}),
        frozenset({22, 3, 4, 28}),
        frozenset({22, 3, 4, 27}),
        frozenset({22, 3, 4, 5}),
        # dummy group 1 dual rotor constraints {39, 40, 41}
        frozenset({39, 22, 3, 2}),
        frozenset({39, 22, 3, 21}),
        frozenset({40, 22, 3, 2}),
        frozenset({40, 22, 3, 21}),
        frozenset({41, 22, 3, 2}),
        frozenset({41, 22, 3, 21}),
        # dummy group 2 dual rotor constraints {42, 37, 38}
        frozenset({37, 21, 3, 2}),
        frozenset({38, 21, 3, 2}),
        frozenset({42, 21, 3, 2}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=shp2_099_1_ex7_to_ex9,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=pathlib.Path("shp2_ex7_to_ex9_corrections"),
    )


def test_higher_order_corrections(sulfur_hexafluoride_to_tetrafluoride_mapping):
    """Test finding higher order corrections using hexafluoride to tetrafluoride as an example."""
    from openff.toolkit import ForceField

    # patch the force field with a dummy parameter for sulfur hexafluoride
    ff = ForceField("openff-2.0.0.offxml")
    angle_handler = ff.get_parameter_handler("Angles")
    angle_handler.add_parameter(
        {
            "smirks": "[#9:1]-[#16:2]-[#9:3]",
            "angle": 90 * offunit.degree,
            "k": 100 * offunit.kilocalorie_per_mole / offunit.radian**2,
        },
        before=0,
    )
    corrections = _derive_uniform_valence_pruning(
        sulfur_hexafluoride_to_tetrafluoride_mapping, ff.to_string()
    )
    # check the tetrafluoride case first this should be blank
    state_0 = corrections["lambda_0"]
    assert not state_0.removed_impropers
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_angles
    assert not state_0.softened_angles
    assert not state_0.removed_dihedrals

    # now check the sulfur hexafluoride end which should be the higher order junction
    state_1 = corrections["lambda_1"]
    # check the blank terms first
    assert not state_1.removed_impropers
    assert not state_1.stiffened_angles
    assert not state_1.stiffened_dihedrals
    assert not state_1.removed_dihedrals
    # check that there are 2 removed redundant angles done by the higher order function
    # should involve the dummy atoms {3, 6}
    assert state_1.removed_angles == {
        frozenset({1, 3, 5}),  # dummy atom 1
        frozenset({1, 5, 6}),  # dummy atom 2
    }
    # there should be 6 softened angles related to the 2 dummy atoms done by the triple correction
    assert state_1.softened_angles == {
        frozenset({1, 2, 3}),  # dummy atom 1
        frozenset({0, 1, 3}),
        frozenset({1, 3, 4}),
        frozenset({1, 2, 6}),  # dummy atom 2
        frozenset({0, 1, 6}),
        frozenset({1, 4, 6}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=sulfur_hexafluoride_to_tetrafluoride_mapping,
        corrections=corrections,
        force_field=ff.to_string(),
        output_dir=pathlib.Path("sulfur_hexafluoride_corrections"),
    )


# def test_higher_order_group_corrections(pentafluorosulfanylbenzene_to_sulfur_hexafluoride_mapping):
#     from openff.toolkit import ForceField
#     # patch the force field with a dummy parameter for sulfur hexafluoride
#     ff = ForceField("openff-2.0.0.offxml")
#     angle_handler = ff.get_parameter_handler("Angles")
#     angle_handler.add_parameter(
#         {
#             # fake angle for the F-S-F/C angles
#             "smirks": "[#9:1]-[#16:2]-[*:3]",
#             "angle": 90 * offunit.degree,
#             "k": 100 * offunit.kilocalorie_per_mole / offunit.radian ** 2,
#         },
#         before=0
#     )
#     proper_handler = ff.get_parameter_handler("ProperTorsions")
#     proper_handler.add_parameter(
#         {
#             # fake proper torsion for the F-S-C-C torsions
#             "smirks": "[#9:1]-[#16:2]-[*:3]-[*:4]",
#             "phase": [0.0 * offunit.degree],
#             "periodicity":  [1,],
#             "k": [1 * offunit.kilocalorie_per_mole],
#         },
#         before=0
#     )
#     corrections = _derive_dummy_junction_corrections(pentafluorosulfanylbenzene_to_sulfur_hexafluoride_mapping, ff.to_string())
#     print(corrections)


def test_find_dummy_junctions_no_dummies(chloroethane_to_fluoroethane_mapping):
    """Make sure no dummy groups are found when we have a 1:1 mapping."""
    # check each end state
    for mol in [
        chloroethane_to_fluoroethane_mapping.componentA,
        chloroethane_to_fluoroethane_mapping.componentB,
    ]:
        dummy_junctions = _find_dummy_junctions_rdkit(mol.to_rdkit(), set())
        assert not dummy_junctions


def test_find_dummy_junctions_single_dummy_atom(chloroethane_to_ethane_mapping):
    """Make sure we can find the single dummy atom junctions for both cases."""
    for mol, core in zip(
        [
            chloroethane_to_ethane_mapping.componentA.to_rdkit(),
            chloroethane_to_ethane_mapping.componentB.to_rdkit(),
        ],
        [
            chloroethane_to_ethane_mapping.componentA_to_componentB.keys(),
            chloroethane_to_ethane_mapping.componentA_to_componentB.values(),
        ],
    ):
        dummy_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core}
        dummy_junctions = _find_dummy_junctions_rdkit(mol, dummy_atoms)
        assert len(dummy_junctions) == 1
        junction = dummy_junctions[0]
        assert junction.dummies == {0}
        assert junction.junction_atom == 1
        assert junction.physical == {2, 3, 4}


def test_find_dummy_junctions_dummy_group(toluene_to_pyridine_mapping):
    """Make sure we can find the single dummy atom junction when we have a dummy group."""
    # only toluene should need corrections
    mol = toluene_to_pyridine_mapping.componentA.to_rdkit()
    core = set(toluene_to_pyridine_mapping.componentA_to_componentB.keys())
    dummy_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core}
    dummy_junctions = _find_dummy_junctions_rdkit(mol, dummy_atoms)
    assert len(dummy_junctions) == 1
    junction = dummy_junctions[0]
    # dual anchor with a single branch
    assert junction.dummies == {0}
    assert junction.junction_atom == 1
    assert junction.physical == {2, 6}


def test_find_dummy_junctions_same_ring(isoquinoline):
    """Make sure that dummy atoms in the same ring as the anchor are not tagged - this avoids bond breaking cases."""
    # make a fake mapping which only maps one ring of the isoquinoline
    mol = isoquinoline.to_rdkit()
    # the 4 carbons and their hydrogens from the broken ring
    dummy_atoms = {6, 7, 8, 9, 13, 14, 15, 16}
    dummy_junctions = _find_dummy_junctions_rdkit(mol, dummy_atoms)
    assert not dummy_junctions


def test_find_dummy_junctions_multiple_dummy_branches(jmc_28_to_jmc_30_mapping):
    """Make sure we can find the multiple dummy junctions connected to a single core anchor"""
    mol = jmc_28_to_jmc_30_mapping.componentA.to_rdkit()
    core = set(jmc_28_to_jmc_30_mapping.componentA_to_componentB.keys())
    dummy_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core}
    dummy_junctions = _find_dummy_junctions_rdkit(mol, dummy_atoms)
    assert len(dummy_junctions) == 1
    junction = dummy_junctions[0]
    # terminal anchor triple branch case
    assert junction.dummies == {36, 37, 38}
    assert junction.junction_atom == 35
    assert junction.physical == {21}


def test_find_dummy_junctions_multiple_dummy_groups(ejm_50_to_ejm_55_mapping):
    """Make sure we can find the multiple dummy groups with different core anchors in a single molecule ejm_50."""
    mol = ejm_50_to_ejm_55_mapping.componentA.to_rdkit()
    core = set(ejm_50_to_ejm_55_mapping.componentA_to_componentB.keys())
    dummy_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in core}
    dummy_junctions = _find_dummy_junctions_rdkit(mol, dummy_atoms)
    assert len(dummy_junctions) == 2
    junctions_by_anchor = {j.junction_atom: j for j in dummy_junctions}
    # grab the junction for atom 19 - dual anchor double branch
    junction = junctions_by_anchor[19]
    assert junction.dummies == {29, 30}
    assert junction.physical == {17, 31}
    # check the other terminal junction
    junction = junctions_by_anchor[31]
    assert junction.dummies == {32}
    assert junction.physical == {19}


@pytest.mark.parametrize(
    "smiles, groups",
    [
        pytest.param(
            "[H:3][C:1]([H:4])([H:5])[C:2]([H:6])([H:7])[H:8]",
            {2, 3, 4, 5, 6, 7},
            id="Ethane",
        ),
        pytest.param("[H:3][C:1]([H:4])([H:5])[O:2][H:6]", {2, 3, 4, 5}, id="Methanol"),
        pytest.param(
            "[H:3][C:1]([H:4])([H:5])[S:2][H:6]", {2, 3, 4, 5}, id="Methanethiol"
        ),
        pytest.param(
            "[H:3][C:1]([H:4])([H:5])[N:2]([H:6])[H:7]",
            {2, 3, 4, 5, 6},
            id="Methylamine",
        ),
    ],
)
def test_find_free_rotors(smiles, groups):
    """Make sure we can find the free rotors in the given molecule smiles."""
    mol = Molecule.from_mapped_smiles(smiles)
    free_rotors = _find_free_rotors(mol.to_rdkit())
    assert free_rotors == groups


def test_find_improper_corrections_no_improper(chloroethane_to_ethane_mapping):
    """Make sure no improper corrections are found when non are present"""
    ff = ForceField("openff-2.0.0.offxml")
    mol = chloroethane_to_ethane_mapping.componentA.to_openff()
    labels = ff.label_molecules(mol.to_topology())[0]
    # make the junction for this transformation, Cl is atom 0 which is the dummy
    junction = JunctionData(
        junction_atom=1,
        dummies={
            0,
        },
        physical={2, 3, 4},
    )
    corrections = _find_improper_corrections(junction, labels)
    assert not corrections.removed_impropers


def test_find_improper_corrections_single_improper(toluene_to_pyridine_mapping):
    """Make sure the single improper correction is found which spans a dummy and physical anchor atoms."""
    ff = ForceField("openff-2.0.0.offxml")
    mol = toluene_to_pyridine_mapping.componentA.to_openff()
    labels = ff.label_molecules(mol.to_topology())[0]
    # make the junction for this transformation, Cl is atom 0 which is the dummy
    junction = JunctionData(
        junction_atom=1,
        dummies={
            0,
        },
        physical={2, 6},
    )
    corrections = _find_improper_corrections(junction, labels)
    # the removed improper should cover all junction atoms
    assert corrections.removed_impropers == {frozenset({0, 1, 2, 6})}


def test_find_improper_corrections_double_improper(ejm_31_to_ejm_49_mapping):
    """Make sure only a single improper of two available are removed which spans a dummy and physical anchor atoms."""
    ff = ForceField("openff-2.0.0.offxml")
    mol = ejm_31_to_ejm_49_mapping.componentB.to_openff()
    labels = ff.label_molecules(mol.to_topology())[0]
    junction = JunctionData(
        junction_atom=17,
        dummies={
            31,
        },
        physical={16, 18},
    )
    corrections = _find_improper_corrections(junction, labels)
    assert corrections.removed_impropers == {frozenset({16, 17, 18, 31})}


def test_prune_dihedrals_single_dihedral():
    """Make sure no dihedrals are pruned when only a single possible dihedral is available."""
    corrections = CorrectionData(
        removed_dihedrals={
            frozenset({0, 1, 2, 3}),
        }
    )
    # some fake target dihedrals
    target_dihedrals = {(0, 1, 2, 3), (0, 1, 2, 4)}
    corrections = _prune_multiple_path_dihedrals(target_dihedrals, corrections)
    assert corrections.removed_dihedrals == {
        frozenset({0, 1, 2, 3}),
    }


def test_prune_dihedrals_multiple_dihedrals():
    """Make sure a single dihedral is added to the removed list when two are available."""
    corrections = CorrectionData()
    # this is the case from the docstring of the prune method
    target_dihedrals = {(0, 1, 2, 4), (0, 1, 3, 4)}
    corrections = _prune_multiple_path_dihedrals(target_dihedrals, corrections)
    assert corrections.removed_dihedrals == {
        # the dihedral with the larger sum of indices should be removed
        frozenset({0, 1, 3, 4}),
    }


def test_get_heaviest_anchor_dihedral_no_ring(ejm_42_to_ejm_54_mapping):
    """Make sure we can find the expected heaviest anchor dihedral when the options are not in a ring."""
    # take the ejm_42 end state and one of the dual anchor junctions
    mol = ejm_42_to_ejm_54_mapping.componentA.to_rdkit()
    junction_atom = 19
    # junction atom, dummies and physical atoms are excluded
    excluded_atoms = {19, 29, 30, 31, 17}
    dummy_junction_dihedrals = {
        (29, 19, 31, 34),  # dummy group 1
        (29, 19, 31, 32),
        (29, 19, 17, 18),
        (29, 19, 17, 16),
        (30, 19, 31, 34),  # dummy group 2
        (30, 19, 31, 32),
        (30, 19, 17, 18),
        (30, 19, 17, 16),
    }
    heavy_atom = _get_heaviest_dihedral_anchor(
        dummy_junction_dihedrals, excluded_atoms, mol, junction_atom
    )
    assert heavy_atom == 18


def test_get_heaviest_anchor_in_ring(chlorobenzene):
    """Make sure that an anchor in a ring is picked if the core junction is the same ring as the options even if the non-ring
    atom is heavier"""
    mol = chlorobenzene.to_rdkit()
    # pick a hydrogen next to the CL as the dummy atom (7)
    junction_atom = 2
    excluded_atoms = {2, 7, 1, 3}
    dummy_junction_dihedrals = {(7, 2, 1, 0), (7, 2, 1, 6), (7, 2, 3, 8), (7, 2, 3, 4)}
    heavy_atom = _get_heaviest_dihedral_anchor(
        dummy_junction_dihedrals, excluded_atoms, mol, junction_atom
    )
    # make sure its atom 6 a carbon despite the CL (0) being the heaviest option
    assert heavy_atom == 6


#################################################################################
# Test applying the corrections to the HTF
#################################################################################


def test_toluene_pyridine_angle_corrections_htf(toluene_to_pyridine_mapping):
    """Make sure that the angle corrections are correctly applied while making the HTF."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    corrections = _derive_dummy_junction_corrections(
        toluene_to_pyridine_mapping,
        settings.forcefield_settings.small_molecule_forcefield,
    )

    htf = make_htf(toluene_to_pyridine_mapping, settings, corrections=corrections)
    toluene = toluene_to_pyridine_mapping.componentA.to_openff()
    pyridine = toluene_to_pyridine_mapping.componentB.to_openff()
    ff = ForceField(settings.forcefield_settings.small_molecule_forcefield + ".offxml")
    toluene_labels = ff.label_molecules(toluene.to_topology())[0]
    pyridine_labels = ff.label_molecules(pyridine.to_topology())[0]
    mapping = toluene_to_pyridine_mapping.componentA_to_componentB

    corrected_system = htf.hybrid_system
    # extract the forces
    forces = {force.getName(): force for force in corrected_system.getForces()}
    # only the angle and bond forces should be affected
    # angles not interpolated should be in the standard force
    standard_angle_force = forces["HarmonicAngleForce"]
    # there should only be 7 angles not interpolated in the simulation:
    # 1 anchor angle as 1 removed
    # 3 angles in the dummy group
    # 3 angles anchoring each methyl H to the core junction atom
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 7

    for i in range(num_angles):
        p1, p2, p3, angle_eq, k = standard_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # make sure there are dummy atoms in all angles
        assert set(angle).intersection(htf._atom_classes["unique_old_atoms"])

        # now compare the parameters to the toluene labels
        toluene_angle = toluene_labels["Angles"][angle]
        # angle_eq
        assert angle_eq == toluene_angle.angle.m_as(offunit.radian) * ommunit.radian
        # k
        assert (
            k
            == toluene_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            * ommunit.kilojoule_per_mole
            / ommunit.radian**2
        )

    # there should be 17 interpolated angle terms covering fully mapped and scaled angles
    custom_angle_force = forces["CustomAngleForce"]
    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 17
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"

    # track angles turned off
    removed_angles = set()
    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        # p1, p2, p3 are the index in toluene/pyridine get the expected parameters from the labels
        if (
            p1 in htf._atom_classes["unique_old_atoms"]
            or p3 in htf._atom_classes["unique_old_atoms"]
        ):
            # this angle involves at least one old dummy atom from toluene
            toluene_angle = toluene_labels["Angles"][(p1, p2, p3)]
            # lambda_0 angle
            assert params[0] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == toluene_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )
            # lambda_1 angle stays the same
            assert params[2] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_1 k is removed
            assert params[3] == 0 * offunit.kilojoule_per_mole / offunit.radian**2
            removed_angles.add(frozenset({p1, p2, p3}))
        # there are no dummy atoms in pyridine so all others must be fully mapped
        else:
            # fully mapped angle
            toluene_angle = toluene_labels["Angles"][(p1, p2, p3)]
            e1 = mapping[p1]
            e2 = mapping[p2]
            e3 = mapping[p3]
            pyridine_angle = pyridine_labels["Angles"][(e1, e2, e3)]
            # lambda_0 angle should be the toluene angle
            assert params[0] == toluene_angle.angle.m_as(offunit.radian)
            # lambda_0 k should be the toluene k
            assert params[1] == toluene_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )
            # lambda_1 angle should be the pyridine angle
            assert params[2] == pyridine_angle.angle.m_as(offunit.radian)
            # lambda_1 k should be the pyridine k
            assert params[3] == pyridine_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )

    # make sure that the corrections extracted from the htf match what we expected
    assert removed_angles == corrections["lambda_1"].removed_angles


def test_toluene_pyridine_torsion_corrections_htf(toluene_to_pyridine_mapping):
    """Make sure that the torsion corrections are correctly applied while making the HTF."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    corrections = _derive_dummy_junction_corrections(
        toluene_to_pyridine_mapping,
        settings.forcefield_settings.small_molecule_forcefield,
    )

    htf = make_htf(toluene_to_pyridine_mapping, settings, corrections=corrections)
    toluene = toluene_to_pyridine_mapping.componentA.to_openff()
    pyridine = toluene_to_pyridine_mapping.componentB.to_openff()
    ff = ForceField(settings.forcefield_settings.small_molecule_forcefield + ".offxml")
    toluene_labels = ff.label_molecules(toluene.to_topology())[0]
    pyridine_labels = ff.label_molecules(pyridine.to_topology())[0]
    mapping = toluene_to_pyridine_mapping.componentA_to_componentB

    corrected_system = htf.hybrid_system
    # extract the forces
    forces = {force.getName(): force for force in corrected_system.getForces()}

    # check the torsions
    standard_torsion_force = forces["PeriodicTorsionForce"]
    num_standard_torsions = standard_torsion_force.getNumTorsions()
    # Note the 5 improper torsions should be conserved which should give us 15 improper potentials in total but
    # due to the degenerate order in which the improper torsions match only 2 (6 terms total) are stored in this force the rest
    # are in the interpolated force
    assert num_standard_torsions == 26
    for i in range(num_standard_torsions):
        p1, p2, p3, p4, periodicity, phase, k = (
            standard_torsion_force.getTorsionParameters(i)
        )
        torsion = (p1, p2, p3, p4)
        # if its fully mapped make sure the parameters match
        if len(htf._atom_classes["unique_old_atoms"].intersection(torsion)) == 0:
            # now compare the parameters to the toluene and pyridine labels they should be the same at both end states
            # check if we have a proper or improper torsion
            if torsion in toluene_labels["ProperTorsions"]:
                toluene_torsion = toluene_labels["ProperTorsions"][torsion]
                e1 = mapping[p1]
                e2 = mapping[p2]
                e3 = mapping[p3]
                e4 = mapping[p4]
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
                e1 = mapping[p1]
                e2 = mapping[p2]
                e3 = mapping[p3]
                e4 = mapping[p4]
                pyridine_torsion = pyridine_labels["ImproperTorsions"][(e2, e1, e3, e4)]
            # check against toluene parameters
            assert periodicity in toluene_torsion.periodicity
            term_index = toluene_torsion.periodicity.index(periodicity)
            assert (
                phase
                == toluene_torsion.phase[term_index].m_as(offunit.radian)
                * ommunit.radian
            )
            assert (
                k
                == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                * ommunit.kilojoule_per_mole
                * improper_scale
            )
            # check against pyridine parameters
            assert periodicity in pyridine_torsion.periodicity
            term_index = pyridine_torsion.periodicity.index(periodicity)
            assert (
                phase
                == pyridine_torsion.phase[term_index].m_as(offunit.radian)
                * ommunit.radian
            )
            assert (
                k
                == pyridine_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                * ommunit.kilojoule_per_mole
                * improper_scale
            )

        # this could also be an anchor torsion which is kept and should use the toluene only values
        elif len(htf._atom_classes["unique_old_atoms"].intersection(torsion)) > 1:
            # we don't check improper torsions as non should be kept
            toluene_torsion = toluene_labels["ProperTorsions"][torsion]
            assert periodicity in toluene_torsion.periodicity
            term_index = toluene_torsion.periodicity.index(periodicity)
            assert (
                phase
                == toluene_torsion.phase[term_index].m_as(offunit.radian)
                * ommunit.radian
            )
            assert (
                k
                == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                * ommunit.kilojoule_per_mole
            )

    # check the interpolated terms - this force has the scaled torsions and impropers
    # Note some impropers are incorrectly scaled here due to the degenerate ordering mentioned above
    # this should have no effect on the energy as the same value is used for both end states
    custom_torsion_force = forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"
    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 35

    # track the removed torsions and impropers
    removed_torsions = set()
    removed_improper_torsions = set()
    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        # p1, p2, p3, p4 are the index in toluene/pyridine get the expected parameters from the labels
        if htf._atom_classes["unique_old_atoms"].intersection({p1, p2, p3, p4}):
            # this is a torsion involving at least one old dummy atom from toluene and should be removed
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
            assert (
                params[2]
                == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                * improper_scale
            )
            # lambda_1 periodicity stays the same
            assert params[3] in toluene_torsion.periodicity
            # lambda_1 phase stays the same
            assert params[4] == toluene_torsion.phase[term_index].m_as(offunit.radian)
            assert params[5] == 0.0

            if improper_scale == 1.0:
                removed_torsions.add(frozenset({p1, p2, p3, p4}))
            else:
                removed_improper_torsions.add(frozenset({p1, p2, p3, p4}))

        else:
            # this is a fully mapped torsion which is changing between toluene and pyridine
            # however the HTF only allows scaling to or from zero potentials not between potentials so the check
            # is complicated
            if params[0] == 0.0 and params[1] == 0.0 and params[2] == 0.0:
                # this is a pyridine torsion which is zeroed at lambda_0
                # map the torsion
                e1 = mapping[p1]
                e2 = mapping[p2]
                e3 = mapping[p3]
                e4 = mapping[p4]
                if (e1, e2, e3, e4) in pyridine_labels["ProperTorsions"]:
                    pyridine_torsion = pyridine_labels["ProperTorsions"][
                        (e1, e2, e3, e4)
                    ]
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
                assert params[4] == pyridine_torsion.phase[term_index].m_as(
                    offunit.radian
                )
                # lambda_1 k
                assert (
                    params[5]
                    == pyridine_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                    * improper_scale
                )
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
                assert params[1] == toluene_torsion.phase[term_index].m_as(
                    offunit.radian
                )
                # lambda_0 k
                assert (
                    params[2]
                    == toluene_torsion.k[term_index].m_as(offunit.kilojoule_per_mole)
                    * improper_scale
                )

    assert removed_torsions == corrections["lambda_1"].removed_dihedrals
    assert removed_improper_torsions == corrections["lambda_1"].removed_impropers


def test_propane_dimethyl_ether_angle_corrections_htf(
    propane_to_dimethyl_ether_mapping,
):
    """Make sure that angle corrections are applied correctly for the dual anchor double branch case."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    corrections = _derive_dummy_junction_corrections(
        propane_to_dimethyl_ether_mapping,
        settings.forcefield_settings.small_molecule_forcefield,
    )
    htf = make_htf(propane_to_dimethyl_ether_mapping, settings, corrections=corrections)
    propane = propane_to_dimethyl_ether_mapping.componentA.to_openff()
    dimethyl_ether = propane_to_dimethyl_ether_mapping.componentB.to_openff()
    ff = ForceField(settings.forcefield_settings.small_molecule_forcefield + ".offxml")
    propane_labels = ff.label_molecules(propane.to_topology())[0]
    dimethyl_ether_labels = ff.label_molecules(dimethyl_ether.to_topology())[0]
    mapping = propane_to_dimethyl_ether_mapping.componentA_to_componentB

    corrected_system = htf.hybrid_system
    # extract the forces
    forces = {force.getName(): force for force in corrected_system.getForces()}

    # angles not interpolated should be in the standard force
    standard_angle_force = forces["HarmonicAngleForce"]
    # there should only be 3 angles not interpolated in the simulation:
    # 1 anchor angle for each dummy (2 total)
    # 1 dummy spanning angle
    num_angles = standard_angle_force.getNumAngles()
    assert num_angles == 3

    for i in range(num_angles):
        p1, p2, p3, angle_eq, k = standard_angle_force.getAngleParameters(i)
        angle = (p1, p2, p3)
        # make sure there are dummy atoms in all angles
        assert set(angle).intersection(htf._atom_classes["unique_old_atoms"])

        # now compare the parameters to the propane labels
        propane_angle = propane_labels["Angles"][angle]
        # angle_eq
        assert angle_eq == propane_angle.angle.m_as(offunit.radian) * ommunit.radian
        # k
        assert (
            k
            == propane_angle.k.m_as(offunit.kilojoule_per_mole / offunit.radian**2)
            * ommunit.kilojoule_per_mole
            / ommunit.radian**2
        )

    # there should be 15 interpolated angle terms covering fully mapped and scaled angles
    custom_angle_force = forces["CustomAngleForce"]
    num_angles = custom_angle_force.getNumAngles()
    assert num_angles == 15
    # there should be a single global parameter for lambda
    assert custom_angle_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_angle_force.getGlobalParameterName(0) == "lambda_angles"

    # track angles turned off
    removed_angles = set()
    for i in range(num_angles):
        p1, p2, p3, params = custom_angle_force.getAngleParameters(i)
        if (
            p1 in htf._atom_classes["unique_old_atoms"]
            or p3 in htf._atom_classes["unique_old_atoms"]
        ):
            # this angle involves at least one old dummy atom from propane
            propane_angle = propane_labels["Angles"][(p1, p2, p3)]
            # lambda_0 angle
            assert params[0] == propane_angle.angle.m_as(offunit.radian)
            # lambda_0 k
            assert params[1] == propane_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )
            # lambda_1 angle stays the same
            assert params[2] == propane_angle.angle.m_as(offunit.radian)
            # lambda_1 k is removed
            assert params[3] == 0 * offunit.kilojoule_per_mole / offunit.radian**2
            removed_angles.add(frozenset({p1, p2, p3}))
        # there are no dummy atoms in dimethyl-ether so all others must be fully mapped
        else:
            # fully mapped angle
            propane_angle = propane_labels["Angles"][(p1, p2, p3)]
            e1 = mapping[p1]
            e2 = mapping[p2]
            e3 = mapping[p3]
            dimethyl_ether_angle = dimethyl_ether_labels["Angles"][(e1, e2, e3)]
            # lambda_0 angle should be the propane angle
            assert params[0] == propane_angle.angle.m_as(offunit.radian)
            # lambda_0 k should be the propane k
            assert params[1] == propane_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )
            # lambda_1 angle should be the dimethyl_ether angle
            assert params[2] == dimethyl_ether_angle.angle.m_as(offunit.radian)
            # lambda_1 k should be the dimethyl_ether k
            assert params[3] == dimethyl_ether_angle.k.m_as(
                offunit.kilojoule_per_mole / offunit.radian**2
            )

    # make sure that the corrections extracted from the htf match what we expected
    assert removed_angles == corrections["lambda_1"].removed_angles


def test_propane_dimethyl_ether_torsion_corrections_htf(
    propane_to_dimethyl_ether_mapping,
):
    """Make sure that torsion corrections are applied correctly for the dual anchor double branch case with a stiffening
    for terminal free rotor anchor groups."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    corrections = _derive_dummy_junction_corrections(
        propane_to_dimethyl_ether_mapping,
        settings.forcefield_settings.small_molecule_forcefield,
    )
    htf = make_htf(propane_to_dimethyl_ether_mapping, settings, corrections=corrections)
    propane = propane_to_dimethyl_ether_mapping.componentA.to_openff()
    dimethyl_ether = propane_to_dimethyl_ether_mapping.componentB.to_openff()
    ff = ForceField(settings.forcefield_settings.small_molecule_forcefield + ".offxml")
    propane_labels = ff.label_molecules(propane.to_topology())[0]
    dimethyl_ether_labels = ff.label_molecules(dimethyl_ether.to_topology())[0]
    mapping = propane_to_dimethyl_ether_mapping.componentA_to_componentB

    corrected_system = htf.hybrid_system
    # extract the forces
    forces = {force.getName(): force for force in corrected_system.getForces()}

    # check the torsions
    standard_torsion_force = forces["PeriodicTorsionForce"]
    num_standard_torsions = standard_torsion_force.getNumTorsions()
    # no torsions should be preserved in this case
    assert num_standard_torsions == 0

    # check the interpolated terms - this force has the scaled torsions
    custom_torsion_force = forces["CustomTorsionForce"]
    # there should be a single global parameter for lambda
    assert custom_torsion_force.getNumGlobalParameters() == 1
    # make sure it has the correct name
    assert custom_torsion_force.getGlobalParameterName(0) == "lambda_torsions"
    num_torsions = custom_torsion_force.getNumTorsions()
    assert num_torsions == 26

    # track the removed torsions
    removed_torsions = set()
    stiffened_torsions = set()

    stiffened_k = htf._stiffened_k_correction.value_in_unit(ommunit.kilojoule_per_mole)

    for i in range(num_torsions):
        p1, p2, p3, p4, params = custom_torsion_force.getTorsionParameters(i)
        if htf._atom_classes["unique_old_atoms"].intersection({p1, p2, p3, p4}):
            # this is a torsion involving at least one old dummy atom from propane could be removed
            # or stiffened
            propane_torsion = propane_labels["ProperTorsions"][(p1, p2, p3, p4)]

            # check the interpolated k value to work out what the torsion is doing
            new_k = params[5]
            if new_k == 0.0:
                removed_torsions.add(frozenset({p1, p2, p3, p4}))
                # we can do standard checks
                # lambda_0 periodicity
                assert params[0] in propane_torsion.periodicity
                term_index = propane_torsion.periodicity.index(params[0])
                # lambda_0 phase
                assert params[1] == propane_torsion.phase[term_index].m_as(
                    offunit.radian
                )
                # lambda_0 k
                assert params[2] == propane_torsion.k[term_index].m_as(
                    offunit.kilojoule_per_mole
                )
                # lambda_1 periodicity stays the same
                assert params[3] in propane_torsion.periodicity
                # lambda_1 phase stays the same
                assert params[4] == propane_torsion.phase[term_index].m_as(
                    offunit.radian
                )

            elif new_k == stiffened_k:
                stiffened_torsions.add(frozenset({p1, p2, p3, p4}))
                # the torsion should use the internal stiffening values
                assert params[0] == params[3] == 1.0  # stiffened periodicity is 1
                assert params[1] == params[4] == 0.0  # stiffened phase
                assert params[2] == 0.0  # initial k

        else:
            # this is a fully mapped torsion which is changing between propane and dimethyl-ether
            # however the HTF only allows scaling to or from zero potentials not between potentials so the check
            # is complicated
            if params[0] == 0.0 and params[1] == 0.0 and params[2] == 0.0:
                # this is a dimethyl-ether torsion which is zeroed at lambda_0
                # map the torsion
                e1 = mapping[p1]
                e2 = mapping[p2]
                e3 = mapping[p3]
                e4 = mapping[p4]
                dimethyl_ether_torsion = dimethyl_ether_labels["ProperTorsions"][
                    (e1, e2, e3, e4)
                ]
                # lambda_1 periodicity
                assert params[3] in dimethyl_ether_torsion.periodicity
                term_index = dimethyl_ether_torsion.periodicity.index(params[3])
                # lambda_1 phase
                assert params[4] == dimethyl_ether_torsion.phase[term_index].m_as(
                    offunit.radian
                )
                # lambda_1 k
                assert (
                    params[5]
                    == dimethyl_ether_torsion.k[term_index].m_as(
                        offunit.kilojoule_per_mole
                    )
                    * 1
                    / dimethyl_ether_torsion.idivf[
                        term_index
                    ]  # this term has a non 1 idivf
                )
            elif params[3] == 0.0 and params[4] == 0.0 and params[5] == 0.0:
                # this is a propane torsion which is zeroed at lambda_1
                propane_torsion = propane_labels["ProperTorsions"][(p1, p2, p3, p4)]
                # lambda_0 periodicity
                assert params[0] in propane_torsion.periodicity
                term_index = propane_torsion.periodicity.index(params[0])
                # lambda_0 phase
                assert params[1] == propane_torsion.phase[term_index].m_as(
                    offunit.radian
                )
                # lambda_0 k
                assert params[2] == propane_torsion.k[term_index].m_as(
                    offunit.kilojoule_per_mole
                )
    # check stiffened first
    assert stiffened_torsions == corrections["lambda_1"].stiffened_dihedrals
    # due to the stiffened torsions also being removed first due to the way the HTF is set up we need to add them
    # to the ref before we check
    assert (
        removed_torsions
        == corrections["lambda_1"].removed_dihedrals
        | corrections["lambda_1"].stiffened_dihedrals
    )


def test_toluene_pyridine_corrections_htf_solvent(toluene_to_pyridine_mapping):
    """Make sure the same angle corrections are applied in a solvated system, in this case there should be no difference."""
    settings = RelativeHybridTopologyProtocol.default_settings()
    corrections = _derive_dummy_junction_corrections(
        toluene_to_pyridine_mapping,
        settings.forcefield_settings.small_molecule_forcefield,
    )

    # make a water solvent to add to the system
    solvent = SolventComponent()
    # limit the number of waters added
    settings.solvation_settings.number_of_solvent_molecules = 100
    settings.solvation_settings.solvent_padding = None
    htf = make_htf(
        toluene_to_pyridine_mapping, settings, solvent=solvent, corrections=corrections
    )

    hybrid_topology = htf.omm_hybrid_topology
    # check how many waters were added
    num_waters = sum((1 for r in hybrid_topology.residues() if r.name == "HOH"))
    assert num_waters == 100

    # the corrections should not change in this case
    assert corrections == htf._valence_corrections_terms
