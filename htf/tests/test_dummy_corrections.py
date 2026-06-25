import pathlib
from dataclasses import fields

from htf.utils import _derive_dummy_junction_corrections, make_htf, _draw_dummy_corrections
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol


def test_single_terminal_corrections(ejm_50_to_ejm_42_mapping, tmp_path):
    """Test single dummy group terminal corrections for a single and triple branch type."""
    corrections = _derive_dummy_junction_corrections(ejm_50_to_ejm_42_mapping, "openff-2.0.0.offxml")
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
        frozenset({32, 19, 30, 31})
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
        frozenset({32, 19, 29, 31}), # dummy group 1 redundant dihedrals
        frozenset({32, 19, 30, 31}),
        frozenset({33, 19, 29, 31}), # dummy group 2 redundant dihedrals
        frozenset({33, 19, 30, 31}),
        frozenset({34, 19, 29, 31}), # dummy group 3 redundant dihedrals
        frozenset({34, 19, 30, 31})
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_50_to_ejm_42_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_50_to_ejm_42_corrections"
    )


def test_tyk2_dummy_group_interaction_corrections(ejm_50_to_ejm_55_mapping, tmp_path):
    """Test corrections which should leave interactions between two dummy groups not on the same core junction atom (ejm_50)."""
    corrections = _derive_dummy_junction_corrections(ejm_50_to_ejm_55_mapping, "openff-2.0.0.offxml")
    # check the simple non-corrected end state first ejm_55
    # we need no corrections for this triple terminal junction due to the geometry of the junction and lack of dummy groups
    state_0 = corrections["lambda_0"]
    # make sure all corrections are blank
    assert not state_0.softened_angles
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    assert not state_0.removed_impropers
    assert not state_0.removed_dihedrals
    assert not state_0.removed_angles

    # check the more complicated interacting dummy group case ejm_50
    # this has a terminal and dual anchor double branch dummy group which are connected
    state_1 = corrections["lambda_1"]
    # first check the blank corrections
    assert not state_1.stiffened_dihedrals
    assert not state_1.stiffened_angles
    assert not state_1.softened_angles
    assert not state_1.removed_impropers
    # make sure 2 redundant angles are removed around the dual anchor junction involving dummy atoms {29, 30}
    assert state_1.removed_angles == {
        frozenset({19, 29, 31}), # dummy atom 1 on dual anchor junction
        frozenset({19, 30, 31}) # dummy atom 2 on dual anchor junction
    }
    # make sure 2 redundant dihedrals are removed which terminate in the dummy atoms {29, 30}
    assert state_1.removed_dihedrals == {
        frozenset({16, 17, 19, 29}), # redundant dihedral on dummy atom 1
        frozenset({16, 17, 19, 30}) # redundant dihedral on dummy atom 2
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_50_to_ejm_55_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=pathlib.Path("ejm_50_to_ejm_55_corrections")
    )


def test_tyk2_dual_junction_large_single_branch_corrections(ejm_31_to_ejm_49_mapping, tmp_path):
    """Test the dual junction single branch corrections for a large dummy group which has an improper in the junction (ejm_49).
    Here we want to make sure the improper spanning the core and the dummy junction atom is removed while the improper
    spanning the dummy group and the core junction atom is kept.
    """
    corrections = _derive_dummy_junction_corrections(ejm_31_to_ejm_49_mapping, "openff-2.0.0.offxml")
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
    assert state_1.removed_angles == {
        frozenset({17, 18, 19})
    }
    # make sure the 4 redundant dihedrals are removed
    assert state_1.removed_dihedrals == {
        frozenset({16, 17, 19, 28}), # dummy junction redundant dihedral
        frozenset({17, 18, 19, 29}), # dummy group redundant dihedrals
        frozenset({17, 18, 19, 30}),
        frozenset({17, 18, 19, 31})
    }

    # check the dual anchor single branch large group ejm_49
    state_0 = corrections["lambda_0"]
    # check the blank corrections first
    assert not state_0.softened_angles
    assert not state_0.stiffened_angles
    assert not state_0.stiffened_dihedrals
    # make sure the redundant angle is removed linking the dummy junction atom 31
    assert state_0.removed_angles == {
        frozenset({17, 18, 31})
    }
    # make sure only 1 improper is removed which spanns the core and dummy junction
    assert state_0.removed_impropers == {
        frozenset({16, 17, 18, 31})
    }
    # make sure the 3 redundant dihedrals are removed
    assert state_0.removed_dihedrals == {
        frozenset({16, 17, 27, 31}), # dummy junction redundant dihedral
        frozenset({17, 18, 30, 31}), # dummy group redundant dihedrals
        frozenset({32, 17, 18, 31})
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_31_to_ejm_49_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_31_to_ejm_49_corrections"
    )


def test_find_toluene_to_pyridine_corrections(toluene_to_pyridine_mapping, tmp_path):
    """Simple test case with corrections at one end (toluene), dual anchor single branch correction type."""
    corrections = _derive_dummy_junction_corrections(toluene_to_pyridine_mapping, "openff-2.0.0.offxml")
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
        frozenset({0, 1, 2, 3}), # dummy anchor constraint
        frozenset({0, 1, 2, 10}), # dummy anchor constraint
        frozenset({0, 1, 6, 14}), # dummy anchor constraint
        frozenset({0, 1, 2, 8}), # dummy group dual anchor constraint
        frozenset({0, 1, 2, 9}), # dummy group dual anchor constraint
        frozenset({0, 1, 2, 7}), # dummy group dual anchor constraint
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=toluene_to_pyridine_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "toluene_to_pyridine_corrections"
    )


def test_propane_to_dimethyl_ether_corrections(propane_to_dimethyl_ether_mapping, tmp_path):
    """Free rotor test case with corrections at one end (dimethyl ether), dual anchor two branch correction type."""
    corrections = _derive_dummy_junction_corrections(propane_to_dimethyl_ether_mapping, "openff-2.0.0.offxml")
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
        frozenset({0, 1, 9}), # extra dummy angle constraint
        frozenset({0, 1, 10}), # extra other dummy angle constraint
    }
    # check we have 10 removed in total 5 for each group
    assert state_1.removed_dihedrals == {
        frozenset({0, 9, 5, 7}), # dummy group 1 redundant dihedrals
        frozenset({0, 1, 2, 9}),
        frozenset({0, 1, 3, 9}),
        frozenset({0, 9, 5, 6}),
        frozenset({0, 1, 4, 9}),
        frozenset({0, 10, 5, 6}), # dummy group 2 redundant dihedrals
        frozenset({0, 1, 2, 10}),
        frozenset({0, 1, 10, 4}),
        frozenset({0, 10, 5, 7}),
        frozenset({0, 1, 10, 3})
    }
    # check we have two stiffened dihedrals as the anchor is a free rotor
    # they should also terminate in the same atom
    assert state_1.stiffened_dihedrals == {
        # note both dihedrals terminate in atom 8
        frozenset({8, 0, 5, 9}),  # dummy group 1
        frozenset({8, 0, 10, 5}) # dummy group 2
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=propane_to_dimethyl_ether_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "propane_to_dimethyl_ether_corrections"
    )


def test_triple_corrections_tyk2(ejm_31_to_ejm_42_mapping, tmp_path):
    """Triple non-planar junction correction type at both ends, end stateB (lambda_0 corrections) is a dummy group with pruned group anchor dihedral constraints."""
    corrections = _derive_dummy_junction_corrections(ejm_31_to_ejm_42_mapping, "openff-2.0.0.offxml")
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
        frozenset({17, 19, 30})
    }
    # check that the 2 possible dihedrals coupling the dummy are removed
    assert state_1.removed_dihedrals == {
        # both dihedrals terminate in 30
        frozenset({17, 18, 19, 30}),
        frozenset({16, 17, 19, 30})
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
        frozenset({19, 29, 31})
    }
    # check the removed dihedrals
    assert state_0.removed_dihedrals == {
        frozenset({16, 17, 19, 31}), # dummy junction redundant dihedrals
        frozenset({17, 18, 19, 31}),
        frozenset({32, 19, 30, 31}), # dummy group redundant dihedrals anchor 1
        frozenset({33, 19, 30, 31}),
        frozenset({34, 19, 30, 31}),
        frozenset({32, 19, 29, 31}), # dummy group redundant dihedrals anchor 2
        frozenset({33, 19, 29, 31}),
        frozenset({34, 19, 29, 31}),
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=ejm_31_to_ejm_42_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "ejm_31_to_ejm_42_corrections"
    )


def test_terminal_co_linear_tyk2_corrections(jmc_28_to_jmc_30_mapping, tmp_path):
    """Triple terminal junction to a terminal co-linear junction (jmc_30)"""
    corrections = _derive_dummy_junction_corrections(jmc_28_to_jmc_30_mapping, "openff-2.0.0.offxml")
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
        frozenset({35, 34, 19, 21})
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
        frozenset({35, 19, 36, 21}), # dummy atom 36 redundant dihedrals
        frozenset({34, 35, 36, 21}),
        frozenset({37, 35, 19, 21}), # dummy atom 37 redundant dihedrals
        frozenset({37, 34, 35, 21}),
        frozenset({35, 19, 21, 38}), # dummy atom 38 redundant dihedrals
        frozenset({34, 35, 21, 38})
    }
    # draw the corrections
    _draw_dummy_corrections(
        mapping=jmc_28_to_jmc_30_mapping,
        corrections=corrections,
        force_field="openff-2.0.0.offxml",
        output_dir=tmp_path / "jmc_28_to_jmc_30_corrections"
    )