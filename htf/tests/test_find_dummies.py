from htf.utils import _find_dummy_junctions


def test_find_dummies_single_junction(htf_chloro_ethane):
    htf = htf_chloro_ethane["htf"]
    junctions = _find_dummy_junctions(htf)
    expected_junctions = {
        "lambda_0": {0: {"junction_atom": 1, "dummies": [8], "physical": [0, 2, 3, 4]}},
        "lambda_1": {0: {"junction_atom": 1, "dummies": [0], "physical": [2, 3, 4, 8]}},
    }
    assert junctions == expected_junctions

def test_find_dummies_no_dummy(htf_chloro_fluoroethane):
    htf = htf_chloro_fluoroethane["htf"]
    junctions = _find_dummy_junctions(htf)
    expected_junctions = {
        "lambda_0": {},
        "lambda_1": {},
    }
    assert junctions == expected_junctions

def test_find_dummies_triple_junction(htf_propane_chloroethane):
    """Make sure we can find the dummy-core junctions in a triple junction system."""
    htf = htf_propane_chloroethane["htf"]
    junctions = _find_dummy_junctions(htf)
    expected_junctions = {
        # no dummies in lambda_0
        "lambda_0": {},
        # 3 dummy Hs on the terminal carbon in lambda_1
        "lambda_1": {0: {"junction_atom": 5, "dummies": [8, 6, 7], "physical": [0]}},
    }
    assert junctions == expected_junctions

def test_find_dummies_double_junction(htf_propane_dimethyl_ether):
    """Make sure we can find the dummy-core junctions in a double junction system."""
    htf = htf_propane_dimethyl_ether["htf"]
    junctions = _find_dummy_junctions(htf)
    expected_junctions = {
        # no dummies in lambda_0
        "lambda_0": {},
        # 2 dummy Hs on the central carbon in lambda_1
        "lambda_1": {0: {"junction_atom": 0, "dummies": [9, 10], "physical": [1, 5]}},
    }
    assert junctions == expected_junctions

