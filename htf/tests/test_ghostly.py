from htf.utils import _parse_ghostly_output


def test_parse_ghostly_output(ghostly_output_chloroethane_to_ethane):
    corrections = _parse_ghostly_output(ghostly_output_chloroethane_to_ethane)
    expected_corrections = {
        "lambda_0": {
            "bridges": {
                0: {
                    "bridge_atom": 1,
                    "ghosts": [8],
                    "physical": [0, 2, 3, 4],
                    "type": 4,
                }
            },
            "removed_angles": [(2, 1, 8)],
            "removed_dihedrals": [(5, 2, 1, 8), (7, 2, 1, 8), (6, 2, 1, 8)],
            "stiffened_angles": [(4, 1, 8), (3, 1, 8)],
            "softened_angles": {},
        },
        "lambda_1": {
            "bridges": {
                0: {
                    "bridge_atom": 1,
                    "ghosts": [0],
                    "physical": [2, 3, 4, 8],
                    "type": 4,
                }
            },
            "removed_angles": [(0, 1, 2), (0, 1, 3)],
            "removed_dihedrals": [(0, 1, 2, 5), (0, 1, 2, 7), (0, 1, 2, 6)],
            "stiffened_angles": [(0, 1, 4)],
            "softened_angles": {},
        },
    }
    assert corrections == expected_corrections