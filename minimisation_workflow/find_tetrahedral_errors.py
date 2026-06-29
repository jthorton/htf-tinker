"""
Find tetrahedral errors from the minimisation workflow
"""

from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import pathlib
import click
import json
from itertools import combinations


@click.command()
@click.option("--input-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path))
def main(input_dir: pathlib.Path):
    """
    Find tetrahedral errors from the minimisation workflow
    """
    # Get all the SDF files in the input directory
    sdf_files = list(input_dir.glob("*.sdf"))

    # Loop through each SDF file and check for tetrahedral errors
    tet_angles = []
    for sdf_file in sdf_files:
        suppl = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
        for mol in suppl:
            if mol is None:
                continue
            atom_rings = mol.GetRingInfo().AtomRings()
            for atom in mol.GetAtoms():
                if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and atom.GetDegree() == 4 and not any(atom.GetIdx() in ring for ring in atom_rings):
                    # Get the angles between the bonds
                    angles = []
                    neighbors = list(atom.GetNeighbors())
                    for atom_pair in combinations(neighbors, 2):
                        angle = rdMolTransforms.GetAngleDeg(mol.GetConformer(), atom_pair[0].GetIdx(), atom.GetIdx(), atom_pair[1].GetIdx())
                        angles.append({
                            "angle": angle,
                            "atoms": (atom_pair[0].GetIdx(), atom.GetIdx(), atom_pair[1].GetIdx()),
                            "molecule": str(sdf_file)
                        })
                    tet_angles.extend(angles)

    # Check if any angles are outside the range of 109.5 +/- 10 degrees
    tet_errors = [angle["angle"] for angle in tet_angles if angle["angle"] < 99.5 or angle["angle"] > 119.5]
    if tet_errors:
        click.echo(f"Found {len(tet_errors)} tetrahedral errors in {len(sdf_files)} SDF files.")
        click.echo(f"Angles outside the range of 109.5 +/- 10 degrees: {tet_errors}")
    else:
        click.echo(f"No tetrahedral errors found in {len(sdf_files)} SDF files.")

    # save out the raw angles as a json file so we can make a plot across all results later
    with open(input_dir / "tetrahedral_angles.json", "w") as f:
        json.dump(tet_angles, f, indent=4)


if __name__ == "__main__":
    main()