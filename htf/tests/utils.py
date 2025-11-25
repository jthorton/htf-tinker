import numpy as np
import openmm


def _make_system_with_cmap(
        map_sizes: list[int],
        mapped_torsions: list[tuple[int, int, int, int, int, int, int, int, int]] | None = None,
        num_atoms: int = 8
):
    """
    Build an OpenMM System with a CMAP term based on the provided mapping data.
    :param map_sizes: Mapping data:
    :param mapped_torsions:
    :param num_atoms:
    :return:
    """
    assert num_atoms >= 8, "num_atoms must be at least 8 to accommodate mapped torsions"
    system = openmm.System()
    # add dummy forces to avoid errors
    for force in [openmm.NonbondedForce, openmm.HarmonicBondForce, openmm.HarmonicAngleForce, openmm.PeriodicTorsionForce]:
        system.addForce(force())

    for _ in range(num_atoms):
        system.addParticle(12.0)  # Add carbon-like particles

    # create a CMAP force
    cmap_force = openmm.CMAPTorsionForce()

    for map_size in map_sizes:
        # Create a grid for the CMAP
        grid = [0.0] * (map_size * map_size)
        cmap_force.addMap(map_size, grid)

    if mapped_torsions is None:
        # add a single cmap term for all atoms using the first map
        mapped_torsions = [(0, 0, 1, 2, 3, 4, 5, 6, 7)]

    for torsion in mapped_torsions:
        cmap_force.addTorsion(torsion[0], *torsion[1:])

    system.addForce(cmap_force)
    # build a basic topology for the number of atoms bonding each atom to the next
    topology = openmm.app.Topology()
    chain = topology.addChain()
    res = topology.addResidue('RES', chain)
    atoms = []
    for i in range(num_atoms):
        atom = topology.addAtom(f'C{i+1}', openmm.app.element.carbon, res)
        atoms.append(atom)
        if i > 0:
            topology.addBond(atoms[i-1], atoms[i])
    # build a fake set of positions
    positions = openmm.unit.Quantity(np.zeros((num_atoms, 3)), openmm.unit.nanometer)
    return system, topology, positions