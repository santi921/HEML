import numpy as np
from rdkit import Chem

from HEML.utils.data import pdb_to_xyz, filter_other_by_distance, filter_xyz_by_distance
from HEML.utils.xyz2mol import xyz2AC_vdW
from HEML.utils.dictionaries import * 

def shift_and_rotate(
    xyz_list, center=[0, 0, 0], x_axis=[1, 0, 0], y_axis=[0, 1, 0], z_axis=[0, 0, 1]
):
    for i in range(len(xyz_list)):
        xyz_list[i] = xyz_list[i] - center
        xyz_list[i] = np.array(
            [
                np.dot(xyz_list[i], x_axis),
                np.dot(xyz_list[i], y_axis),
                np.dot(xyz_list[i], z_axis),
            ]
        )
    return xyz_list


def ms(x, y, z, radius, resolution=50):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0 : 2 * np.pi : resolution * 2j, 0 : np.pi : resolution * 1j]
    X = radius * np.cos(u) * np.sin(v) + x
    Y = radius * np.sin(u) * np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)


def get_AC(atoms, xyz, covalent_factor=1.4):
    pt = Chem.GetPeriodicTable()
    # xyz to distance matrix
    xyz = np.array(xyz)
    dist_mat = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    num_atoms = xyz.shape[0]

    AC = np.zeros((num_atoms, num_atoms), dtype=int)
    for i in range(num_atoms):
        a_i = atoms[i]
        Rcov_i = pt.GetRcovalent(a_i) * covalent_factor
        for j in range(i + 1, num_atoms):
            a_j = atoms[j]
            Rcov_j = pt.GetRcovalent(a_j) * covalent_factor
            if dist_mat[i, j] <= Rcov_i + Rcov_j:
                AC[i, j] = 1
    return AC


def connectivity_to_list_of_bonds(connectivity_mat):
    bonds = []
    for i in range(len(connectivity_mat)):
        for j in range(i + 1, len(connectivity_mat)):
            if connectivity_mat[i][j] > 0:
                bonds.append([i, j])
    return bonds


def get_nodes_and_edges_from_pdb(
    file="../../data/pdbs_processed/1a4e.pdb", distance_filter=5.0
):

    xyz, charge, atom = pdb_to_xyz(file)
    filtered_xyz = filter_xyz_by_distance(
        xyz, center=[130.581, 41.541, 38.350], distance=distance_filter
    )
    # filtered_charge = filter_other_by_distance(xyz, charge, center = [130.581,  41.541,  38.350], distance = distance_filter)
    filtered_atom = filter_other_by_distance(
        xyz, atom, center=[130.581, 41.541, 38.350], distance=distance_filter
    )
    connectivity_mat, rdkit_mol = xyz2AC_vdW(filtered_atom, filtered_xyz)
    connectivity_mat = get_AC(filtered_atom, filtered_xyz, covalent_factor=1.3)

    bonds = connectivity_to_list_of_bonds(connectivity_mat)
    return filtered_atom, bonds, filtered_xyz

