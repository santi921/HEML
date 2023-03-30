import numpy as np
import networkx as nx
from rdkit import Chem
import plotly.graph_objects as go

from HEML.utils.data import (
    pdb_to_xyz,
    filter_other_by_distance,
    filter_xyz_by_distance,
    filter_by_residue,
    pull_mats_w_label,
)
from HEML.utils.xyz2mol import xyz2AC_vdW
from HEML.utils.fields import split_and_filter, pca
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

def connectivity_filter(filtered_atom, filtered_xyz, connectivity_mat, track=26):
    """
    Filter out atoms that are not connected to the track atom.
    """
    track_index = [i for i in range(len(filtered_atom)) if filtered_atom[i] == track]
    bonds = connectivity_to_list_of_bonds(connectivity_mat)
    graph = nx.from_edgelist(bonds)
    
    if not nx.is_connected(graph):
        largest_cc = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        #print(largest_cc[0].nodes, track_index[0])
        has_track = [i for i in range(len(largest_cc)) if track_index[0] in largest_cc[i].nodes]
        graph_to_use = largest_cc[has_track[0]]
        #print(graph_to_use.nodes)
        filtered_atoms = [filtered_atom[i] for i in graph_to_use.nodes]
        filtered_xyz = [filtered_xyz[i] for i in graph_to_use.nodes]
        filtered_connectivity_mat = np.zeros((len(filtered_atoms), len(filtered_atoms)))
        for i in range(len(filtered_atoms)):
            for j in range(len(filtered_atoms)):
                filtered_connectivity_mat[i][j] = connectivity_mat[filtered_atoms[i]][filtered_atoms[j]]
        track_indices = [i for i in range(len(filtered_atoms)) if filtered_atoms[i] == track]
    return filtered_atoms, filtered_xyz, filtered_connectivity_mat, track_indices

def get_nodes_and_edges_from_pdb(
    file="../../data/pdbs_processed/1a4e.pdb", 
    distance_filter=None, 
    heme_filter=False, 
    filter_connectivity=False
):
    if heme_filter:
        xyz, charge, atom, residues = pdb_to_xyz(file, ret_residues=True)
        
    else:
        xyz, charge, atom = pdb_to_xyz(file)
    
    if distance_filter is not None:
        filtered_xyz = filter_xyz_by_distance(
            xyz, center=[130.581, 41.541, 38.350], distance=distance_filter
        )
        filtered_atom = filter_other_by_distance(
            xyz, atom, center=[130.581, 41.541, 38.350], distance=distance_filter
        )
    
    else: 
        
        filtered_xyz, filtered_atom = filter_by_residue(xyz, atom, residues, "HEM")
        
        if distance_filter is not None:
            filtered_xyz = filter_xyz_by_distance(
                filtered_xyz, center=[130.581, 41.541, 38.350], distance=distance_filter
            )
            filtered_atom = filter_other_by_distance(
                filtered_xyz, filtered_atom, center=[130.581, 41.541, 38.350], distance=distance_filter
            )
        

    connectivity_mat, rdkit_mol = xyz2AC_vdW(filtered_atom, filtered_xyz)
    connectivity_mat = get_AC(filtered_atom, filtered_xyz, covalent_factor=1.3)
    bonds = connectivity_to_list_of_bonds(connectivity_mat)
    
    if filter_connectivity: 
        filtered_atom, filtered_xyz, filtered_connectivity_mat, track_indices = connectivity_filter(filtered_atom, filtered_xyz, connectivity_mat, track=26)
        bonds = connectivity_to_list_of_bonds(filtered_connectivity_mat)
    
    return filtered_atom, bonds, filtered_xyz


def get_cones_viz_from_pca(
    vector_scale=3,
    components=10,
    cutoff=75,
    data_file="../../data/protein_data.csv",
    dir_fields="../../data/cpet/",
    bounds={"x": [-3.0, 3.0], "y": [-3.0, 3.0], "z": [-3.0, 3.0]},
    step_size={"x": 0.3, "y": 0.3, "z": 0.3},
):

    cones = []

    x, _ = pull_mats_w_label(data_file=data_file, dir_fields=dir_fields)
    print("field shape: " + str(x.shape))
    arr_min, arr_max, = np.min(x), np.max(x)
    # getting sign of every element
    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)
    x = (x - arr_min) / np.abs(arr_max - arr_min + 0.1)

    x_untransformed = x
    x_pca, pca_obj = pca(x, verbose=True, pca_comps=components, write=False)
    shape_mat = x.shape

    for ind, pca_comp in enumerate(pca_obj.components_):
        comp_vect_field = pca_comp.reshape(
            shape_mat[1], shape_mat[2], shape_mat[3], shape_mat[4]
        )
        bohr_to_ang = 1.88973
        x, y, z = np.meshgrid(
            np.arange(
                bounds["x"][0] * bohr_to_ang,
                (bounds["x"][1] + step_size["x"]) * bohr_to_ang,
                step_size["x"] * bohr_to_ang,
            ),
            np.arange(
                bounds["y"][0] * bohr_to_ang,
                (bounds["y"][1] + step_size["y"]) * bohr_to_ang,
                step_size["y"] * bohr_to_ang,
            ),
            np.arange(
                bounds["z"][0] * bohr_to_ang,
                (bounds["z"][1] + step_size["z"]) * bohr_to_ang,
                step_size["z"] * bohr_to_ang,
            ),
        )

        u_1, v_1, w_1 = split_and_filter(
            comp_vect_field, cutoff=cutoff, std_mean=True, min_max=False
        )

        cones.append(
            go.Cone(
                x=x.flatten(),
                y=y.flatten(),
                z=z.flatten(),
                u=u_1,
                v=v_1,
                w=w_1,
                sizeref=vector_scale,
                opacity=0.4,
                showscale=False,
                colorscale="Greens",
            )
        )

    return cones


def mat_to_cones(
    mat,
    shape,
    vector_scale=3,
    cutoff=0,
    bounds={"x": [-3.0, 3.0], "y": [-3.0, 3.0], "z": [-3.0, 3.0]},
    step_size={"x": 0.3, "y": 0.3, "z": 0.3},
    bohr_to_ang_conv=False,
    cos_center_scaling=False,
    std_mean=False,
    log1=False,
    unlog1=False,
    min_max=False,
):

    bohr_to_ang = 1
    if bohr_to_ang_conv:
        bohr_to_ang = 1.88973

    comp_vect_field = mat.reshape(shape[1], shape[2], shape[3], shape[4])
    x, y, z = np.meshgrid(
        np.arange(
            bounds["x"][0] * bohr_to_ang,
            (bounds["x"][1] + step_size["x"]) * bohr_to_ang,
            step_size["x"] * bohr_to_ang,
        ),
        np.arange(
            bounds["y"][0] * bohr_to_ang,
            (bounds["y"][1] + step_size["y"]) * bohr_to_ang,
            step_size["y"] * bohr_to_ang,
        ),
        np.arange(
            bounds["z"][0] * bohr_to_ang,
            (bounds["z"][1] + step_size["z"]) * bohr_to_ang,
            step_size["z"] * bohr_to_ang,
        ),
    )
    print(x.shape)
    print(y.shape)
    print(z.shape)
    u_1, v_1, w_1 = split_and_filter(
        comp_vect_field,
        cutoff=cutoff,
        std_mean=std_mean,
        min_max=min_max,
        log1=log1,
        unlog1=unlog1,
        cos_center_scaling=cos_center_scaling,
    )

    return go.Cone(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        u=u_1,
        v=v_1,
        w=w_1,
        sizeref=vector_scale,
        opacity=0.4,
        showscale=False,
        colorscale="Greens",
    )
