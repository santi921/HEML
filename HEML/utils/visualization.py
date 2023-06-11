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
from HEML.utils.data import get_fe_positions, get_N_positions


def check_viz_dict(options):
    """
    sets defaults for viz dict if not present
    """

    if "filter_dict" not in options.keys():
        options["filter_dict"] = {}
    if "cone_dict" not in options.keys():
        options["cone_dict"] = {}
    if "alignment_dict" not in options.keys():
        options["alignment_dict"] = {}

    # defaults for cone dict
    if "cutoff" not in options["cone_dict"].keys():
        print("setting cutoff to default: ", 98)
        options["cone_dict"]["cutoff"] = 98
    if "vector_scale" not in options["cone_dict"].keys():
        print("setting vector_scale to default: ", 5)
        options["cone_dict"]["vector_scale"] = 5
    if "std_mean" not in options["cone_dict"].keys():
        print("setting std_mean to default: ", "False")
        options["cone_dict"]["std_mean"] = False
    if "log1" not in options["cone_dict"].keys():
        print("setting log1 to default: ", "False")
        options["cone_dict"]["log1"] = False
    if "unlog1" not in options["cone_dict"].keys():
        print("setting unlog1 to default: ", "False")
        options["cone_dict"]["unlog1"] = False
    if "min_max" not in options["cone_dict"].keys():
        print("setting min_max to default: ", "False")
        options["cone_dict"]["min_max"] = False
    if "opacity" not in options["cone_dict"].keys():
        print("setting opacity to default: ", 0.5)
        options["cone_dict"]["opacity"] = 0.5

    # alignment dict
    if "alignment_method" not in options["alignment_dict"].keys():
        print("setting alignment method to default: ", "heme")
        options["alignment_dict"]["alignment_method"] = "heme"
    if "center" not in options["alignment_dict"].keys():
        print("setting center to default: ", [0, 0, 0])
        options["alignment_dict"]["center"] = [0, 0, 0]
    if "x_axis" not in options["alignment_dict"].keys():
        print("setting x axis to default: ", [1, 0, 0])
        options["alignment_dict"]["x_axis"] = [1, 0, 0]
    if "y_axis" not in options["alignment_dict"].keys():
        print("setting y axis to default: ", [0, 1, 0])
        options["alignment_dict"]["y_axis"] = [0, 1, 0]

    # filter_dict
    if "filter_connectivity" not in options["filter_dict"].keys():
        print("setting filter connectivity to default: ", True)
        options["filter_dict"]["filter_connectivity"] = True
    if "residue_filter" not in options["filter_dict"].keys():
        print("setting residue filter to default: ", False)
        options["filter_dict"]["residue_filter"] = False
    if "distance_filter" not in options["filter_dict"].keys():
        print("setting distance filter to default: ", 5.0)
        options["filter_dict"]["distance_filter"] = 5.0
    # other defaults
    if "x_axis_range" not in options.keys():
        print("setting x axis range to default: ", [-4, 4])
        options["x_axis_range"] = [-4, 4]
    if "y_axis_range" not in options.keys():
        print("setting y axis range to default: ", [-4, 4])
        options["y_axis_range"] = [-4, 4]
    if "field_dims" not in options.keys():
        print("setting field dims to default: ", (21, 21, 21))
        options["field_dims"] = (21, 21, 21)

    if "save_name" not in options.keys():
        print("setting save name to default: ", "./test")
        options["save_name"] = "./test"
    return options


def shift_and_rotate(
    xyz_list, center=[0, 0, 0], x_axis=[1, 0, 0], y_axis=[0, 1, 0], z_axis=[0, 0, 1]
):
    x_axis = np.array(x_axis) / np.linalg.norm(x_axis)
    y_axis = np.array(y_axis) / np.linalg.norm(y_axis)
    z_axis = np.array(z_axis) / np.linalg.norm(z_axis)
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
    Filter out atoms that are not connected to the track atom. This might need to be retooled for non-heme molecules.
    """
    track_index = [i for i in range(len(filtered_atom)) if filtered_atom[i] == track]
    bonds = connectivity_to_list_of_bonds(connectivity_mat)
    graph = nx.from_edgelist(bonds)

    if not nx.is_connected(graph):
        largest_cc = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]

        has_track = [
            i for i in range(len(largest_cc)) if track_index[0] in largest_cc[i].nodes
        ]
        graph_to_use = largest_cc[has_track[0]]
        # print(graph_to_use.nodes)
        filtered_atoms = [filtered_atom[i] for i in graph_to_use.nodes]
        filtered_xyz = [filtered_xyz[i] for i in graph_to_use.nodes]
        filtered_connectivity_mat = np.zeros((len(filtered_atoms), len(filtered_atoms)))
        for i in range(len(filtered_atoms)):
            for j in range(len(filtered_atoms)):
                filtered_connectivity_mat[i][j] = connectivity_mat[filtered_atoms[i]][
                    filtered_atoms[j]
                ]
        track_indices = [
            i for i in range(len(filtered_atoms)) if filtered_atoms[i] == track
        ]
    return filtered_atoms, filtered_xyz, filtered_connectivity_mat, track_indices


def get_nodes_and_edges_from_pdb(
    filter_dict,
    file="../../data/pdbs_processed/1a4e.pdb",
    center=[130.581, 41.541, 38.350],
):
    """
    Obtains and processes the nodes and edges from a pdb file.
    Takes:
        file: pdb file to be processed
        distance_filter: distance to filter by
        filter_connectivity: whether to filter by connectivity to central molecule
        center: center of molecule to plot
    """
    filtered_atom, filtered_xyz = [], []
    filtered_xyz_dist, filtered_atom_dist = [], []
    xyz, charge, atom, residues = pdb_to_xyz(file, ret_residues=True)
    print("-" * 25 + "Filtering Module Start" + "-" * 25)
    print("number of atoms before filtering: \t", len(atom))
    if filter_dict["distance_filter"] is not False:
        print("distance filter - {} ang".format(filter_dict["distance_filter"]))
        xyz_filter, residues = filter_xyz_by_distance(
            xyz,
            center=center,
            residues=residues,
            distance=filter_dict["distance_filter"],
            ret_residues=True,
        )
        atoms_filter = filter_other_by_distance(
            xyz, atom, center=center, distance=filter_dict["distance_filter"]
        )
        print("# of atoms left after dist filter: \t", len(xyz_filter))
        filtered_xyz_dist.extend(xyz_filter)
        filtered_atom_dist.extend(atoms_filter)

        print("done filtering by distance!")
        xyz = filtered_xyz_dist
        atom = filtered_atom_dist

    if filter_dict["residue_filter"] is not False:
        print("filtering by residue")
        for res in filter_dict["residue_filter"]:
            filtered_xyz_temp, filtered_atom_temp = filter_by_residue(
                xyz, atom, residues, res
            )

            if filtered_atom_temp != []:
                # filtered_atom.extend(filtered_atom_temp)
                # filtered_xyz.extend(filtered_xyz_temp)
                filtered_atom = filtered_atom_temp
                filtered_xyz = filtered_xyz_temp

        assert filtered_atom != [], "No atoms found for residue filter {}".format(
            filter_dict["residue_filter"]
        )
        print("done filtering by residue!")
        print("number of atoms left after res filter: \t", len(filtered_atom))

    else:
        filtered_xyz = xyz
        filtered_atom = atom

    if (
        filter_dict["distance_filter"] is False
        and filter_dict["residue_filter"] is False
    ):
        filtered_xyz_final, charge, filtered_atom_final = pdb_to_xyz(
            file, ret_residues=False
        )

    else:
        filtered_atom_final, filtered_xyz_final = [], []
        for i in range(len(filtered_xyz)):
            if not np.any(np.all(filtered_xyz[i] == filtered_xyz_final)):
                filtered_xyz_final.append(filtered_xyz[i])
                filtered_atom_final.append(filtered_atom[i])

    print("getting connectivity matrix")
    print("number of atoms left: \t\t", len(filtered_atom_final))
    connectivity_mat = get_AC(
        filtered_atom_final, filtered_xyz_final, covalent_factor=1.3
    )
    bonds = connectivity_to_list_of_bonds(connectivity_mat)

    if filter_dict["filter_connectivity"]:
        print("filtering by connectivity")
        (
            filtered_atom_final,
            filtered_xyz_final,
            filtered_connectivity_mat,
            track_indices,
        ) = connectivity_filter(filtered_atom, filtered_xyz, connectivity_mat, track=26)
        bonds = connectivity_to_list_of_bonds(filtered_connectivity_mat)
        print("done filtering by connectivity!")
        print("number of atoms left after conn filter: \t\t", len(filtered_atom_final))
    print("-" * 23 + "Filtering Module Complete" + "-" * 23)
    return filtered_atom_final, bonds, filtered_xyz_final


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
    (
        arr_min,
        arr_max,
    ) = np.min(
        x
    ), np.max(x)
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


def get_molecule_dict(
    alignment_dict, file="../../../data/pdbs_processed/1a4e.pdb", filter_dict=None
):  
    #print(alignment_dict)
    alignment_method = alignment_dict["alignment_method"]

    if alignment_method == "heme":
        print("getting heme alignment from pdb: ".format(file))
        fe_dict = get_fe_positions(file)
        n_dict = get_N_positions(file, fe_ID=fe_dict["id"], fe_xyz=fe_dict["xyz"])

        Fe_pos = fe_dict["xyz"]
        NA_pos = n_dict["N1_xyz"]
        NB_pos = n_dict["N2_xyz"]
        NC_pos = n_dict["N3_xyz"]
        ND_pos = n_dict["N4_xyz"]

        center = np.mean([NA_pos, NB_pos, NC_pos, ND_pos], axis=0)
        x_axis = np.array(NA_pos) - np.array(center)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.array(NB_pos) - np.array(center)
        y_axis = y_axis / np.linalg.norm(y_axis)
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)

    elif alignment_method == "dict":
        print("using dict alignment provided")
        center = alignment_dict["center"]
        x_axis = alignment_dict["x_axis"]
        y_axis = alignment_dict["y_axis"]

        # noramlize axis
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = y_axis / np.linalg.norm(y_axis)
        # z_axis = z_axis / np.linalg.norm(z_axis)
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)

    else:
        if fe_dict == None:
            Fe_pos = [130.581, 41.541, 38.350]

        if n_dict == None:
            NA_pos = [129.775, 39.761, 38.051]
            NB_pos = [130.581, 41.865, 36.409]
            NC_pos = [131.320, 43.348, 38.639]
            ND_pos = [130.469, 41.267, 40.273]

        center = np.mean([NA_pos, NB_pos, NC_pos, ND_pos], axis=0)

    atom_list, bond_list, xyz_list = get_nodes_and_edges_from_pdb(
        file=file, filter_dict=filter_dict, center=center
    )

    ind_fe = [i for i in range(len(atom_list)) if atom_list[i] == 26]
    #print(ind_fe)
    #print(atom_list[ind_fe[0]])
    #print(xyz_list[ind_fe[0]])
    #print("center at helper: " + str(center))
    #print(x_axis, y_axis, z_axis)
    xyz_list = shift_and_rotate(
        xyz_list, center=center, x_axis=x_axis, y_axis=y_axis, z_axis=z_axis
    )
    print("number of atoms:  {}".format(len(atom_list)))
    string_element = "\n{}\n\n".format(len(atom_list))
    for i, atom in enumerate(atom_list):
        string_element += "{} {} {} {}\n".format(
            int_atom_dict[atom], xyz_list[i][0], xyz_list[i][1], xyz_list[i][2]
        )
    dict_input = {"symbols": atom_list, "geometry": xyz_list, "connectivity": bond_list}
    # get index of iron 
    #get element that equals 26 in atom_list
    ind_fe = [i for i in range(len(atom_list)) if atom_list[i] == 26]
    print(ind_fe)
    print(atom_list[ind_fe[0]])
    print(xyz_list[ind_fe[0]])
    return string_element, dict_input


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
    opacity=0.4,
):
    bohr_to_ang = 1
    if bohr_to_ang_conv:
        bohr_to_ang = 1.88973

    comp_vect_field = mat.reshape(shape[1], shape[2], shape[3], shape[4])
    x, y, z = np.meshgrid(
        np.arange(
            bounds["x"][0] * bohr_to_ang,
            (bounds["x"][1]) * bohr_to_ang,
            step_size["x"] * bohr_to_ang,
        ),
        np.arange(
            bounds["y"][0] * bohr_to_ang,
            (bounds["y"][1]) * bohr_to_ang,
            step_size["y"] * bohr_to_ang,
        ),
        np.arange(
            bounds["z"][0] * bohr_to_ang,
            (bounds["z"][1]) * bohr_to_ang,
            step_size["z"] * bohr_to_ang,
        ),
    )
    #print(x.shape)
    #print(y.shape)
    #print(z.shape)
    #print("opacity: {}".format(opacity))
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
        opacity=opacity,
        showscale=False,
        colorscale="Picnic",
    )
