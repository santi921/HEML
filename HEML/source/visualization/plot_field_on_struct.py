import argparse, moly
import numpy as np

from HEML.utils.dictionaries import *
from HEML.utils.visualization import (
    get_nodes_and_edges_from_pdb,
    mat_to_cones, shift_and_rotate,
    get_AC, 
    xyz2AC_vdW,
    connectivity_to_list_of_bonds, 
    connectivity_filter
)

from HEML.utils.data import (
    pdb_to_xyz,
    filter_other_by_distance,
    filter_xyz_by_distance,
    filter_by_residue,
    mat_pull,
    get_fe_positions,
    get_N_positions,
    get_options
)


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

    if filter_dict["distance_filter"] is not False:
        print("distance filter - {} ang".format(filter_dict["distance_filter"]))
        xyz_filter, residues = filter_xyz_by_distance(
            xyz, 
            center=center,
            residues=residues, 
            distance=filter_dict["distance_filter"], 
            ret_residues=True
        )
        filtered_xyz_dist.extend(xyz_filter)
        filtered_atom_dist.extend(filter_other_by_distance(
            xyz, atom, center=center, distance=filter_dict["distance_filter"]
        ))

        if filter_dict["residue_filter"] is not False: 
            xyz = filtered_xyz_dist
            atom = filtered_atom_dist
    
    
    if filter_dict["residue_filter"] is not False:
        #print(residues)
        for res in filter_dict["residue_filter"]: 
            filtered_xyz_temp, filtered_atom_temp = filter_by_residue(xyz, atom, residues, res)
            #print(filtered_atom_temp)
            if filtered_atom_temp != []:
                filtered_atom.extend(filtered_atom_temp)
                filtered_xyz.extend(filtered_xyz_temp) 
        assert filtered_atom != [], "No atoms found for residue filter {}".format(filter_dict["residue_filter"])
    else: 
        filtered_xyz = xyz
        filtered_atom = atom

    if filter_dict["distance_filter"] is False and filter_dict["residue_filter"] is False:
        filtered_xyz_final, charge, filtered_atom_final = pdb_to_xyz(file)
    else: 
        filtered_atom_final, filtered_xyz_final = [], []
        for i in range(len(filtered_xyz)):
            if not np.any(np.all(filtered_xyz[i] == filtered_xyz_final)):
                filtered_xyz_final.append(filtered_xyz[i])
                filtered_atom_final.append(filtered_atom[i])
    
    connectivity_mat = get_AC(filtered_atom_final, filtered_xyz_final, covalent_factor=1.3)
    bonds = connectivity_to_list_of_bonds(connectivity_mat)
    
    if filter_dict["filter_connectivity"] : 
        filtered_atom_final, filtered_xyz_final, filtered_connectivity_mat, track_indices = connectivity_filter(filtered_atom, filtered_xyz, connectivity_mat, track=26)
        bonds = connectivity_to_list_of_bonds(filtered_connectivity_mat)
    
    return filtered_atom_final, bonds, filtered_xyz_final


def get_molecule_dict(
        file = "../../../data/pdbs_processed/1a4e.pdb",
        alignment_dict=None, 
        filter_dict=None):
    
    alignment_method = alignment_dict["alignment_method"]

    if alignment_method == "heme": 
        print("getting heme alignment from pdb: ".format(file))
        fe_dict = get_fe_positions(file)
        n_dict = get_N_positions(
            file, 
            fe_ID=fe_dict["id"], 
            fe_xyz=fe_dict["xyz"]
        )

        Fe_pos = fe_dict["xyz"]
        NA_pos = n_dict["N1_xyz"]
        NB_pos = n_dict["N2_xyz"]
        NC_pos = n_dict["N3_xyz"]
        ND_pos = n_dict["N4_xyz"]
    
            
        center = np.mean([NA_pos, NB_pos, NC_pos, ND_pos], axis = 0)
        x_axis = np.array(NA_pos) - np.array(center)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.array(NB_pos) - np.array(center)
        y_axis = y_axis / np.linalg.norm(y_axis)
        z_axis = np.cross(y_axis, x_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)

    elif alignment_method == "dict": 
        center = alignment_dict["center"]
        x_axis = alignment_dict["x_axis"]
        y_axis = alignment_dict["y_axis"]
        
        # noramlize axis 
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = y_axis / np.linalg.norm(y_axis)
        #z_axis = z_axis / np.linalg.norm(z_axis)
        z_axis = np.cross(x_axis, y_axis)
        z_axis = z_axis / np.linalg.norm(z_axis)
    
    else:
        if fe_dict == None:
            Fe_pos = [130.581,  41.541,  38.350]
        
        if n_dict==None:        
            NA_pos = [129.775,  39.761,  38.051]
            NB_pos = [130.581,  41.865,  36.409]
            NC_pos = [131.320,  43.348,  38.639]
            ND_pos = [130.469,  41.267,  40.273]

        center = np.mean([NA_pos, NB_pos, NC_pos, ND_pos], axis = 0)


    atom_list, bond_list, xyz_list = get_nodes_and_edges_from_pdb(
        file=file, 
        filter_dict=filter_dict,
        center=center)

    xyz_list = shift_and_rotate(
        xyz_list, 
        center = center, 
        x_axis = x_axis,
        y_axis = y_axis,
        z_axis = z_axis
    )
    #diag_dict = {"Fe": [], "N": []}
    print("number of atoms: ".format(len(atom_list)))
    string_element = "\n{}\n\n".format(len(atom_list))
    for i, atom in enumerate(atom_list):
        string_element +="{} {} {} {}\n".format(
            int_atom_dict[atom], xyz_list[i][0], xyz_list[i][1], xyz_list[i][2])
    
    dict_input = {"symbols": atom_list, "geometry": xyz_list, "connectivity": bond_list}
    return string_element, dict_input


def plot_field(
    field_file, 
    cone_options, 
    alignment_options,
    filter_options,
    molecule_file="../../../data/pdbs_processed/1a4e.pdb",
    dimensions=(21, 21, 21),
    show=False, 
    save=True,
    x_axis_range=[-4, 4],
    y_axis_range=[-4, 4],
    ):

    
    string_element, dict_input = get_molecule_dict(
        molecule_file, 
        alignment_options,
        filter_options, 
    )
    
    component = mat_to_cones(
        mat_pull(field_file), 
        (1, dimensions[0], dimensions[1], dimensions[2], 3),
        bohr_to_ang_conv=True, 
        vector_scale=cone_options["vector_scale"],
        std_mean=cone_options["std_mean"], 
        log1=cone_options["log1"], 
        unlog1=cone_options["unlog1"],
        min_max=cone_options["min_max"],
        cutoff=cone_options["cutoff"],)

    fig = moly.Figure()
    molecule = moly.Molecule.from_data(string_element, dtype="string")
    fig.add_molecule("heme", molecule)
    fig.add_trace(component)
    # add cube to show the center of the molecule from -3 to 3
    fig.fig.update_layout(yaxis_range=y_axis_range, xaxis_range=x_axis_range)

    config = {
    'toImageButtonOptions': {
        'format': 'svg', # one of png, svg, jpeg, webp
        'filename': 'custom_image',
        'height': 1000,
        'width': 1000,
        'scale':6 # Multiply title/legend/axis/canvas sizes by this factor
    }
    }
    fig.fig.update_scenes(
        xaxis_visible=False, 
        yaxis_visible=False,
        zaxis_visible=False)
    
    fig.fig.update(
        layout_showlegend=False)
    
    camera = dict(
        up=dict(x=0, y=0, z=0.7),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0.25, y=0.25, z=0.3)
    )

    fig.fig.update_layout(scene_camera=camera, dragmode='orbit')
    
    if show: 
        fig.fig.show(config=config)
    if save: 
        fig.fig.write_html("./{}.html".format(field_file.split("/")[-1].split(".")[0]), config=config)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-options", help="location of options file", default="./plot_options.json"
    )

    parser.add_argument(
        "-pdb_file", help="location of options file", default="./options/options.json"
    )

    parser.add_argument(
        "-field_file", help="number of magnitudes to calculate", required=True
    )

    options_loc = str(parser.parse_args().options)
    field_file = str(parser.parse_args().field_file)
    pdb_file = str(parser.parse_args().pdb_file)    
    options = get_options(options_loc)
    
    if "filter_dict" not in options.keys():
        options["filter_dict"] = {}
    if "cone_dict" not in options.keys():
        options["cone_dict"] = {}
    if "alignment_dict" not in options.keys():
        options["alignment_dict"] = {}
    
    # defaults for cone dict
    if "cutoff" not in options["cone_dict"].keys():
        options["cone_dict"]["cutoff"] = 98
    if "vector_scale" not in options["cone_dict"].keys():
        options["cone_dict"]["vector_scale"] = 5
    if "std_mean" not in options["cone_dict"].keys():
        options["cone_dict"]["std_mean"] = False
    if "log1" not in options["cone_dict"].keys():
        options["cone_dict"]["log1"] = False
    if "unlog1" not in options["cone_dict"].keys():
        options["cone_dict"]["unlog1"] = False
    if "min_max" not in options["cone_dict"].keys():

        options["cone_dict"]["min_max"] = False
    
    # alignment dict
    if "alignment_method" not in options["alignment_dict"].keys():
        print("setting alignment method to default")
        options["alignment_dict"]["alignment_method"] = "heme"
    if "center" not in options["alignment_dict"].keys():
        print("setting center to default")
        options["alignment_dict"]["center"] = [0, 0, 0]
    if "x_axis" not in options["alignment_dict"].keys():
        print("setting x axis to default")
        options["alignment_dict"]["x_axis"] = [1, 0, 0]
    if "y_axis" not in options["alignment_dict"].keys():
        print("setting y axis to default")
        options["alignment_dict"]["y_axis"] = [0, 1, 0]
    
    # filter_dict 
    if "filter_connectivity" not in options["filter_dict"].keys():
        print("setting filter connectivity to default")
        options["filter_dict"]["filter_connectivity"] = True
    if "distance_to_filter" not in options["filter_dict"].keys():
        print("setting distance to filter to default")
        options["filter_dict"]["distance_to_filter"] = 5.0
    if "residue_filter" not in options["filter_dict"].keys():
        print("setting residue filter to default")
        options["filter_dict"]["filter_residues"] = None
    if "distance_filter" not in options["filter_dict"].keys():
        print("setting distance filter to default")
        options["filter_dict"]["distance_filter"] = True

    # other defaults
    if "x_axis_range" not in options.keys():
        print("setting x axis range to default")
        options["x_axis_range"] = [-4, 4]
    if "y_axis_range" not in options.keys():
        print("setting y axis range to default")
        options["y_axis_range"] = [-4, 4]
    if "field_dims" not in options.keys():
        print("setting field dims to default")
        options["field_dims"] = (21, 21, 21)
    

    plot_field(
        field_file,
        molecule_file=pdb_file,
        dimensions=options["field_dims"],
        show=options["show"],
        save=options["save"],
        filter_options=options["filter_dict"],
        cone_options=options["cone_dict"],
        alignment_options=options["alignment_dict"],
    )

main()