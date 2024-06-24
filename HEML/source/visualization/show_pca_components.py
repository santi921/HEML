from plotly.subplots import make_subplots
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go


import argparse, moly
import numpy as np


# from HEML.utils.data import *
from HEML.utils.data import pull_mats_from_MD_folder
from HEML.utils.fields import pca, unwrap_pca
from HEML.utils.visualization import mat_to_cones, get_molecule_dict, check_viz_dict
from HEML.utils.data import get_options


def plot_field(
    cone_options,
    alignment_options,
    filter_options,
    x, 
    pca_obj,
    meta,
    molecule_file="../../../data/pdbs_processed/3wxo.pdb",
    dimensions=(21, 21, 21),
    show=False,
    save=True,
    x_axis_range=[-4, 4],
    y_axis_range=[-4, 4],
    pca_comp_to_show=0,
    file_save="pca_component.html",
):
    print(alignment_options)
    string_element, dict_input = get_molecule_dict(
        file=molecule_file,
        alignment_dict=alignment_options,
        filter_dict=filter_options,
    )

    #x_pca, _ = pca(
    #    x, verbose=True, pca_comps=int(pca_comp_to_show + 1), whitening=True,
    #    pca=pca_obj
    #)
    shape_mat = x.shape
    print("dimensions: ", dimensions)
    print("shape_mat: ", shape_mat)
    pca_comp = pca_obj.components_[pca_comp_to_show].reshape(
        1, shape_mat[1], shape_mat[2], shape_mat[3], shape_mat[4]
    )

    component = mat_to_cones(
        pca_comp,
        (1, shape_mat[1], shape_mat[2], shape_mat[3], shape_mat[4]),
        bohr_to_ang_conv=True,
        bounds={"x": meta["bounds_x"], "y": meta["bounds_y"], "z": meta["bounds_z"]},
        step_size={
            "x": meta["step_size_x"],
            "y": meta["step_size_y"],
            "z": meta["step_size_z"],
        },
        vector_scale=cone_options["vector_scale"],
        std_mean=cone_options["std_mean"],
        log1=cone_options["log1"],
        unlog1=cone_options["unlog1"],
        min_max=cone_options["min_max"],
        cutoff=cone_options["cutoff"],
        opacity=cone_options["opacity"],
        sparsify=cone_options["sparsify"],
        sparse_factor=cone_options["sparsify_factor"],
    )

    fig = moly.Figure()
    molecule = moly.Molecule.from_data(string_element, dtype="string")
    fig.add_molecule("molecule", molecule, style="tubes")
    fig.add_trace(component)
    fig.fig.update_layout(yaxis_range=y_axis_range, xaxis_range=x_axis_range)

    config = {
        "toImageButtonOptions": {
            "format": "png",  # one of png, svg, jpeg, webp
            "filename": "custom_image",
            "height": 1000,
            "width": 1000,
            "scale": 6,  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)
    fig.fig.update(layout_showlegend=False)

    camera = dict(
        up=dict(x=0, y=0, z=0.2),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0.6, y=0.3, z=0.15),
    )
    fig.fig.update_layout(scene_camera=camera, dragmode="orbit")

    if show:
        fig.fig.show(config=config)

    if save:
        fig.fig.write_html("{}".format(file_save), config=config)


def main():
    """
    Base implementation of overlaying field on molecular structure some of the options are:
    filter_dict: dictionary of atoms to filter out of the pdb file
        - filter_connectivity (bool): whether to filter out atoms that are not connected to the biggest graph
        - filter_distance (float or false): whether to filter out atoms that are not within a certain distance of the center
        - residue_filter(list or false): filter residues that are not in the list
    cone_dict: dictionary of options for the cones
        - vector_scale: scale of the vectors
        - std_mean: whether to use the std or mean for the color
        - log1: whether to log the data
        - unlog1: whether to unlog the data (keep as false, this is for some other function)
        - min_max: whether to use the min max for the color
        - cutoff: cutoff for the color

    alignment_dict: dictionary of options for the alignment of the molecular center
        - center: center of the molecule
        - x_axis: x axis of the molecule
        - y_axis: y axis of the molecule
        - alignment_method: heme (finds center based on heme logic), dict (uses the center and axis from the dict), none (guesses)
    field_dims: dimensions of the field
    x_axis_range: range of the x axis
    y_axis_range: range of the y axis
    show: whether to show the plot
    save: whether to save the plot as html


    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-options", help="location of options file", default="./plot_options.json"
    )

    parser.add_argument(
        "-pdb_file",
        help="location of options file",
        default="../../../data/pdbs_processed/1a4e.pdb",
    )

    parser.add_argument(
        "-data_root", help="root of fields to compute pca on", required=True
    )
    parser.add_argument("-pca_components", help="which pca to show", default=0)

    pca_comps_to_show = int(parser.parse_args().pca_components)
    options_loc = str(parser.parse_args().options)
    pdb_file = str(parser.parse_args().pdb_file)
    data_root = str(parser.parse_args().data_root)
    options = get_options(options_loc)
    options = check_viz_dict(options)

    # change this to get the data from the crystal folder
    x, _, _ = pull_mats_from_MD_folder(
        data_file="../../../data/protein_data.csv",
        root_dir=data_root,
    )

    #x_test, _, _ = pull_mats_from_MD_folder(
    #    data_file="../../../data/protein_data.csv",
    #    root_dir="../../../data/fields_test/",
    #)

    # concat x, x_test
    #x = np.concatenate((x, x_test), axis=0)
    
    meta = {
        "bounds_x": [-4, 4],
        "bounds_y": [-4, 4],
        "bounds_z": [-4, 4],
        "step_size_x": 0.4,
        "step_size_y": 0.4,
        "step_size_z": 0.4,
    }

    arr_min = np.min(x)
    arr_max = np.max(x)

    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)

    mat_transform, pca_obj = pca(x, verbose=True, pca_comps=10, whitening=True)
    

    for i in range(pca_comps_to_show):
        plot_field(
            molecule_file=pdb_file,
            dimensions=options["field_dims"],
            x=x,
            pca_obj=pca_obj,
            meta=meta, 
            show=options["show"],
            save=options["save"],
            x_axis_range=options["x_axis_range"],
            y_axis_range=options["y_axis_range"],
            filter_options=options["filter_dict"],
            cone_options=options["cone_dict"],
            alignment_options=options["alignment_dict"],
            pca_comp_to_show=i,
            file_save=options["save_name"] + str(i) + ".html",
        )


main()
