import argparse, moly
import numpy as np

from HEML.utils.dictionaries import *
from HEML.utils.visualization import mat_to_cones, check_viz_dict, get_molecule_dict
from HEML.utils.data import mat_pull, get_options


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
    save_name="field.html",
):
    string_element, dict_input = get_molecule_dict(
        file=molecule_file,
        alignment_dict=alignment_options,
        filter_dict=filter_options,
    )

    mat = mat_pull(field_file, verbose=True)
    meta = mat_pull(field_file, meta_data=True)
    component = mat_to_cones(
        mat,
        (1, dimensions[0], dimensions[1], dimensions[2], 3),
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
        sparse_factor=cone_options["sparse_factor"],
    )

    fig = moly.Figure(figsize=(800, 800))
    molecule = moly.Molecule.from_data(string_element, dtype="string")
    fig.add_molecule("molecule", molecule, style="tubes")
    fig.add_trace(component)
    fig.fig.update_layout(yaxis_range=y_axis_range, xaxis_range=x_axis_range)

    config = {
        "toImageButtonOptions": {
            "format": "svg",  # one of png, svg, jpeg, webp
            "filename": "custom_image",
            "height": 1000,
            "width": 1000,
            "scale": 6,  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.fig.update_scenes(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False)

    fig.fig.update(layout_showlegend=False)

    camera = dict(
        up=dict(x=0, y=0, z=0.7),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=0.25, y=0.25, z=0.3),
    )

    fig.fig.update_layout(scene_camera=camera, dragmode="orbit")

    if show:
        fig.fig.show(config=config)
    if save:
        fig.fig.write_html("{}".format(save_name), config=config)


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
            - sparsify: whether to sparsify the field
            - sparse_factor: factor to sparsify the field by
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
        "-pdb_file", help="location of options file", default="./options/options.json"
    )

    parser.add_argument(
        "-field_file", help="number of magnitudes to calculate", required=True
    )

    print("overlays field on molecular structure")
    print("field file: {}".format(parser.parse_args().field_file))
    print("pdb file: {}".format(parser.parse_args().pdb_file))
    options_loc = str(parser.parse_args().options)
    field_file = str(parser.parse_args().field_file)
    pdb_file = str(parser.parse_args().pdb_file)
    options = get_options(options_loc)
    options = check_viz_dict(options)

    plot_field(
        field_file,
        molecule_file=pdb_file,
        dimensions=options["field_dims"],
        show=options["show"],
        save=options["save"],
        save_name=options["save_name"],
        x_axis_range=options["x_axis_range"],
        y_axis_range=options["y_axis_range"],
        filter_options=options["filter_dict"],
        cone_options=options["cone_dict"],
        alignment_options=options["alignment_dict"],
    )


main()
