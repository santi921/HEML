import plotly.graph_objects as go
import moviepy.editor as mpy
import io
from PIL import Image
import numpy as np
import moly
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from glob import glob
from tqdm import tqdm
from HEML.utils.data import mat_pull
from HEML.utils.visualization import (
    mat_to_cones,
    get_nodes_and_edges_from_pdb,
    shift_and_rotate,
)
from HEML.utils.dictionaries import *


def fetch_heme_as_string(molecule_file="../../../data/pdbs_processed/1a4e.pdb"):
    # THIS IS HARD CODED FOR HEMES CHANGE IF YOU WANT TO USE THIS FOR SOMETHING ELSE
    atom_list, bond_list, xyz_list = get_nodes_and_edges_from_pdb(
        molecule_file, distance_filter=5.5
    )

    NA_pos = [129.775, 39.761, 38.051]
    NB_pos = [130.581, 41.865, 36.409]
    NC_pos = [131.320, 43.348, 38.639]
    ND_pos = [130.469, 41.267, 40.273]
    Fe_pos = [130.581, 41.541, 38.350]
    center = np.mean([NA_pos, NB_pos, NC_pos, ND_pos], axis=0)
    x_axis = np.array(NA_pos) - np.array(center)
    x_axis = x_axis / np.linalg.norm(x_axis)
    y_axis = np.array(NB_pos) - np.array(center)
    y_axis = y_axis / np.linalg.norm(y_axis)
    z_axis = np.cross(y_axis, x_axis)
    z_axis = z_axis / np.linalg.norm(z_axis)
    xyz_list = shift_and_rotate(
        xyz_list, center=center, x_axis=x_axis, y_axis=y_axis, z_axis=z_axis
    )
    diag_dict = {"Fe": [], "N": []}
    string_element = "\n{}\n\n".format(len(atom_list))
    for i, atom in enumerate(atom_list):
        # check if nitrogen
        if atom == 7:
            diag_dict["N"].append(xyz_list[i])
        if atom == 26:
            diag_dict["Fe"].append(xyz_list[i])
        string_element += "{} {} {} {}\n".format(
            int_atom_dict[atom], xyz_list[i][0], xyz_list[i][1], xyz_list[i][2]
        )

    return string_element


def plotly_fig2array(fig):
    # convert Plotly fig to  an array
    fig_bytes = fig.to_image(format="png", scale=2)
    buf = io.BytesIO(fig_bytes)
    img = Image.open(buf)
    return np.asarray(img)


def animate_fields(molecule_file="../../../data/pdbs_processed/1a4e.pdb"):
    folder_plot = "../../../data/1u5u_amber_md/"

    frames_in = []
    # get all the files in the folder ending in .dat
    files = glob(folder_plot + "/*.dat")
    shape = mat_pull(files[0]).shape
    shape = (1, shape[0], shape[1], shape[2], shape[3])
    meta_data = mat_pull(files[0], meta_data=True)
    step_size = {
        "x": meta_data["step_size_x"],
        "y": meta_data["step_size_y"],
        "z": meta_data["step_size_z"],
    }

    bounds = {
        "x": meta_data["bounds_x"],
        "y": meta_data["bounds_y"],
        "z": meta_data["bounds_z"],
    }
    scalar = 1.8897259885789
    meta_data_bohr = {
        "step_size_x": meta_data["step_size_x"] * scalar,
        "step_size_y": meta_data["step_size_y"] * scalar,
        "step_size_z": meta_data["step_size_z"] * scalar,
        "bounds_x": [i * scalar for i in meta_data["bounds_x"]],
        "bounds_y": [i * scalar for i in meta_data["bounds_y"]],
        "bounds_z": [i * scalar for i in meta_data["bounds_z"]],
    }

    hemes_string = fetch_heme_as_string(molecule_file)
    fig_mol = moly.Figure()
    molecule = moly.Molecule.from_data(hemes_string, dtype="string")
    # fig_mol.add_molecule("heme", molecule)
    fig = fig_mol.fig

    # Layout
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                nticks=10,
                range=[1 * i for i in meta_data_bohr["bounds_x"]],
            ),
            yaxis=dict(
                nticks=10,
                range=[1 * i for i in meta_data_bohr["bounds_y"]],
            ),
            zaxis=dict(
                nticks=10,
                range=[1 * i for i in meta_data_bohr["bounds_z"]],
            ),
        ),
        title="Slices in volumetric data",
        width=2000,
        height=2000,
    )
    """    fig.update_layout(
            scene_camera=dict(
                up=dict(x=0, y=0, z=0.7),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=10.5, y=10.5, z=10.3)
            )
        )"""

    def make_frame(t):
        cones = mat_to_cones(
            mat_pull(files[int(fps * t)]),
            shape,
            vector_scale=2,
            cutoff=70,
            step_size=step_size,
            bounds=bounds,
            bohr_to_ang_conv=True,
        )

        fig.data = []
        fig_mol.add_molecule("heme", molecule)
        fig.add_traces(cones)
        fig.update_layout(
            scene_camera=dict(
                up=dict(x=0, y=0, z=0.7),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=0.5, y=0.5, z=0.5),
            )
        )
        return plotly_fig2array(fig)

    fps = 2
    animation = mpy.VideoClip(make_frame, duration=len(files[::10]) / fps)
    animation.write_gif("test.gif", fps=fps)


animate_fields()
