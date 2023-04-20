import moly
import numpy as np
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


def animate_fields(folder_plot, molecule_file="../../../data/pdbs_processed/1a4e.pdb"):
    hemes_string = fetch_heme_as_string(molecule_file)
    fig_mol = moly.Figure()
    molecule = moly.Molecule.from_data(hemes_string, dtype="string")
    fig_mol.add_molecule("heme", molecule)

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

    print("Gathering frames....")
    for i in tqdm(range(int(len(files) / 10))):
        if i == 0:
            frame_init = mat_to_cones(
                mat_pull(files[i]),
                shape,
                vector_scale=1,
                step_size=step_size,
                bounds=bounds,
                bohr_to_ang_conv=True,
            )
        frames_in.append(
            mat_to_cones(
                mat_pull(files[i]),
                shape,
                vector_scale=1,
                step_size=step_size,
                bounds=bounds,
                bohr_to_ang_conv=True,
            )
        )

    fig_mol.fig.add_trace(frame_init)
    fig = fig_mol.fig
    # go.Figure(data=frame_init)
    fig.frames = [
        go.Frame(data=[frames_in[i]], name=str(i)) for i in range(len(frames_in))
    ]

    sliders = [
        {
            "pad": {"b": 10, "t": 60},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [[f.name], frame_args(0)],
                    "label": str(k),
                    "method": "animate",
                }
                for k, f in enumerate(fig.frames)
            ],
        }
    ]

    # Layout
    fig.update_layout(
        scene=dict(
            xaxis=dict(
                nticks=10,
                range=meta_data_bohr["bounds_x"],
            ),
            yaxis=dict(
                nticks=10,
                range=meta_data_bohr["bounds_y"],
            ),
            zaxis=dict(
                nticks=10,
                range=meta_data_bohr["bounds_z"],
            ),
        ),
        title="Slices in volumetric data",
        width=1000,
        height=1000,
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None, frame_args(50)],
                        "label": "&#9654;",  # play symbol
                        "method": "animate",
                    },
                    {
                        "args": [[None], frame_args(0)],
                        "label": "&#9724;",  # pause symbol
                        "method": "animate",
                    },
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 70},
                "type": "buttons",
                "x": 0.1,
                "y": 0,
            }
        ],
        scene_camera=dict(
            up=dict(x=0, y=0, z=0.7),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=0.5, y=0.5, z=0.3),
        ),
        sliders=sliders,
    )
    config = {
        "toImageButtonOptions": {
            "format": "svg",  # one of png, svg, jpeg, webp
            "filename": "custom_image",
            "height": 2000,
            "width": 2000,
            "scale": 6,  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }
    fig.show(config=config)


def frame_args(duration):
    return {
        "frame": {"duration": duration},
        "mode": "immediate",
        "fromcurrent": True,
        "transition": {"duration": duration, "easing": "linear"},
    }
