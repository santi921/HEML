import moly
from HEML.utils.visualization import get_molecule_dict, mat_to_cones
from HEML.utils.data import mat_pull


def test_atom_filter():
    file = "./test_data/test.pdb"
    keys_check = ["symbols", "geometry", "connectivity"]
    align_dict = {
        "alignment_method": "heme",
        "center": [130.581, 41.541, 38.350],
        "x_axis": [1, 0, 0],
        "y_axis": [0, 1, 0],
    }

    filter_dict = {
        "distance_filter": 10.0,
        "residue_filter": False,
        "filter_connectivity": False,
    }

    string, filter_distance_dict = get_molecule_dict(
        file=file, alignment_dict=align_dict, filter_dict=filter_dict
    )

    assert all(key in filter_distance_dict.keys() for key in keys_check)
    assert len(filter_distance_dict["symbols"]) == 443

    filter_dict["distance_filter"] = 6.0
    string, filter_distance_dict = get_molecule_dict(
        file=file, alignment_dict=align_dict, filter_dict=filter_dict
    )

    assert all(key in filter_distance_dict.keys() for key in keys_check)
    assert len(filter_distance_dict["symbols"]) == 98

    filter_dict["residue_filter"] = [
        "TYR",
        "AZI",
        "HEM",
        "HIS",
    ]
    string, filter_distance_dict = get_molecule_dict(
        file=file, alignment_dict=align_dict, filter_dict=filter_dict
    )

    assert all(key in filter_distance_dict.keys() for key in keys_check)
    assert len(filter_distance_dict["symbols"]) == 68

    filter_dict["filter_connectivity"] = True
    string, filter_distance_dict = get_molecule_dict(
        file=file, alignment_dict=align_dict, filter_dict=filter_dict
    )
    assert all(key in filter_distance_dict.keys() for key in keys_check)
    assert len(filter_distance_dict["symbols"]) == 61


def test_moly():
    file = "./test_data/test.pdb"
    test_mat = "./test_data/test.dat"
    bound = 5
    step = 0.4

    cones = mat_to_cones(
        mat=mat_pull(test_mat),
        shape=[1, 21, 21, 21, 3],
        cutoff=99,
        vector_scale=5,
        bounds={"x": [-bound, bound], "y": [-bound, bound], "z": [-bound, bound]},
        step_size={"x": step, "y": step, "z": step},
        bohr_to_ang_conv=True,
    )

    align_dict = {
        "alignment_method": "heme",
        "center": [130.581, 41.541, 38.350],
        "x_axis": [1, 0, 0],
        "y_axis": [0, 1, 0],
    }

    filter_dict = {
        "distance_filter": 6.0,
        "residue_filter": [
            "TYR",
            "AZI",
            "HEM",
            "HIS",
        ],
        "filter_connectivity": True,
    }

    string, filter_distance_dict = get_molecule_dict(
        file=file, alignment_dict=align_dict, filter_dict=filter_dict
    )

    fig = moly.Figure()
    molecule = moly.Molecule.from_data(string, dtype="string")
    fig.add_molecule("heme", molecule)
    fig.fig.update_layout(yaxis_range=[-5, 5], xaxis_range=[-5, 5])
    fig.add_trace(cones)
