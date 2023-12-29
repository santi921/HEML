from HEML.utils.data import (
    fetch_charges_dict,
    pull_mats_w_label,
    mat_pull,
)


def test_mat_pull():
    """Test the mat_pull function"""
    test_mat = "./test_data/test.dat"
    mat = mat_pull(test_mat)
    # print(mat.shape)
    assert mat.shape == (21, 21, 21, 3), "mat_pull failed"
    meta_data = mat_pull(test_mat, meta_data=True)
    check_keys = [
        "first_line",
        "steps_x",
        "steps_y",
        "steps_z",
        "step_size_x",
        "step_size_y",
        "step_size_z",
        "bounds_x",
        "bounds_y",
        "bounds_z",
    ]
    for key in check_keys:
        assert key in meta_data.keys(), "mat_pull failed"


def test_dataset():
    x, y, meta = pull_mats_w_label(
        dir_fields="./test_data/dataset/",
        data_file="./test_data/protein_data.csv",
        meta_data=True,
    )
    assert x.shape == (187, 21, 21, 21, 3), "pull_mats_w_label failed"
    assert y.shape == (187, 3), "pull_mats_w_label failed"
    meta_keys = [
        "first_line",
        "steps_x",
        "steps_y",
        "steps_z",
        "step_size_x",
        "step_size_y",
        "step_size_z",
        "bounds_x",
        "bounds_y",
        "bounds_z",
    ]
    for key in meta_keys:
        assert key in meta.keys(), "pull_mats_w_label failed"


def test_charge_pull():
    charge_list = fetch_charges_dict("./test_data/test.pqr")
    len(charge_list) == 29460, "fetch_charges_dict failed. Wrong number of charges"
    keys = ["position", "charge", "radius"]
    for charge_dict in charge_list:
        for key, value in charge_dict.items():
            assert key in keys, "fetch_charges_dict failed. Wrong keys"
            assert (
                type(value) == float or type(value) == list
            ), "fetch_charges_dict failed. Wrong values"
    probe_list = [charge_list[0], charge_list[100], charge_list[-1]]
    assert probe_list[0]["position"] == [37.122, -33.316, 31.089]
    assert probe_list[0]["charge"] == -0.133
    assert probe_list[0]["radius"] == 1.55

    assert probe_list[1]["position"] == [36.751, -19.603, 24.233]
    assert probe_list[1]["charge"] == 0.096
    assert probe_list[1]["radius"] == 1.2

    assert probe_list[2]["position"] == [-8.58, -9.655, -15.042]
    assert probe_list[2]["charge"] == 0.102
    assert probe_list[2]["radius"] == 1.2
