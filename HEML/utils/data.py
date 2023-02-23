import os, json
from rdkit import Chem
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from sklearn.cluster import AffinityPropagation
from HEML.utils.xyz2mol import xyz2AC_vdW
from copy import deepcopy


atom_int_dict = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Br": 35,
    "Fe": 26,
    "FE": 26,
    "I": 53,
}

atom_mass = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984,
    "P": 30.9738,
    "S": 32.065,
    "Cl": 35.453,
    "Br": 79.904,
    "Fe": 55.845,
    "I": 126.904,
}


int_atom_dict = {
    1: "H",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    15: "P",
    16: "S",
    17: "Cl",
    35: "Br",
    26: "Fe",
    53: "I",
}


atomic_size = {
    "H": 0.5,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "Fe": 1.80,
    "I": 1.98,
}


atom_colors = {
    "H": "white",
    "C": "black",
    "N": "blue",
    "O": "red",
    "F": "orange",
    "P": "green",
    "S": "yellow",
    "Cl": "green",
    "Br": "brown",
    "Fe": "orange",
    "I": "purple",
}


def get_options(options_file="./options.json"):
    """
    Get options from options.json file and create folders if they don't exist.
    Takes
        options_file: path to options.json file
    Returns
        options: dictionary of options
    """
    with open(options_file) as f:
        options = json.load(f)
    for key in options:
        if "folder" in key:
            if not os.path.exists(options[key]):
                os.makedirs(options[key])

    return options


def create_folders(folder_name):
    """
    Creates the folders for the turbomole calculations.
    Takes:
        folder_name: the folder where the folders should be created
    """

    if not os.path.exists("{}/no_charges".format(folder_name)):
        os.makedirs("{}/no_charges".format(folder_name))
    if not os.path.exists("{}/embedding".format(folder_name)):
        os.makedirs("{}/embedding".format(folder_name))

    if not os.path.exists("{}/embedding/o".format(folder_name)):
        os.makedirs("{}/embedding/o".format(folder_name))
    if not os.path.exists("{}/embedding/oh".format(folder_name)):
        os.makedirs("{}/embedding/oh".format(folder_name))
    if not os.path.exists("{}/embedding/normal".format(folder_name)):
        os.makedirs("{}/embedding/normal".format(folder_name))

    if not os.path.exists("{}/no_charges/o".format(folder_name)):
        os.makedirs("{}/no_charges/o".format(folder_name))
    if not os.path.exists("{}/no_charges/oh".format(folder_name)):
        os.makedirs("{}/no_charges/oh".format(folder_name))
    if not os.path.exists("{}/no_charges/normal".format(folder_name)):
        os.makedirs("{}/no_charges/normal".format(folder_name))


def check_if_file_is_empty(file):
    if os.stat(file).st_size == 0:
        return True
    else:
        return False


def check_if_dict_has_None(dict):
    for key, value in dict.items():
        if value is None:
            return True

    if dict is {}:
        return True
    else:
        return False


def break_up_line(str_process):
    split_str = str_process.split("-")
    if split_str[0] == "":
        ret_1 = "-" + split_str[1]
    else:
        ret_1 = split_str[0]
    return ret_1, "-" + split_str[-1]


def spacefinder(List_String):
    try:
        len(List_String) == 11
    except:
        print("Does your pdb contain charges? List of strings should have 11 values.")

    slen1 = len(List_String[1])
    slen2 = len(List_String[2])
    slen5 = len(List_String[5])
    slen6 = len(List_String[6])
    slen7 = len(List_String[7])
    slen8 = len(List_String[8])
    slen9 = len(List_String[9])
    slen10 = len(List_String[10])

    if List_String[0] == "HETATM":
        backlen1 = 4 - (slen1 - 1)
    else:
        backlen1 = 6 - (slen1 - 1)
    if slen2 > 3:
        backlen2 = 2 - (slen2 - 3)
        backlen3 = 3 - (slen2 - 2)
    else:
        backlen2 = 2
        backlen3 = 3 - (slen2 - 1)
    backlen5 = 3 - (slen5 - 1)
    backlen6 = 11 - (slen6 - 1)
    backlen7 = 7 - (slen7 - 1)
    backlen8 = 7 - (slen8 - 1)
    backlen9 = 6 - (slen9 - 1)
    backlen10 = 6 - (slen10 - 1)

    outstring = (
        List_String[0]
        + (" " * backlen1)
        + List_String[1]
        + (" " * backlen2)
        + List_String[2]
        + (" " * backlen3)
        + List_String[3]
        + " "
        + List_String[4]
        + (" " * backlen5)
        + List_String[5]
        + (" " * backlen6)
        + List_String[6]
        + (" " * backlen7)
        + List_String[7]
        + (" " * backlen8)
        + List_String[8]
        + (" " * backlen9)
        + List_String[9]
        + (" " * backlen10)
        + List_String[10]
    )
    return outstring


def pdb_to_xyz(file):
    with open(file, "r") as f:
        lines = f.readlines()
    xyz = []
    charge = []
    atom = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            # xyz.append([float(i) for i in line.split()[6:9]])
            charge.append(float(line.split()[-2]))
            atom.append(atom_int_dict[line.split()[-1]])
    return xyz, charge, atom


def filter_xyz_by_distance(xyz, center=[0, 0, 0], distance=5):
    xyz = np.array(xyz, dtype=float)
    center = np.array(center, dtype=float)
    return xyz[np.linalg.norm(xyz - center, axis=1) < distance]


def filter_other_by_distance(xyz, other, center=[0, 0, 0], distance=5):
    xyz = np.array(xyz, dtype=float)
    center = np.array(center, dtype=float)
    mask = np.linalg.norm(xyz - center, axis=1) < distance
    mask = [i for i in range(len(mask)) if mask[i]]
    return [other[i] for i in mask]


def get_N_positions(file, fe_ID, fe_xyz):
    print(file)
    N_ID, N_ID2, N_ID3, N_ID4 = None, None, None, None
    N1_xyz, N2_xyz, N3_xyz, N4_xyz = None, None, None, None

    with open(file, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()

        shift = 0
        if len(line[0]) > 6:
            shift = -1

        if (
            "HETATM" in line[0]
            and ("NA" in line[2 + shift] and "HEM" or "HEC" in line[3 + shift])
            and line[4 + shift] == fe_ID.split(":")[0]
        ):
            N_ID = f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}"
            try:
                N1_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
            except:
                N1_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
        if (
            "HETATM" in line[0]
            and ("NB" in line[2 + shift] and "HEM" or "HEC" in line[3 + shift])
            and line[4 + shift] == fe_ID.split(":")[0]
        ):
            N_ID2 = f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}"
            try:
                N2_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
            except:
                N2_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
        if (
            "HETATM" in line[0]
            and ("NC" in line[2 + shift] and "HEM" or "HEC" in line[3 + shift])
            and line[4 + shift] == fe_ID.split(":")[0]
        ):
            N_ID3 = f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}"
            try:
                N3_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
            except:
                N3_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
        if (
            "HETATM" in line[0]
            and ("ND" in line[2 + shift] and "HEM" or "HEC" in line[3 + shift])
            and line[4 + shift] == fe_ID.split(":")[0]
        ):
            N_ID4 = f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}"
            try:
                N4_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
            except:
                N4_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]

    # if there are missing N atoms, find all of the nitrogens and assign them to the closest N atom

    full_N_dict = {}
    if N_ID == None or N_ID2 == None or N_ID3 == None or N_ID4 == None:
        count = 0
        for j in readfile:
            line = j.split()

            shift = 0
            if len(line[0]) > 6:
                shift = -1

            heme_id = fe_ID.split(":")[0]

            if (
                "HETATM" in line[0]
                and ("N" in line[2 + shift] and "HEM" or "HEC" in line[3 + shift])
                and line[4 + shift] == heme_id
            ):
                try:
                    full_N_dict[count] = {
                        "id": f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}",
                        "xyz": [float(j[31:38]), float(j[38:45]), float(j[46:54])],
                        "distance_to_iron": np.linalg.norm(
                            np.array(
                                [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                            )
                            - np.array(fe_xyz)
                        ),
                    }
                except:
                    full_N_dict[count] = {
                        "id": f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}",
                        "xyz": [float(line[-5]), float(line[-4]), float(line[-3])],
                        "distance_to_iron": np.linalg.norm(
                            np.array(
                                [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                            )
                            - np.array(fe_xyz)
                        ),
                    }
                count += 1

        # sort the dictionary by distance to iron
        full_N_dict = {
            k: v
            for k, v in sorted(
                full_N_dict.items(), key=lambda item: item[1]["distance_to_iron"]
            )
        }
        # relabel the keys as 0,1,2,3
        full_N_dict = {i: full_N_dict[j] for i, j in enumerate(full_N_dict.keys())}
        # get first 4 values of the dictionary
        full_N_dict = dict(list(full_N_dict.items())[0:4])
        # assign the N_IDs
        N_ID = full_N_dict[0]["id"]
        N_ID2 = full_N_dict[1]["id"]
        N_ID3 = full_N_dict[2]["id"]
        N_ID4 = full_N_dict[3]["id"]
        # assign the N_xyz
        N1_xyz = full_N_dict[0]["xyz"]
        N2_xyz = full_N_dict[1]["xyz"]
        N3_xyz = full_N_dict[2]["xyz"]
        N4_xyz = full_N_dict[3]["xyz"]

    assert N_ID != None, "Nitrogens 1 were not found"
    assert N_ID2 != None, "Nitrogens 2 were not found"
    assert N_ID3 != None, "Nitrogens 3 were not found"
    assert N_ID4 != None, "Nitrogens 4 were not found"

    mean_N_xyz = np.mean(np.array([N1_xyz, N2_xyz, N3_xyz, N4_xyz]), axis=0)

    nitrogen_dict = {
        "mean_N_xyz": mean_N_xyz,
        "N_ID1": N_ID,
        "N_ID2": N_ID2,
        "N_ID3": N_ID3,
        "N_ID4": N_ID4,
        "N1_xyz": N1_xyz,
        "N2_xyz": N2_xyz,
        "N3_xyz": N3_xyz,
        "N4_xyz": N4_xyz,
    }

    return nitrogen_dict


def get_fe_positions(file):
    fe_ID, fe_xyz = None, None
    with open(file, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if "HETATM" in line[0] and ("FE" in line[2] or "FE" in line[1]):
            shift = 0
            if "FE" in line[1]:
                shift = -1
            fe_ID = f"{line[4+shift]}:{line[5+shift]}:{line[2+shift]}"
            fe_xyz = [line[6 + shift], line[7 + shift], line[8 + shift]]
            fe_xyz = [float(x) for x in fe_xyz]
            fe_xyz = np.array(fe_xyz)
            break

    return {"id": fe_ID, "xyz": fe_xyz}


def get_ligand_info(file, fe_xyz):
    best_crit_dist = 10.0
    fe_crit_dist = 10.0
    best_crit = None

    with open(file, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        sg_cond = "ATOM" in line[0] and "SG" in line[2] and "CYS" in line[3]
        oh_cond = "ATOM" in line[0] and "OH" in line[2] and "TYR" in line[3]
        nend_cond = (
            "ATOM" in line[0]
            and (("NE2" in line[2]) or ("ND1") in line[2])
            and ("HIS" in line[3] or "HD1" in line[3])
        )

        if sg_cond or oh_cond or nend_cond:
            crit_ID = f"{line[4]}:{line[5]}:{line[2]}"
            # catch any clumped elements
            x = line[6]
            y = line[7]
            z = line[8]

            if len(line[6]) > 8:
                x, y = break_up_line(x)
                z = line[7]

            if len(line[7]) > 8:
                y, z = break_up_line(y)

            crit_xyz = [x, y, z]
            crit_xyz = [float(x) for x in crit_xyz]
            crit_xyz = np.array(crit_xyz)
            crit_dist = np.linalg.norm(fe_xyz - crit_xyz)

            if crit_dist < fe_crit_dist and crit_dist < best_crit_dist:
                best_crit_dist = crit_dist
                best_crit = crit_ID

    ligand_dict = {
        "best_crit_dist": best_crit_dist,
        "best_crit": best_crit,
        "crit_xyz": crit_xyz,
    }

    return ligand_dict


def split_and_filter(mat, cutoff=95, min_max=True, std_mean=False):

    arr_mean, arr_std, arr_min, arr_max = (
        np.mean(mat),
        np.std(mat),
        np.min(mat),
        np.max(mat),
    )
    if min_max:
        mat = (mat - arr_min) / (arr_max - arr_min + 10e-10)

    if std_mean:
        mat = (mat - arr_mean) / (arr_std)

    try:
        u = mat[0][:, :, :, 0].flatten()
        v = mat[0][:, :, :, 1].flatten()
        w = mat[0][:, :, :, 2].flatten()
    except:
        u = mat[:, :, :, 0].flatten()
        v = mat[:, :, :, 1].flatten()
        w = mat[:, :, :, 2].flatten()

    component_distro = [
        np.sqrt(u[ind] ** 2 + v[ind] ** 2 + w[ind] ** 2) for ind in range(len(u))
    ]
    cutoff = np.percentile(component_distro, cutoff)

    for ind, i in enumerate(component_distro):
        if i < cutoff:
            u[ind], v[ind], w[ind] = 0, 0, 0

    u = np.around(u, decimals=2)
    v = np.around(v, decimals=2)
    w = np.around(w, decimals=2)

    return u, v, w


def mat_pull(file):

    with open(file) as f:
        lines = f.readlines()

    steps_x = 2 * int(lines[0].split()[2]) + 1
    steps_y = 2 * int(lines[0].split()[3]) + 1
    steps_z = 2 * int(lines[0].split()[4][:-1]) + 1
    mat = np.zeros((steps_x, steps_y, steps_z, 3))

    # gap_x = round(np.abs(float(lines[steps_x*steps_y + 7].split()[0]) - float(lines[7].split()[0])), 4)
    # gap_y = round(np.abs(float(lines[steps_x+8].split()[1]) - float(lines[7].split()[1])), 4)
    # gap_z = round(np.abs(float(lines[8].split()[2]) - float(lines[7].split()[2])), 4)

    for ind, i in enumerate(lines[7:]):
        line_split = i.split()
        # print(i)
        mat[
            int(ind / (steps_z * steps_y)),
            int(ind / steps_z % steps_y),
            ind % steps_z,
            0,
        ] = float(line_split[-3])
        mat[
            int(ind / (steps_z * steps_y)),
            int(ind / steps_z % steps_y),
            ind % steps_z,
            1,
        ] = float(line_split[-2])
        mat[
            int(ind / (steps_z * steps_y)),
            int(ind / steps_z % steps_y),
            ind % steps_z,
            2,
        ] = float(line_split[-1])
    return mat


def pull_mats_w_label(
    dir_data="../../../data/protein_data.csv", dir_fields="../../../data/cpet/"
):

    x, y = [], []
    df = pd.read_csv(dir_data)
    y_count, h_count, c_count = 0, 0, 0
    for row in df.iterrows():
        # print(row[1]['name'])
        cpet_name = dir_fields + "efield_cox_" + row[1]["name"] + ".dat"

        if os.path.exists(cpet_name):
            x.append(mat_pull(cpet_name))
            if row[1]["label"] == "Y":
                y.append([1, 0, 0])
                y_count += 1
            elif row[1]["label"] == "H":
                y.append([0, 1, 0])
                h_count += 1
            else:
                y.append([0, 0, 1])
                c_count += 1
    print(y_count, h_count, c_count)
    return np.array(x), np.array(y)


def aug_all(mat, target, xy=True, z=False, mut=False):
    full_aug = []
    full_aug_target = []

    for ind, i in enumerate(mat):
        x_aug, y_aug = augment_mat_field(i, target[ind])
        [full_aug.append(j) for j in x_aug]
        [full_aug_target.append(j) for j in y_aug]

    return np.array(full_aug), np.array(full_aug_target)


def augment_mat_field(mat, target, xy=True, z=False):
    aug_target = []
    aug_mat = []

    if xy:
        x_flip = np.array(np.flip(mat, axis=0), dtype=float)
        y_flip = np.array(np.flip(mat, axis=1), dtype=float)
        xy_flip = np.array(np.flip(np.flip(mat, axis=1), axis=0), dtype=float)

        x_flip[:, :, :, 0] = -1 * x_flip[:, :, :, 0]
        y_flip[:, :, :, 1] = -1 * y_flip[:, :, :, 1]
        xy_flip[:, :, :, 0] = -1 * xy_flip[:, :, :, 0]
        xy_flip[:, :, :, 1] = -1 * xy_flip[:, :, :, 1]

        aug_mat.append(mat)
        aug_mat.append(x_flip)
        aug_mat.append(y_flip)
        aug_mat.append(xy_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)

    if z:
        z_flip = np.array(np.flip(mat, axis=2), dtype=float)
        xz_flip = np.array(np.flip(np.flip(mat, axis=0), axis=2), dtype=float)
        yz_flip = np.array(np.flip(np.flip(mat, axis=1), axis=2), dtype=float)
        xyz_flip = np.array(
            np.flip(np.flip(np.flip(mat, axis=2), axis=1), axis=0), dtype=float
        )

        z_flip[:, :, :, 0] = -1 * z_flip[:, :, :, 0]
        xz_flip[:, :, :, 0] = -1 * xz_flip[:, :, :, 0]
        xz_flip[:, :, :, 2] = -1 * xz_flip[:, :, :, 2]
        yz_flip[:, :, :, 1] = -1 * yz_flip[:, :, :, 1]
        yz_flip[:, :, :, 2] = -1 * yz_flip[:, :, :, 2]
        xyz_flip[:, :, :, 0] = -1 * xyz_flip[:, :, :, 0]
        xyz_flip[:, :, :, 1] = -1 * xyz_flip[:, :, :, 1]
        xyz_flip[:, :, :, 2] = -1 * xyz_flip[:, :, :, 2]

        aug_mat.append(z_flip)
        aug_mat.append(xz_flip)
        aug_mat.append(yz_flip)
        aug_mat.append(xyz_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)

    return aug_mat, aug_target


def pca(mat, pca=None, pca_comps=10, verbose=False, write=False):

    mat_transform = mat.reshape(
        mat.shape[0], mat.shape[1] * mat.shape[2] * mat.shape[3] * mat.shape[4]
    )
    if pca == None:
        pca = PCA(n_components=pca_comps)
        mat_transform = pca.fit_transform(mat_transform)

    else:
        mat_transform = pca.transform(mat_transform)

    cum_explained_var = []
    for i in range(0, len(pca.explained_variance_ratio_)):
        if i == 0:
            cum_explained_var.append(pca.explained_variance_ratio_[i])
        else:
            cum_explained_var.append(
                pca.explained_variance_ratio_[i] + cum_explained_var[i - 1]
            )

    pc0 = pca.components_[0]
    # print(np.shape(pc0))
    pc0 = pc0.reshape(1, mat.shape[1], mat.shape[2], mat.shape[3], mat.shape[4])

    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "cone"}]])
    x, y, z = np.meshgrid(
        np.arange(-3.0, 3.3, 0.3), np.arange(-3.0, 3.3, 0.3), np.arange(-3.0, 3.3, 0.3)
    )
    u = pc0[0][:, :, :, 0].flatten()
    v = pc0[0][:, :, :, 1].flatten()
    w = pc0[0][:, :, :, 2].flatten()

    comp_vect_field = pc0.reshape(
        mat.shape[1], mat.shape[2], mat.shape[3], mat.shape[4]
    )

    u_1, v_1, w_1 = split_and_filter(
        comp_vect_field, cutoff=98, std_mean=True, min_max=False
    )

    vector_scale = 3
    fig.add_trace(
        go.Cone(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            u=u_1,
            v=v_1,
            w=w_1,
            sizeref=vector_scale,
        ),
        row=1,
        col=1,
    )

    if write:
        fig.write_html("./out_pca.html")

    if verbose:
        print(pca.explained_variance_ratio_)

    return mat_transform, pca


def unwrap_pca(mat, pca, shape):
    mat = pca.inverse_transform(mat)
    mat = mat.reshape(len(mat), shape[1], shape[2], shape[3], shape[4])
    return mat


def helmholtz_hodge_decomp_approx(
    file="../../data/cpet/efield_cox_1sj21.dat", show=False
):
    Vf = mat_pull(file)
    NX, NY, NZ = Vf[:, :, :, 1].shape
    # print(NX, NY, NZ)
    Vfx = Vf[:, :, :, 0]
    Vfy = Vf[:, :, :, 1]
    Vfz = Vf[:, :, :, 2]

    vx_f = np.fft.fftn(Vfx)
    vy_f = np.fft.fftn(Vfy)
    vz_f = np.fft.fftn(Vfz)

    kx = np.fft.fftfreq(NX).reshape(NX, 1, 1)
    ky = np.fft.fftfreq(NY).reshape(NY, 1)
    kz = np.fft.fftfreq(NZ)
    k2 = kx**2 + ky**2 + kz**2
    k2[0, 0, 0] = 1.0  # to avoid inf. we do not care about the k=0 component

    div_Vf_f = vx_f * kx + vy_f * ky + vz_f * kz  # * 1j
    V_compressive_overk = div_Vf_f / k2
    V_compressive_x = np.fft.ifftn(
        V_compressive_overk * kx
    )  # [:,np.newaxis,np.newaxis])
    V_compressive_y = np.fft.ifftn(V_compressive_overk * ky)
    V_compressive_z = np.fft.ifftn(V_compressive_overk * kz)

    V_solenoidal_x = Vfx - V_compressive_x
    V_solenoidal_y = Vfy - V_compressive_y
    V_solenoidal_z = Vfz - V_compressive_z

    # check if the solenoidal part really divergence-free
    divVs = np.fft.ifftn(
        (
            np.fft.fftn(V_solenoidal_x) * kx
            + np.fft.fftn(V_solenoidal_y) * ky
            + np.fft.fftn(V_solenoidal_z) * kz
        )
        * 1j
        * 2.0
        * np.pi
    )

    # print('div_solenoidal max:', abs(divVs).max())
    # check the power in solenoidal and compressive components
    # print('variance:')
    # print( 'solenoidal x,y,z:', V_solenoidal_x.var(), V_solenoidal_y.var(), V_solenoidal_z.var())
    # print('compressive x,y,z:', V_compressive_x.var(), V_compressive_y.var(), V_compressive_z.var())

    if show:
        # plot one slice of the decomposed field on X-Y plane
        X, Y = np.meshgrid(range(NY), range(NX))
        scale = 1
        plt.figure()
        plt.quiver(
            X,
            Y,
            V_solenoidal_x[:, :, 0] / scale,
            V_solenoidal_y[:, :, 0] / scale,
            norm=True,
        )
        plt.figure()
        plt.quiver(
            X,
            Y,
            V_compressive_x[:, :, 0] / scale,
            V_compressive_y[:, :, 0] / scale,
            norm=True,
        )
        plt.ion()
        plt.show()

    solenoidal = {
        "x": V_solenoidal_x.real,
        "y": V_solenoidal_y.real,
        "z": V_solenoidal_z.real,
    }
    compressize = {
        "x": V_compressive_x.real,
        "y": V_compressive_y.real,
        "z": V_compressive_z.real,
    }
    return solenoidal, compressize


def compress(distance_matrix):
    compressed_dictionary = {}
    affinity = AffinityPropagation(affinity="precomputed")
    affinity.fit(distance_matrix)
    cluster_centers_indices = affinity.cluster_centers_indices_
    labels = list(affinity.labels_)
    n_clusters_ = len(cluster_centers_indices)

    print(f"Estimated number of clusters: {n_clusters_}")
    # get count of a value in a list
    for i in range(n_clusters_):
        compressed_dictionary[i] = {
            "count": str(labels.count(i)),
            "index_center": str(cluster_centers_indices[i]),
        }
    return compressed_dictionary


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


def fetch_charges_dict(file_name="test.pqr"):
    """
    Given a list of dictionaries with element and position, traverse a pqr file and get the charges from the file
    EXCLUDING ELEMENTS IN THE LIST OF DICTIONARIES
    Takes:
        list of dictionaries with element and position
    Returns:
        list of dictionaries with element, position and charge
    """

    pqr_dict = []
    # get the lines of the pqr file
    with open(file_name, "r") as f:
        lines = f.readlines()

    for line in lines:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        charge = float(line[54:61].strip())
        radius = float(line[62:68].strip())
        if np.abs(charge) >= 0.01:
            pqr_dict.append({"position": [x, y, z], "charge": charge, "radius": radius})

    return pqr_dict


def get_elements(file_name):
    """
    Open xyz and get all the elements in the file
    Takes:
        file_name: the name of the xyz file
    Returns:
        elements: a list of elements
    """
    elements = []
    with open(file_name, "r") as f:
        lines = f.readlines()
    for line in lines:
        if line[0:2].strip().isalpha():
            elements.append(line.split()[0])
    return elements


def put_charges_in_turbo_files(folder_name, charges_dict):
    """
    Traverses subdirectories in embedding folder and puts charges in the turbomole file

    Takes
        folder_name: the folder where the turbomole files are located
        charges_dict: a dictionary with the charges
    Returns: Nothing

    """
    # find folder named embedding and go into all subfolders
    for root, dirs, files in os.walk(folder_name):
        for file in files:
            if file.endswith("control"):
                print("editing control file with charges from pqr dictionary")
                # remove last line of file - the $end
                with open(os.path.join(root, file), "r") as f:
                    lines = f.readlines()
                with open(os.path.join(root, file), "w") as f:
                    for line in lines[:-1]:
                        f.write(line)

                # append dictionary to end of file
                with open(os.path.join(root, file), "a") as f:
                    f.write("$point_charges\n")
                    for charge in charges_dict:
                        f.write(
                            "\t{} {} {} {}\n".format(
                                charge["position"][0],
                                charge["position"][1],
                                charge["position"][2],
                                charge["charge"],
                            )
                        )
                        # f.write("CHARGE " + str(charge["charge"]) + " " + str(charge["radius"]) + " " + str(charge["position"][0]) + " " + str(charge["position"][1]) + " " + str(charge["position"][2]) + "")
                    f.write("$end\n")


def get_frozen_atoms(file_name):
    """
    get the two carbons most out of the plane to freeze
    Takes:
        file_name: the name of the xyz file
    Returns:
        frozen_atoms: a binary list of frozen atoms
    """

    carbon_xyz, ind_carbons = get_carbon_xyz_from_file(file_name)

    cross = get_cross_vector(file_name)
    fe_dict = get_fe_positions(file_name)
    n_dict = get_N_positions(file_name, fe_dict["xyz"])

    mean_xyz = n_dict["mean_N_xyz"]
    dot_list = [np.dot(i[0] - n_dict["mean_N_xyz"], cross) for i in carbon_xyz]
    dot_list = np.array(dot_list) / np.linalg.norm(cross)

    # filter for coplanar carbons
    carbon_planar_ind = []  # the indices of the coplanar carbons
    for ind, i in enumerate(carbon_xyz):
        if dot_list[ind] < 0.5:
            carbon_planar_ind.append(ind)

    # get the four furthest, in plane carbons
    carbon_planar_xyz = np.array(carbon_xyz)[carbon_planar_ind]
    diff = carbon_planar_xyz - mean_xyz
    distances = np.apply_along_axis(np.linalg.norm, 1, diff)

    furthest_ind = np.argsort(distances)[::-1][:4]
    most_out_of_plane_ind = np.argsort(dot_list)[::-1][:2]
    # combine the two lists
    frozen_atom_ind = np.concatenate((furthest_ind, most_out_of_plane_ind))
    return_list = [ind_carbons[i] for i in frozen_atom_ind]

    return return_list


def get_carbon_xyz_from_file(file_name):
    """
    get the xyz coordinates of the carbons in the file
    Takes:
        file_name: the name of the xyz file
    Returns:
        carbon_xyz: a list of the xyz coordinates of the carbons
    """
    carbon_xyz, ind = [], []
    with open(file_name, "r") as f:
        lines = f.readlines()
        # go through all the lines and find nitrogens
        for line_ind, line in enumerate(lines):
            if line[0:2].strip().isalpha():
                if line.split()[0] == "C" or line.split()[0] == "c":
                    carbon_xyz.append(
                        [
                            float(line.split()[1]),
                            float(line.split()[2]),
                            float(line.split()[3]),
                        ]
                    )
                    ind.append(line_ind - 1)
    return carbon_xyz, ind


def get_cross_vector(file_name):

    # find the four nitrogens closest to the iron
    fe_info = get_fe_positions(file_name)
    fe_xyz = fe_info["xyz"]
    fe_ID = fe_info["id"]
    nitrogen_info = get_N_positions(file_name, fe_ID, fe_xyz)

    direction_1 = nitrogen_info["N1_xyz"] - nitrogen_info["mean_N_xyz"]
    direction_2 = nitrogen_info["N2_xyz"] - nitrogen_info["mean_N_xyz"]
    direction_3 = nitrogen_info["N3_xyz"] - nitrogen_info["mean_N_xyz"]
    # compute cross and take the two most orthogonal directions
    dot_12 = np.dot(direction_1, direction_2)
    dot_13 = np.dot(direction_1, direction_3)
    dot_23 = np.dot(direction_2, direction_3)

    if dot_23 > dot_13 and dot_23 > dot_12:
        direction_1 = direction_3
    if dot_13 > dot_12 and dot_13 > dot_23:
        direction_2 = direction_3

    direction_1 /= np.linalg.norm(direction_1)
    direction_2 /= np.linalg.norm(direction_2)
    cross = np.cross(direction_1, direction_2)
    cross = cross / np.linalg.norm(cross)
    return cross
