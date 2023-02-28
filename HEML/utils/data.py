import os, json
import pandas as pd
import numpy as np
from copy import deepcopy
from HEML.utils.dictionaries import * 


def get_options(options_file="./options.json", create_folders=True):
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
        sg_cond = (
            "ATOM" in line[0] 
            and ("SG" in line[2])
            and ("CY1" in line[3] or "CYS" in line[3]))
        
        oh_cond = (
            "ATOM" in line[0] 
            and ("OH" in line[2])
            and ("TYR" in line[3] or "TR1" in line[3]))
        
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


def mat_pull(file, meta_data=False):

    with open(file) as f:
        lines = f.readlines()

    if meta_data:
        
        steps_x = 2 * int(lines[0].split()[2]) + 1
        steps_y = 2 * int(lines[0].split()[3]) + 1
        steps_z = 2 * int(lines[0].split()[4][:-1]) + 1
        x_size = float(lines[0].split()[-3])
        y_size = float(lines[0].split()[-2])
        z_size = float(lines[0].split()[-1])

        meta_dict = {
            "steps_x": steps_x,
            "steps_y": steps_y,
            "steps_z": steps_z,
            "step_size_x": np.round(x_size / float(lines[0].split()[2]), 4),
            "step_size_y": np.round(y_size / float(lines[0].split()[3]), 4),
            "step_size_z": np.round(z_size / float(lines[0].split()[4][:-1]), 4),
            "first_line": lines[0]
        }

        return meta_dict
    else: 
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
    data_file="../../../data/protein_data.csv", dir_fields="../../../data/cpet/"
):

    x, y = [], []
    df = pd.read_csv(data_file)
    print(df.shape)
    y_count, h_count, c_count = 0, 0, 0
    for row in df.iterrows():
        #print(row[1]['name'])
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
