import numpy as np
import pandas as pd
import os
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist, squareform
from io import StringIO

import os, re, argparse
import numpy as np
from glob import glob
from HEML.utils.data import (
    get_options,
    check_if_file_is_empty,
    get_fe_positions,
    get_c1_positions,
    get_ligand_info,
    get_N_positions,
    check_if_dict_has_None,
)
from HEML.utils.mol2topqr import mol2_to_pqr_folder


def find_closest_atoms(df_re, df_all):
    # Create a result list
    closest_atoms = []
    # Extract coordinates from the dataframes
    coords_re = df_re[["x", "y", "z"]].values
    coords_all = df_all[["x", "y", "z"]].values
    # Total number of Re atoms
    total_re = len(coords_re)
    # Iterate over each Re atom
    for i, re_atom in enumerate(coords_re):
        # Print progress
        # print(f"Processing Re atom {i+1} of {total_re}...")
        # Calculate the distances from the current Re atom to all other atoms
        distances = cdist([re_atom], coords_all)[0]
        # Add the distances to df_all as a new column
        df_all["distance"] = distances
        # Sort df_all by distance
        df_sorted = df_all.sort_values("distance")
        # Get the three closest carbons and closest chloride
        carbons = df_sorted[df_sorted["Atom"] == "C"][:3]
        chloride = df_sorted[df_sorted["Atom"] == "Cl"].iloc[0]
        # Combine Re atom, closest carbons, and closest chloride
        df_atoms = pd.concat([df_re.iloc[[i]], carbons, pd.DataFrame([chloride])])
        # Append it to the result list
        closest_atoms.append(df_atoms)
        # Drop the distance column for the next iteration
        df_all.drop("distance", axis=1, inplace=True)
    return closest_atoms


def remove_furthest_carbon(df_nearest_C_Cl):
    # Create a result list
    modified_dfs = []
    # Iterate over each dataframe in the input list
    for df in df_nearest_C_Cl:
        # Get the coordinates of the atoms
        coords = df[["x", "y", "z"]].values
        # Calculate the pairwise distances between atoms
        distances = squareform(pdist(coords))
        # Get the distance from each carbon to the chloride
        # Assuming the chloride is the last row in the dataframe
        carbon_distances = distances[-1, 1:4]
        # Get the index of the carbon that is furthest from the chloride
        furthest_carbon_index = np.argmax(carbon_distances) + 1
        # Get the label of the furthest carbon from the DataFrame's index
        furthest_carbon_label = df.index[furthest_carbon_index]
        # print(furthest_carbon_label)
        # Drop the furthest carbon
        df_modified = df.drop(furthest_carbon_label)
        # Append the modified dataframe to the result list
        modified_dfs.append(df_modified)
    return modified_dfs


def calculate_vectors_and_cross_product(df_list):
    results = []

    # Iterate over each dataframe in the input list
    for df in df_list:
        dict_info = {}
        # Get the Re atom (assumed to be the first row in the dataframe)
        re_atom = df.iloc[0]
        re_coords = re_atom[["x", "y", "z"]].values.astype(float)
        # Get the Cl atom (assumed to be the last row in the dataframe)
        cl_atom = df.iloc[-1]
        cl_coords = cl_atom[["x", "y", "z"]].values.astype(float)
        # Get the carbon atoms (assumed to be the second and third rows in the dataframe)
        c_atoms = df.iloc[1:3]
        c_coords = c_atoms[["x", "y", "z"]].values.astype(float)
        # Calculate vectors from carbons to Re
        c_vectors = c_coords - re_coords
        # Calculate vector from Cl to Re
        cl_vector = cl_coords - re_coords
        # Compute the cross products of the carbon vectors in both orders
        cross_product_1 = np.cross(c_vectors[0], c_vectors[1])
        cross_product_2 = np.cross(c_vectors[1], c_vectors[0])
        # normalize
        cross_product_1 = cross_product_1 / np.linalg.norm(cross_product_1)
        cross_product_2 = cross_product_2 / np.linalg.norm(cross_product_2)
        # Decide which carbon to label as "carbon x" and which to label as "carbon y"
        if np.linalg.norm(cross_product_1 - cl_vector) > np.linalg.norm(
            cross_product_2 - cl_vector
        ):
            c_atoms.index = ["carbon x", "carbon y"]

        else:
            c_atoms.index = ["carbon y", "carbon x"]

        # Combine the atoms to a single dataframe and append to results
        # results.append(df_result)
        dict_info["center"] = re_coords
        dict_info["axis_1"] = c_atoms.loc["carbon x"][["x", "y", "z"]].values
        dict_info["axis_2"] = c_atoms.loc["carbon y"][["x", "y", "z"]].values
        dict_info["cross_product_1"] = cross_product_1
        results.append(dict_info)
    return results


def get_rhe_dictionaries(pdb):
    with open(pdb, "r") as f:
        lines = f.readlines()
    filtered_lines = [line for line in lines if line.startswith(("HETATM"))]
    filtered_data = "\n".join(filtered_lines)
    data = StringIO(filtered_data)
    print(data)

    column_names = [
        "Name",
        "num",
        "Atom_raw",
        "Res",
        "Chain",
        "res_num",
        "x",
        "y",
        "z",
        "radius",
        "charge",
    ]

    df = pd.read_csv(data, delimiter=" ", skipinitialspace=True, names=column_names)
    # print(df.iloc[0])
    # strip numbers from atom_raw
    df["Atom"] = df["Atom_raw"].str.replace("\d+", "")
    df = df.loc[~df["Atom"].str.contains("H")]
    df = df.loc[~df["Atom"].str.contains("N")]
    # cols_to_drop = ["Atom_full", "Res", "a", "b", "Res #"]
    # df = df.drop(cols_to_drop, axis=1)
    re_df = df[df["Atom"] == "RE"].copy()
    # print(re_df.iloc[0])
    re_df = re_df.reset_index().drop(columns="index")
    df = df.reset_index().drop(columns="index")

    test = find_closest_atoms(re_df, df)
    modified_list = remove_furthest_carbon(test)
    final_result = calculate_vectors_and_cross_product(modified_list)

    return final_result


def write_charges_processed(
    output,
    pqr_file,
    center_xyz,
    zero_active=False,
    zero_radius=False,
    ligands_to_zero_radius=3.0,
):
    openfile = open(pqr_file)
    readfile = openfile.readlines()
    openfile.close
    with open(output, "w") as outfile:
        for j in readfile:
            line_split = re.split(r"(\s+)", j)
            # strip whitespace from list
            line_split = [x.strip() for x in line_split]
            # remove empty strings from list
            line_split = list(filter(None, line_split))

            if zero_active and ("HETATM" in line_split[0]):
                temp_write = j[:56] + "0.000" + j[61:]
                outfile.write(temp_write)

            elif zero_radius:
                xyz_str_x = j[30:38]
                xyz_str_y = j[38:46]
                xyz_str_z = j[46:54]
                distance = np.sqrt(
                    (float(xyz_str_x) - float(center_xyz[0])) ** 2
                    + (float(xyz_str_y) - float(center_xyz[1])) ** 2
                    + (float(xyz_str_z) - float(center_xyz[2])) ** 2
                )

                if distance < ligands_to_zero_radius:
                    lig_str = j[17:21].strip()
                    print("zeroing distance:{} w/ dist {}".format(lig_str, distance))
                    # print(xyz_str_x, xyz_str_y, xyz_str_z)
                    # print(j)
                    temp_write = j[:56] + "0.000" + j[61:]
                    outfile.write(temp_write)
                else:
                    outfile.write(j)
            else:
                outfile.write(j)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    # parser.add_argument("--zero_active", help="zero active site", default=True)
    # parser.add_argument(
    #    "--zero_everything_charged", help="zero everything charged", default=False
    # )
    # store t/f if --box is used
    parser.add_argument("--box", help="box", action="store_true")
    parser.add_argument("--box_size", help="box size", default=4.0)
    parser.add_argument("--density", help="density", default=10)
    parser.add_argument("--samples", help="samples", default=10000)
    parser.add_argument("--bins", help="bins", default=25)
    parser.add_argument("--step_size", help="step size", default=0.001)
    parser.add_argument("--zero_radius", help="tf zero by radius", action="store_true")

    pdb = "/home/santiagovargas/Downloads/movie_phenanthroline_wt_ryan_CO_tetramer_prepped_1_minenergy.pdb_renumbered.pdb"

    options_loc = parser.parse_args().options
    # zero_active = parser.parse_args().zero_active
    # zero_everything_charged = parser.parse_args().zero_everything_charged
    box = bool(parser.parse_args().box)
    box_size = float(parser.parse_args().box_size)
    density = int(parser.parse_args().density)
    samples = int(parser.parse_args().samples)
    bins = int(parser.parse_args().bins)
    step_size = float(parser.parse_args().step_size)
    zero_radius = bool(parser.parse_args().zero_radius)

    options = get_options(options_loc)
    outdir = options["processed_charges_folder"]
    outdir_cpet = options["cpet_folder"]
    charges_directory = options["charges_folder"]
    if zero_radius:
        ligands_to_zero_radius = options["zero_radius"]
        print("zeroing active site radius from iron: {}".format(ligands_to_zero_radius))

    fail = 0
    filelist = glob(charges_directory + "*pqr")
    # check if there are no files in the directory with the correct extension
    if len(filelist) == 0:
        print(
            "No files in the directory with the correct extension, checking mol2 and converting to pqr"
        )
        filelist = glob(charges_directory + "*mol2")
        if len(filelist) == 0:
            print("No files in the directory with the correct extension, exiting")
            exit()
        else:
            mol2_to_pqr_folder(charges_directory)
            filelist = glob(charges_directory + "*pqr")

    for i in filelist:
        # new filename
        filename = os.path.basename(i)
        listname = filename.split(".")

        # check that output isn't already there in the processed directory
        output = f"{outdir}{listname[0]}.pqr"
        if os.path.exists(output):
            print("output file already exists")

        # checks that input file is not empty
        if check_if_file_is_empty(i):
            print("input file is empty")
            pass

        else:
            fail_cond = False

            try:
                rhe_list = get_rhe_dictionaries(i)
                assert len(rhe_list) == 4, "rhe list is not 4 long"

            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1
            # print(rhe_list)
            filename = os.path.basename(i)
            listname = filename.split(".")

            if not fail_cond:
                print("not failed")
                for ind, rhe_dict in enumerate(rhe_list):
                    center_xyz = rhe_dict["center"]
                    axis_1 = rhe_dict["axis_1"]
                    axis_2 = rhe_dict["axis_2"]

                    output = f"{outdir}{listname[0]}_{ind}.pqr"

                    if os.path.exists(output):
                        print("output file already exists")

                    print("writing to {}".format(output))
                    # print("center: {}".format(center_xyz))
                    # print("axis_1: {}".format(axis_1))
                    # print("axis_2: {}".format(axis_2))
                    write_charges_processed(
                        output,
                        i,
                        center_xyz,
                        zero_active=False,
                        zero_radius=zero_radius,
                        ligands_to_zero_radius=ligands_to_zero_radius,
                    )
                    file_name = i.split("/")[-1].split(".")[0]

                    if box:
                        options = open(
                            f"{outdir_cpet}options_field_{file_name}_{ind}.txt", "w+"
                        )

                        options.write(
                            f"align {center_xyz[0]}:{center_xyz[1]}:{center_xyz[2]} {axis_1[0]}:{axis_1[1]}:{axis_1[2]} {axis_2[0]}:{axis_1[1]}:{axis_1[2]}\n"
                        )
                        options.write(f"%plot3d \n")
                        options.write(f"    show false \n")
                        options.write(
                            "    volume box {} {} {} \n".format(
                                box_size, box_size, box_size
                            )
                        )
                        options.write(
                            "    density {} {} {} \n".format(density, density, density)
                        )
                        options.write(
                            "output {}efield_cox_{}_{}.dat \n".format(
                                outdir_cpet, file_name, ind
                            )
                        )
                        options.write(f"end \n")
                        options.close()
                    else:
                        options = open(
                            f"{outdir_cpet}options_topol_{file_name}_{ind}.txt", "w+"
                        )

                        options.write(
                            f"align {center_xyz[0]}:{center_xyz[1]}:{center_xyz[2]} {axis_1[0]}:{axis_1[1]}:{axis_1[2]} {axis_2[0]}:{axis_1[1]}:{axis_1[2]}\n"
                        )
                        options.write(f"%topology \n")
                        options.write(
                            "    volume box {} {} {} \n".format(
                                box_size, box_size, box_size
                            )
                        )
                        options.write("    stepSize {} \n".format(step_size))
                        options.write("    samples {} \n".format(samples))
                        options.write(
                            "    sampleOutput efield_topo_{}_{} \n".format(
                                file_name, ind
                            )
                        )
                        options.write("    bins {} \n".format(bins))
                        options.write(f"end \n")
                        options.close()

    print("fail count: {}".format(fail))
