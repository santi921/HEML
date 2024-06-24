import os, re, argparse
from glob import glob
import numpy as np
from HEML.utils.data import (
    get_options,
    check_if_file_is_empty,
    get_c1_positions,
    get_fe_positions,
)

from HEML.utils.mol2topqr import mol2_to_pqr_folder

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    parser.add_argument("--zero_active", help="zero active site", default=True)
    parser.add_argument("--zero_radius", help="tf zero by radius", action="store_true")

    options_loc = parser.parse_args().options
    zero_active = parser.parse_args().zero_active
    zero_radius = parser.parse_args().zero_radius

    options = get_options(options_loc)
    outdir = options["processed_charges_folder"]
    outdir_cpet = options["cpet_folder"]
    charges_directory = options["charges_folder"]

    if zero_active:
        ligands_to_zero = options["ligands_to_zero"]
        print("zeroing active site for: {}".format(ligands_to_zero))

    if zero_radius:
        ligands_to_zero_radius = options["zero_radius"]
        print("zeroing active site radius from iron: {}".format(ligands_to_zero_radius))

    print(outdir)
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
            openfile = open(i)
            readfile = openfile.readlines()
            openfile.close
            filename = i.split(".")[0]

            # Look for FE in each line, and get chain:res:atom ID and xyz coords if exists.
            fail_cond = True

            try:
                c1_dict = get_c1_positions(i)
                fe_dict = get_fe_positions(i)
                assert c1_dict["id"] != None
                print(c1_dict["id"], c1_dict["id"])
                fail_cond = False

            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1

            filename = os.path.basename(i)
            listname = filename.split(".")
            if not fail_cond:
                with open(output, "w") as outfile:
                    print("processing from file: {}".format(i))

                    for j in readfile:

                        tf_zero = False
                        # grab from column 17 to 21 inclusive
                        lig_str = j[17:21].strip()
                        
                        if zero_active and lig_str in ligands_to_zero:
                            tf_zero = True
                            print("zeroing ligand", lig_str)

                        if not tf_zero and zero_radius:
                            xyz_str_x = j[30:38]
                            xyz_str_y = j[38:46]
                            xyz_str_z = j[46:54]
                            distance = np.sqrt(
                                (float(xyz_str_x) - float(fe_dict["xyz"][0])) ** 2
                                + (float(xyz_str_y) - float(fe_dict["xyz"][1])) ** 2
                                + (float(xyz_str_z) - float(fe_dict["xyz"][2])) ** 2
                            )

                            if distance < ligands_to_zero_radius:
                                tf_zero = True
                                print(
                                    "zeroing distance:{} w/ dist {}".format(
                                        lig_str, distance
                                    )
                                )

                            if tf_zero:
                                temp_write = j[:56] + "0.000" + j[61:]
                                outfile.write(temp_write)
                            else:
                                outfile.write(j)
                        else:
                            outfile.write(j)

                file_name = i.split("/")[-1].split(".")[0]
                options = open(f"{outdir_cpet}options_magni_{file_name}.txt", "w+")
                options.write(f"%field \n")
                c1_dict = get_c1_positions(i)
                options.write(
                    "    locations {}:{}:{} \n".format(
                        c1_dict["xyz"][0], c1_dict["xyz"][1], c1_dict["xyz"][2]
                    )
                )
                options.write(f"output {outdir_cpet}{file_name}_mag.dat \n")
                options.write(f"end \n")
                options.close()

    print("fail count: {}".format(fail))
