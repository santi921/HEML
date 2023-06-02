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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    parser.add_argument("--zero_active", help="zero active site", default=True)
    parser.add_argument(
        "--zero_everything_charged", help="zero everything charged", default=False
    )
    # store t/f if --box is used
    parser.add_argument("--box", help="box", action="store_true")
    parser.add_argument("--box_size", help="box size", default=4.0)
    parser.add_argument("--density", help="density", default=10)
    parser.add_argument("--samples", help="samples", default=10000)
    parser.add_argument("--bins", help="bins", default=25)
    parser.add_argument("--step_size", help="step size", default=0.001)
    parser.add_argument("--zero_radius", help="tf zero by radius", action="store_true")
    parser.add_argument(
        "--carbene",
        help="if computing on carbene attached to heme, this will center the calc at the mean of the fe-c bond",
        action="store_true",
    )

    options_loc = parser.parse_args().options
    zero_active = parser.parse_args().zero_active
    zero_everything_charged = parser.parse_args().zero_everything_charged
    box = bool(parser.parse_args().box)
    box_size = float(parser.parse_args().box_size)
    density = int(parser.parse_args().density)
    samples = int(parser.parse_args().samples)
    bins = int(parser.parse_args().bins)
    step_size = float(parser.parse_args().step_size)
    zero_radius = bool(parser.parse_args().zero_radius)
    carbene_tf = bool(parser.parse_args().carbene)

    options = get_options(options_loc)
    outdir = options["processed_charges_folder"]
    outdir_cpet = options["cpet_folder"]
    charges_directory = options["charges_folder"]
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

            ligand_dict = {}
            nitrogen_dict = {}

            try:
                fe_dict = get_fe_positions(i)
                assert fe_dict["id"] != None
                print(fe_dict["id"], fe_dict["xyz"])
                ligand_dict = get_ligand_info(i, fe_dict["xyz"])
                nitrogen_dict = get_N_positions(i, fe_dict["id"], fe_dict["xyz"])
                nitro_none = check_if_dict_has_None(nitrogen_dict)
                ligand_none = check_if_dict_has_None(ligand_dict)
                if carbene_tf:
                    c1_dict = get_c1_positions(i)
                    carbene_none = check_if_dict_has_None(c1_dict)
                    
                    if carbene_none:
                        print("carbene none")
                        fail_cond = True
                    else:
                        print(fe_dict["xyz"], c1_dict["xyz"])
                        print(np.mean([fe_dict["xyz"], c1_dict["xyz"]], axis=0))
                        print(np.array(fe_dict["xyz"]) - np.array(c1_dict["xyz"]))
                        mean_xyz = np.mean([fe_dict["xyz"], c1_dict["xyz"]], axis=0)
                        fail_cond = False

                if not nitro_none and not ligand_none:
                    fail_cond = False

                    if ligand_dict["best_crit_dist"] > 4.0:
                        print(ligand_dict["best_crit_dist"])
                        print(
                            f"ERROR: No cysteine/tyrosine/histine ligand found for {i}.\n"
                        )
                        fail += 1
                        continue

                    else:
                        print(f'Nitro 1 {nitrogen_dict["N_ID1"]}')
                        print(f'Nitro 2 {nitrogen_dict["N_ID2"]}')
                        print(f'Nitro 3 {nitrogen_dict["N_ID3"]}')
                        print(f'Nitro 4 {nitrogen_dict["N_ID4"]}')
                        print(f'ligand of note {ligand_dict["best_crit"]}\n')
                        center = [float(i) for i in nitrogen_dict["mean_N_xyz"]]

            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1

            filename = os.path.basename(i)
            listname = filename.split(".")

            if not fail_cond:
                ligand_identifier = ligand_dict["best_crit"].split(":")

                with open(output, "w") as outfile:
                    for j in readfile:
                        line_split = re.split(r"(\s+)", j)
                        # remove items in list that are empty strings
                        #line_split = [x for x in line_split if x != ""]
                        # strip whitespace from list
                        line_split = [x.strip() for x in line_split]
                        # remove empty strings from list
                        line_split = list(filter(None, line_split))
                        #print(line_split)

                        cond = (
                            line_split[8] == ligand_identifier[0]
                            and line_split[10] == ligand_identifier[1]
                        )
                        if zero_active and ("HETATM" in line_split[0] or cond):
                            temp_write = j[:56] + "0.000" + j[61:]
                            outfile.write(temp_write)

                        elif carbene_tf and line_split[3] in ["CB1"]:
                            temp_write = j[:56] + "0.000" + j[61:]
                            outfile.write(temp_write)                        

                        elif zero_everything_charged and line_split[3] in [
                            "ASP",
                            "GLU",
                            "LYS",
                            "ARG",
                            "HIS",
                        ]:
                            temp_write = j[:56] + "0.000" + j[61:]
                            outfile.write(temp_write)

                        elif zero_radius:
                            xyz_str_x = j[30:38]
                            xyz_str_y = j[38:46]
                            xyz_str_z = j[46:54]
                            if carbene_tf:
                                distance = np.sqrt(
                                    (float(xyz_str_x) - float(mean_xyz[0])) ** 2
                                    + (float(xyz_str_y) - float(mean_xyz[1])) ** 2
                                    + (float(xyz_str_z) - float(mean_xyz[2])) ** 2
                                )

                            else:
                                distance = np.sqrt(
                                    (float(xyz_str_x) - float(fe_dict["xyz"][0])) ** 2
                                    + (float(xyz_str_y) - float(fe_dict["xyz"][1])) ** 2
                                    + (float(xyz_str_z) - float(fe_dict["xyz"][2])) ** 2
                                )

                            if distance < ligands_to_zero_radius:
                                lig_str = j[17:21].strip()
                                print(
                                    "zeroing distance:{} w/ dist {}".format(
                                        lig_str, distance
                                    )
                                )
                                temp_write = j[:56] + "0.000" + j[61:]
                                outfile.write(temp_write)
                            else:
                                outfile.write(j)
                        else:
                            outfile.write(j)

                file_name = i.split("/")[-1].split(".")[0]
                if box:
                    options = open(f"{outdir_cpet}options_field_{file_name}.txt", "w+")
                    if carbene_tf:
                        nitro_axis_1 = [nitrogen_dict["N1_xyz"][0]-nitrogen_dict["mean_N_xyz"][0]+mean_xyz[0], nitrogen_dict["N1_xyz"][1]-nitrogen_dict["mean_N_xyz"][1]+mean_xyz[1], nitrogen_dict["N1_xyz"][2]-nitrogen_dict["mean_N_xyz"][2]+mean_xyz[2]]
                        nitro_axis_2 = [nitrogen_dict["N2_xyz"][0]-nitrogen_dict["mean_N_xyz"][0]+mean_xyz[0], nitrogen_dict["N2_xyz"][1]-nitrogen_dict["mean_N_xyz"][1]+mean_xyz[1], nitrogen_dict["N2_xyz"][2]-nitrogen_dict["mean_N_xyz"][2]+mean_xyz[2]]
                        
                        options.write(
                            f'align {mean_xyz[0]}:{mean_xyz[1]}:{mean_xyz[2]} {nitro_axis_1[0]}:{nitro_axis_1[1]}:{nitro_axis_1[2]} {nitro_axis_2[0]}:{nitro_axis_2[1]}:{nitro_axis_2[2]}\n'
                        )
                    else:
                        options.write(
                            f'align {nitrogen_dict["mean_N_xyz"][0]}:{nitrogen_dict["mean_N_xyz"][1]}:{nitrogen_dict["mean_N_xyz"][2]} {nitrogen_dict["N1_xyz"][0]}:{nitrogen_dict["N1_xyz"][1]}:{nitrogen_dict["N1_xyz"][2]} {nitrogen_dict["N2_xyz"][0]}:{nitrogen_dict["N2_xyz"][1]}:{nitrogen_dict["N2_xyz"][2]}\n'
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
                    options.write(f"output {outdir_cpet}efield_cox_{file_name}.dat \n")
                    options.write(f"end \n")
                    options.close()
                else:

                    options = open(f"{outdir_cpet}options_topol_{file_name}.txt", "w+")
                    if carbene_tf:
                        nitro_axis_1 = [nitrogen_dict["N1_xyz"][0]-nitrogen_dict["mean_N_xyz"][0]+mean_xyz[0], nitrogen_dict["N1_xyz"][1]-nitrogen_dict["mean_N_xyz"][1]+mean_xyz[1], nitrogen_dict["N1_xyz"][2]-nitrogen_dict["mean_N_xyz"][2]+mean_xyz[2]]
                        nitro_axis_2 = [nitrogen_dict["N2_xyz"][0]-nitrogen_dict["mean_N_xyz"][0]+mean_xyz[0], nitrogen_dict["N2_xyz"][1]-nitrogen_dict["mean_N_xyz"][1]+mean_xyz[1], nitrogen_dict["N2_xyz"][2]-nitrogen_dict["mean_N_xyz"][2]+mean_xyz[2]]
                        options.write(
                            f'align {mean_xyz[0]}:{mean_xyz[1]}:{mean_xyz[2]} {nitro_axis_1[0]}:{nitro_axis_1[1]}:{nitro_axis_1[2]} {nitro_axis_2[0]}:{nitro_axis_2[1]}:{nitro_axis_2[2]}\n'
                        )
                    else:
                        options.write(
                            f'align {nitrogen_dict["mean_N_xyz"][0]}:{nitrogen_dict["mean_N_xyz"][1]}:{nitrogen_dict["mean_N_xyz"][2]} {nitrogen_dict["N1_xyz"][0]}:{nitrogen_dict["N1_xyz"][1]}:{nitrogen_dict["N1_xyz"][2]} {nitrogen_dict["N2_xyz"][0]}:{nitrogen_dict["N2_xyz"][1]}:{nitrogen_dict["N2_xyz"][2]}\n'
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
                        "    sampleOutput efield_topo_{file_name} \n \n"
                    )
                    options.write("    bins {} \n".format(bins))
                    options.write(f"end \n")
                    options.close()

    print("fail count: {}".format(fail))
