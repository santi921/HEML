import os, re, argparse
from glob import glob
from HEML.utils.data import (
    get_options, 
    check_if_file_is_empty,
    get_fe_positions,
    get_ligand_info, 
    get_N_positions,
    check_if_dict_has_None
)
from HEML.utils.mol2topqr import mol2_to_pqr_folder

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    parser.add_argument(
        "--zero_active", help="zero active site", default=True
    )
    parser.add_argument(
        "--zero_everything_charged", help="zero everything charged", default=False
    )
    # store t/f if --box is used
    parser.add_argument(
        "--box", help="box", action="store_true"
        )

    parser.add_argument(
        "--box_size", help="box size", default=4.0
    )

    parser.add_argument(
        "--density", help="density", default=10
    )
    parser.add_argument(
        "--samples", help="samples", default=3000
    )
    parser.add_argument(
        "--bins", help="bins", default=25
    )
    parser.add_argument(
        "--step_size", help="step size", default=0.001
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

    options = get_options(options_loc)
    outdir = options["processed_charges_folder"]
    outdir_cpet = options["cpet_folder"]
    charges_directory = options["charges_folder"]

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

            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1

            filename = os.path.basename(i)
            listname = filename.split(".")
            center = [float(i) for i in nitrogen_dict["mean_N_xyz"]]

            if not fail_cond:
                ligand_identifier = ligand_dict["best_crit"].split(":")

                with open(output, "w") as outfile:
                    for j in readfile:
                        line_split = re.split(r"(\s+)", j)
                        cond = (
                            line_split[8] == ligand_identifier[0]
                            and line_split[10] == ligand_identifier[1]
                        )
                        
                        #x, y, z = float(j[30:38]), float(j[38:46]), float(j[46:54])
                        #box_conditional = x > center[0] + box_size or x < center[0] - box_size or y > center[1] + box_size or y < center[1] - box_size or z > center[2] + box_size or z < center[2] - box_size
                        
                        #if box_conditional:
                        #    temp_write = j[:56] + "0.000" + j[61:]
                        #    outfile.write(temp_write)

                        if zero_active and ("HETATM" in line_split[0] or cond):
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

                        else:
                            outfile.write(j)

                file_name = i.split("/")[-1].split(".")[0]
                if box:

                    options = open(f"{outdir_cpet}options_field_{file_name}.txt", "w+")
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

                    options = open(
                        f"{outdir_cpet}options_topol_{file_name}.txt", "w+"
                    )
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
                    options.write("    sampleOutput {} \n".format(file_name))
                    options.write("    bins {} \n".format(bins))
                    options.write(f"end \n")
                    options.close()

    print("fail count: {}".format(fail))
