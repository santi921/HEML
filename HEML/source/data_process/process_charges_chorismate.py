import os, re, argparse
import numpy as np
from glob import glob
from HEML.utils.data import (
    get_options,
    check_if_file_is_empty,
)
from HEML.utils.mol2topqr import mol2_to_pqr_folder

def get_CHR_positions(i,atom_list):
    with open(i, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if atom_list[0] in line[2] and "CHR" in line[3]:
            C1_xyz = [line[6], line[7], line[8]]
            C1_xyz = [float(x) for x in C1_xyz]
            C1_xyz = np.array(C1_xyz)
            break

    with open(i, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if atom_list[1] in line[2] and "CHR" in line[3]:
            C3_xyz = [line[6], line[7], line[8]]
            C3_xyz = [float(x) for x in C1_xyz]
            C3_xyz = np.array(C1_xyz)
            break

    with open(i, "r") as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if atom_list[2] in line[2] and "CHR" in line[3]:
            C5_xyz = [line[6], line[7], line[8]]
            C5_xyz = [float(x) for x in C1_xyz]
            C5_xyz = np.array(C1_xyz)
            break
    return C1_xyz, C3_xyz, C5_xyz

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    # store t/f if --box is used
    parser.add_argument("--box", help="box", action="store_true")
    parser.add_argument("--box_size", help="box size", default=2.0)
    parser.add_argument("--density", help="density", default=10)
    parser.add_argument("--samples", help="samples", default=10000)
    parser.add_argument("--bins", help="bins", default=25)
    parser.add_argument("--step_size", help="step size", default=0.001)

    options_loc = parser.parse_args().options
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

            fail_cond = True

            ligand_dict = {}
            nitrogen_dict = {}

            try:
                #Extract C1, C3, and C5 of CHR residue
                atom_list=["C1","C3","C5"]
                C1_xyz, C3_xyz, C5_xyz = get_CHR_positions(i,atom_list)
                assert C1_xyz != None
                assert C3_xyz != None
                assert C5_xyz != None
                print(C1_xyz,C3_xyz,C5_xyz)
                C1_C5_mean=np.mean([C1_xyz,C5_xyz],axis=0)

            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1

            filename = os.path.basename(i)
            listname = filename.split(".")

            if not fail_cond:

                with open(output, "w") as outfile:
                    for line in readfile:
                        line_split = re.split(r"(\s+)", line)
                        # remove items in list that are empty strings
                        # line_split = [x for x in line_split if x != ""]
                        # strip whitespace from list
                        line_split = [x.strip() for x in line_split]
                        # remove empty strings from list
                        line_split = list(filter(None, line_split))
                        # print(line_split)
                        
                        # Check if the line represents an atom
                        if len(line_split) > 0 and line_split[0] == "ATOM":
                            if line_split[3] == "CHR":
                                if not first_chr_found:
                                # Zero the charge and set flags
                                    temp_write = line[:56] + "0.000" + line[61:]
                                    outfile.write(temp_write)
                                    processing_chr = True
                                elif processing_chr:
                                    # Continue processing the rest of CHR atoms
                                    temp_write = line[:56] + "0.000" + line[61:]
                                    outfile.write(temp_write)
                                else:
                                    # Write the line as is for CHR atoms after the first instance
                                    outfile.write(line)
                            else:
                                if processing_chr:
                                    # We've reached the end of the first CHR instance
                                    first_chr_found = True
                                    processing_chr = False
                                    outfile.write(line)
                        else:
                            # Write non-atom lines as is
                            outfile.write(line)


                file_name = i.split("/")[-1].split(".")[0]
                
                if box:
                    options = open(f"{outdir_cpet}options_field_{file_name}.txt", "w+")

                    options.write(
                        f"align {C1_C5_mean[0]}:{C1_C5_mean[1]}:{C1_C5_mean[2]} {C1_xyz[0]}:{C1_xyz[1]}:{C1_xyz[2]} {C3_xyz[0]}:{C3_xyz[1]}:{C3_xyz[2]}\n"
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
                    options.write(f"output {outdir_cpet}efield_box_{file_name}.dat \n")
                    options.write(f"end \n")
                    options.close()

                else:
                    options = open(f"{outdir_cpet}options_topol_{file_name}.txt", "w+")
                    options.write(
                        f"align {C1_C5_mean[0]}:{C1_C5_mean[1]}:{C1_C5_mean[2]} {C1_xyz[0]}:{C1_xyz[1]}:{C1_xyz[2]} {C3_xyz[0]}:{C3_xyz[1]}:{C3_xyz[2]}\n"
                    )
                    options.write(f"%topology \n")
                    options.write(
                        "    volume box {} {} {} \n".format(
                            box_size, box_size, box_size
                        )
                    )
                    options.write("    stepSize {} \n".format(step_size))
                    options.write("    samples {} \n".format(samples))
                    options.write(f"    sampleOutput efield_topo_{file_name} \n \n")
                    options.write("    bins {} \n".format(bins))
                    options.write(f"end \n")
                    options.close()

    print("fail count: {}".format(fail))
