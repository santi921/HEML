import os, re, argparse
from glob import glob
from HEML.utils.data import *
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
    parser.add_argument(
        "--box", help="box", default=False
    )
    parser.add_argument(
        "--box_size", help="box size", default=4.0
    )
    options_loc = parser.parse_args().options
    zero_active = parser.parse_args().zero_active
    zero_everything_charged = parser.parse_args().zero_everything_charged
    box = bool(parser.parse_args().box)
    box_size = parser.parse_args().box_size

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

            try:
                fe_dict = get_fe_positions(i)
                assert fe_dict["id"] != None
                print(fe_dict["id"], fe_dict["xyz"])
                fail_cond = False
                print(fe_dict["id"], fe_dict["xyz"])
            except:
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1


            filename = os.path.basename(i)
            listname = filename.split(".")

            if not fail_cond:
                
                with open(output, "w") as outfile:
                    for j in readfile:

                        line_split = re.split(r"(\s+)", j)

                        if zero_active and ("HETATM" in line_split[0]):
                            temp_write = j[:56] + "0.000" + j[61:]
                            outfile.write(temp_write)
                        else:
                            outfile.write(j)

                file_name = i.split("/")[-1].split(".")[0]
                options = open(f"{outdir_cpet}options_mag_{file_name}.txt", "w+")
                options.write(f"%field \n")
                options.write(f"    locations {fe_dict[0]}:{fe_dict[1]}:{fe_dict[2]} \n")
                options.write(f"output {outdir_cpet}efield_mag_FE_{file_name}.dat \n")
                options.write(f"end \n")
                options.close()


    print("fail count: {}".format(fail))
