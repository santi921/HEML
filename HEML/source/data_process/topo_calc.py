import os, re, argparse
import numpy as np
from glob import glob
from HEML.utils.data import get_options, get_fe_positions, get_N_positions
from HEML.utils.mol2topqr import mol2_to_pqr_folder
from CPET.source.topo_calc import Topo_calc
from CPET.utils.parser import parse_pqr, filter_pqr_radius, filter_pqr_residue


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    parser.add_argument("--concur_slip", help="concurrent slip", default=8)

    options_loc = parser.parse_args().options
    concur_slip = int(parser.parse_args().concur_slip)
    options = get_options(options_loc)
    options["concur_slip"] = int(concur_slip)
    outdir_cpet = options["cpet_folder"]
    charges_directory = options["charges_folder"]

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

    print(filelist)
    for i in range(len(filelist)):
        ind_select = np.random.randint(0, len(filelist))
        file_pqr = filelist[ind_select]
        file_pqr_short = file_pqr.split("/")[-1].split(".")[0]

        if not os.path.isfile(outdir_cpet + file_pqr_short + ".txt"):
            # try:
            options["path_to_pqr"] = file_pqr
            fe_dict = get_fe_positions(file_pqr)
            n_dict = get_N_positions(
                file_pqr, fe_ID=fe_dict["id"], fe_xyz=fe_dict["xyz"]
            )

            options["center"] = list(fe_dict["xyz"])
            options["x"] = list(n_dict["N1_xyz"])
            options["y"] = list(n_dict["N2_xyz"])
            topo = Topo_calc(options)
            print("Computing topo for file: ", file_pqr_short)
            hist = topo.compute_topo()
            np.savetxt(outdir_cpet + file_pqr_short + ".txt", hist)
            # except:
            #    print("Failed on file: ", file_pqr_short)
            #    fail += 1
        else:
            print("File already exists, skipping")
