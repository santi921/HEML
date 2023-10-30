from HEML.utils.cpet import run_box_calcs
from HEML.utils.data import get_options
import os, argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    options_loc = parser.parse_args().options
    options = get_options(options_loc)

    cpet_folder = options["cpet_folder"]
    processed_charges_folder = options["processed_charges_folder"]
    cpet_loc = options["cpet_loc"]
    print(os.listdir(cpet_folder))
    print(os.listdir(processed_charges_folder))
    run_box_calcs(cpet_loc, cpet_folder, processed_charges_folder)


main()
