import argparse
from HEML.utils.cpet import run_mag_calcs
from HEML.utils.data import get_options


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )

    parser.add_argument(
        "--num", help="number of magnitudes to calculate", default=10000
    )

    options_loc = parser.parse_args().options
    num = int(parser.parse_args().num)

    options = get_options(options_loc)

    cpet_folder = options["cpet_folder"]
    cpet_loc = options["cpet_loc"]
    processed_charges_folder = options["processed_charges_folder"]
    run_mag_calcs(cpet_loc, cpet_folder, processed_charges_folder)


main()
