import argparse
from HEML.utils.cpet import run_sweep_topology_calcs
from HEML.utils.data import get_options


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    parser.add_argument(
        "--num", help="number of topologies to calculate", default=10000
    )
    parser.add_argument("--threads", help="number of threads to use", default=16)
    options_loc = parser.parse_args().options
    num = int(parser.parse_args().num)
    threads = int(parser.parse_args().threads)
    options = get_options(options_loc)

    # cpet_folder = options["cpet_folder"]
    cpet_loc = options["cpet_loc"]
    outdir_sweep = options["sweep_folder"]
    # sweep_config = options["sweep_parameters"]
    # sweep_folders, config_list = sweep_config_to_folders_and_base_confs(sweep_config)
    processed_charges_folder = options["processed_charges_folder"]
    run_sweep_topology_calcs(
        cpet_loc, outdir_sweep, processed_charges_folder, num, threads
    )


main()
