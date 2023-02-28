import subprocess, os, random, argparse
from glob import glob

from HEML.utils.data import get_options


def main():
    # setup parser to get location of options file

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    options_loc = parser.parse_args().options
    options = get_options(options_loc)

    folder_processed = options["processed_pdb_folder"]
    out_folder = options["charges_folder"]
    chargefw2_loc = options["chargefw2_loc"]

    files = glob(folder_processed + "*.pdb")

    for ind in range(len(files)*10):
        i = random.choice(files)
        bool_exists = os.path.exists(out_folder + i.split("/")[-1] + ".pqr")
        print(out_folder + i.split("/")[-1] + ".pqr")
        print(bool_exists)
        if not bool_exists:
            result = subprocess.run(
                [
                    chargefw2_loc,
                    "--mode",
                    "charges",
                    "--input-file",
                    i,
                    "--chg-out-dir",
                    out_folder,
                    "--read-hetatm",
                    "TRUE",
                ]
            )
        files.remove(i)


main()
