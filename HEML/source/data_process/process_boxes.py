# get charges - preprocess
import subprocess, os, random
import numpy as np
import matplotlib.pyplot as plt


def main():
    files = os.listdir("../../data/charge_processed/")
    for ind in range(len(files)):
        i = random.choice(files).split("_")[0][:-4]
        print(i)
        bool_exists = os.path.exists("./cpet/efield_cox_" + i + ".dat")
        if not bool_exists:
            result = subprocess.run(
                [
                    "/ocean/projects/che160019p/shared/CPET/bin/cpet",
                    "-o",
                    "./cpet/options_" + i + ".txt",
                    "-p",
                    "./charge_processed/" + i + ".pqr",
                    "-v",
                ]
            )

        print(
            "/ocean/projects/che160019p/shared/CPET/bin/cpet",
            " -o "
            + "./cpet/options_"
            + i
            + ".txt"
            + " -p "
            + "./charge_processed/"
            + i
            + ".pqr",
        )


main()

#'/ocean/projects/che160019p/shared/CPET/bin/cpet' -o ./cpet/options_2wm41.txt -p ./charge_processed/2wm41.pqr
