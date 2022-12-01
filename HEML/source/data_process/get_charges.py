import subprocess, os, random, json
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from HEML.utils.data import get_options


def main():

    options = get_options("./options.json")
    folder_processed = options["processed_pdb_folder"]
    out_folder = options["charges_folder"]
    chargefw2_loc = options["chargefw2_loc"]

    files = glob(folder_processed)

    for ind in range(len(files)):
        i = random.choice(files)
        bool_exists = os.path.exists(out_folder + i.split("/")[-1] + ".pqr")
        print(out_folder + i.split("/")[-1] + ".pqr")
        print(bool_exists)
        if(not bool_exists):
            result = subprocess.run([
            chargefw2_loc, 
            '--mode', 
            'charges',
            '--input-file', 
            i,
            '--chg-out-dir',
            out_folder,
            '--read-hetatm',
            'TRUE' 
            ])
            
main()

