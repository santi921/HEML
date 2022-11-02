#get charges - preprocess
import subprocess, os, random
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

# hpc files
#folder_processed = "/ocean/projects/che160019p/santi92/pdbs_processed_heme/*"
#out_folder = "/ocean/projects/che160019p/santi92/heme_charges/"

# local
folder_processed = "../../../data/pdbs_processed/*"
out_folder = "../../../data/charges/"

files = glob(folder_processed)

for ind in range(len(files)):
    i = random.choice(files)
    bool_exists = os.path.exists(out_folder + i.split("/")[-1] + ".pqr")
    print(out_folder + i.split("/")[-1] + ".pqr")
    print(bool_exists)
    if(not bool_exists):
        result = subprocess.run([
        '/ocean/projects/che160019p/shared/ChargeFW2/bin/chargefw2', 
        '--mode', 
        'charges',
        '--input-file', 
        i,
        '--chg-out-dir',
        out_folder,
        '--read-hetatm',
        'TRUE' 
        ])
        


# rename all files ending in *pdb1 to *pdb 
# for i in glob(out_folder + "*pdb1"):
#     os.rename(i, i[:-1])
