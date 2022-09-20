#get charges - preprocess
import subprocess, os, random
import numpy as np
import matplotlib.pyplot as plt


files = os.listdir("../../data/pdbs_processed/")
for ind in range(len(files)):
    i = random.choice(files)
    bool_exists = os.path.exists("../../data/charges/" + i)
    if(not bool_exists):
        result = subprocess.run([
        '/ocean/projects/che160019p/shared/ChargeFW2/bin/chargefw2', 
        '--mode', 
        'charges',
        '--input-file', 
        '../../data/pdbs_processed/' + i,
        '--chg-out-dir',
        '../../data/charges/',
        '--read-hetatm',
        'TRUE' 
        ])
        
