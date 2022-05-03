import pandas as pd 
import numpy as np 
from glob import glob 
import os 



# ./cpet/efield_cox_5o4k.dat
def mat_pull(file):
    with open(file) as f: 
        lines = f.readlines()

    steps_x = 2 * int(lines[0].split()[2]) 
    steps_y = 2 * int(lines[0].split()[3]) 
    steps_z = 2 * int(lines[0].split()[4][:-1]) 
    gap_x = round(np.abs(float(lines[steps_x*steps_y + 7].split()[0]) - float(lines[7].split()[0])), 4)
    gap_y = round(np.abs(float(lines[steps_x+8].split()[1]) - float(lines[7].split()[1])), 4)
    gap_z = round(np.abs(float(lines[8].split()[2]) - float(lines[7].split()[2])), 4)
    mat = np.zeros((steps_x, steps_y, steps_z, 3))
    
    for ind, i in enumerate(lines[7:]):
        line_split = i.split()
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 0] = float(line_split[-3])
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 1] = float(line_split[-2])
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 2] = float(line_split[-1])
    return mat  


def pull_mats_w_label(dir_dat):

    x, y = [], []
    df = pd.read_csv("protein_data.csv")
    y_count, h_count, c_count = 0, 0, 0
    for row in df.iterrows():
        cpet_name = "./cpet/efield_cox_" + row[1]['name'] + ".dat"
        if(os.path.exists(cpet_name)):
            
            x.append(mat_pull(cpet_name))

            if(row[1]['label'] == 'Y'):
                y.append([1,0,0])
                y_count += 1
            elif(row[1]['label'] == 'H'):
                y.append([0,1,0])
                h_count += 1
            else:
                y.append([0,0,1])
                c_count += 1
    print(y_count, h_count, c_count)
    return np.array(x), np.array(y)

