import pandas as pd 
import numpy as np 




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


def pull_mats(dir_dat):
    from glob import glob 
    dir_data = glob("./cpet/*dat")
    mat_arr = []
    for i in dir_data:
        mat_arr.append(mat_pull(i))


pull_mats('./cpet')

