from posixpath import split
from tkinter import E
import numpy as np
import glob, os, re
    
# zero heteroatom / iron 
# other options include zeroing non active site residues or all charged residues

def break_up_line(str_process):
    split_str = str_process.split('-')
    if (split_str[0] == ''):
        ret_1 = "-" + split_str[1]
    else: 
        ret_1 = split_str[0]
    return ret_1, "-" + split_str[-1]

def spacefinder(List_String):
    try:
        len(List_String) == 11
    except:
        print("Does your pdb contain charges? List of strings should have 11 values.")

    slen1 = len(List_String[1])
    slen2 = len(List_String[2])
    slen5 = len(List_String[5])
    slen6 = len(List_String[6])
    slen7 = len(List_String[7])
    slen8 = len(List_String[8])
    slen9 = len(List_String[9])
    slen10 = len(List_String[10])

    if List_String[0] == "HETATM":
        backlen1 = 4 - (slen1-1)
    else:
        backlen1 = 6 - (slen1-1)
    if slen2 > 3:
        backlen2 = 2 - (slen2-3)
        backlen3 = 3 - (slen2-2)
    else:
        backlen2 = 2
        backlen3 = 3 - (slen2-1)
    backlen5 = 3 - (slen5-1)
    backlen6 = 11 - (slen6-1)
    backlen7 = 7 - (slen7-1)
    backlen8 = 7 - (slen8-1)
    backlen9 = 6 - (slen9-1)
    backlen10 = 6 - (slen10-1)

    outstring = List_String[0] + (" " * backlen1) + List_String[1] + (" " * backlen2) + List_String[2] + (" " * backlen3) \
        + List_String[3] + " " + List_String[4] + (" " * backlen5) + List_String[5] + (" " * backlen6) + List_String[6] \
             + (" " * backlen7) + List_String[7]  + (" " * backlen8) + List_String[8]  + (" " * backlen9) + List_String[9] \
                  + (" " * backlen10) + List_String[10]
    return outstring


if __name__ == "__main__" :
    zero_active = True
    zero_everything_charged = False
    
    n1_xyz_list, n2_xyz_list, n3_xyz_list, n4_xyz_list = [], [], [], []
    crit_xyz_list = []

    fail = 0
    filelist = glob.glob('../../data/charges/*.pqr')

    for i in filelist:
        print(i)
        #output.write(f'Opening {i}...\n')
        openfile = open(i)
        readfile = openfile.readlines()
        openfile.close
        best_crit = ''
        best_N = ''
        fe_crit_dist = 10.0
        fe_N_dist = 10.0
        best_crit_dist = 10.0
        best_N_dist = 10.0
        filename = i.split('.')[0]

        # Look for FE in each line, and get chain:res:atom ID and xyz coords if exists.

        for j in readfile:
            line = j.split()
            if 'HETATM' in line[0] and ('FE' in line[2] or "FE" in line[1]):
                shift = 0 
                if("FE" in line[1]):
                    shift = -1
                fe_ID = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}'
                fe_xyz = [line[6+shift], line[7+shift], line[8+shift]]
                fe_xyz = [float(x) for x in fe_xyz]
                fe_xyz = np.array(fe_xyz)
                break

                

        # Look for CYS ligand in each line, get ID and xyz coordinates, find nearest sulfur to FE.
        best_crit_dist = 100

        for j in readfile:
            line = j.split()
            
            sg_cond = 'ATOM' in line[0] and 'SG' in line[2] and 'CYS' in line[3]
            oh_cond = 'ATOM' in line[0] and 'OH' in line[2] and 'TYR' in line[3]
            nend_cond = 'ATOM' in line[0] and (('NE2' in line[2]) or ("ND1") in line[2]) and 'HIS' in line[3]

            shift = 0 
            if(len(line[0]) > 6):
                shift = -1

            if 'HETATM' in line[0] and ('NA' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                try:
                    N1_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                except: 
                    N1_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]

            if 'HETATM' in line[0] and ('NB' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID2 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                try:
                    N2_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                except: 
                    N2_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]

            if 'HETATM' in line[0] and ('NC' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID3 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                try:
                    N3_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                except: 
                    N3_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
            
            if 'HETATM' in line[0] and ('ND' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID4 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                try:
                    N4_xyz = [float(j[31:38]), float(j[38:45]), float(j[46:54])]
                except: 
                    N4_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
            
            
                mean_N_xyz =np.mean(np.array([N1_xyz, N2_xyz, N3_xyz, N4_xyz]), axis=0)
            
            if (sg_cond or oh_cond or nend_cond):
                crit_ID = f'{line[4]}:{line[5]}:{line[2]}'
                # catch any clumped elements
                x = line[6]
                y = line[7]
                z = line[8]

                if(len(line[6]) > 8):
                    x, y = break_up_line(x)
                    z = line[7]

                if(len(line[7]) > 8):
                    y, z = break_up_line(y)

                crit_xyz = [x, y, z]
                crit_xyz = [float(x) for x in crit_xyz]
                crit_xyz = np.array(crit_xyz)
                crit_dist = np.linalg.norm(fe_xyz - crit_xyz)

                if (crit_dist < fe_crit_dist and crit_dist < best_crit_dist):
                    best_crit_dist = crit_dist
                    best_crit = crit_ID
                    
                    
        if best_crit_dist > 4.0:
            print(crit_dist)
            print(best_crit_dist)
            print(f'ERROR: No cysteine/tyrosine/histine ligand found for {i}.\n')
            fail += 1
            continue
        
        #print(f'Cysteine ligand distance from {best_crit} is {best_crit_dist} angstrom.')
        print(f'Nitro 1 {N_ID}')
        print(f'Nitro 2 {N_ID2}')
        print(f'Nitro 3 {N_ID3}')
        print(f'Nitro 4 {N_ID4}')
        print(f'ligand of note {best_crit}\n')

        n1_xyz_list.append(N1_xyz - fe_xyz)
        n2_xyz_list.append(N2_xyz - fe_xyz)
        n3_xyz_list.append(N3_xyz - fe_xyz)
        n4_xyz_list.append(N4_xyz - fe_xyz)
        crit_xyz_list.append(crit_xyz - fe_xyz)
        
        # get info for each
        ligand_identifier = best_crit.split(':')
        


    print(np.mean(np.array(n1_xyz_list), axis=0))
    print(np.mean(np.array(n2_xyz_list ), axis=0))
    print(np.mean(np.array(n3_xyz_list) , axis=0))
    print(np.mean(np.array(n4_xyz_list), axis=0))
    print(np.mean(np.array(crit_xyz_list), axis=0))
    print("fail count: {}".format(fail))

