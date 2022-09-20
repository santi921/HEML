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
                N1_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
            if 'HETATM' in line[0] and ('NB' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID2 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                N2_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
            if 'HETATM' in line[0] and ('NC' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID3 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
                N3_xyz = [float(line[-5]), float(line[-4]), float(line[-3])]
            if 'HETATM' in line[0] and ('ND' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
                N_ID4 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
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

        # get info for each
        ligand_identifier = best_crit.split(':')
        
        # new filename
        filename = os.path.basename(i)
        listname = filename.split('.')
        output = f'../../data/charge_processed/{listname[0]}.pqr'
        
        with open(output, 'w') as outfile:
            for j in readfile:
                line_split = re.split(r'(\s+)', j)
                cond = ((line_split[8] == ligand_identifier[0] and line_split[10] == ligand_identifier[1]) )
                
                if (zero_active and ('HETATM' in line_split[0] or cond)):
                    temp_write = j[:56] + '0.000' + j[61:]
                    outfile.write(temp_write)
                if(zero_everything_charged and line_split[3] in ["ASP", "GLU", "LYS", "ARG", "HIS"]):
                    temp_write = j[:56] + '0.000' + j[61:]
                    outfile.write(temp_write)
                else: 
                    outfile.write(j)


        file_name = i.split("charges")[-1][1:].split('.')[0]
        options = open(f'../../data/cpet/options_{file_name}.txt', 'w+')
        options.write(f'align {mean_N_xyz} {N_ID} {N_ID2}\n')
        options.write(f'%plot3d \n')
        options.write(f'    show false \n')
        options.write(f'    volume box 3.0 3.0 3.0 \n')
        options.write(f'    density 10 10 10 \n')
        options.write(f'output ../../data/cpet/efield_cox_{file_name}.dat \n')
        options.write(f'end \n')                   
        options.close()

    print("fail count: {}".format(fail))

