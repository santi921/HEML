from posixpath import split
import numpy as np
import glob, os, re
    
# zero heteroatom / iron 

def break_up_line(str_process):
    split_str = str_process.split('-')
    if (split_str[0] == ''):
        ret_1 = "-" + split_str[1]
    else: 
        ret_1 = split_str[0]
    return ret_1, "-" + split_str[-1]


fail = 0
filelist = glob.glob('./charges/*.pqr')

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
        if 'HETATM' in line[0] and ('NB' in line[2+shift] and 'HEM' or 'HEC' in line[3+shift]) and line[4+shift] == fe_ID.split(":")[0]:
            N_ID2 = f'{line[4+shift]}:{line[5+shift]}:{line[2+shift]}' 
        
        
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
    print(f'ligand of note {best_crit}\n')

    # get info for each
    ligand_identifier = best_crit.split(':')
    
    # new filename
    filename = os.path.basename(i)
    listname = filename.split('.')
    output = f'./charge_processed/{listname[0]}.pqr'
    
    with open(output, 'w') as outfile:
        for j in readfile:
            line_split = re.split(r'(\s+)', j)
            cond = ((line_split[8] == ligand_identifier[0] and line_split[10] == ligand_identifier[1]) )
            if ('HETATM' in line_split[0] or cond):
                #print(j[:55])
                temp_write = j[:56] + '0.000' + j[61:]
                
                outfile.write(temp_write)
            else: 
                outfile.write(j)


    file_name = i.split("charges")[-1][1:].split('.')[0]
    options = open(f'./cpet/options_{file_name}.txt', 'w+')
    options.write(f'align {fe_ID} {best_crit} {N_ID}\n')
    options.write(f'%plot3d \n')
    options.write(f'    show false \n')
    options.write(f'    volume box 3.0 3.0 3.0 \n')
    options.write(f'    density 4 4 4 \n')
    options.write(f'output ./cpet/efield_cox_{file_name}.dat \n')
    options.write(f'end \n')                   
    options.close()


