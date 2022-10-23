from posixpath import split
from tkinter import E
import numpy as np
import os, re
from glob import glob
from HEML.utils.data_utils import *

if __name__ == "__main__" :
    zero_active = True
    zero_everything_charged = False
    
    box = False
    box_size = 3.0

    fail = 0
    outdir = "/ocean/projects/che160019p/santi92/processed_charges/"
    outdir_cpet = "/ocean/projects/che160019p/santi92/cpet/"
    charges_directory = "/ocean/projects/che160019p/santi92/heme_charges/*pqr" 

    #outdir = "../../data/charge_processed/"
    #outdir_cpet = "../../data/cpet/"
    #charges_directory = "../../data/charges/*pqr" 
    
    filelist = glob(charges_directory)
    #print(filelist)
    for i in filelist:


        # new filename
        filename = os.path.basename(i)
        listname = filename.split('.')

        #check that output isn't already there in the processed directory
        output = f'{outdir}{listname[0]}.pqr'
        if os.path.exists(output):
            print("output file already exists")
            pass

        #checks that input file is not empty
        elif(check_if_file_is_empty(i)): 
            print("input file is empty")
            pass

        else:        
            openfile = open(i)
            readfile = openfile.readlines()
            openfile.close
            filename = i.split('.')[0]

            # Look for FE in each line, and get chain:res:atom ID and xyz coords if exists.

            fail_cond = True

            ligand_dict = {}
            nitrogen_dict = {}

            try: 
                
                fe_id, fe_xyz = get_fe_positions(i)
                assert fe_id != None
                print(fe_id, fe_xyz)
                ligand_dict = get_ligand_info(i, fe_xyz)
                nitrogen_dict = get_N_positions(i, fe_id, fe_xyz)
                nitro_none = check_if_dict_has_None(nitrogen_dict)
                ligand_none = check_if_dict_has_None(ligand_dict)
                if(not nitro_none and not ligand_none):
                    fail_cond = False

            
                    if ligand_dict["best_crit_dist"] > 4.0:
                        print(ligand_dict["best_crit_dist"])
                        print(f'ERROR: No cysteine/tyrosine/histine ligand found for {i}.\n')
                        fail += 1
                        continue
                    else: 
                        print(f'Nitro 1 {nitrogen_dict["N_ID1"]}')
                        print(f'Nitro 2 {nitrogen_dict["N_ID2"]}')
                        print(f'Nitro 3 {nitrogen_dict["N_ID3"]}')
                        print(f'Nitro 4 {nitrogen_dict["N_ID4"]}')
                        print(f'ligand of note {ligand_dict["best_crit"]}\n')

            except: 
                fail_cond = True
                print("Failed File: ".format(i))
                fail += 1 




            filename = os.path.basename(i)
            listname = filename.split('.')

            if(not fail_cond):  
                ligand_identifier = ligand_dict["best_crit"].split(':')

                with open(output, 'w') as outfile:
                    for j in readfile:
                        line_split = re.split(r'(\s+)', j)
                        cond = ((line_split[8] == ligand_identifier[0] and line_split[10] == ligand_identifier[1]) )
                        
                        if (zero_active and ('HETATM' in line_split[0] or cond)):
                            temp_write = j[:56] + '0.000' + j[61:]
                            outfile.write(temp_write)
                        elif(zero_everything_charged and line_split[3] in ["ASP", "GLU", "LYS", "ARG", "HIS"]):
                            temp_write = j[:56] + '0.000' + j[61:]
                            outfile.write(temp_write)
                        else: 
                            outfile.write(j)

                file_name = i.split("charges")[-1][1:].split('.')[0]

                if(box):
                    density = 10            
                    options = open(f'{outdir_cpet}options_field_{file_name}.txt', 'w+')
                    options.write(f'align {nitrogen_dict["mean_N_xyz"][0]}:{nitrogen_dict["mean_N_xyz"][1]}:{nitrogen_dict["mean_N_xyz"][2]} {nitrogen_dict["N1_xyz"][0]}:{nitrogen_dict["N1_xyz"][1]}:{nitrogen_dict["N1_xyz"][2]} {nitrogen_dict["N2_xyz"][0]}:{nitrogen_dict["N2_xyz"][1]}:{nitrogen_dict["N2_xyz"][2]}\n')
                    options.write(f'%plot3d \n')
                    options.write(f'    show false \n')
                    options.write('    volume box {} {} {} \n'.format(box_size,box_size,box_size))
                    options.write('    density {} {} {} \n'.format(density, density, density))
                    options.write(f'output {outdir}efield_cox_{file_name}.dat \n')
                    options.write(f'end \n')                   
                    options.close()
                else: 
                    samples = 1000
                    bins = 20
                    step_size = 0.001
                    options = open(f'{outdir_cpet}options_topology_{file_name}.txt', 'w+')
                    options.write(f'align {nitrogen_dict["mean_N_xyz"][0]}:{nitrogen_dict["mean_N_xyz"][1]}:{nitrogen_dict["mean_N_xyz"][2]} {nitrogen_dict["N1_xyz"][0]}:{nitrogen_dict["N1_xyz"][1]}:{nitrogen_dict["N1_xyz"][2]} {nitrogen_dict["N2_xyz"][0]}:{nitrogen_dict["N2_xyz"][1]}:{nitrogen_dict["N2_xyz"][2]}\n')
                    options.write(f'%topology \n')
                    options.write('    volume box {} {} {} \n'.format(box_size,box_size,box_size))
                    options.write('    stepSize {} \n'.format(step_size))
                    options.write('    samples {} \n'.format(samples))
                    options.write('    sampleOutput {} \n'.format(file_name))
                    options.write('    bins {} \n'.format(bins))
                    options.write(f'end \n')                   
                    options.close()
        

    print("fail count: {}".format(fail))

