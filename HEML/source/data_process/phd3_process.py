import os 
from HEML.utils.setup_phd3 import addh, write_dict_to_xyz, extract_heme_and_ligand_from_pdb
import numpy as np 



def main():
    
    # get subfolders in current directory
    for protein_name in os.listdir():
        if(os.path.isdir(protein_name)):
            print(protein_name)

            folder_name = "./" + protein_name
            #check if pdb is in folder
            if os.path.exists(os.path.join(folder_name, protein_name + ".pdb")):
                pdb_file = protein_name + ".pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_movie.pdb")):
                pdb_file = protein_name + "_movie.pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_todo_process.pdb")):
                pdb_file = protein_name + "_todo_process.pdb"

            # pull hemes and ligands from pdb - dict of xyz subset
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file)
            # write xyz to file - write xyz
            xyz_file_name = write_dict_to_xyz(folder_name, protein_name, dict_xyz)
            
            #convert xyz to pdb
            os.system("obabel -i xyz {} -o pdb -O ./{}/{}_heme.pdb". format(xyz_file_name, protein_name, protein_name))

    # add h to pdb 
    #addh("./{folder_name}/{folder_name}_heme.pdb")
    
    # convert pdb back to xyz
    #os.system("obabel -i pdb ./{folder_name}/{folder_name}_heme.pdb -o xyz -O ./{folder_name}/{folder_name}_heme.xyz")
    
    # convert xyz to coord 
    #os.system("./x2t ./{folder_name}/{folder_name}_heme.xyz > ./{folder_name}/coord")

    # run setupphd3.py 
    #os.system('cd {}'.format(folder_name))
    #os.system('setupphd3.py')
    #os.system('cd ..')

main()