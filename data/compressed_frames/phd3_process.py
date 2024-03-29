import os 
from HEML.utils.setup_phd3 import addh, write_dict_to_xyz, extract_heme_and_ligand_from_pdb, get_N_positions
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

            # write with =O
            # pull hemes and ligands from pdb - dict of xyz subset
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_oh = False, add_o = True)
            # write xyz to file - write xyz
            xyz_file_name = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_oh = False, add_o = True)
            #convert xyz to pdb
            os.system("obabel -i xyz {} -o pdb -O ./{}/{}_heme.pdb". format(xyz_file_name, protein_name, protein_name))


            # write with =OH
            # pull hemes and ligands from pdb - dict of xyz subset
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_o = False, add_oh = True)
            # write xyz to file - write xyz
            xyz_file_name = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_o = False, add_oh = True)
            #convert xyz to pdb
            os.system("obabel -i xyz {} -o pdb -O ./{}/{}_heme.pdb". format(xyz_file_name, protein_name, protein_name + "_oh"))

            # write vanilla
            # pull hemes and ligands from pdb - dict of xyz subset
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_o = False, add_oh = False)
            # write xyz to file - write xyz
            xyz_file_name = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_oh = False, add_o = False)            
            #convert xyz to pdb
            os.system("obabel -i xyz {} -o pdb -O ./{}/{}_heme.pdb". format(xyz_file_name, protein_name, protein_name + "_o"))

main()
