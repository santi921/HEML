import os 
from HEML.utils.setup_phd3 import extract_heme_and_ligand_from_pdb, xtb_sanitize_and_save
from HEML.utils.data import get_options
import numpy as np 



def main():
    options = get_options("./options.json")
    root = options["compressed_proteins_folder"]

    # get subfolders in current directory
    for protein_name in os.listdir(root):
        if(os.path.isdir(os.path.join(root,protein_name))):
            print(protein_name)
            folder_name = os.path.join(root, protein_name)
            #check if pdb is in folder
            if os.path.exists(os.path.join(folder_name, protein_name + ".pdb")):
                pdb_file = protein_name + ".pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_movie.pdb")):
                pdb_file = protein_name + "_movie.pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_todo_process.pdb")):
                pdb_file = protein_name + "_todo_process.pdb"

            # write with =O
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_oh = False, add_o = True, freeze=True)
            xyz_file_name_1 = xtb_sanitize_and_save(folder_name, protein_name,dict_xyz, add_oh = False, add_o = True)
            # write xyz to file - write xyz
            #xyz_file_name_1 = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_oh = False, add_o = True)
            
            # write with -OH
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_o = False, add_oh = True, freeze = False)
            xyz_file_name_2 = xtb_sanitize_and_save(folder_name, protein_name,dict_xyz, add_o = False, add_oh = True)
            # write xyz to file - write xyz
            #xyz_file_name_2 = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_o = False, add_oh = True)            

            # write vanilla =O
            dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file, add_o = True, add_oh = False, freeze = False)
            xyz_file_name_3 = xtb_sanitize_and_save(folder_name, protein_name,dict_xyz, add_o = True, add_oh = False)
            # write xyz to file - write xyz
            #xyz_file_name_3 = write_dict_to_xyz(folder_name, protein_name, dict_xyz, add_o = True, add_oh = False)            
            
            #convert xyz to pdb
            os.system("obabel -i xyz {} -o pdb -O {}/{}_heme.pdb". format(os.path.join(folder_name, xyz_file_name_1), folder_name, protein_name))
            os.system("obabel -i xyz {} -o pdb -O {}/{}_heme.pdb". format(os.path.join(folder_name, xyz_file_name_2), folder_name, protein_name+ "_oh"))
            os.system("obabel -i xyz {} -o pdb -O {}/{}_heme.pdb". format(os.path.join(folder_name, xyz_file_name_3), folder_name, protein_name+ "_o"))

main()