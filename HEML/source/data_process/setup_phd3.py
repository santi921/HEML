import os 
from HEML.utils.setup_phd3 import write_dict_to_xyz, extract_heme_and_ligand_from_pdb

def main():
    
    folder_name = "./108_3hb6"
    pdb_file = "108_3hb6.pdb"
    
    # pull hemes and ligands from pdb - dict of xyz subset
    dict_xyz = extract_heme_and_ligand_from_pdb(folder_name, pdb_file)
    
    # write xyz to file - write xyz
    xyz_file_name = write_dict_to_xyz(folder_name, folder_name[2:], dict_xyz)
    
    #convert xyz to pdb
    os.system("obabel -i xyz ./{folder_name}/{xyz_file_name} -o pdb -O ./{folder_name}/{folder_name}_heme.pdb")

    # add h to pdb 
    addh("./{folder_name}/{folder_name}_heme.pdb")
    
    # convert pdb back to xyz
    os.system("obabel -i pdb ./{folder_name}/{folder_name}_heme.pdb -o xyz -O ./{folder_name}/{folder_name}_heme.xyz")
    
    # convert xyz to coord 
    os.system("./x2t ./{folder_name}/{folder_name}_heme.xyz > ./{folder_name}/coord")

    # run setupphd3.py 
    os.system('cd {}'.format(folder_name))
    os.system('setupphd3.py')
    os.system('cd ..')

