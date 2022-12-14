import os, random
from HEML.utils.data import *
from glob import glob

def create_folders():
    """Create a folder for each protein in the current directory and move the
        protein into the folder.    
    """
    # Get the current directory
    current_dir = os.getcwd()
    # Get the list of files in the current directory
    files = os.listdir(current_dir)
    # Loop over the files
    for file in files:
        if file[-3:] == "pdb":
            # Get the name of the protein
            protein = file.split(".")[0].split("_")[1]
            frame = file.split(".")[0].split("_")[0]
            # Create a folder for the protein
            os.mkdir(frame + "_" + protein)
            # Move the protein into the folder
            os.rename(file, frame + "_" + protein + "/" + file)


def move_charges_into_folder(charges_root, compressed_frame_folder):
    """Move the charges into the folder of the frame.
    Takes: 
        charges_root: the root directory of the charges
        compressed_frame_folder: the root directory of the compressed frames
    """
    #get list of frame number and protein name from compressed_frame_folder
    # get name of every folder in compressed_frame_folder
    folders = [f for f in os.listdir(compressed_frame_folder) if os.path.isdir(os.path.join(compressed_frame_folder, f))]
    for i in folders: 
        frame = i.split("_")[0]
        protein = i.split("_")[1]
        # check if the charges file exists
        if os.path.exists(os.path.join(charges_root, frame + "_" + protein + "_movie.pqr")):
            # move the charges file into the folder
            os.system("cp " + os.path.join(charges_root, frame + "_" + protein + "_movie.pqr") + " " + os.path.join(compressed_frame_folder, i))


def get_element_and_xyz(line):
    
    line_split = line.split()
    shift = 0 
    if(len(line_split[0]) > 6):
        shift = -1
    x = line_split[6 + shift]
    y = line_split[7 + shift]
    z = line_split[8 + shift]

    if(len(line_split[6]) > 8):
        x, y = break_up_line(x)
        z = line_split[7]
    if(len(line_split[7]) > 8):
        y, z = break_up_line(y)
    
    xyz = [x, y, z]
    xyz = [float(x) for x in xyz]
    xyz = np.array(xyz)
    element = line.split()[-1]

    return {"element":element, "xyz": xyz, "line": line}
    

def extract_heme_and_ligand_from_pdb(root, file): 
    """
    Extract the heme from the pdb files and save them in a new folder.
    Takes: 
        root: the root directory of the pdb files
    """
    file_folder = os.path.join(root,file)
    out_list = []
    fe_dict = get_fe_positions(file_folder)
    assert fe_dict["id"] != None
    ligand_dict = get_ligand_info(file_folder, fe_dict["xyz"])
    ligand_none = check_if_dict_has_None(ligand_dict)
    
    if(not ligand_none):
        fail_cond = False

        if ligand_dict["best_crit_dist"] > 4.0:
            print(ligand_dict["best_crit_dist"])
            print(f'ERROR: No cysteine/tyrosine/histine ligand found for {file_folder}.\n')
            fail += 1
            #continue
    
    #iterate lines of pdb file
    with open(file_folder, "r") as f:
        for line in f:
            
            if 'HETATM' in line.split()[0]:
                
                shift = 0 
                heme_cond = line[17:20] == "HEM"
                heme_chain_cond = line[21] == fe_dict["id"].split(":")[0]
                heme_id_cond = line[22:26].strip() == fe_dict["id"].split(":")[1].strip()
                hetero_cond = 'HETATM' in line.split()[0]

                if(heme_cond and heme_id_cond and heme_chain_cond and hetero_cond):
                    out_list.append(get_element_and_xyz(line))
            
            line_split = line.split()
            
            if(len(line_split) > 3):
                sg_cond = 'CYS' in line_split[3]
                oh_cond = 'TYR' in line_split[3]
                nend_cond = 'HIS' in line_split[3]

                if (sg_cond or oh_cond or nend_cond):
                    
                    ligand_chain_cond = line[21] == ligand_dict["best_crit"].split(":")[0]
                    ligand_id_cond = line[22:26].strip() == ligand_dict["best_crit"].split(":")[1].strip()
                    anisou_cond = 'ANISOU' in line_split[0]
                    if(ligand_chain_cond and ligand_id_cond and not anisou_cond):
                        out_list.append(get_element_and_xyz(line))

    return out_list            
                    

def addh(pdb_file): 
    """
    Add hydrogens to the pdb files and overwrite the olds files.
    """
    from chimera import runCommand as rc
    from os import system as run
    
    print(pdb_file)
    print("open " + pdb_file)

    rc("open " + pdb_file)
    rc("delete solvent")
    rc("addh")
    rc("write #0 " + pdb_file)
    print("write #0 " + pdb_file)
    rc("close session")


def write_dict_to_xyz(folder, name, dict_xyz):
    # write the xyz file
    xyz_file_name = os.path.join(folder, name)
    xyz_file_name += ".xyz"

    with open(xyz_file_name, "w") as f:
        # iterate over dictionary 
        f.write(str(len(dict_xyz)) + "\n\n")

        for i in dict_xyz: 
            f.write("{} {}  {}  {}\n".format(i["element"], i["xyz"][0], i["xyz"][1], i["xyz"][2]))
            
    return xyz_file_name