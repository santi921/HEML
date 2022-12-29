import os, json, ast
from chimera import runCommand as rc
from os import system as run
import numpy as np 

def addh(pdb_file): 
    """
    Add hydrogens to the pdb files and overwrite the olds files.
    """
    print(pdb_file)
    print("open " + pdb_file)

    rc("open " + pdb_file)
    rc("delete solvent")
    rc("addh")
    rc("write #0 " + pdb_file[:-4] + "_h.pdb")
    print("write #0 " + pdb_file[:-4] + "_h.pdb")
    rc("close session")


def create_folders(folder_name):
    """
    Creates the folders for the turbomole calculations.
    Takes:
        folder_name: the folder where the folders should be created
    """

    if not os.path.exists("{}/no_charges".format(folder_name)):
        os.makedirs("{}/no_charges".format(folder_name))
    if not os.path.exists("{}/embedding".format(folder_name)):
        os.makedirs("{}/embedding".format(folder_name))

    if not os.path.exists("{}/embedding/o".format(folder_name)):
        os.makedirs("{}/embedding/o".format(folder_name))
    if not os.path.exists("{}/embedding/oh".format(folder_name)):
        os.makedirs("{}/embedding/oh".format(folder_name))
    if not os.path.exists("{}/embedding/normal".format(folder_name)):
        os.makedirs("{}/embedding/normal".format(folder_name))

    if not os.path.exists("{}/no_charges/o".format(folder_name)):
        os.makedirs("{}/no_charges/o".format(folder_name))
    if not os.path.exists("{}/no_charges/oh".format(folder_name)):
        os.makedirs("{}/no_charges/oh".format(folder_name))
    if not os.path.exists("{}/no_charges/normal".format(folder_name)):
        os.makedirs("{}/no_charges/normal".format(folder_name))


def setup_turbomole(folder_name):
    """
    Sets up the turbomole calculation.
    Takes:
        folder: the folder where the calculation should be set up
    """
    os.chdir(folder_name)
    os.system('setupturbomole.py -t')
    os.chdir("../../..")


def write_json(folder, frozen_atoms = [], atoms_present = [], charge = 0): 
    """
    Writes a json file for the turbomole calculation. 
    Takes: 
        folder: the folder where the json file should be written
        frozen_atoms: a list of atoms that should be frozen in the calculation
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    basic_dict = {
        "geometry": { 
            "cartesians": True,
            "idef": { "idef_on": False },
            "freeze_stretch": ["4,5"],
            "ired": False,
            "iaut": { 
                "iaut_on": False, 
                "bonds": ["4,5"]    
            }
        },
        "dft": {
            "dft_on": True,
            "func": "tpss",
            "grid": "m4"           
        },
        "scf": {
            "iter": 300,
            "conv": 5
        },
        "basis": { "all": "def2-SVP" },
        "stp": {
            "itvc": 0,
            "trad": 0.1
        },
        "open_shell": {
            "open_shell_on": False,
        },
        "cosmo": 4,
        "freeze_atoms": [],
        "calculation": "geo",
        "geo_iterations": 200,
        "weight": False,
        "gcart": None,
        "denconv": None,
        "rij": True,
        "marij": True,
        "dsp": True,
        "charge": 0
    }

    if atoms_present == []:
        basic_dict["basis"] = {
            "all": "def2-SVP",
            "fe": "def2-TZVP",
            "n": "def2-TZVP",
            "s": "def2-TZVP",
            "o": "def2-TZVP"
        }
    else:
        basic_dict["basis"]["all"] = "def2-SVP"
        for atom in atoms_present:
            if atom == "Fe" or atom == "fe":
                basic_dict["basis"]["fe"] = "def2-TZVP"
            elif atom == "N" or atom == "n":
                basic_dict["basis"]["n"] = "def2-TZVP"
            elif atom == "S" or atom == "s":
                basic_dict["basis"]["s"] = "def2-TZVP"
            elif atom == "O" or atom == "o":
                basic_dict["basis"]["o"] = "def2-TZVP"
            else: 
                basic_dict["basis"][atom.lower()] = "def2-SVP"
    if charge != 0:
        basic_dict["charge"] = charge
        if charge == -2: 
            basic_dict["open_shell"]["open_shell_on"] = True
            basic_dict["open_shell"]["unpaired"] = 1
        
    if frozen_atoms != []:
        basic_dict['freeze_atoms'] = frozen_atoms
    
    
    with open(folder + "definput.json", "w") as outfile:
        json.dump(basic_dict, outfile, indent = 4)


def get_options(options_file = "./options.json"):
    """
    Get options from options.json file and create folders if they don't exist.
    Takes 
        options_file: path to options.json file
    Returns
        options: dictionary of options
    """
    with open(options_file) as f:
        options = json.load(f)
    for key in options:
        if "folder" in key:
            if not os.path.exists(options[key]):
                os.makedirs(options[key])
    
    return options


def fetch_charges_dict(file_name = 'test.pqr'):
    """
    Given a list of dictionaries with element and position, traverse a pqr file and get the charges from the file
    EXCLUDING ELEMENTS IN THE LIST OF DICTIONARIES
    Takes: 
        list of dictionaries with element and position
    Returns:
        list of dictionaries with element, position and charge
    """
    
    pqr_dict = []
    # get the lines of the pqr file
    with open(file_name, 'r') as f:
        lines = f.readlines()
    
    for line in lines: 
        x       = float(line[30:38].strip())
        y       = float(line[39:46].strip())
        z       = float(line[47:54].strip())
        charge  = float(line[55:61].strip())
        radius  = float(line[62:68].strip())    
        if np.abs(charge) >= 0.01:
            pqr_dict.append({"position": [x,y,z], "charge": charge, "radius": radius})

    return pqr_dict


def get_elements(file_name): 
    """
    Open xyz and get all the elements in the file
    Takes:
        file_name: the name of the xyz file
    Returns:
        elements: a list of elements
    """
    elements = []
    with open(file_name, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0:2].strip().isalpha():
            elements.append(line.split()[0])
    return elements


def put_charges_in_turbo_files(folder_name, charges_dict): 
    """
    Traverses subdirectories in embedding folder and puts charges in the turbomole file
    
    Takes   
        folder_name: the folder where the turbomole files are located
        charges_dict: a dictionary with the charges
    Returns: Nothing 
    
    """
    # find folder named embedding and go into all subfolders
    for root, dirs, files in os.walk(folder_name):
        for file in files:
            if file.endswith("control"):
                print("editing control file with charges from pqr dictionary")
                # remove last line of file - the $end
                with open(os.path.join(root, file), 'r') as f:
                    lines = f.readlines()
                with open(os.path.join(root, file), 'w') as f:
                    for line in lines[:-1]:
                        f.write(line)

                # append dictionary to end of file
                with open(os.path.join(root, file), 'a') as f:
                    f.write("$point_charges\n")
                    for charge in charges_dict:
                        f.write("\t{} {} {} {}\n".format(charge["position"][0], charge["position"][1], charge["position"][2], charge["charge"]))
                        #f.write("CHARGE " + str(charge["charge"]) + " " + str(charge["radius"]) + " " + str(charge["position"][0]) + " " + str(charge["position"][1]) + " " + str(charge["position"][2]) + "")
                    f.write("$end\n")


def get_frozen_atoms(file_name):
    """
    get the two carbons most out of the plane to freeze
    Takes:
        file_name: the name of the xyz file
    Returns:
        frozen_atoms: a binary list of frozen atoms
    """
    frozen_atoms, nitrogen_xyz, carbon_xyz = [], [], []

    with open(file_name, 'r') as f:
        lines = f.readlines()
    # go through all the lines and find nitrogens
    for line in lines:
        if line[0:2].strip().isalpha():
            if line.split()[0] == "N":
                nitrogen_xyz.append(
                    [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
                    )
            if line.split()[0] == "Fe":
                fe_xyz = [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
            if line.split()[0] == "C":
                carbon_xyz.append(
                    [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
                    )
            
    # find the four nitrogens closest to the iron
    nitrogen_xyz = np.array(nitrogen_xyz)
    fe_xyz = np.array(fe_xyz)
    distances = np.linalg.norm(nitrogen_xyz - fe_xyz, axis = 1)
    closest_nitrogens = np.argsort(distances)[:4]
    closest_nitrogens_xyz = nitrogen_xyz[closest_nitrogens]
    mean_nitrogen_xyz = np.mean(closest_nitrogens_xyz, axis = 0)
    # use the mean nitrogen to find the two most out of plane
    
    mean_xyz = mean_nitrogen_xyz
    direction_1 = closest_nitrogens[0] - mean_nitrogen_xyz
    direction_2 = closest_nitrogens[1] - mean_nitrogen_xyz
    direction_3 = closest_nitrogens[2] - mean_nitrogen_xyz
    # compute cross and take the two most orthogonal directions
    dot_12 = np.dot(direction_1, direction_2)
    dot_13 = np.dot(direction_1, direction_3)
    dot_23 = np.dot(direction_2, direction_3)

    if(dot_23 > dot_13 and dot_23 > dot_12):
        direction_1 = direction_3
    if(dot_13 > dot_12 and dot_13 > dot_23):
        direction_2 = direction_3

        
    direction_1 /= np.linalg.norm(direction_1)
    direction_2 /= np.linalg.norm(direction_2)
    cross = np.cross(direction_1, direction_2)

    # find the two most out of plane carbons
    carbon_xyz = np.array(carbon_xyz)
    
    # find the two most out of plane carbons
    dot_list = [np.dot(i[0] - mean_xyz, cross) for i in carbon_xyz]
    dot_list = np.array(dot_list)
    most_out_of_plane = np.argsort(dot_list)[::-1][:2]
    # get distances to iron
    distances = np.linalg.norm(carbon_xyz - fe_xyz, axis = 1)
    # get four carbons furthest from iron
    furthest_carbons = np.argsort(distances)[::-1][:4]


    # add the two most out of plane carbons to the frozen atoms list
    for i in range(len(lines)):
        if lines[i][0:2].strip().isalpha():
            if lines[i].split()[0] == "C":
                if i in most_out_of_plane or i in furthest_carbons:
                    frozen_atoms.append(True)
                else:
                    frozen_atoms.append(False)
            else:
                frozen_atoms.append(False)

    return frozen_atoms
    


def main():

    options = get_options("./options.json")
    root = options["compressed_proteins_folder"]
    x2t_loc = options["x2t_loc"]

    for protein_name in os.listdir(root):
        if(os.path.isdir(os.path.join(root,protein_name))):
            print(protein_name)
            folder_name = root + protein_name

            # add h to pdb 
            addh("{}/{}_heme.pdb".format(folder_name, protein_name))
            addh("{}/{}_oh_heme.pdb".format(folder_name, protein_name))
            addh("{}/{}_o_heme.pdb".format(folder_name, protein_name))

            # convert pdb back to xyz
            os.system("obabel -i pdb {}/{}_heme_h.pdb -o xyz -O {}/{}_heme_h.xyz".format(folder_name, protein_name, folder_name, protein_name))
            os.system("obabel -i pdb {}/{}_oh_heme_h.pdb -o xyz -O {}/{}_oh_heme_h.xyz".format(folder_name, protein_name, folder_name, protein_name))
            os.system("obabel -i pdb {}/{}_o_heme_h.pdb -o xyz -O {}/{}_o_heme_h.xyz".format(folder_name, protein_name, folder_name, protein_name))

            # make three folders for o, oh, and normal heme
            create_folders(folder_name)
            
            # convert xyz to coord 
            os.system("{} {}/{}_heme_h.xyz > {}/no_charges/normal/coord".format(x2t_loc, folder_name, protein_name, folder_name))
            os.system("{} {}/{}_o_heme_h.xyz > {}/no_charges/o/coord".format(x2t_loc, folder_name, protein_name, folder_name))
            os.system("{} {}/{}_oh_heme_h.xyz > {}/no_charges/oh/coord".format(x2t_loc, folder_name, protein_name, folder_name))
            os.system("{} {}/{}_heme_h.xyz > {}/embedding/normal/coord".format(x2t_loc, folder_name, protein_name, folder_name))
            os.system("{} {}/{}_o_heme_h.xyz > {}/embedding/o/coord".format(x2t_loc, folder_name, protein_name, folder_name))
            os.system("{} {}/{}_oh_heme_h.xyz > {}/embedding/oh/coord".format(x2t_loc, folder_name, protein_name, folder_name))

            # get element types from xyz file
            frozen_atoms_oh = get_frozen_atoms("{}/{}_o_heme_h.xyz".format(folder_name, protein_name))
            frozen_atoms_o = get_frozen_atoms("{}/{}_oh_heme_h.xyz".format(folder_name, protein_name))
            
            elements = get_elements("{}/{}_heme_h.xyz".format(folder_name, protein_name))
            
            # write json file for turbomole 
            write_json("{}/no_charges/o/".format(folder_name), frozen_atoms = [frozen_atoms_o], charge=-3, atoms_present=elements)
            write_json("{}/no_charges/oh/".format(folder_name), frozen_atoms = [frozen_atoms_oh], charge=-3, atoms_present=elements)
            write_json("{}/no_charges/normal/".format(folder_name), frozen_atoms = [frozen_atoms_o], charge=-2, atoms_present=elements)
            write_json("{}/embedding/o/".format(folder_name), frozen_atoms = [frozen_atoms_o], charge=-3, atoms_present=elements)
            write_json("{}/embedding/oh/".format(folder_name), frozen_atoms = [frozen_atoms_oh], charge=-3, atoms_present=elements)
            write_json("{}/embedding/normal/".format(folder_name), frozen_atoms = [frozen_atoms_o], charge=-2, atoms_present=elements)
        
            setup_turbomole("{}/no_charges/o/".format(folder_name))
            setup_turbomole("{}/no_charges/oh/".format(folder_name))
            setup_turbomole("{}/no_charges/normal/".format(folder_name))
            setup_turbomole("{}/embedding/o/".format(folder_name))
            setup_turbomole("{}/embedding/oh/".format(folder_name))
            setup_turbomole("{}/embedding/normal/".format(folder_name))

            # add charges to the embedding folders
            # find pqr file in folder
            pqr_file = [f for f in os.listdir(folder_name) if f.endswith(".pqr")][0]
            charges_dict = fetch_charges_dict(os.path.join(folder_name, pqr_file))
            print("-"*20 + "charges fetched" + "-"*20)
            put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/oh/"), charges_dict)
            put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/o/"), charges_dict)
            put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/normal/"), charges_dict)

main()
