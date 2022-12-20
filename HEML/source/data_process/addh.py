import os, json, ast
from chimera import runCommand as rc
from os import system as run

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


def write_json(folder, frozen_atoms = []): 
    """
    Writes a json file for the turbomole calculation. 
    Takes: 
        folder: the folder where the json file should be written
        frozen_atoms: a list of atoms that should be frozen in the calculation
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    basic_dict = {}
    basic_dict["geometry"] = {
        "cartesians": True,
    }
    basic_dict["geometry"]["idef"] = {
        "idef_on": False,

    }
    basic_dict["geometry"]["ired"] = False
    basic_dict["geometry"]["iaut"] = {
        "iaut_on": False,
    }

    basic_dict["dft"] = {
        "dft_on": False,
        "func": "b3lyp",
        "grid": "m4"
    }

    basic_dict["scf"] = {
        "iter": 300,
        "conv": 5
    }
    basic_dict["stp"] = {
        "itvc": 0,
        "trad": 0.1
    }

    basic_dict["open_shell"] = {
        "open_shell_on": False,
        "unpaired": 0
    }

    basic_dict["basis"] = {
        "all": "def2-SVP",
        "fe": "def2-TZVP",
        "n": "def2-TZVP",
        "s": "def2-TZVP",
        "o": "def2-TZVP"
    }

    basic_dict["cosmo"] = 4
    basic_dict["freeze_atoms"] = []
    basic_dict["calculation"] = "geo"
    basic_dict["geo_iterations"] = 200
    basic_dict["weight"] = False
    basic_dict["gcart"] = None 
    basic_dict["denconv"] = None
    basic_dict["rij"] = True
    basic_dict["marij"] = True
    basic_dict["dsp"] = True
    basic_dict["charges"] = -2

    basic_dict_json = json.loads(basic_dict)

    if frozen_atoms != []:
        basic_dict_json['freeze_atoms'] = frozen_atoms
    # remove unicode from json
    dict_no_unicode = {}
    for key, value in basic_dict_json.items():
        if value is dict :
            dict_no_unicode[key] = {}
            for key2, value2 in value.items():
                dict_no_unicode[key][key2] = str(value2)
        else:
            dict_no_unicode[key] = str(value)

    dict = dict_no_unicode
    with open(folder + "definput.json", "w") as outfile:
        json.dump(dict, outfile)


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


def main():

    options = get_options("./options.json")
    root = options["compressed_proteins_folder"]
    x2t_loc = options["x2t_loc"]

    for protein_name in os.listdir(root):
        if(os.path.isdir(protein_name)):

            print(protein_name)
            folder_name = root + protein_name
            #check if pdb is in folder
            if os.path.exists(os.path.join(folder_name, protein_name + ".pdb")):
                pdb_file = protein_name + ".pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_movie.pdb")):
                pdb_file = protein_name + "_movie.pdb"
            if os.path.exists(os.path.join(folder_name, protein_name + "_todo_process.pdb")):
                pdb_file = protein_name + "_todo_process.pdb"

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

            # write json file for turbomole 
            write_json("{}/no_charges/o/".format(folder_name), frozen_atoms = [])
            write_json("{}/no_charges/oh/".format(folder_name), frozen_atoms = [])
            write_json("{}/no_charges/normal/".format(folder_name), frozen_atoms = [])
            write_json("{}/embedding/o/".format(folder_name), frozen_atoms = [])
            write_json("{}/embedding/oh/".format(folder_name), frozen_atoms = [])
            write_json("{}/embedding/normal/".format(folder_name), frozen_atoms = [])
        
            setup_turbomole("{}/no_charges/o/".format(folder_name))
            setup_turbomole("{}/no_charges/oh/".format(folder_name))
            setup_turbomole("{}/no_charges/normal/".format(folder_name))
            setup_turbomole("{}/embedding/o/".format(folder_name))
            setup_turbomole("{}/embedding/oh/".format(folder_name))
            setup_turbomole("{}/embedding/normal/".format(folder_name))
main()