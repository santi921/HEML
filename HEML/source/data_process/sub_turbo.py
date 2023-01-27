import os, json
import numpy as np 

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

from turbomoleio.input.define import DefineRunner

from HEML.utils.data import (
    create_folders, 
    get_options, 
    fetch_charges_dict, 
    get_elements,
    put_charges_in_turbo_files, 
    get_frozen_atoms
)


def setup_turbomole(folder_name):
    """
    Sets up the turbomole calculation.
    Takes:
        folder: the folder where the calculation should be set up
    """
    os.chdir(folder_name)
    os.system('setupturbomole.py -t')
    os.chdir("../../..")


def submit_turbomole(folder_name, check_if_done = True):
    os.chdir(folder_name)   
    # check if file called GEO_OPT_CONVERGED exists
    
    if check_if_done:
        if os.path.exists("GEO_OPT_CONVERGED"):
            os.chdir("../../..")
            return

    os.system("sbatch ./submit.sh")
    os.chdir("../../..")


def define_and_submit_turbomoleio(
    folder_name, 
    n = 4, 
    t = 24, 
    submit = False, 
    check_if_done = True, 
    frozen_atoms = [],
    atoms_present = [],
    charge = 0
    ):
    
    if check_if_done:
        if os.path.exists("GEO_OPT_CONVERGED"):
            os.chdir("../../..")
            return

    dp = get_dictionary(frozen_atoms, atoms_present, charge)
    dr = DefineRunner(parameters=dp)
    dr.run_full()
    os.system(f'sed -i "s/scforbitalshift  closedshell=.05/scforbitalshift  closedshell=.3 /" {folder_name}/control')
    if submit: submit_turbomole(folder_name, n, t, check_if_done)


def get_dictionary(frozen_atoms = [], atoms_present = [], charge = 0):

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
        "iter": 500,
        "conv": 4
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
    "geo_iterations": 600,
    "weight": False,
    "gcart": None,
    "denconv": None,
    "rij": False,
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

    return basic_dict


def write_json(folder, frozen_atoms = [], atoms_present = [], charge = 0): 
    """
    Writes a json file for the turbomole calculation. 
    Takes: 
        folder: the folder where the json file should be written
        frozen_atoms: a list of atoms that should be frozen in the calculation
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    basic_dict = get_dictionary(frozen_atoms, atoms_present, charge)
    
    with open(folder + "definput.json", "w") as outfile:
        json.dump(basic_dict, outfile, indent = 4)

   
def check_submitted(folder):
    """
    Checks if the protein has been submitted to sbatch
    Takes:
        folder: the folder of the protein
    Returns:
        True if the protein has been submitted, False otherwise
    
    """

    # check if there's a .sh file and a slurm* file 
    sh_file = False
    slurm_file = False

    for file in os.listdir(folder):
        if file.endswith(".sh"):
            sh_file = True
        if file.startswith("slurm-"):
            slurm_file = True
    
    if sh_file and slurm_file:
        return True
    return False


def clean_up(folder, filter="GEO_OPT_FAILED"):
    # check if there's a file with name filter in the folder
    # if there is, remove all the files in the folder
    for file in os.listdir(folder):
        if file.endswith(filter):        
            # remove every file that isn't a .sh file or a slurm* file or coord file
            for file in os.listdir(folder):
                if not file.endswith(".sh") and not file.startswith("slurm-") and not file.endswith(".coord"):
                    os.remove(os.path.join(folder, file))


def write_launch_sbatch(folder, time = 24, cpus=4, submit_tf = False):
    """
    Writes a launch sbatch file for the protein
    Takes:
        folder: the folder of the protein
        submit_tf: if True, submits the sbatch file
    """
    # write the sbatch file 
    with open(folder + "launch.sh", "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH -q RM-shared\n")
        outfile.write("#SBATCH -t {}:00:00\n".format(time))
        outfile.write("#SBATCH -J o\n")
        outfile.write("#SBATCH --no-requeue\n")
        outfile.write("#SBATCH -N 1\n")
        outfile.write("#SBATCH --ntasks-per-node {}\n".format(cpus))
        outfile.write("#SBATCH --mail-user=santi92\n")
        outfile.write("#SBATCH --mail-type=FAIL\n")
        outfile.write("echo \"\n\tSettings MKL_DEBUG_CPU_TYPE=5\n\"\n")
        outfile.write("export MKL_DEBUG_CPU_TYPE=5\n")
        outfile.write("qqdir=$SLURM_SUBMIT_DIR")
        outfile.write("cd $qqdir\n")
        outfile.write("echo \"\n\tSUBMIT DIRECTORY = $qqdir\n\tCORES = {} \n\tTIME = {} \n\tSCRATCH = $LOCAL \n\tDATE = `date`\n\"\n".format(cpus, time))
        outfile.write("echo \"\n\tAdding OpenBabel to Path\n\"\n")
        outfile.write("export PATH=/ocean/projects/che160019p/shared/openbabel-2.4.0/bin/:$PATH\n")
        outfile.write("unalis -a\n")
        outfile.write("export PYTHONUNBUFFERED=1\n")
        outfile.write("time /ocean/projects/che160019p/santi92/phd3/phd3/bin/runturbomole.py -n {} -t {}\n".format(cpus, time))
        outfile.write("exit 0\n")

    if submit_tf:
        os.system("sbatch " + folder + "launch.sh")
        

def main():
    submit_tf = False
    only_submit = False
    cleanup_tf = True
    embedd_tf = False

    options = get_options("./options.json")
    root = options["compressed_proteins_folder"]
    x2t_loc = options["x2t_loc"]

    for ind, protein_name in enumerate(os.listdir(root)):
        if(os.path.isdir(os.path.join(root,protein_name))):
            print(protein_name)
            try: 
                folder_name = root + protein_name
                # check if the protein has been submitted to sbatch 
                if not check_submitted(folder_name):
                    if cleanup_tf: 
                        clean_up(folder_name, filter="GEO_OPT_FAILED")
                        if not only_submit: 
                            # add h to pdb 
                          
                            # make three folders for o, oh, and normal heme
                            create_folders(folder_name)
                            
                            # convert xyz to coord 
                            os.system("{} {}/{}_heme_h.xyz > {}/no_charges/normal/coord".format(x2t_loc, folder_name, protein_name, folder_name))
                            os.system("{} {}/{}_o_heme_h.xyz > {}/no_charges/o/coord".format(x2t_loc, folder_name, protein_name, folder_name))
                            os.system("{} {}/{}_oh_heme_h.xyz > {}/no_charges/oh/coord".format(x2t_loc, folder_name, protein_name, folder_name))
                            os.system("{} {}/{}_heme_h.xyz > {}/embedding/normal/coord".format(x2t_loc, folder_name, protein_name, folder_name))
                            os.system("{} {}/{}_o_heme_h.xyz > {}/embedding/o/coord".format(x2t_loc, folder_name, protein_name, folder_name))
                            os.system("{} {}/{}_oh_heme_h.xyz > {}/embedding/oh/coord".format(x2t_loc, folder_name, protein_name, folder_name))

                            # get some info from xyz
                            frozen_atoms_oh = get_frozen_atoms("{}/{}_oh_heme_h.xyz".format(folder_name, protein_name))
                            frozen_atoms_o = get_frozen_atoms("{}/{}_o_heme_h.xyz".format(folder_name, protein_name))
                            frozen_atoms_heme = get_frozen_atoms("{}/{}_heme_h.xyz".format(folder_name, protein_name))

                            elements = get_elements("{}/{}_heme_h.xyz".format(folder_name, protein_name))
                            
                            if embedd_tf:

                                define_and_submit_turbomoleio("{}/embedding/oh/".format(folder_name), frozen_atoms_oh, elements, submit=False, charge=-3)
                                define_and_submit_turbomoleio("{}/embedding/o/".format(folder_name), frozen_atoms_o, elements, submit=False, charge=-3)
                                define_and_submit_turbomoleio("{}/embedding/normal/".format(folder_name), frozen_atoms_heme, elements, submit=False, charge=-2)

                                # find pqr file in folder
                                pqr_file = [f for f in os.listdir(folder_name) if f.endswith(".pqr")][0]
                                charges_dict = fetch_charges_dict(os.path.join(folder_name, pqr_file))
                                print("-"*20 + "charges fetched" + "-"*20)
                                
                            
                                put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/oh/"), charges_dict)
                                put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/o/"), charges_dict)
                                put_charges_in_turbo_files(os.path.join(folder_name, "/embedding/normal/"), charges_dict)
                            
                            else:
                                define_and_submit_turbomoleio("{}/embedding/oh/".format(folder_name), frozen_atoms_oh, elements, submit=True, charge=-3)
                                define_and_submit_turbomoleio("{}/embedding/o/".format(folder_name), frozen_atoms_o, elements, submit=True, charge=-3)
                                define_and_submit_turbomoleio("{}/embedding/normal/".format(folder_name), frozen_atoms_heme, elements, submit=True, charge=-2)
                        
                        if submit_tf:
                            submit_turbomole("{}/embedding/o/".format(folder_name), t = 24, n = 4)
                            submit_turbomole("{}/embedding/oh/".format(folder_name), t = 24, n = 4)
                            submit_turbomole("{}/embedding/normal/".format(folder_name), t = 24, n = 4)
                        
                        else: 
                            print("not submitting calculations")

                        print("done with {} of {}".format(ind, len(os.listdir(root))))
                
                if only_submit: 
                    submit_turbomole("{}/embedding/o/".format(folder_name), t = 24, n = 4)
                    submit_turbomole("{}/embedding/oh/".format(folder_name), t = 24, n = 4)
                    submit_turbomole("{}/embedding/normal/".format(folder_name), t = 24, n = 4)
                    
            except:
                print("error with {}".format(protein_name))
                continue


main()
