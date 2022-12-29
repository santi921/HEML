import os, random, json
from chimera import runCommand as rc
from os import system as run
from glob import glob


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
    targets = []

    with open("../../../data/het_list.txt") as f:
        del_list = f.readlines()
        del_list = [i[0:-1] for i in del_list]
    del_list = del_list[0:3] 

    # read json file for folder locations  

    options = get_options("./options.json")
    pdb_folder = str(options["pdb_folder"])
    output_folder = str(options["processed_pdb_folder"])
    charges_folder = str(options["charges_folder"])

    #files = os.listdir(pdb_folder)
    files = glob(pdb_folder + "*.pdb*")
    files_out = glob(output_folder+"*")
    files_out_charges = os.listdir(charges_folder)

    #print(random.sample(files, 5000))
    for i in random.sample(files,100):
        print("-"* 40)
        # check if file is already processed
        if not os.path.exists(output_folder + i.split("/")[-1]):
            print("pro_" + i.split("/")[-1] + ".pqr" in files_out_charges)
            if("pro_" + i.split("/")[-1] + ".pqr" not in files_out_charges): 
                rc("open " + pdb_folder + i.split("/")[-1])
                rc("delete solvent")
                rc("addh")
                [rc("delete :" + j) for j in del_list]
                rc("write #0 " + pdb_folder + "pro_" + i.split("/")[-1])
                print("write #0 " + pdb_folder + "pro_" + i.split("/")[-1])
                rc("close session")
                print("processed h + solvent")
                run("chimera --nogui " +  pdb_folder[:-1] + "pro_" + i.split("/")[-1] + " incompleteSideChains.py")
                run("mv ./temp.pdb " + output_folder + i.split("/")[-1])
                os.remove(pdb_folder + "pro_" + i.split("/")[-1])


main()
