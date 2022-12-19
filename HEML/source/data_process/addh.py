import os, json
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
    #run("mv ./temp.pdb " + pdb_file)

def write_json(folder, frozen_atoms = []): 
    """
    Writes a json file for the turbomole calculation. 
    Takes: 
        folder: the folder where the json file should be written
        frozen_atoms: a list of atoms that should be frozen in the calculation
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    basic_dict = '''
    { 
        "geometry": {
            "cartesian": false,
            "idef": {"idef_on": false},
            "ired" : false,
            "iaut": {"iaut_on": false}
        },
        "dft": {
            "dft_on" : true,
            "func" : "tpss",
            "grid" : "m4"
        },
        "scf": {
            "iter": 300,
            "conv": 5
        },
        "stp": {
            "itvc": 0,
            "trad": 0.1
        },
        "open_shell" : {
            "open_shell_on" : true,
            "unpaired" : 1
        },
        "basis" : {
            "all" : "def2-SVP",
            "fe" : "def2-TZVP",
            "n" : "def2-TZVP",
            "s" : "def2-TZVP",
            "o" : "def2-TZVP"
        },
        "cosmo": 4,
        "freeze_atoms" : [],
        "calculation" : "geo",
        "geo_iterations" : 200,
        "weight" : false,
        "gcart" : null,
        "denconv" : null,
        "rij": true,
        "marij": true,
        "dsp": true,
        "charges": -2
    }'''
    
    basic_dict_json = json.loads(basic_dict)

    if frozen_atoms != []:
        basic_dict_json['freeze_atoms'] = frozen_atoms

    json.dump(basic_dict_json, folder + "definput.json", indent = 4)


def main():

    for protein_name in os.listdir("./"):
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


            # add h to pdb 
            addh("{}/{}_heme.pdb".format(folder_name, protein_name))
            addh("{}/{}_oh_heme.pdb".format(folder_name, protein_name))
            addh("{}/{}_o_heme.pdb".format(folder_name, protein_name))

            # convert pdb back to xyz
            os.system("obabel -i pdb ./{}/{}_heme_h.pdb -o xyz -O ./{}/{}_heme_h.xyz".format(folder_name, folder_name, folder_name, folder_name))
            os.system("obabel -i pdb ./{}/{}_oh_heme_h.pdb -o xyz -O ./{}/{}_oh_heme_h.xyz".format(folder_name, folder_name, folder_name, folder_name))
            os.system("obabel -i pdb ./{}/{}_o_heme_h.pdb -o xyz -O ./{}/{}_o_heme_h.xyz".format(folder_name, folder_name, folder_name, folder_name))

            # make three folders for o, oh, and normal heme

            if not os.path.exists("{}/o".format(folder_name)):
                os.makedirs("{}/o".format(folder_name))
            if not os.path.exists("{}/oh".format(folder_name)):
                os.makedirs("{}/oh".format(folder_name))
            if not os.path.exists("{}/normal".format(folder_name)):
                os.makedirs("{}/normal".format(folder_name))

            # convert xyz to coord 
            os.system("./x2t ./{}/{}_heme_h.xyz > ./{}/normal/coord".format(folder_name, folder_name, folder_name))
            os.system("./x2t ./{}/{}_o_heme_h.xyz > ./{}/o/coord".format(folder_name, folder_name, folder_name))
            os.system("./x2t ./{}/{}_oh_heme_h.xyz > ./{}/oh/coord".format(folder_name, folder_name, folder_name))

            # write json file for turbomole 
            write_json("{}/o/".format(folder_name), frozen_atoms = [])
            write_json("{}/oh/".format(folder_name), frozen_atoms = [])
            write_json("{}/normal/".format(folder_name), frozen_atoms = [])
            
            # run setupphd3.py 
            os.chdir("{}".format(folder_name))
            os.system('setupturbomole.py')
            os.chdir("..")

            os.chdir("{}/normal/".format(folder_name))
            os.system('setupturbomole.py')
            os.chdir("../..")

            os.chdir("{}/o/".format(folder_name))
            os.system('setupturbomole.py')
            os.chdir("../..")
            

            os.chdir("{}/oh/".format(folder_name))
            os.system('setupturbomole.py')
            os.chdir("../..")
            

main()