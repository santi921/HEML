import os 
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

            # convert pdb back to xyz
            os.system("obabel -i pdb ./{}/{}_heme_h.pdb -o xyz -O ./{}/{}_heme_h.xyz".format(folder_name, folder_name, folder_name, folder_name))

            # convert xyz to coord 
            os.system("./x2t ./{}/{}_heme_h.xyz > ./{}/coord".format(folder_name, folder_name, folder_name))

            # run setupphd3.py 
            os.system('cd {}'.format(folder_name))
            os.system('setupphd3.py')
            os.system('cd ..')

main()