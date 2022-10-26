

import os 
from glob import glob
 
def find_files_that_start_with_0_and_delete(folder):
    files = os.listdir(folder)
    for i in files:
        if(i[0] == "0"):
            os.remove(folder + i)
            print("removed {}".format(i))


def find_files_that_dont_start_with_pro_and_delete(folder):
    files = os.listdir(folder)
    for i in files:
        if(i[0:3] != "pro"):
            os.remove(folder + i)
            print("removed {}".format(i))


def move_proteins_to_folders(top_files_folder = "/ocean/projects/che160019p/santi92/cpet/", protein_loc = 3):
    
    top_files = glob(top_files_folder + "*.top")
    top_files_sans_folder = [i.split("/")[-1] for i in top_files]
    protein_set = []
    for i in top_files_sans_folder: 
        protein_id = i.split("_")[protein_loc]
        if (protein_id not in protein_set):
            protein_set.append(protein_id)

    # make folders for each protein
    for i in protein_set:
        #check if folder exists
        if not os.path.exists(top_files_folder + i):
            os.mkdir(top_files_folder + i)
            print("made folder {}".format(i))
    
    # copy files to folders
    for ind, i in enumerate(top_files):
        protein_id = i.split("/")[-1].split("_")[protein_loc]
        if(ind%1000==0): 
            print("copied {} to {}".format(i, top_files_folder + protein_id))

        os.system("cp {} {}".format(i, top_files_folder + protein_id))
        