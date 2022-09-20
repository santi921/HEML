import os, random
from chimera import runCommand as rc
from os import system as run
#from Dockprep import prep

# todo:
targets = []

with open("../../data/het_list.txt") as f:
    del_list = f.readlines()
    del_list = [i[0:-1] for i in del_list]

#del_list = del_list[0:3] 

files = os.listdir("../../data/pdbs/")


for i in files:
    rc("open ../../data/pdbs/" + i)
    rc("delete solvent")
    rc("addh")
    [rc("delete :" + j) for j in del_list]
    print(i)
    rc("write #0 ../../data/pdbs_processed/" + i)
    rc("close session")
    run("chimera --nogui ../../data/pdbs_processed/" + i + " incompleteSideChains.py")
    run("mv ./temp.pdb ../../data/pdbs_processed/" + i)
rc("stop now")



#files_processed = os.listdir("../../data/pdbs_processed/")
#for i in files_processed:
    