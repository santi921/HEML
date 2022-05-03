import os, random
from chimera import runCommand as rc
from os import system as run
#from Dockprep import prep

# todo:
targets = []

with open("het_list.txt") as f:
    del_list = f.readlines()
    del_list = [i[0:-1] for i in del_list]

#del_list = del_list[0:3] 

files = os.listdir("./pdbs/")


for i in files:
    rc("open ./pdbs/" + i)
    rc("delete solvent")
    rc("addh")
    [rc("delete :" + j) for j in del_list]
    print(i)
    rc("write #0 ./pdbs_processed/" + i)
    rc("close session")
    run("chimera --nogui ./pdbs_processed/" + i + " incompleteSideChains.py")
    run("mv ./temp.pdb ./pdbs_processed/" + i)
rc("stop now")
