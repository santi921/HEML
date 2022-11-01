import os, random
from chimera import runCommand as rc
from os import system as run
from glob import glob
#from Dockprep import prep

# todo:
targets = []

with open("../../data/het_list.txt") as f:
    del_list = f.readlines()
    del_list = [i[0:-1] for i in del_list]

#del_list = del_list[0:3] 
pdb_folder = "/ocean/projects/che160019p/santi92/heme_traj/*"
out_folder = "/ocean/projects/che160019p/santi92/pdbs_processed_heme/"

#pdb_folder = "../../data/pdbs/"
#files = os.listdir(pdb_folder)
files = glob(pdb_folder)
files_out = glob(out_folder+"*")
files_out_charges = os.listdir("/ocean/projects/che160019p/santi92/heme_charges/")

#print(random.sample(files, 5000))
for i in random.sample(files,100000):
    print("-"* 40)
    try:
        print("pro_" + i.split("/")[-1] + ".pqr" in files_out_charges)
        #print("pro_" + i.split("/")[-1] + ".pqr")
        #if(out_folder + "pro_" + i.split("/")[-1] not in files_out):
        if("pro_" + i.split("/")[-1] + ".pqr" not in files_out_charges):        
            rc("open " + pdb_folder[:-1] + i.split("/")[-1])
            rc("delete solvent")
            rc("addh")
            [rc("delete :" + j) for j in del_list]
            #print(i)
            rc("write #0 " + pdb_folder[:-1] + "pro_" + i.split("/")[-1])
            print("write #0 " + pdb_folder[:-1] + "pro_" + i.split("/")[-1])
            rc("close session")
            print("processed h + solvent")
            run("chimera --nogui " +  pdb_folder[:-1] + "pro_" + i.split("/")[-1] + " incompleteSideChains.py")
            run("mv ./temp.pdb " + out_folder + i.split("/")[-1])
            os.remove(pdb_folder[:-1] + "pro_" + i.split("/")[-1])
    except: 
        passrc("stop now")



#files_processed = os.listdir("../../data/pdbs_processed/")
#for i in files_processed:
    
