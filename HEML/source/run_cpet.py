import os 
from glob import glob
from random import sample, choice

def main():
    
    #files_target = glob("/ocean/projects/che160019p/santi92/cpet/options_topology*.txt")
    #files_done = os.listdir("/ocean/projects/che160019p/santi92/cpet/")
    #charges_folder = os.listdir("/ocean/projects/che160019p/santi92/charges_processed/")

    files_target = glob('../../data/cpet/options_topology*.txt')
    files_done = os.listdir("../../data/cpet/")
    charges_folder = "../../data/charge_processed/"

    for i in range(100):
        file = choice(files_target)
        protein=file.split("_")[-1].split(".")[0]#.split("_")[-1]
        print("protein file: {}".format(protein))
        if(protein+".top" not in files_done):
            
            launch_str = "./cpet -p {} -t 16 -o {} ".format('{}.pqr'.format(charges_folder+protein), file)
            print(launch_str)    
            os.system(launch_str)

        print("cpet done running")
        #os.system("mv {}_0.top ../../data/cpet/{}.top".format(protein, protein))
    print("done running cpet")
main()
