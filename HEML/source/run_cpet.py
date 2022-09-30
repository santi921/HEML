import os 
from glob import glob
from random import sample, choice

def main():
    
    #files_target = glob('../../data/cpet/options_field*.txt')
    files_target = glob('../../data/cpet/options_topology*.txt')
    files_done = os.listdir("../../data/cpet/")
    
    for i in range(100):
        file = choice(files_target)
        protein=file.split("_")[-1].split(".")[0]#.split("_")[-1]
        print("protein file: {}".format(protein))
        if(protein+".top" not in files_done):
            #print("deez")    
            os.system("./cpet -p {} -t 16 -o {} ".format('../../data/charge_processed/{}.pqr'.format(protein), file))
        print("cpet done running")
        #os.system("mv {}_0.top ../../data/cpet/{}.top".format(protein, protein))
    print("done running cpet")
main()
