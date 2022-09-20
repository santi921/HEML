import os 
from glob import glob


def main():
    
    #files_target = glob('../../data/cpet/options_field*.txt')
    files_target = glob('../../data/cpet/options_topology*.txt')
    for file in files_target:
        protein=file.split("_")[-1].split(".")[0]#.split("_")[-1]
        os.system("./cpet -p {} -o {} ".format('../../data/charge_processed/{}.pqr'.format(protein), file))
        os.system("mv {}_0.top ../../data/cpet/{}.top".format(protein, protein))
    print("done running cpet")
main()