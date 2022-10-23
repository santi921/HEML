import os 
from glob import glob
from random import sample, choice

def main():
    cpet_dir = "/ocean/projects/che160019p/santi92/cpet/"
    files_target = glob(cpet_dir + "options_topology*.txt")
    files_done = os.listdir(cpet_dir)
    charges_folder = "/ocean/projects/che160019p/santi92/processed_charges/"

    #files_target = glob('../../data/cpet/options_topology*.txt')
    #files_done = os.listdir("../../data/cpet/")
    #charges_folder = "../../data/charge_processed/"

    for i in range(10000):
        file = choice(files_target)
        protein=file.split("/")[-1][17:]
        #.split("_")[-1].split(".")[0]#.split("_")[-1]
        print("protein file: {}".format(protein))

        if(protein+".top" not in files_done):            
            launch_str = "./cpet -p {} -t 16 -o {} ".format('{}.pqr'.format(charges_folder+protein[:-4]), file)
            print(launch_str)    
            os.system(launch_str)
        print("cpet done running")

        os.system("mv {}_0.top {}{}.top".format(protein[:-4], cpet_dir, protein[:-4]))
        #pro_1262_movie_2wm4_addh_nohet_fixmer_long_run_2.txt_0.top'
    print("done running cpet")
main()
