import os 

def main(): 
    # list of directories to copy from
    dirs = ["run{}".format(i+1) for i in range(5)]
    # copy all files from these directories ending in *top to target directory
    target = "../../merged_topologies/WT"
    for ind, dir in enumerate(dirs):
        for file in os.listdir(dir):
            if file.endswith("top"):
                print("{}/{} --> {}/{}_{}".format(dir, file, target, ind, file))
                os.system("cp {}/{} {}/{}_{}".format(dir, file, target, ind, file))
main()