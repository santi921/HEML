import os, json
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from HEML.utils.data import compress, get_options
from random import choice

def main():
    options = get_options("./options.json")    
    root = options["cpet_folder"]
    
    # create a new folder with all compressed dictionaries
    folders = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]
    merge_dict = {}
    compressed_folder = root + "/compressed/"
    if not os.path.exists(compressed_folder):
        os.mkdir(compressed_folder)

    for i in range(len(folders)):
        folder = folders[i]
        if os.path.exists(os.path.join(root, folder) + "/compressed_dictionary.json"):
            with open(os.path.join(root, folder) + "/compressed_dictionary.json", 'r') as inputfile:
                compressed_dictionary = json.load(inputfile)

            for k, v in compressed_dictionary.items():
                name_center = v["name_center"]
                if not os.path.exists(compressed_folder + name_center):
                    # get name of center from full path 
                    os.system("cp " + name_center + " " + compressed_folder + name_center.split("/")[-1])
            os.system("cp " + os.path.join(root, folder) + "/compressed_dictionary.json " + compressed_folder + folder + "_compressed_dictionary.json")
            #folder_name = folder.split("_")[0]
            #merge_dict[]
        
    topo_files = [compressed_folder+f for f in os.listdir(compressed_folder) if f.endswith(".top")]
    #topo_files = 
    # sorts the files in some way
    topo_files.sort(key=lambda i: i.split("_")[0])
    with open(compressed_folder + 'topo_file_list.txt', 'w') as file_list:
        for i in topo_files:
            file_list.write(f'{i} \n')

    histograms = make_histograms(topo_files)
    distance_matrix = construct_distance_matrix(histograms)
    
    with open(compressed_folder + "distance_matrix.dat", 'w') as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f'{col} ')
            outputfile.write("\n")

main()