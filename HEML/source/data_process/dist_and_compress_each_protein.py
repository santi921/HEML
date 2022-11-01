import os, json
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from HEML.utils.data import compress
from random import choice

def main():
    root = "/ocean/projects/che160019p/santi92/cpet/"
    # get list of folders in directory specified by user
    folders = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

    # for each folder, run the distance mat calculations
    for i in range(len(folders)):
        folder = choice(folders)
        print(folder)
        if os.path.exists(os.path.join(root, folder) + "/compressed_distance_matrix.json"):
            print("compressed already exists")
            
        else: 
            topo_files = [os.path.join(root, folder)+"/"+f for f in os.listdir(os.path.join(root, folder)) if f.endswith(".top")]
            # sorts the files in some way
            topo_files.sort(key=lambda i: i.split("_")[0])

            with open(os.path.join(root, folder) + '/topo_file_list.txt', 'w') as file_list:
                for i in topo_files:
                    file_list.write(f'{i} \n')

            histograms = make_histograms(topo_files)
            distance_matrix = construct_distance_matrix(histograms)
            
            with open(os.path.join(root, folder) + "/distance_matrix.dat", 'w') as outputfile:
                for row in distance_matrix:
                    for col in row:
                        outputfile.write(f'{col} ')
                    outputfile.write("\n")

            compress_dictionary = compress(distance_matrix)

            for k, v in compress_dictionary.items():
                compress_dictionary[k]["name_center"] = topo_files[int(v["index_center"])]

            with open(os.path.join(root, folder) + "/compressed_dictionary.json", 'w') as outputfile:
                json.dump(compress_dictionary, outputfile)

main()