import os, json, argparse
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from HEML.utils.data import get_options
from HEML.utils.fields import compress
from random import choice


def main():
    """
    Usage: select a folder with MD runs in sub folder and this script will compute the distance matrix and the generate the compressed 
    dictionary for each run. It will also copy the compressed topologies to a compressed folder in the root directory.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--options", help="location of options file", default="./options/options.json"
    )
    
    options_loc = parser.parse_args().options
    options = get_options(options_loc, create_folders=False)
    
    root = options["target_folder"]
    output_folder = options["compressed_output_folder"]

    # get list of folders in directory specified by user
    folders = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]
    
    for i in ["output", "compressed_out"]:
        if i in folders: folders.remove(i)
    # sort the folders by name 
    
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    #remove folder output if it exists from list of folders
    
    # for each folder, run the distance mat calculations and put them in the output folder
    for ind in range(len(folders)):
        folder = folders[ind]
        
        # merge root, folder and cpet
        path_target = os.path.join(root, folder, "cpet/")
        print(path_target)

        topo_files = [
            path_target + f
            for f in os.listdir(path_target)
            if f.endswith(".top")
        ] # get the list of topologies in the folder
        print("number of topologies in folder \"{}\": {}".format(folder, len(topo_files)))
        
        # sorts the files in some way
        topo_files.sort(key=lambda i: i.split("_")[0])

        with open(
            output_folder + "{}_topo_file_list.txt".format(ind), "w"
        ) as file_list:
            for i in topo_files:
                file_list.write(f"{i} \n")

        histograms = make_histograms(topo_files)
        distance_matrix = construct_distance_matrix(histograms)

        with open(
            output_folder + "{}_distance_matrix.dat".format(ind), "w"
        ) as outputfile:
            for row in distance_matrix:
                for col in row:
                    outputfile.write(f"{col} ")
                outputfile.write("\n")
        print("constucted distance matrix for folder \"{}\"".format(folder))
        compress_dictionary = compress(distance_matrix)

        print("moving central topologies to compressed folder...")

        for k, v in compress_dictionary.items():
            compress_dictionary[k]["name_center"] = topo_files[
                int(v["index_center"])
            ]
        print("saving compressed dictionary...")
        with open(
            output_folder + "{}_compressed_dictionary.json".format(ind), "w"
        ) as outputfile:
            json.dump(compress_dictionary, outputfile)

        for k, v in compress_dictionary.items():         
            name_center = v["name_center"]
            #print(name_center)
            #print("{}".format(name_center))
            #print("{}/{}_{}".format(output_folder, ind, name_center.split("/")[-1]))
            if not os.path.exists(output_folder + name_center):
                # get name of center from full path
                os.system(
                    "cp {} {}{}_{}".format( 
                        name_center, 
                        output_folder, 
                        ind, 
                        name_center.split("/")[-1])
                )

    # make a list of all the compressed topologies
    topo_files = [
        output_folder + f
        for f in os.listdir(output_folder)
        if f.endswith(".top")
    ]
    topo_files.sort(key=lambda i: i.split("_")[0])
    with open(output_folder + "topo_file_list.txt", "w") as file_list:
        for i in topo_files:
            file_list.write(f"{i} \n")

    histograms = make_histograms(topo_files)
    distance_matrix = construct_distance_matrix(histograms)

    with open(output_folder + "distance_matrix.dat", "w") as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f"{col} ")
            outputfile.write("\n")


main()
