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

    #root = "./"
    output_folder = "./"
    path_target = './'


    topo_files = [
        path_target + f for f in os.listdir(path_target) if f.endswith(".top")
    ]  # get the list of topologies in the folder
    print('number of topologies in folder "{}": {}'.format(path_target, len(topo_files)))

    # sorts the files in some way
    topo_files.sort(key=lambda i: i.split("_")[0])

    with open(
        output_folder + "topo_file_list.txt", "w"
    ) as file_list:
        for i in topo_files:
            file_list.write(f"{i} \n")

    histograms = make_histograms(topo_files)
    distance_matrix = construct_distance_matrix(histograms)

    with open(
        output_folder + "distance_matrix.dat", "w"
    ) as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f"{col} ")
            outputfile.write("\n")
    print('constucted distance matrix for folder "{}"'.format(path_target))


main()
