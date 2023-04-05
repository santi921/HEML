import os, json, argparse
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from HEML.utils.data import get_options
from HEML.utils.fields import compress
from random import choice
import argparse

def main():
    """
    Usage: select a folder with MD runs in sub folder and this script will compute the distance matrix and the generate the compressed 
    dictionary for each run. It will also copy the compressed topologies to a compressed folder in the root directory.
    """

    #root = "./"
    output_folder = "./"
    path_target = './'

    argparser = argparse.ArgumentParser()
    argparser.add_argument('dampening', type=float, help='dampening factor', default=0.5)
    argparser.add_argument('maxits', type=int, help='maximum number of iterations', default=1000)
    argparser.add_argument('--compress', action='store_true', help='compress the topologies')

    args = argparser.parse_args()

    damping = args.dampening
    maxits = args.maxits
    compress = bool(args.compress)

    topo_files = [
        path_target + f for f in os.listdir(path_target) if f.endswith(".top")
    ]  # get the list of topologies in the folder
    print('number of topologies in folder "{}": {}'.format(path_target, len(topo_files)))

    # sorts the files in some way
    #topo_files.sort(key=lambda i: i.split("_")[0])
    # sort list 
    topo_files.sort()
    
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

    if compress: 
        compress_dictionary = compress(distance_matrix, damping=damping, maxits=maxits)

main()
