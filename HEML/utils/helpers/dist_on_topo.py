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

    # root = "./"
    output_folder = "./compressed_out/"
    path_target = "./"

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-dampening", type=float, help="dampening factor", default=0.5
    )
    argparser.add_argument(
        "-max_iter", type=int, help="maximum number of iterations", default=1000
    )
    argparser.add_argument(
        "--compress", action="store_true", help="compress the topologies"
    )

    args = argparser.parse_args()

    damping = float(args.dampening)
    max_iter = float(args.max_iter)
    compress_tf = bool(args.compress)
    # if outputfolder doesnt exist, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    topo_files = [
        path_target + f for f in os.listdir(path_target) if f.endswith(".top")
    ]  # get the list of topologies in the folder
    print(
        'number of topologies in folder "{}": {}'.format(path_target, len(topo_files))
    )

    # sorts the files in some way
    # topo_files.sort(key=lambda i: i.split("_")[0])
    # sort list
    topo_files.sort()

    with open(output_folder + "topo_file_list.txt", "w") as file_list:
        for i in topo_files:
            file_list.write(f"{i} \n")

    histograms = make_histograms(topo_files)
    distance_matrix = construct_distance_matrix(histograms)

    # construct distance matrix
    with open(output_folder + "distance_matrix.dat", "w") as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f"{col} ")
            outputfile.write("\n")

    print('constucted distance matrix for folder "{}"'.format(path_target))

    if compress_tf:
        compress_dictionary = compress(
            distance_matrix, damping=damping, max_iter=int(max_iter), names=topo_files
        )
        with open(output_folder + "loc_compressed_dictionary.json", "w") as outputfile:
            json.dump(compress_dictionary, outputfile, indent=4)

        for k, v in compress_dictionary.items():
            if (
                k != "total_count"
                and k != "silhouette"
                and k != "labels"
                and k != "n_clusters"
            ):
                print(v)
                name_center = v["name"]
                if not os.path.exists(output_folder + name_center):
                    # get name of center from full path
                    os.system(
                        "cp {} {}/{}".format(name_center, output_folder, name_center)
                    )


main()
