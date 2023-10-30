import numpy as np
import os, json, argparse
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from HEML.utils.data import get_options
from HEML.utils.fields import compress
from HEML.utils.analysis import simple_resonance_analysis
from random import choice
import pandas as pd
import pickle


def nice_output(results_dict, output_folder):
    list_means = []
    list_stds = []
    list_n_clusters = []
    list_cluster_distribution = []
    list_cluster_centers = []
    list_sweep_name = []
    list_cluster_silhouette = []

    for k, v in results_dict.items():
        list_sweep_name.append(k)
        list_means.append(v["dist_mean"])
        list_stds.append(v["dist_std"])
        list_n_clusters.append(v["compress_dictionary"]["n_clusters"])
        list_cluster_silhouette.append(v["compress_dictionary"]["silhouette"])
        list_cluster_centers_temp = []
        list_cluster_distribution_temp = []
        for sub_k, sub_v in v["compress_dictionary"].items():
            if sub_k.isnumeric():
                list_cluster_centers_temp.append(sub_v["index_center"])
                percentage_raw = float(sub_v["percentage"])
                # percentage_raw = "{:.0f}".format(percentage_raw)
                list_cluster_distribution_temp.append(percentage_raw)

        # get indices that would sort list_cluster_distribution_temp
        in_sort = np.argsort(list_cluster_distribution_temp)[::-1]
        # sort list_cluster_distribution_temp
        list_cluster_distribution_temp = np.array(list_cluster_distribution_temp)[
            in_sort
        ]
        # sort list_cluster_centers_temp
        list_cluster_centers_temp = np.array(list_cluster_centers_temp)[in_sort]
        list_cluster_centers.append(list_cluster_centers_temp)
        list_cluster_distribution.append(list_cluster_distribution_temp)

    # print the above lists in a nice table
    print(
        "sweep_name\t\t\t\t mean \t std \t n_clusters \t silhou. \t cluster_distribution \t cluster_centers"
    )

    for i in range(len(list_sweep_name)):
        raw_distr = list_cluster_distribution[i]
        # only show up to top 3
        raw_distr = raw_distr[:3]
        raw_distr = np.round(raw_distr, 1)
        raw_index = list_cluster_centers[i]
        raw_index = raw_index[:3]

        print(
            "{}\t {:.2f} \t {:.2f} \t {} \t\t {:.2f} \t\t {} \t{}".format(
                list_sweep_name[i],
                list_means[i],
                list_stds[i],
                list_n_clusters[i],
                list_cluster_silhouette[i],
                raw_distr,
                raw_index,
            )
        )
    # save the above table to a file
    with open(output_folder + "summary.txt", "w") as outfile:
        outfile.write(
            "sweep_name\t\t\t\t mean \t std \t n_clusters \t silhou. \t cluster_distribution \t cluster_centers \n"
        )
        for i in range(len(list_sweep_name)):
            raw_distr = list_cluster_distribution[i]
            # only show up to top 3
            raw_distr = raw_distr[:3]
            raw_distr = np.round(raw_distr, 1)
            raw_index = list_cluster_centers[i]
            raw_index = raw_index[:3]

            outfile.write(
                "{}\t {:.2f} \t {:.2f} \t {} \t\t {:.2f} \t\t {} \t{} \n".format(
                    list_sweep_name[i],
                    list_means[i],
                    list_stds[i],
                    list_n_clusters[i],
                    list_cluster_silhouette[i],
                    raw_distr,
                    raw_index,
                )
            )


def elementwise_compare(results_dict, output_folder, sweep_folders):
    distance_of_distances = np.zeros((len(sweep_folders), len(sweep_folders)))

    # elementwise comparision of distance matrices
    for ind_i, i in enumerate(sweep_folders):
        for ind_j, j in enumerate(sweep_folders):
            matrix1 = results_dict[i]["distance_matrix"]
            matrix2 = results_dict[j]["distance_matrix"]

            if np.array_equal(matrix1, matrix2):
                distance = 0
            else:
                distance = np.linalg.norm(matrix1 - matrix2)
            distance_of_distances[ind_i, ind_j] = distance
    # make distance_of_distances a dataframe

    df = pd.DataFrame(distance_of_distances, index=sweep_folders, columns=sweep_folders)
    df.to_csv(output_folder + "distance_of_distances.csv")


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
    print(options)
    sweep_root = options["sweep_folder"]
    output_folder = sweep_root + "sweep_out/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    reload = bool(options["sweep_parameters"]["reload"])
    damping = float(options["damping"]) if "damping" in options else 0.5
    max_iter = int(options["max_iter"]) if "max_iter" in options else 1000
    sweep_folders = os.listdir(sweep_root)
    sweep_folders = [
        f for f in sweep_folders if os.path.isdir(os.path.join(sweep_root, f))
    ]
    for i in ["output", "compressed_out", "sweep_out"]:
        if i in sweep_folders:
            sweep_folders.remove(i)
    # sort the folders by name
    sweep_folders.sort()
    sweep_folders
    results_dict = {}

    n_topos_per_folder = -1
    for ind, i in enumerate(sweep_folders):
        folder = sweep_folders[ind]
        results_dict[folder] = {}
        # merge root, folder and cpet
        path_target = os.path.join(sweep_root, folder)
        topo_files = [
            f for f in os.listdir(path_target) if f.endswith(".top")
        ]  # get the list of topologies in the folder
        print('number of topologies in folder "{}": {}'.format(folder, len(topo_files)))
        if n_topos_per_folder == -1:
            n_topos_per_folder = len(topo_files)
        assert n_topos_per_folder == len(
            topo_files
        ), "number of topologies in each folder must be the same"
        topo_files.sort()
        topo_file_name = path_target + "/topo_file_list.txt"
        with open(topo_file_name, "w") as file_list:
            for i in topo_files:
                file_list.write(f"{i} \n")
        topo_files_with_path = [path_target + "/" + f for f in topo_files]
        histograms = make_histograms(topo_files_with_path)
        if reload:
            if os.path.exists(path_target + "/distance_matrix.dat"):
                print("loading distance matrix from file")
                distance_matrix = np.loadtxt(path_target + "/distance_matrix.dat")

        distance_matrix = construct_distance_matrix(histograms)
        with open(path_target + "/distance_matrix.dat", "w") as outputfile:
            for row in distance_matrix:
                for col in row:
                    outputfile.write(f"{col} ")
                outputfile.write("\n")
        print('constructed distance matrix for folder "{}"'.format(folder))
        compress_dictionary = compress(
            distance_matrix, damping=damping, max_iter=max_iter
        )

        # add names to dictionary of files in each cluster
        labels = compress_dictionary["labels"]
        topo_files = [i.strip() for i in topo_files]
        for i in range(len(labels)):
            if "files" not in compress_dictionary[str(labels[i])]:
                compress_dictionary[str(labels[i])]["files"] = []
            compress_dictionary[str(labels[i])]["files"].append(topo_files[i])

        # compute simple resonance analysis
        # compress_dictionary = simple_resonance_analysis(
        #    compress_dictionary, run_key="run"
        # )

        # get mean and std of distance matrix
        mean = distance_matrix.mean()
        std = distance_matrix.std()
        results_dict[folder]["compress_dictionary"] = compress_dictionary
        results_dict[folder]["dist_mean"] = mean
        results_dict[folder]["dist_std"] = std
        results_dict[folder]["distance_matrix"] = distance_matrix

    # with open(output_folder + "raw_results_dict.json", "w") as outfile:
    #    json.dump(results_dict, outfile)
    # save as a pickle file
    with open(output_folder + "raw_results_dict.pkl", "wb") as outfile:
        pickle.dump(results_dict, outfile)
    nice_output(results_dict, output_folder)
    elementwise_compare(results_dict, output_folder, sweep_folders)


main()
