import os, json, argparse
from HEML.utils.cpet import make_histograms, construct_distance_matrix, read_distance_matrix
from HEML.utils.data import get_options
from HEML.utils.fields import compress
from HEML.utils.analysis import simple_resonance_analysis
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
    # check if the key damping exists in the options file, if not set it to 0.5
    damping = float(options["damping"]) if "damping" in options else 0.5
    max_iter = int(options["max_iter"]) if "max_iter" in options else 1000
    reload_distance_matrix = bool(options["reload_dist_mat"]) if "reload_dist_mat" in options else False

    # get list of folders in directory specified by user
    folders = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

    for i in ["output", "compressed_out"]:
        if i in folders:
            folders.remove(i)
    # sort the folders by name
    folders.sort()

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    # remove folder output if it exists from list of folders

    # for each folder, run the distance mat calculations and put them in the output folder
    for ind in range(len(folders)):
        folder = folders[ind]

        # merge root, folder and cpet
        path_target = os.path.join(root, folder, "cpet/")

        topo_files = [
            path_target + f for f in os.listdir(path_target) if f.endswith(".top")
        ]  # get the list of topologies in the folder
        print('number of topologies in folder "{}": {}'.format(folder, len(topo_files)))

        # sorts the files in some way
        # topo_files.sort(key=lambda i: i.split("_")[0])
        topo_files.sort()

        print(topo_files)
        with open(
            output_folder + "{}_topo_file_list.txt".format(ind), "w"
        ) as file_list:
            for i in topo_files:
                file_list.write(f"{i} \n")

        if reload_distance_matrix:
            try:
                distance_matrix = read_distance_matrix(output_folder + "{}_distance_matrix.dat".format(ind))
            except:
                print("Pre-existing distance matrix not found")
        else:
            histograms = make_histograms(topo_files)
            distance_matrix = construct_distance_matrix(histograms)
        
            with open(
                output_folder + "{}_distance_matrix.dat".format(ind), "w"
            ) as outputfile:
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

        print("moving central topologies to compressed folder...")
        for k, v in compress_dictionary.items():
            if (
                k != "total_count"
                and k != "silhouette"
                and k != "labels"
                and k != "n_clusters"
            ):
            if k.isnumeric():
                compress_dictionary[k]["name_center"] = topo_files[
                    int(v["index_center"])
                ]

        print("saving compressed dictionary...")
        with open(
            output_folder + "{}_compressed_dictionary.json".format(ind), "w"
        ) as outputfile:
            json.dump(compress_dictionary, outputfile)

        for k, v in compress_dictionary.items():
            if (
                k != "total_count"
                and k != "silhouette"
                and k != "labels"
                and k != "n_clusters"
            ):
            if k.isnumeric():
                name_center = v["name_center"]
                if not os.path.exists(output_folder + name_center):
                    # get name of center from full path
                    os.system(
                        "cp {} {}{}_{}".format(
                            name_center, output_folder, ind, name_center.split("/")[-1]
                        )
                    )

    # make a list of all the compressed topologies
    topo_files = [
        output_folder + f for f in os.listdir(output_folder) if f.endswith(".top")
    ]
    # topo_files.sort(key=lambda i: i.split("_")[0])
    topo_files.sort()

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
