import os
from glob import glob
from random import choice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


matplotlib.rcParams.update(
    {  # Use mathtext, not LaTeX
        "text.usetex": False,
        # Use the Computer modern font
        "font.family": "serif",
        "font.serif": "cmr10",
        "mathtext.fontset": "cm",
        # Use ASCII minus
        "axes.unicode_minus": False,
        "font.size": 16,
    }
)


def run_box_calcs(cpet_path, target_path, charges_dir):
    files_target = glob(target_path + "options_field*.txt")
    files_done = os.listdir(target_path)
    charges_dir = charges_dir

    for i in range(20000):
        file = choice(files_target)
        files_target.remove(file)
        protein = file.split("/")[-1][14:]  # works for protein movies

        print("protein file: {}".format(protein))

        if protein + ".top" not in files_done:
            launch_str = "{} -p {} -t 16 -o {} ".format(
                cpet_path, "{}.pqr".format(charges_dir + protein[:-4]), file
            )
            print(launch_str)
            os.system(launch_str)
        print("cpet done running")

    print("done running cpet")


def run_mag_calcs(cpet_path, target_path, charges_dir):
    files_target = glob(target_path + "options_mag*.txt")
    files_done = os.listdir(target_path)

    # filter files that dont end in .dat
    files_done = [i for i in files_done if i[-4:] == ".dat"]
    charges_dir = charges_dir

    for i in range(20000):
        file = choice(files_target)
        files_target.remove(file)
        protein = file.split("/")[-1][14:]  # works for protein movies

        # .split("_")[-1].split(".")[0]#.split("_")[-1]
        print("protein file: {}".format(protein))

        if protein + ".mag.dat" not in files_done:
            launch_str = "{} -p {} -o {} ".format(
                cpet_path, "{}.pqr".format(charges_dir + protein[:-4]), file
            )
            print(launch_str)
            os.system(launch_str)


def run_topology_calcs(cpet_path, target_path, charges_dir, num=10000, threads=16):
    files_target = glob(target_path + "options_topol*.txt")
    files_done = os.listdir(target_path)
    charges_dir = charges_dir

    for i in range(num):
        file = choice(files_target)
        files_target.remove(file)
        protein = file.split("/")[-1][14:]

        print("protein file: {}".format(protein.split(".")[0]))

        if protein.split(".")[0] + ".top" not in files_done:
            launch_str = "{} -p {} -t {} -o {} ".format(
                cpet_path, "{}.pqr".format(charges_dir + protein[:-4]), threads, file
            )
            print(launch_str)
            os.system(launch_str)
        print("cpet done running")

        os.system(
            "mv efield_topo_{}_0.top {}{}.top".format(
                protein[:-4], target_path, protein[:-4]
            )
        )
    print("done running cpet")


def run_sweep_topology_calcs(cpet_path, sweep_root, charges_dir, num=10000, threads=16):
    charges_dir = charges_dir
    # in sweep root get all of the folders in that directory
    sweep_folders = os.listdir(sweep_root)
    # for each folder in the sweep root
    files_target = glob(sweep_root + sweep_folders[0] + "/options_topol*.txt")
    files_done = glob(sweep_root + sweep_folders[0] + "/*.top")
    print(files_done)
    # files_done = os.listdir(target_path)

    for i in range(num):
        file = choice(files_target)
        files_target.remove(file)

        protein = file.split("/")[-1][14:]  # this is currently hard coded
        print("protein file: {}".format(protein.split(".")[0]))

        # check last folder to make sure all are done
        if (
            sweep_root + sweep_folders[-1] + protein.split(".")[0] + ".top"
            not in files_done
        ):
            # inner loop to iterate over all folders
            for folder_single in sweep_folders:
                cpet_options_file_single = (
                    sweep_root + folder_single + "/" + file.split("/")[-1]
                )
                full_path_topo = sweep_root + folder_single
                launch_str = "{} -p {} -t {} -o {} ".format(
                    cpet_path,
                    "{}.pqr".format(charges_dir + protein[:-4]),
                    threads,
                    cpet_options_file_single,
                )
                print(launch_str)
                os.system(launch_str)
                os.system(
                    "mv topo_{}_0.top {}/{}.top".format(
                        protein[:-4], full_path_topo, protein[:-4]
                    )
                )
                print("cpet done running")
        else:
            print("protein {} already done".format(protein.split(".")[0]))

    print("reached terminate condition")


def make_histograms(topo_files, plot=False):
    histograms = []

    for topo_file in topo_files:
        curvatures, distances = [], []

        with open(topo_file) as topology_data:
            for line in topology_data:
                if line.startswith("#"):
                    continue

                line = line.split(",")
                distances.append(float(line[0]))
                curvatures.append(float(line[1]))

        # bins is number of histograms bins in x and y direction (so below is 100x100 bins)
        # range gives xrange, yrange for the histogram
        a, b, c, q = plt.hist2d(
            distances,
            curvatures,
            bins=200,
            range=[[0, 10], [0, 30]],
            norm=matplotlib.colors.LogNorm(),
            density=True,
            cmap="jet",
        )

        NormConstant = 0
        for j in a:
            for m in j:
                NormConstant += m

        actual = []
        for j in a:
            actual.append([m / NormConstant for m in j])

        actual = np.array(actual)
        histograms.append(actual.flatten())
        if plot:
            plt.show()

    return np.array(histograms)


def distance_numpy(hist1, hist2):
    a = (hist1 - hist2) ** 2
    b = hist1 + hist2
    return np.sum(np.divide(a, b, out=np.zeros_like(a), where=b != 0)) / 2.0


def construct_distance_matrix(histograms):
    matrix = np.diag(np.zeros(len(histograms)))
    for i, hist1 in enumerate(histograms):
        for j, hist2 in enumerate(histograms[i + 1 :]):
            j += i + 1
            matrix[i][j] = distance_numpy(hist1, hist2)
            matrix[j][i] = matrix[i][j]

    return matrix


def config_to_folder(single_sweep_config):
    sweep_string = "hist_{}_step_{}_samples_{}_box_{}".format(
        int(single_sweep_config["hist_bins"]),
        str(single_sweep_config["step_size"]).split(".")[-1],
        int(single_sweep_config["samples"]),
        str(single_sweep_config["box_size"]).replace(".", "_"),
    )
    return sweep_string


def config_to_folder(single_sweep_config):
    sweep_string = "hist_{}_step_{}_samples_{}_box_{}".format(
        int(single_sweep_config["hist_bins"]),
        str(single_sweep_config["step_size"]).split(".")[-1],
        int(single_sweep_config["samples"]),
        str(single_sweep_config["box_size"]).replace(".", "_"),
    )
    return sweep_string


def sweep_config_to_folders_and_base_confs(sweep_parameters):
    base_config = sweep_parameters["base_config"]
    sweep_configs = []
    for sample in sweep_parameters["samples"]:
        single_sweep_config = base_config.copy()
        single_sweep_config["samples"] = sample
        sweep_configs.append(single_sweep_config)

    for step in sweep_parameters["step_size"]:
        single_sweep_config = base_config.copy()
        single_sweep_config["step_size"] = step
        sweep_configs.append(single_sweep_config)

    for hist in sweep_parameters["hist_bins"]:
        single_sweep_config = base_config.copy()
        single_sweep_config["hist_bins"] = hist
        sweep_configs.append(single_sweep_config)

    for box in sweep_parameters["box_size"]:
        single_sweep_config = base_config.copy()
        single_sweep_config["box_size"] = box
        sweep_configs.append(single_sweep_config)

    sweep_folders = [config_to_folder(config) for config in sweep_configs]

    return sweep_folders, sweep_configs
