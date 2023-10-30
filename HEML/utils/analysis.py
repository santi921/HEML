import numpy as np
from copy import deepcopy


def get_run_number(file, run_key="run"):
    # get the run number
    run_number = file.split("-")
    # get the string with run_key in it
    run_number = [i for i in run_number if run_key in i]
    # remove the run_key
    run_number = run_number[0].replace(run_key, "")
    # convert to integer
    run_number = int(run_number)
    return run_number


def compute_resonance_times(
    files, full_files, run_key="run", split_by_run=False, summary_stats=True
):
    """compute the resonance times for a given cluster, separated by run_key.
    Then compute the mean and number of samples in each run
    Takes:
        files: list of files in a cluster
        full_files: list of all files
        run_key: string that separates the run number from the rest of the file name
        split_by_run: if True, returns a dictionary with the same keys as res_dict,
            but with the values being lists of resonance times for each run
        summary_stats: if True, returns means, max, and number of samples in each run. If False,
            returns the full list of resonance times for each run
    Returns:
        res_dict: dictionary with mean, number of samples, and max resonance time
        if split_by_run is True, then returns a dictionary with the same keys as res_dict,
    """

    if len(files) <= 1:  # esp useful with filtered_resonance_analysis
        return {"mean": 1, "n_entries": 0, "max": 1}
    res_dict = {}
    # print("files: ", len(files))
    res_time_temp = 0
    for ind_file_cluster, file in enumerate(files[:-1]):
        run_number = get_run_number(file, run_key=run_key)
        next_file = files[ind_file_cluster + 1]
        next_file_run_number = get_run_number(next_file, run_key=run_key)
        # get the index of the next_file in full_files
        ind_next_file = full_files.index(next_file)
        ind_current_file = full_files.index(file)

        res_time_temp += 1
        # check if run_number is in res_dict
        if run_number not in res_dict.keys():
            res_dict[run_number] = []

        # check if the next file is the next file in full_files
        # print(ind_next_file - ind_current_file)
        if next_file_run_number != run_number:
            # if not, append the res_time_temp to the list
            res_dict[run_number].append(res_time_temp)
            # reset the res_time_temp
            res_time_temp = 0

        elif ind_next_file - ind_current_file != 1:
            # if so, append the res_time_temp to the list
            res_dict[run_number].append(res_time_temp)
            # reset the res_time_temp
            res_time_temp = 0

        elif ind_file_cluster == len(files) - 2:
            # if so, append the res_time_temp to the list
            res_dict[run_number].append(res_time_temp)

    if summary_stats:
        # compute mean and number of samples in each run
        for k, v in res_dict.items():
            res_dict[k] = {
                "mean": float(np.mean(v)),
                "n_entries": int(len(v)),
                "max": int(np.max(v)),
            }

        if split_by_run:
            return res_dict
        else:
            res_dict_aggregate = {"mean": [], "n_entries": [], "max": []}
            for k, v in res_dict.items():
                res_dict_aggregate["mean"].append(float(v["mean"]))
                res_dict_aggregate["n_entries"].append(int(v["n_entries"]))
                res_dict_aggregate["max"].append(int(v["max"]))

            # print("res dict pre: ", res_dict_aggregate)
            res_dict_aggregate["mean"] = float(
                len(files) / np.sum(res_dict_aggregate["n_entries"])
            )

            res_dict_aggregate["n_entries"] = int(
                np.sum(res_dict_aggregate["n_entries"])
            )
            res_dict_aggregate["max"] = int(np.max(res_dict_aggregate["max"]))
            return res_dict_aggregate
    else:
        res_time = []
        # each res_dict merge the lists
        for k, v in res_dict.items():
            res_time.extend(v)
        return res_time


def simple_resonance_analysis(compress_dictionary, run_key="run", raw=False):
    """
    Computes simple statistics on mean and std of resonance times in each cluster
    """
    # first go through each cluster and collect the names of the files and order them
    full_files = []

    # summary_stats = not raw
    for k, v in compress_dictionary.items():
        # check if k is a number
        if k.isnumeric():
            files_single = deepcopy(v["files"])
            full_files.extend(files_single)
    # remove duplicates
    full_files = list(set(full_files))
    # sort the files
    full_files.sort()
    # get the number of files
    n_files = len(full_files)
    print("number of topo files in compressed dict: {}".format(n_files))

    # go through each cluster
    for k, v in compress_dictionary.items():
        # check if k is a number
        if k.isnumeric():
            if raw:
                compress_dictionary[k]["res_times"] = compute_resonance_times(
                    v["files"], full_files, run_key=run_key, summary_stats=False
                )
            else:
                compress_dictionary[k][
                    "resonance_info_simple"
                ] = compute_resonance_times(
                    v["files"],
                    full_files,
                    run_key=run_key,
                    split_by_run=False,
                    summary_stats=True,
                )
                # print the results for this cluster
                print(
                    "cluster {}: {}".format(
                        k, compress_dictionary[k]["resonance_info_simple"]
                    )
                )

    return compress_dictionary


def filtered_resonance_analysis(compress_dictionary, run_key="run", raw=False):
    """
    Computes simple statistics on mean resonance times in each cluster
    """
    assert (
        "boundary_file_names" in compress_dictionary.keys()
    ), "need names of boundary points in dictionary"

    # first go through each cluster and collect the names of the files and order them
    full_files = []

    for k, v in compress_dictionary.items():
        # check if k is a number
        if k.isnumeric():
            files_single = deepcopy(v["files"])
            full_files.extend(files_single)
    # remove duplicates
    full_files = list(set(full_files))
    # sort the files
    full_files.sort()
    # filter out the boundary points
    full_files = [
        i for i in full_files if i not in compress_dictionary["boundary_file_names"]
    ]
    # get the number of files
    n_files = len(full_files)
    print("number of topo files in compressed dict: {}".format(n_files))

    # go through each cluster
    for k, v in compress_dictionary.items():
        # check if k is a number
        if k.isnumeric():
            if raw:
                compress_dictionary[k]["res_times_filtered"] = compute_resonance_times(
                    [
                        i
                        for i in v["files"]
                        if i not in compress_dictionary["boundary_file_names"]
                    ],
                    full_files,
                    run_key=run_key,
                    summary_stats=False,
                )

            else:
                compress_dictionary[k][
                    "resonance_info_filtered"
                ] = compute_resonance_times(
                    [
                        i
                        for i in v["files"]
                        if i not in compress_dictionary["boundary_file_names"]
                    ],
                    full_files,
                    run_key=run_key,
                    split_by_run=False,
                    summary_stats=True,
                )
                # print the results for this cluster
                print(
                    "cluster {}: {}".format(
                        k, compress_dictionary[k]["resonance_info_filtered"]
                    )
                )

    return compress_dictionary
