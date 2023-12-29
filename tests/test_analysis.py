import os
from HEML.utils.fields import compress
from HEML.utils.cpet import make_histograms, construct_distance_matrix


# Test clustering analysis
def test_dist_and_compress():
    # test distance matrix
    # test compress dictionary
    dampening = 0.5
    max_iter = 1000
    compress_tf = True
    output_folder = "./test_data/dataset/compressed_out/"
    path_target = "./test_data/dataset/topo/"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    topo_files = [
        path_target + f for f in os.listdir(path_target) if f.endswith(".top")
    ]
    print(
        'number of topologies in folder "{}": {}'.format(path_target, len(topo_files))
    )
    topo_files.sort()

    with open(output_folder + "topo_file_list.txt", "w") as file_list:
        for i in topo_files:
            file_list.write(f"{i} \n")

    histograms = make_histograms(topo_files)
    distance_matrix = construct_distance_matrix(histograms)
    print(distance_matrix.shape)

    assert distance_matrix.shape == (34, 34)
    # construct distance matrix
    with open(output_folder + "distance_matrix.dat", "w") as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f"{col} ")
            outputfile.write("\n")

    compress_dictionary = compress(
        distance_matrix, damping=dampening, max_iter=int(max_iter), names=topo_files
    )
    print(compress_dictionary.keys())
    keys_nominal = [
        "boundary_inds",
        "silhouette",
        "labels",
        "n_clusters",
        "total_count",
    ]
    keys_cluster = ["0", "1"]

    for k in keys_nominal:
        assert k in compress_dictionary.keys(), "nominal keys not found correctly"

    for k in keys_cluster:
        assert k in compress_dictionary.keys(), "clusters not found correctly"
