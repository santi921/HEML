import os, json
from turtle import distance 
from HEML.utils.cpet import make_histograms, construct_distance_matrix
from sklearn.cluster import AffinityPropagation

def compress(distance_matrix):


    compressed_dictionary = {}
    affinity = AffinityPropagation(affinity='precomputed')
    affinity.fit(distance_matrix)
    cluster_centers_indices = affinity.cluster_centers_indices_
    labels = affinity.labels_
    n_clusters_ = len(cluster_centers_indices)

    print(f'Estimated number of clusters: {n_clusters_}')
    #get count of a value in a list 
    for i in range(n_clusters_):
        compressed_dictionary[i] = {"count":labels.count(i), "index_center" : cluster_centers_indices[i]}
    return compressed_dictionary


def main():
    root = "/ocean/projects/che160019p/santi92/cpet/"
    # get list of folders in directory specified by user
    folders = [f for f in os.listdir(root) if os.path.isdir(os.path.join(root, f))]

    # for each folder, run the distance mat calculations
    for folder in folders:
        topo_files = [f for f in os.listdir(os.path.join(root, folder)) if f.endswith(".top")]
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
            compress_dictionary[k]["name_center"] = topo_files[v["index_center"]]
        
        # save dinctionary as json file
        with open(os.path.join(root, folder) + "/compressed_dictionary.json", 'w') as outputfile:
            json.dump(compress_dictionary, outputfile)

main()