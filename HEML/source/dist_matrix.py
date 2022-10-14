
#!/usr/bin/env python3

import os
import matplotlib
import timeit
import matplotlib.pyplot as plt
import numpy as np

# Customize matplotlib
matplotlib.rcParams.update({  # Use mathtext, not LaTeX
    'text.usetex': False,
    # Use the Computer modern font
    'font.family': 'serif',
    'font.serif': 'cmr10',
    'mathtext.fontset': 'cm',
    # Use ASCII minus
    'axes.unicode_minus': False,
    'font.size': 16
})


def make_histograms(topo_files):
    histograms = []

    for topo_file in topo_files:

        distances = []
        curvatures = []

        with open(topo_file) as topology_data:
            for line in topology_data:
                if line.startswith("#"):
                    continue

                line = line.split(",")
                distances.append(float(line[0]))
                curvatures.append(float(line[1]))

        # bins is number of histograms bins in x and y direction (so below is 100x100 bins)
        # range gives xrange, yrange for the histogram
        a, b, c, q = plt.hist2d(distances, curvatures, bins=200, range=[[0, 10], [0, 30]], norm=matplotlib.colors.LogNorm(), density=True, cmap='jet')

        NormConstant = 0
        for j in a:
            for m in j:
                NormConstant += m

        actual = []
        for j in a:
            actual.append([m/NormConstant for m in j])

        actual = np.array(actual)
        histograms.append(actual.flatten()) 
        #plt.show()

    return np.array(histograms)

def distance_numpy(hist1, hist2):
    a = (hist1-hist2)**2
    b = hist1+hist2
    return np.sum(np.divide(a, b, out=np.zeros_like(a), where=b!=0))/2.0

def construct_distance_matrix(histograms):
    matrix = np.diag(np.zeros(len(histograms)))
    for i, hist1 in enumerate(histograms):
        for j,hist2 in enumerate(histograms[i+1:]):
            j += i+1
            matrix[i][j] = distance_numpy(hist1, hist2)
            matrix[j][i] = matrix[i][j]

    return matrix

def main():
    topo_files = [f for f in os.listdir(".") if f.endswith(".top")]
    # sorts the files in some way
    topo_files.sort(key=lambda i: i.split("_")[0])
    with open('topo_file_list.txt', 'w') as file_list:
        for i in topo_files:
            file_list.write(f'{i} \n')
    histograms = make_histograms(topo_files)


    distance_matrix = construct_distance_matrix(histograms)
    print(distance_matrix)
    plt.matshow(distance_matrix, cmap='gray', vmin=0, vmax=0.25)
    plt.show()
    
    with open("distance_matrix.dat", 'w') as outputfile:
        for row in distance_matrix:
            for col in row:
                outputfile.write(f'{col} ')
            outputfile.write("\n")

if __name__ == "__main__":
    main()
