
#!/usr/bin/env python3

import os
from HEML.utils.cpet import make_histograms, construct_distance_matrix
import matplotlib.pyplot as plt
import matplotlib 

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
