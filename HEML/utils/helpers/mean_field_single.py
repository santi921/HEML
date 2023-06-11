import argparse
import numpy as np
from glob import glob
from HEML.utils.data import mat_pull
from HEML.utils.fields import save_numpy_as_dat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
       "--target_folder", help="location of folder with fields", default="./"
    )
    target_folder = parser.parse_args().target_folder
    # for every *dat in folder compute the mean field 
    dat_files = glob("{}*.dat".format(target_folder))
    fields = []
    for field_file in dat_files:
        # find all files ending in *pdb
        # loop over all files
        dict_meta_data = mat_pull(field_file, meta_data=True)
        handle = field_file.split("/")[-1].split(".")[0]
        # read in the file
        field = mat_pull(field_file)
        # append the field to the list of fields
        fields.append(field)
        mean_field= np.mean(field, axis=(0,1,2))
        # save the average field in the target folder
        np.save("{}/average_{}.npy".format(target_folder, handle), mean_field)
        # save mean data in the target folder in original format
    fields = np.array(fields)
    average_field = np.mean(fields, axis=(0, 1, 2, 3))
    
    # now reiterae of fields and compare to average field
    for ind, field in enumerate(fields): 
        difference = np.mean(field, axis=(0,1,2, 3)) - average_field        
        header = dat_files[ind].split("/")[-1].split(".")[0]
        print(header, difference)
main()

