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
    # find all files ending in *pdb
    files = glob("{}*.dat".format(target_folder))
    fields = []
    # loop over all files
    for ind, file in enumerate(files):
        if ind == 0:
            dict_meta_data = mat_pull(file, meta_data=True)
        # read in the file
        field = mat_pull(file)
        
        # append the field to the list of fields
        fields.append(field)


    fields = np.array(fields)
    print(fields.shape)
    # average the fields over axis 1, 2, and 3
    average_field = np.mean(fields, axis=0)
    # save the average field in the target folder
    np.save("{}average_field.npy".format(target_folder), average_field)
    # save mean data in the target folder in original format
    save_numpy_as_dat(dict_meta_data, average_field, "{}average_field.dat".format(target_folder))
    test_save_mat = mat_pull("{}average_field.dat".format(target_folder))
    # print rmse between the two
    print(np.sqrt(np.mean((average_field - test_save_mat)**2)))
    # assert that the saved field is the same as the average field
    assert np.allclose(average_field, test_save_mat,atol=1.e-4)

main()