import os, json
from chimera import runCommand as rc
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)


def addh(pdb_file):
    """
    Add hydrogens to the pdb files and overwrite the olds files.
    """
    print(pdb_file)
    print("open " + pdb_file)

    rc("open " + pdb_file)
    rc("delete solvent")
    rc("addh")
    rc("write #0 " + pdb_file[:-4] + "_h.pdb")
    print("write #0 " + pdb_file[:-4] + "_h.pdb")
    rc("close session")


def get_options(options_file="./options/options.json"):
    """
    Get options from options.json file and create folders if they don't exist.
    Takes
        options_file: path to options.json file
    Returns
        options: dictionary of options
    """
    with open(options_file) as f:
        options = json.load(f)
    for key in options:
        if "folder" in key:
            if not os.path.exists(options[key]):
                os.makedirs(options[key])

    return options


def main():

    options = get_options("./options/options.json")
    root = options["compressed_proteins_folder"]

    for ind, protein_name in enumerate(os.listdir(root)):
        if os.path.isdir(os.path.join(root, protein_name)):
            print(protein_name)
            try:
                folder_name = root + protein_name
                # check if the protein has been submitted to sbatch

                # add h to pdb
                addh("{}/{}_heme.pdb".format(folder_name, protein_name))
                addh("{}/{}_oh_heme.pdb".format(folder_name, protein_name))
                addh("{}/{}_o_heme.pdb".format(folder_name, protein_name))

                # convert pdb back to xyz
                os.system(
                    "obabel -i pdb {}/{}_heme_h.pdb -o xyz -O {}/{}_heme_h.xyz".format(
                        folder_name, protein_name, folder_name, protein_name
                    )
                )
                os.system(
                    "obabel -i pdb {}/{}_oh_heme_h.pdb -o xyz -O {}/{}_oh_heme_h.xyz".format(
                        folder_name, protein_name, folder_name, protein_name
                    )
                )
                os.system(
                    "obabel -i pdb {}/{}_o_heme_h.pdb -o xyz -O {}/{}_o_heme_h.xyz".format(
                        folder_name, protein_name, folder_name, protein_name
                    )
                )

                print("done with {} of {}".format(ind, len(os.listdir(root))))

            except:
                print("error with {}".format(protein_name))
                continue


main()
