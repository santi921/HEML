from glob import glob


def main():
    # find all files ending in *pdb
    files = glob("*.pdb")
    
    # loop over all files
    for file in files:
        index_cut = -1
        # open the file
        with open(file, "r") as f:
            # read the file
            lines = f.readlines()
            # close the file
            f.close()

        for i, line in enumerate(lines):
            if line.startswith("TER"):
                index_cut = int(i)
                break

        # overwrite the file with only the lines before the first "TER"
        if index_cut > 0:
            with open(file, "w") as f:
                f.writelines(lines[:index_cut])
                f.close()


main()
