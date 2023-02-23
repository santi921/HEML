from glob import glob


def main():
    files = glob("./*pdb")

    for file in files:
        ind = 1
        print(file)
        f = open(str(0) + "_" + file.split("/")[-1], "w")

        with open(file) as f_old:
            lines = f_old.readlines()
        for line in lines:
            if "MODEL" in line:
                f.write("MODEL     1\n")
                print(line)
            elif "ENDMDL" in line:
                print("split")
                ind += 1
                f.close()
                f = open(str(ind) + "_" + file.split("/")[-1], "w")
            else:
                f.write(line)


main()
