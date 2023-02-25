import glob


def main():

    charge_files = glob.glob("../../../data/pdbs/*")
    list = []
    diag_list = []

    for i in charge_files:
        with open(i) as f:
            lines = f.readlines()
            for line in lines:
                if line.split()[0] == "HETATM":

                    item = line.split()[3]
                    len_res = len(item)

                    if len_res > 4:
                        # print(item)
                        item_temp = item
                        item = line.split()[2][-4:]
                        len_res = len(item)

                    if len_res == 1 and len(line.split()[2]) > 6:
                        item = line.split()[2][-4:]
                        len_res = len(item)

                    if len_res != 3 and len_res != 4:  # debug
                        diag_list.append(item)

                    # if(len(line.split()) == 11):
                    #    print(line.split())
                    #    print(line[0:30])

                    inclusion = item in list

                    if not inclusion:
                        list.append(item)

    list.sort()
    diag_list.sort()
    list = sorted(set(list))
    print(sorted(set(diag_list)))

    standard_res = [
        "HEM",
        "AHEM",
        "BHEM",
        "CHEM",
        "DHEM",
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "PYL",
        "SEC",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "DAL",
        "DAR",
        "DSG",
        "DAS",
        "DCY",
        "DGN",
        "DGL",
        "DHI",
        "DIL",
        "DLE",
        "DLY",
        "MED",
        "DPN",
        "DPR",
        "DSN",
        "DTH",
        "DTR",
        "DTY",
        "DVA",
    ]

    for thing in standard_res:
        try:
            list.remove(thing)
        except:
            pass

    with open("../../../data/het_list.txt", "w") as f:
        for item in list:
            f.write("%s\n" % item)


main()
