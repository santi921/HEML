import argparse


# define main
def main():
    names_pqr, charges_pqr, positions_pqr = [], [], []
    names_prmtop, positions_pqr, residue_name_pqr, residue_pointer_pqr  = [], [], [], []

    parser = argparse.ArgumentParser(description="Adds charges from PQR to PRMTOP file")
    parser.add_argument("-p", "--PRMTOP_file", help="PRMTOP file", required=True)
    parser.add_argument("-q", "--PQR_file", help="PQR file", required=True)
    parser.add_argument(
        "-o",
        "--output_file",
        help="Output PRMTOP file, if none set, will write new prmtop w/name of pqr file",
        required=False,
        default=None,
    )
    args = parser.parse_args()

    file_prmtop = args.PRMTOP_file
    file_pqr = args.PQR_file

    # read pqr file
    with open(file_pqr, "r") as f:
        lines = f.readlines()

    # read prmtop file
    with open(file_prmtop, "r") as f:
        lines_prmtop = f.readlines()

    # get charges, names, positions from pqr file
    for line_ind, line in enumerate(lines):
        assert line.startswith("ATOM") or line.startswith(
            "HETATM"
        ), "PQR file is not valid there is a line with a weird start value, check line {}".format(
            line_ind
        )
        charges_pqr.append(float(line[54:61].strip()))
        names_pqr.append(line[12:16].strip())
        #x_pos = line[].strip() # todo
        #y_pos = line[].strip() # todo 
        #z_pos = line[].strip() # todo
        #positions_pqr.append([x_pos, y_pos, z_pos])

    # sanitize -0.0 charges to 0.0
    for charge_index, charge in enumerate(charges_pqr):
        if charges_pqr[charge_index] == 0:
            charges_pqr[charge_index] = 0.0

    # get charges, names, positions from prmtop file
    charge_block, atom_block, position_block, residue_pointer_block, residue_name_block = False, False, False, False, False
    
    for line_ind, line in enumerate(lines_prmtop):
        if line.startswith("%FLAG ATOM_NAME"):
            atom_block = True
            num_atom_per_line = len(lines_prmtop[line_ind + 2].split())
            print("Number of atoms per line: {}".format(num_atom_per_line))

        # if lines_prmtop[line_ind-1].startswith("%FLAG CHARGE") and line.startswith("%FORMAT"):
        if line.startswith("%FLAG CHARGE"):
            atom_block = False
            print("Found the charge block in prmtop")
            charge_block = True

        # if lines_prmtop[line_ind-1].startswith("%FLAG ATOMIC_NUMBER") and line.startswith("%FORMAT"):
        if line.startswith("%FLAG ATOMIC_NUMBER"):
            charge_block = False
            print("Found the atomic number block  in prmtop")
            break

        if line.startswith("%FLAG RESIDUE_POINTER"):
            residue_pointer_block = True
            charge_block=
            print("Found residue pointer block in prmtop")
            break
            
        if line.startswith("%FLAG RESIDUE_LABEL"):
            residue_label_block = True
            print("Found residue label block in prmtop")
            break

        if atom_block:
            if not line.startswith("%"):
                # chunk line into strings of 4 characters
                atoms = [line[i : i + 4].strip() for i in range(0, len(line), 4)]
                # remove empty strings
                atoms = list(filter(None, atoms))
                
                if len(atoms) != num_atom_per_line: 
                    print("warning: varying number of atoms per")
                for atom_ind, atom in enumerate(atoms):
                    names_prmtop.append(atom.strip())
                    # print(atom.strip())

        if residue_name_block: 
            if not line.startswith("%"):
                line_res = line.split()
                [residue_name_pqr.append(res_name) for res_name in line_res]
        
        if residue_pointer_block: 
            if not line.startswith("%"):
                line_pointer = line.split()
                [residue_pointer_pqr.append(pointer) for pointer in line_pointer]

    
    # get charges, names from prmtop file
    substitute_charge_index = 0
    full_charge_index = 0
    charge_block = False
    print(
        "Length of residues in prmtop: {}\nLength of Residues in pqr: {}".format(
            len(names_prmtop), len(names_pqr)
        )
    )

    for line_ind, line in enumerate(lines_prmtop):
        if lines_prmtop[line_ind - 1].startswith("%FLAG CHARGE") and line.startswith(
            "%FORMAT"
        ):
            # print("Found the charge block")
            charge_block = True

        if lines_prmtop[line_ind - 1].startswith(
            "%FLAG ATOMIC_NUMBER"
        ) and line.startswith("%FORMAT"):
            charge_block = False
            break

        if charge_block:
            if not line.startswith("%"):
                charges = line.split()
                charges = [float(charge) for charge in charges]

                # substitute charges
                for charge_ind, charge in enumerate(charges):
                    if substitute_charge_index < len(charges_pqr):
                        if names_prmtop[full_charge_index] == "EPW":
                            full_charge_index += 1
                            # print("skipping epw overwrite")]
                            charges[charge_ind] = 0.0
                        else:
                            # print("{} {} {}".format(
                            #    names_prmtop[full_charge_index],
                            #    charges[charge_ind],
                            #    charges_pqr[substitute_charge_index] * 18.2223))

                            charges[charge_ind] = (
                                charges_pqr[substitute_charge_index] * 18.2223
                            )
                            substitute_charge_index += 1
                            full_charge_index += 1

                new_line = ""
                for charge_ind, charge in enumerate(charges):
                    if charge < 0:
                        new_line += " "
                    else:
                        new_line += "  "
                    new_line += "{:5.8E}".format(charge)

                lines_prmtop[line_ind] = new_line + "\n"

    # write new prmtop file
    if args.output_file is None:
        output = file_pqr.split(".")[0] + ".prmtop"
    else:
        output = args.output_file

    with open(output, "w") as f:
        f.writelines(lines_prmtop)

    # print charges in 5E16.8 format
    # [print('{:5.8f}'.format(charge)) for charge in charges_prmtop]


main()
