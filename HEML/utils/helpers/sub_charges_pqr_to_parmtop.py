import argparse

# define main
def main():
    
    names_pqr, charges_pqr = [], []
    names_prmtop, charges_prmtop = [], []

    parser = argparse.ArgumentParser(description='Adds charges from PQR to PRMTOP file')
    parser.add_argument('-p', '--PRMTOP_file', help='PRMTOP file', required=True)
    parser.add_argument('-q', '--PQR_file', help='PQR file', required=True)
    parser.add_argument('-o', '--output_file', help='Output PRMTOP file, if none set, will write new prmtop w/name of pqr file', 
                        required=False,
                        default=None)
    args = parser.parse_args()
    
    file_prmtop = args.PRMTOP_file
    file_pqr = args.PQR_file

    # read pqr file 
    with open(file_pqr, 'r') as f:
        lines = f.readlines()

    # read prmtop file
    with open(file_prmtop, 'r') as f:
        lines_prmtop = f.readlines()

    # get charges, names from pqr file
    for line_ind, line in enumerate(lines):
        assert line.startswith('ATOM') or line.startswith('HETATM'), \
            'PQR file is not valid there is a line with a weird start value, check line {}'. format(line_ind)
        charges_pqr.append(float(line[54:61].strip()))
        names_pqr.append(line[12:16].strip())
    
    
    # get charges, names from prmtop file
    charge_block, atom_block = False, False
    for line_ind, line in enumerate(lines_prmtop):
        if line.startswith('%FLAG ATOM_NAME'):
            atom_block = True
            num_atom_per_line = len(lines_prmtop[line_ind+2].split())
            print("Number of atoms per line: {}".format(num_atom_per_line))

        #if lines_prmtop[line_ind-1].startswith("%FLAG CHARGE") and line.startswith("%FORMAT"):
        if line.startswith("%FLAG CHARGE"):
            atom_block = False
            print("Found the charge block")
            charge_block = True

        #if lines_prmtop[line_ind-1].startswith("%FLAG ATOMIC_NUMBER") and line.startswith("%FORMAT"):
        if line.startswith("%FLAG ATOMIC_NUMBER"):
            charge_block = False
            print("Found the atomic number block")
            break 

        if atom_block:
            if not line.startswith('%'):
                
                # chunk line into strings of 4 characters
                atoms = [line[i:i+4].strip() for i in range(0, len(line), 4)]
                # remove empty strings
                atoms = list(filter(None, atoms))

                if 'EPW' in atoms:
                    #print("EPW line found")
                    atoms.remove('EPW')
                    for atom_ind, atom in enumerate(atoms):
                        names_prmtop.append(atom)

                else:
                    for atom_ind, atom in enumerate(atoms):
                        names_prmtop.append(atom)
            
        if charge_block:
            if not line.startswith('%'):
                charges = line.split()
                for charge_ind, charge in enumerate(charges):
                    #save in 5E16.8
                    charges_prmtop.append(float(charge))
      
    
    names_prmtop = names_prmtop[:len(names_pqr)]
    charges_prmtop = charges_prmtop[:len(charges_pqr)]
    assert names_pqr == names_prmtop, "The atom names in the PQR and PRMTOP files do not match"

    # get charges, names from prmtop file
    
    substitute_charge_index = 0 
    charge_block = False
    
    for line_ind, line in enumerate(lines_prmtop):
        if lines_prmtop[line_ind-1].startswith("%FLAG CHARGE") and line.startswith("%FORMAT"):
            print("Found the charge block")
            charge_block = True

        if lines_prmtop[line_ind-1].startswith("%FLAG ATOMIC_NUMBER") and line.startswith("%FORMAT"):
            charge_block = False
            break 
            
        if charge_block:
            if not line.startswith('%'):
                charges = line.split()
                charges = [float(charge) for charge in charges]
                
                #print(charges)
                # substitute charges
                for charge_ind, charge in enumerate(charges):
                    if substitute_charge_index < len(charges_pqr):
                        if charges_pqr[substitute_charge_index] == 0: 
                            charges_pqr[substitute_charge_index]=0.0
                        charges[charge_ind] = charges_pqr[substitute_charge_index]
                        substitute_charge_index += 1

                new_line = ""
                for charge_ind, charge in enumerate(charges):
                    if charge < 0:
                        new_line += " "
                    else: 
                        new_line += "  "
                    new_line += "{:5.8E}".format(charge)
                
                lines_prmtop[line_ind] = new_line + '\n'
    
    # write new prmtop file
    if args.output_file is None:
        output = file_pqr.split('.')[0] + '.prmtop'
    else: 
        output = args.output_file

    with open(output, 'w') as f:
        f.writelines(lines_prmtop)

    # print charges in 5E16.8 format
    #[print('{:5.8f}'.format(charge)) for charge in charges_prmtop]
    


main()