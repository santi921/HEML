import os 

atomic_radii = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.8,
    'S': 1.8,
    'CL': 1.75,
    'BR': 1.85,
    'I': 1.98,
    'FE': 1.80,
    'CU': 1.40,
    'ZN': 1.39,
    'MG': 1.73,
    'CA': 2.00,
    'NA': 2.27,
    'K': 2.75,
    'LI': 1.82,
    'AL': 1.84,
    'MN': 1.73,
    'CO': 1.70
}


def mol2_to_pqr(file):
    with open(file, 'r') as f:
        lines_mol2 = f.readlines()
    lines_pqr = []
    start_ind = -1
    for ind, line in enumerate(lines_mol2):

        if start_ind > 0: 
            if line.startswith('@<TRIPOS>BOND'):
                break
            else:
                line_split = line.split()
                #print(line)
                atom_id = line_split[0]
                atom_name = line_split[1]
                x = line_split[2]
                y = line_split[3]
                z = line_split[4]
                atom_type = line_split[5]
                subst_id = line_split[6]
                subst_name = line_split[7]
                charge = line_split[8]
                radius = atomic_radii[atom_type]
                chain = 'A'
                
                # just for shobhit's files
                if subst_name == "HM1": 
                    subst_name = "HEM"
                if subst_name == "FE1": 
                    subst_name = "HEM"

                # with element at end 
                #line_pqr = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format('ATOM', int(atom_id), atom_name, ' ', subst_name, chain, int(subst_id), ' ', float(x), float(y), float(z), float(charge), float(radius), atom_type, ' ')
                # without element at end
                line_pqr = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM', int(atom_id), atom_name, ' ', subst_name, chain, int(subst_id), ' ', float(x), float(y), float(z), float(charge), float(radius))
                lines_pqr.append(line_pqr)
        
        if line.startswith('@<TRIPOS>ATOM'):
            start_ind = ind

    with open(file.replace('.mol2', '.pqr'), 'w') as f:
        f.writelines(lines_pqr)


def mol2_to_pqr_folder(dir): 
    for file in os.listdir(dir):
        if file.endswith('.mol2'):
            mol2_to_pqr(os.path.join(dir, file))

