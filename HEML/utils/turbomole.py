import os, json
import numpy as np 
from turbomoleio.input.define import DefineRunner
from turbomoleio.input.utils import get_define_template, validate_parameters
from monty.os import cd


def get_elements(file_name): 
    """
    Open xyz and get all the elements in the file
    Takes:
        file_name: the name of the xyz file
    Returns:
        elements: a list of elements
    """
    elements = []
    with open(file_name, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line[0:2].strip().isalpha():
            elements.append(line.split()[0])
    return elements


def setup_turbomole(folder_name):
    """
    Sets up the turbomole calculation.
    Takes:
        folder: the folder where the calculation should be set up
    """
    os.chdir(folder_name)
    os.system('setupturbomole.py -t')
    os.chdir("../../..")
    write_sbatch(folder_name)


def submit_turbomole(folder_name, check_if_done = True, t = 24, n = 4, job_name = "o"):
    os.chdir(folder_name)   
    # check if file called GEO_OPT_CONVERGED exists
    
    if check_if_done:
        if os.path.exists("GEO_OPT_CONVERGED"):
            os.chdir("../../..")
            print("calculation is complete, not resubmitting")
            return
        
    write_sbatch(
        "./", 
        time = t, 
        cpus=n, 
        steps = 1000, 
        ri = True, 
        rij = True, 
        conv_crit = 5, 
        gcart = 3, 
        submit_tf = True,
        name = job_name, 
        user = "santi92", keep=False)
    
    os.chdir("../../..")


def define_turbomoleio(
    folder_name, 
    check_if_done = True,
    atoms_present = [],
    charge = 0, 
    spin = 0
    ):
    
    if check_if_done:
        if os.path.exists("GEO_OPT_CONVERGED"):
            os.chdir("../../..")
            return

    dp = get_dictionary(atoms_present, charge, spin=spin)
    timeout=15
    log_filepath = folder_name + "turbomoleio.log"
    workdir = folder_name

    print("valid turboio parameters: {}".format(validate_parameters(dp)))

    dr = DefineRunner(
        log_filepath=log_filepath,
        workdir=workdir,
        timeout=timeout,
        parameters=dp
    )
    
    dr._set_metric()

    with cd(dr.workdir):
        with open(dr.log_filepath, "wb") as logfile:
            dr.define = dr._spawn(
                            dr._get_bin_path(), timeout=dr.timeout, logfile=logfile
                )
            #print(os.listdir())
            dr._initialize_control()
            dr._geometry_menu(new_coords=True)
            dr._switch_to_atomic_attribute_menu()
            dr._define_basis_sets()
            dr._switch_to_molecular_orbital_definition_menu()
            dr._extended_hueckel_theory_menu()
            dr._set_dft_options(use_dft=True) # uses functional, gridsize keys 
            dr._set_scf_options()
            dr._set_ri_state()
            dr._quit_general_menu()
            dr._postprocess()
            case = dr._expect(
                                ["define ended normally", "define ended abnormally"],
                                action="check end of define",
                            )

            if case == 0:
                ended_normally = True

    
def get_dictionary( atoms_present = [], charge = 0, spin=0):

    basic_dict = {        
        "ired": False, 
        "desy": True,
        "method": "dft",
        "functional" : "tpss",
        "gridsize": "m4",
        "scfconv": 4,
        "basis": "def2-SV(P)",
        "scfiterlimit": 1000,
        "marij": True,
        "ri": True, 
        "rijk": False,
        "charge": 0,
        "disp": "DFT-D3",
        "use_cosmo": True,
        "epsilon": 4,
        "coord_file": "coord"
    }

    if atoms_present == []:
        basic_dict["basis_atom"] = {
            "fe": "def2-TZVP",
            "n": "def2-TZVP",
            "s": "def2-TZVP",
            "o": "def2-TZVP"
        }
    else:
        atoms_set = list(set(atoms_present))
        basic_dict["basis_atom"] = {}
        for atom in atoms_set:
            if atom == "Fe" or atom == "fe":
                basic_dict["basis_atom"]["fe"] = "def2-TZVP"
            elif atom == "N" or atom == "n":
                basic_dict["basis_atom"]["n"] = "def2-TZVP"
            elif atom == "S" or atom == "s":
                basic_dict["basis_atom"]["s"] = "def2-TZVP"
            elif atom == "O" or atom == "o":
                basic_dict["basis_atom"]["o"] = "def2-TZVP"
            else: 
                basic_dict["basis_atom"][atom.lower()] = "def2-SV(P)"

    if charge != 0:
        basic_dict["charge"] = charge
        basic_dict["unpaired_electrons"] = spin
        #if charge == -2: 
        #    #basic_dict["open_shell"]["open_shell_on"] = True
        #    basic_dict["unpaired_electrons"] = 1
       
    #if frozen_atoms != []:
    #    basic_dict['freeze_atoms'] = frozen_atoms
    #print(basic_dict)
    return basic_dict


def write_json(folder, frozen_atoms = [], atoms_present = [], charge = 0): 
    """
    Writes a json file for the turbomole calculation. 
    Takes: 
        folder: the folder where the json file should be written
        frozen_atoms: a list of atoms that should be frozen in the calculation
    """
    if not os.path.exists(folder):
        os.makedirs(folder)

    basic_dict = get_dictionary(frozen_atoms, atoms_present, charge)
    
    with open(folder + "definput.json", "w") as outfile:
        json.dump(basic_dict, outfile, indent = 4)


def check_submitted(folder):
    """
    Checks if the protein has been submitted to sbatch
    Takes:
        folder: the folder of the protein
    Returns:
        True if the protein has been submitted, False otherwise
    
    """

    # check if there's a .sh file and a slurm* file 
    sh_file = False
    slurm_file = False

    for file in os.listdir(folder):
        if file.endswith(".sh"):
            sh_file = True
        if file.startswith("slurm-"):
            slurm_file = True
    
    if sh_file and slurm_file:
        return True
    return False


def clean_up(folder, filter=None, clear_control_tf = False):
    # check if there's a file with name filter in the folder
    # if there is, remove all the files in the folder
    for file in os.listdir(folder):
        if filter is not None:
            if file.endswith(filter):        
                # remove every file that isn't a .sh file or a slurm* file or coord file 
                for file in os.listdir(folder):
                    if not file.endswith(".sh") and not file.startswith("slurm-") and not file.endswith("coord") and not file.endswith("control"):
                        os.remove(os.path.join(folder, file))
                    if clear_control_tf:
                        if file.endswith("control"):
                            os.remove(os.path.join(folder, file))
        else:
            # remove every file that isn't a .sh file or a slurm* file or coord file 
            for file in os.listdir(folder):
                if not file.endswith(".sh") and not file.startswith("slurm-") and not file.endswith("coord") and not file.endswith("control"):
                    os.remove(os.path.join(folder, file))
                if clear_control_tf:
                    if file.endswith("control"):
                        os.remove(os.path.join(folder, file))


def write_sbatch(
        folder, 
        time = 24, 
        cpus=4, 
        steps = 100, 
        ri = False, 
        rij=False, 
        conv_crit = 6, 
        gcart = 3, 
        submit_tf = False, 
        user = "santi92", 
        name='normal', 
        keep=False):
    """
    Writes a launch sbatch file for the protein
    Takes:
        folder: the folder of the protein
        submit_tf: if True, submits the sbatch file
        time: the time in hours that the protein should be run for
        cpus: the number of cpus that should be used
        steps: the number of steps that should be run - geo opt 
        ri: if True, runs RI
        rij: if True, runs RIJ
        conv_crit: the convergence criterion for the energy
        gcart: the convergence criterion for the geometry
        user: the user that is running the calculation
        name: the name of the job
        keep: if True, keeps the files after the calculation is done
    """

    # write the sbatch file 
    jobex_command = "jobex  -c {} -np {} -c".format(steps, cpus, time)
    if conv_crit:
        jobex_command += " -energy {}".format(conv_crit)
    if ri:
        jobex_command += " -ri"
    if rij:
        jobex_command += " -rij"
    if keep: 
        jobex_command += " -keep"
    if gcart:
        jobex_command += " -gcart {}".format(gcart)
    
    jobex_command += "\n"

    with open(folder + "launch.sh", "w") as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("#SBATCH -p RM-shared\n")
        outfile.write("#SBATCH --mem=10000\n")
        outfile.write("#SBATCH -t {}:00:00\n".format(time))
        outfile.write("#SBATCH -J {}\n".format(name))
        outfile.write("#SBATCH --no-requeue\n")
        outfile.write("#SBATCH -N 1\n")
        outfile.write("#SBATCH --ntasks-per-node {}\n".format(cpus))
        outfile.write("#SBATCH --mail-user={}\n".format(user))
        outfile.write("#SBATCH --mail-type=FAIL\n")
        outfile.write("echo \"\n\tSettings MKL_DEBUG_CPU_TYPE=5\n\"\n")
        outfile.write("module load mkl\n")
        outfile.write("export MKL_DEBUG_CPU_TYPE=5\n")
        outfile.write("qqdir=$SLURM_SUBMIT_DIR")
        outfile.write("cd $qqdir\n")
        outfile.write("echo \"\n\tSUBMIT DIRECTORY = $qqdir\n\tCORES = {} \n\tTIME = {} \n\tSCRATCH = $LOCAL \n\tDATE = `date`\n\"\n".format(cpus, time))
        outfile.write("echo \"\n\tAdding OpenBabel to Path\n\"\n")
        outfile.write("export PATH=/ocean/projects/che160019p/shared/openbabel-2.4.0/bin/:$PATH\n")
        #outfile.write("unalis -a\n")
        outfile.write("export PYTHONUNBUFFERED=1\n")
        # increase memory limit - 16 gb
        #outfile.write("ulimit -v 4000000\n")
        #outfile.write("time /ocean/projects/che160019p/santi92/phd3/phd3/bin/runturbomole.py -n {} -t {}\n".format(cpus, time))
        outfile.write(jobex_command)
        outfile.write("exit 0\n")

    if submit_tf:
        os.system("sbatch " + folder + "launch.sh")
        

def add_frozen_atoms(folder, frozen_atoms):
    """
    Adds the frozen atoms to the coord file at corresponding positions 
    Takes:
        folder: the folder of the protein
        frozen_atoms: a list of frozen atoms
    
    """  
    # iterate through atom and add " f" to end of line corresponding to frozen atom line 
    with open(folder + "coord", "r") as infile:
        lines = infile.readlines()
        for atom in frozen_atoms:
            lines[atom] = lines[atom].rstrip("\n") + " f\n"
            
    # write the new coord file
    with open(folder + "coord", "w") as outfile:
        outfile.writelines(lines)


def get_fe_positions(file):
    fe_ID, fe_xyz = None, None
    with open(file, 'r') as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if ('FE' in line or "Fe" in line):
            fe_ID = "0"
            fe_xyz = [line[1], line[2], line[3]]
            fe_xyz = [float(x) for x in fe_xyz]
            fe_xyz = np.array(fe_xyz)
            break

    return {"id": fe_ID, "xyz": fe_xyz}


def get_cross_vector(file_name): 

    # find the four nitrogens closest to the iron
    fe_info = get_fe_positions(file_name)
    fe_xyz = fe_info["xyz"]
    fe_ID = fe_info["id"]
    nitrogen_info = get_N_positions(file_name, fe_xyz)


    
    direction_1 = nitrogen_info["N1_xyz"] - nitrogen_info["mean_N_xyz"]
    direction_2 = nitrogen_info["N2_xyz"] - nitrogen_info["mean_N_xyz"]
    direction_3 = nitrogen_info["N3_xyz"] - nitrogen_info["mean_N_xyz"]
    # compute cross and take the two most orthogonal directions
    dot_12 = np.dot(direction_1, direction_2)
    dot_13 = np.dot(direction_1, direction_3)
    dot_23 = np.dot(direction_2, direction_3)

    if(dot_23 > dot_13 and dot_23 > dot_12):
        direction_1 = direction_3
    if(dot_13 > dot_12 and dot_13 > dot_23):
        direction_2 = direction_3


    direction_1 /= np.linalg.norm(direction_1)
    direction_2 /= np.linalg.norm(direction_2)
    cross = np.cross(direction_1, direction_2)
    cross = cross / np.linalg.norm(cross)
    return cross 


def get_carbon_xyz_from_file(file_name):
    """
    get the xyz coordinates of the carbons in the file
    Takes:
        file_name: the name of the xyz file
    Returns:
        carbon_xyz: a list of the xyz coordinates of the carbons
    """
    carbon_xyz, ind = [], []
    with open(file_name, 'r') as f:
        lines = f.readlines()
        # go through all the lines and find nitrogens
        for line_ind, line in enumerate(lines):
            if line[0:2].strip().isalpha():
                if line.split()[0] == "C" or line.split()[0] == "c":
                    carbon_xyz.append(
                        [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
                        )
                    ind.append(line_ind-1)
    return carbon_xyz, ind


def get_N_positions(file, fe_xyz):

    N_xyz_list = []
    with open(file, 'r') as f:
        readfile = f.readlines()

    for j in readfile:
        line = j.split()
        if ('N' in line or "n" in line):
            N_xyz_list.append([float(line[1]), float(line[2]), float(line[3])])

        
    # get distances of N atoms to the Fe atom
    N_dist = [np.linalg.norm(np.array(i) - np.array(fe_xyz)) for i in N_xyz_list]
    # get the indices of the N atoms closest to the Fe atom
    N_indices = np.argsort(N_dist)[:4]
    # get the xyz coordinates of the N atoms closest to the Fe atom
    N_xyz = [N_xyz_list[i] for i in N_indices]
    # get distance of N atoms to the Fe atom
    N_dist = [N_dist[i] for i in N_indices]

    
    # assign the N_IDs
    mean_N_xyz = np.mean(np.array(N_xyz), axis=0)

    nitrogen_dict = {
        "mean_N_xyz": mean_N_xyz,
        "N_ID1": N_indices[0],
        "N_ID2": N_indices[1],
        "N_ID3": N_indices[2],
        "N_ID4": N_indices[3],
        "N1_xyz": N_xyz[0],
        "N2_xyz": N_xyz[1],
        "N3_xyz": N_xyz[2],
        "N4_xyz": N_xyz[3]
    }

    return nitrogen_dict 


def get_frozen_atoms(file_name):
    """
    get the two carbons most out of the plane to freeze
    Takes:
        file_name: the name of the xyz file
    Returns:
        frozen_atoms: a binary list of frozen atoms
    """

    
    carbon_xyz, ind_carbons = get_carbon_xyz_from_file(file_name)
    
    cross = get_cross_vector(file_name)
    fe_dict = get_fe_positions(file_name)
    n_dict = get_N_positions(file_name, fe_dict["xyz"])
    
    mean_xyz = n_dict["mean_N_xyz"]
    dot_list = [np.dot(i[0] - n_dict["mean_N_xyz"], cross) for i in carbon_xyz]
    dot_list = np.array(dot_list) / np.linalg.norm(cross)

    # filter for coplanar carbons    
    carbon_planar_ind = [] # the indices of the coplanar carbons
    for ind, i in enumerate(carbon_xyz):
        if dot_list[ind] < 0.5:
            carbon_planar_ind.append(ind)

    # get the four furthest, in plane carbons
    carbon_planar_xyz = np.array(carbon_xyz)[carbon_planar_ind]
    distances = np.linalg.norm(carbon_planar_xyz - mean_xyz, axis = 1)
    
    furthest_ind = np.argsort(distances)[::-1][:4]
    most_out_of_plane_ind = np.argsort(dot_list)[::-1][:2]
    # combine the two lists
    frozen_atom_ind = np.concatenate((furthest_ind, most_out_of_plane_ind))
    return_list = [ind_carbons[i] for i in frozen_atom_ind]

    return return_list


def fetch_charges_dict(file_name = 'test.pqr'):
    """
    Given a list of dictionaries with element and position, traverse a pqr file and get the charges from the file
    EXCLUDING ELEMENTS IN THE LIST OF DICTIONARIES
    Takes: 
        list of dictionaries with element and position
    Returns:
        list of dictionaries with element, position and charge
    """
    
    pqr_dict = []
    # get the lines of the pqr file
    with open(file_name, 'r') as f:
        lines = f.readlines()
    
    for line in lines: 
        x       = float(line[30:38].strip())
        y       = float(line[39:46].strip())
        z       = float(line[47:54].strip())
        charge  = float(line[55:61].strip())
        radius  = float(line[62:68].strip())    
        if np.abs(charge) >= 0.01:
            pqr_dict.append({"position": [x,y,z], "charge": charge, "radius": radius})

    return pqr_dict
