import argparse
import os
import subprocess
import time
import re 
import psutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inp', help='path to input file')
    parser.add_argument('-s', dest='step', help='value of the step')
    parser.add_argument('-a', dest='A', help='value of A')
    parser.add_argument('-b', dest='B', help='value of B')
    parser.add_argument('-t', dest='transition', help='value of transition')
    parser.add_argument('-p', dest='product', help='value of product')
    parser.add_argument('-email', dest='email', help='email address')
    args = parser.parse_args()

    # Set variables
    inp = args.inp
    step = args.step
    print("running step {}".format(step))
    A = args.A
    B = args.B
    transition = args.transition
    product = args.product
    email = args.email
    substrate = None

    # Read source file
    source_file = open(inp, 'r')
    source_contents = source_file.read()
    source_file.close()

    # iterate over the lines in the source file
    for line in source_contents.splitlines():
        if line.startswith('parsefile'):
            parse_amber_file = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('system'):
            system = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('parm'):
            parm = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('frame'):
            frame = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('trajin'):
            trajin = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('resname'):
            resname = line.split('=')[1].split(" ")[0]
            #remove everything outside the quotes
            resname = resname.split('"')[1]
        elif line.startswith('numberofres'):
            numberofres = line.split('=')[1].split("#")[0].strip()
        elif line.startswith('basis'):
            basis = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('charge'):
            charge = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('unp'):
            unp = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('nodes'):
            nodes = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('tleapinput'):
            tleapinput = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith('substrate'):
            substrate = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith("qmmm_root"):
            qmmm_root = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith("scfiterlimit"):
            scfiterlimit = line.split('=')[1].split(" ")[0].strip()
        elif line.startswith("pdb"):
            pdb=line.split('=')[1].split(" ")[0].strip()
    
    print("donne pasing input file")
    user = os.getenv('USER')
    host = os.uname()[1]
    file_name = tleapinput

    if step == "0":
        
        print("Running Step 0")
        if not os.path.exists('scratch'):
            os.mkdir('scratch')
        
        os.chdir('scratch')
        
        with open('parsefile', 'w') as f:
            f.write(source_contents)

        # Read the lines containing atom types
        with open(file_name) as f:
            atom_type_lines = [line for line in f if '"sp3"' in line]
        
        # Initialize a dictionary to store atom type values
        atom_type_values = {}
        
        # Iterate over the atom type lines and extract the values
        for line in atom_type_lines:
            atom_type = line.split('"')[1]
            value = line.split('"')[3]
            atom_type_values[atom_type] = value
        
        # Print the extracted values
        for atom_type, value in atom_type_values.items():
            print(f"{atom_type}: {value}")
        
        # Read the contents of the parse_amber.tcl file
        #parse_amber_file = "./parse_amber.tcl"
        with open(parse_amber_file) as f:
            content = f.read()
        
        # Make the required substitutions
        for atom_type, atom_value in atom_type_values.items():
            if atom_value == "Fe":
                #print("print Fe")
                content = content.replace("{ return 26 }", f"- {atom_type} {{ return 26 }}")
            elif atom_value == "N":
                #print("sub N")
                content = content.replace("{ return 7 }", f"- {atom_type} {{ return 7 }}")
            elif atom_value == "O":
                #print("sub O")
                content = content.replace("{ return 8 }", f"- {atom_type} {{ return 8 }}")
            elif atom_value == "S":
                #print("sub S")
                content = content.replace("{ return 16 }", f"- {atom_type} {{ return 16 }}")
            elif atom_value == "C":
                #print("sub C")
                content = content.replace("{ return 6 }", f"- {atom_type} {{ return 6 }}")
            else:
                print(f"Atom type {atom_type} not substituted in parse_amber.tcl file.")
                exit(1)
        
        # Write the modified contents back to the file
        with open(parse_amber_file, 'w') as f:
            print("writing to parse_amber.tcl")
            f.write(content)
            #filename = args.pdb
            filename=pdb
            # Use the split method to extract the desired part
            parts = filename.split("_")
            system = parts[2]
            run = parts[3]
            frame = parts[4]
            
            if "." in frame: frame = frame.split(".")[0] 
            
            print(f"system: {system}")
            print(f"run: {run}")
            print(f"frame: {frame}")


        with open(f"{system}_{run}_{frame}.in", "w") as f:
            print("writing to frame.in")
            f.write(f"""parm {parm}
                trajin {trajin} {frame} {frame}
                autoimage
                trajout {system}_{run}_{frame}.inpcrd restart
                trajout {system}_{run}_{frame}.pdb 
                run
                exit
                """)

        with open(f"{system}_{run}_{frame}.in") as f:
            print(f.read())

        with open(f"strip_{system}_{run}_{frame}.in", "w") as f:
            print("writing to water_strip.in")
            f.write(f"""parm {parm}
            trajin {system}_{run}_{frame}.inpcrd
            reference {system}_{run}_{frame}.inpcrd
            strip :Na+,Cl-
            strip !(:{numberofres}<:10.0) outprefix stripped10
            trajout stripped10.{system}_{run}_{frame}.inpcrd restart
            run
            exit
            """)

        with open(f"strip_{system}_{run}_{frame}.in") as f:
            print(f.read())


        with open(f"rc_{system}_{run}_{frame}.in", "w") as f:
            print("writing to rc.in")
            f.write(f"""parm stripped10.*.prmtop
            trajin stripped10.{system}_{run}_{frame}.inpcrd
            trajout rc.pdb
            trajout rc.rst restart
            run
            exit
            """)

        with open(f"rc_{system}_{run}_{frame}.in") as f:
            print(f.read())

        #CREATING RC COMPLEX FILES
        print("running cpptraj to generate pdb, water stripped pdb, and rc files")

        cpptraj_command = f"cpptraj -i {system}_{run}_{frame}.in > {system}_{run}_{frame}.out"
        process = subprocess.Popen(cpptraj_command, shell=True)
        process.wait()
        if os.path.isfile(f"{system}_{run}_{frame}.out") and \
            "Error" not in open(f"{system}_{run}_{frame}.out").read():
            
            print("Generated Frame PDB")
        else:
            print("Cpptraj Error")
            exit()

        cpptraj_command = f"cpptraj -i strip_{system}_{run}_{frame}.in > strip_{system}_{run}_{frame}.out"
        process = subprocess.Popen(cpptraj_command, shell=True)
        process.wait()
        if os.path.isfile(f"strip_{system}_{run}_{frame}.out") and \
              "Error" not in open(f"strip_{system}_{run}_{frame}.out").read():
            print("Generated WaterStripped PDB and prmtop")
        else:
            print("Cpptraj Error")
            exit()

        cpptraj_command = f"cpptraj -i rc_{system}_{run}_{frame}.in > rc_{system}_{run}_{frame}.out"
        process = subprocess.Popen(cpptraj_command, shell=True)
        process.wait()
        if os.path.isfile(f"rc_{system}_{run}_{frame}.out") and "Error" not in open(f"rc_{system}_{run}_{frame}.out").read():
            print("Generated RC files for RC_OPT")
        else:
            print("Cpptraj Error")
            exit()

        print("Copying prmtop file: stripped10.{base}.prmtop to rc.prmtop")
        os.system(f"cp stripped10.*.prmtop rc.prmtop")
        os.system("sed -i '9s/1/0/' rc.prmtop")

        # ------------------- CREATING RC COMPLEX FILES ------------------- 
        os.system(f"nohup cpptraj -i {system}_{run}_{frame}.in > {system}_{run}_{frame}.out &")
        process = os.getpid()
        while process in os.popen("pgrep cpptraj").read().split():
            time.sleep(1)
        if "Error" in open(f"{system}_{run}_{frame}.out").read():
            print("Cpptraj Error")
            exit()
        else:
            print("Generated Frame PDB")
        os.system(f"nohup cpptraj -i strip_{system}_{run}_{frame}.in > strip_{system}_{run}_{frame}.out &")
        process = os.getpid()
        while process in os.popen("pgrep cpptraj").read().split():
            time.sleep(1)
        if "Error" in open(f"strip_{system}_{run}_{frame}.out").read():
            print("Cpptraj Error")
            exit()
        else:
            print("Generated WaterStripped PDB and prmtop")
        os.system(f"nohup cpptraj -i rc_{system}_{run}_{frame}.in > rc_{system}_{run}_{frame}.out &")
        process = os.getpid()
        while process in os.popen("pgrep cpptraj").read().split():
            time.sleep(1)
        if "Error" in open(f"rc_{system}_{run}_{frame}.out").read():
            print("Cpptraj Error")
            exit()
        else:
            print("Generated RC files for RC_OPT")
        os.system(f"cp stripped10.*.prmtop rc.prmtop")
        with open("rc.prmtop", "r+") as f:
            lines = f.readlines()
            lines[8] = lines[8].replace("1", "0")
            f.seek(0)
            f.write("".join(lines))
            f.truncate()


        # ------------------- MAKING QMMM MODEL -------------------
        print("Creating QMMM Model")
        with open(f"QM_MM_{system}_{run}_{frame}.tcl", "w") as f:
            f.write(f"mol load pdb rc.pdb\n")
            if substrate != None:
                f.write(f'atomselect top "same residue as (within 8 of (resname {resname} {substrate}))"\n')
            else: 
                f.write(f'atomselect top "same residue as (within 8 of (resname {resname}))"\n')
            
            f.write(f'atomselect0 num\n')
            f.write(f'atomselect0 writepdb MM_{system}_{run}_{frame}.pdb\n')
            f.write(f'set myfile [open mm_{system}_{run}_{frame}.txt w]\n')
            f.write(f'puts $myfile [atomselect0 list]\n')
            f.write(f'close $myfile\n')
            if substrate != None:
                f.write(f'atomselect top "(resname {resname} and not backbone and not type HA H) or (resname {substrate})"\n')
            else:
                f.write(f'atomselect top "(resname {resname} and not backbone and not type HA H)"\n')
            f.write(f'atomselect1 num\n')
            f.write(f'atomselect1 writepdb QM_{system}_{run}_{frame}.pdb\n')
            f.write(f'atomselect1 writexyz QM_{system}_{run}_{frame}.xyz\n')
            f.write(f'set myfile1 [open qm_{system}_{run}_{frame}.txt w]\n')
            f.write(f'puts $myfile1 [atomselect1 list]\n')
            f.write(f'close $myfile1\n')
            if substrate != None:
                f.write(f'atomselect top "resname {resname} and type HA or resname FE1 or resname HM1 and name NA or resname CB1 and name C1 or resname {substrate} and name C3"\n')
            else: 
                f.write(f'atomselect top "resname {resname} and type HA or resname FE1 or resname HM1 and name NA or resname CB1 and name C1"\n')
            
            f.write(f'set resid [atomselect1 get resid]\n')
            f.write(f'foreach elementid $resid {{dict set tmp $elementid 1}}\n')
            f.write(f'set id [dict keys $tmp]\n')
            f.write(f'set resname [atomselect1 get resname]\n')
            f.write(f'foreach elementname $resname {{dict set tmp2 $elementname 1}}\n')
            f.write(f'set name [dict keys $tmp2]\n')
            f.write(f'set myresidues [open qm_mm_{system}_{run}_{frame}.sh w]\n')
            f.write(f'puts $myresidues "resid=($id)"\n')
            f.write(f'puts $myresidues "resname=($name)"\n')
            f.write(f'puts $myresidues "myresidues=()"\n')
            f.write(f'puts $myresidues "n=${{#resname[@]}}"\n')
            f.write(f'puts $myresidues "for i in $(seq 1 $n);"\n')
            f.write(f'do\n')
            f.write(f'puts $myresidues "myresidues+=${{resname[i-1]}}${{resid[i-1]}}"\n')
            f.write(f'done\n')
            f.write(f'puts $myresidues "cat > myresidues_{system}_{run}_{frame}.dat <<ENDOFFILE"\n')
            f.write(f'puts $myresidues "set res [ pdb_to_res \\"rc.pdb\\"]"\n')
            myresidues_line = "set myresidues  \\[ inlist function=combine residues= \\\\\\\$res sets= {\\\${myresidues\\[*]}} target=QM ]"
            f.write(f'puts $myresidues "{myresidues_line}\n')
            f.write(f'puts $myresidues "ENDOFFILE')
            f.write(f'close $myresidues\n')
            f.write(f'exit\n')


        # print messages
        print(f"Using rc.pdb")
        print(f"Using Residues:{resname}")
        if substrate != None:
            print(f"Using Substrate:{substrate}")

        # run VMD command
        os.system(f"vmd -dispdev text -e QM_MM_{system}_{run}_{frame}.tcl")

        # change file permissions and execute shell script
        os.system(f"chmod +x qm_mm_{system}_{run}_{frame}.sh")
        os.system(f"./qm_mm_{system}_{run}_{frame}.sh")

        # AWK script to increment values by 1
        addone_script = """
        BEGIN{
        RS = " "
        }

        {
        a = $1
        ++a
        printf( "%d " , a )
        }
        """

        # execute AWK script on input file and redirect output to new file
        with open(f"mm_{system}_{run}_{frame}.txt", "r") as f:
            input_file_contents = f.read()
            output_file_contents = os.popen(f"echo '{input_file_contents}' | awk '{addone_script}'").read()
            with open(f"MM_{system}_{run}_{frame}.dat", "w") as f_out:
                f_out.write("set active { " + output_file_contents.strip() + " }\n")


        # execute AWK script on input file and redirect output to new file
        with open(f"mm_{system}_{run}_{frame}.txt", "r") as f:
            input_file_contents = f.read()
            output_file_contents = os.popen(f"echo '{input_file_contents}' | awk '{addone_script}'").read()
            with open(f"MM_{system}_{run}_{frame}.dat", "w") as f_out:
                f_out.write("set active { " + output_file_contents.strip() + " }\n")

        # execute another AWK script on input file and redirect output to new file
        with open(f"qm_{system}_{run}_{frame}.txt", "r") as f:
            input_file_contents = f.read()
            output_file_contents = os.popen(f"echo '{input_file_contents}' | awk '{addone_script}'").read()
            with open(f"QM_{system}_{run}_{frame}.dat", "w") as f_out:
                f_out.write("set qm_atoms { " + output_file_contents.strip() + " }\n")


        # Define the path for the system, run, and frame
        qmmm_path = os.path.join(qmmm_root, f"{system}_{run}_{frame}")
        # Define the path for the RC optimization directory
        rc_path = os.path.join(qmmm_path, "1-rc-opt")

        # create directories
        def check_create_warn(folder_path):
            if not os.path.exists(folder_path):
                os.mkdir(folder_path)
            elif not os.path.isdir(folder_path):
                print(f"{folder_path} already exists but is not a directory")
        
        for folder in [qmmm_root, qmmm_path, rc_path]:
            check_create_warn(folder)

        print("Copying files to the 1-rc-opt directory")
        # Copy files to the 1-rc-opt directory
        os.system(f"cp rc.pdb {rc_path}/.")
        os.system(f"cp rc.rst {rc_path}/.")
        os.system(f"cp rc.prmtop {rc_path}/.")
        os.system(f"cp QM_{system}_{run}_{frame}.dat {rc_path}/QM.dat")
        os.system(f"cp MM_{system}_{run}_{frame}.dat {rc_path}/MM.dat")
        os.system(f"cp myresidues_{system}_{run}_{frame}.dat {rc_path}/myresidues.dat")
        os.system(f"cp parse_amber.tcl {rc_path}/.")
        os.system(f"cp input.in {rc_path}/.")


        # Change the working directory
        os.chdir(f'{qmmm_root}/{system}_{run}_{frame}/1-rc-opt')

        print("Creating RC Optimization Input File")

        with open('RC_dlfind.chm', 'w') as f:
            f.write("# adenine - Amber example with polarisation turned off\n")
            f.write("# hybrid with electrostaic embedding\n")
            f.write("global sys_name_id\n")
            f.write("source parse_amber.tcl\n")
            f.write("source MM.dat\n")
            f.write("source QM.dat\n")
            f.write("source myresidues.dat\n")
            f.write("set sys_name_id rc\n")
            f.write("set prmtop rc.prmtop\n")
            f.write("set inpcrd rc.rst\n")
            f.write("load_amber_coords inpcrd=$inpcrd prmtop=$prmtop coords=rc.c\n")
            f.write("# # for the time being we have to calculate an energy to be able to call list_amber_atom_charges\n")
            f.write("energy energy=e coords=rc.c theory=dl_poly  : [ list \\\n")
            f.write("amber_prmtop_file=$prmtop \\\n")
            f.write("scale14 = [ list [ expr 1 / 1.2 ] 0.5  ] \\\n")
            f.write("mxexcl=2000  \\\n")
            f.write("mxlist=40000 \\\n")
            f.write("cutoff=1000 \\\n")
            f.write("use_pairlist = no \\\n")
            f.write("save_dl_poly_files = yes \\\n")
            f.write("exact_srf=yes \\\n")
            f.write("list_option=none ]\n")
            f.write("set atom_charges [ list_amber_atom_charges ]\n")
            f.write("# optimize geometry with distance A-B fixed\n")
            f.write("dl-find coords=rc.c maxcycle=999 active_atoms=$active residues=$myresidues list_option=full result=${sys_name_id}.opt.c \\\n")
            f.write("theory=hybrid : [ list \\\n")
            f.write("coupling= shift \\\n")
            f.write("qm_region= $qm_atoms \\\n")
            f.write("atom_charges= $atom_charges \\\n")
            f.write("qm_theory= turbomole : [list   \\\n")
            f.write("read_control= yes \\\n")
            f.write(str("scratchdir=ocean/projects/che160019p/") + str(user) + str("/temp \\\n"))
            f.write("hamiltonian= b3-lyp \\\n")
            f.write("scftype= uhf  ]  \\\n")
            f.write("mm_theory= dl_poly  : [ list \\\n")
            f.write("amber_prmtop_file= $prmtop \\\n")
            f.write("exact_srf=yes \\\n")
            f.write("use_pairlist=no \\\n")
            f.write("mxlist=40000 \\\n")
            f.write("cutoff=1000 \\\n")
            f.write("mxexcl=2000  \\\n")
            f.write("debug_memory=no \\\n")
            f.write("scale14 = [ list [ expr 1 / 1.2 ] 0.5  ] \\\n")
            f.write("conn= rc.c \\\n")
            f.write("save_dl_poly_files = yes \\\n")
            f.write("list_option=none ]]\n\n")
            f.write('# save structure\n')
            f.write("read_pdb  file= ${sys_name_id}.pdb  coords=hybrid.dl_poly.coords \\\n")
            f.write("write_pdb file= ${sys_name_id}.opt.pdb coords= ${sys_name_id}.opt.c \\\n")
            f.write("write_xyz file= ${sys_name_id}.QMregion.opt.xyz coords=hybrid.${qm_theory}.coords \\\n\n")
            f.write(" exit\n")



    if step == "1":
        print("Running Step 1")

        job = os.getcwd()
        jobname = "RC-Optimization"
        subprocess.Popen(["nohup", "chemsh", "rc_dlfind.chm"], stdout=open("rc_dlfind.log", "w"), stderr=subprocess.STDOUT)
        process = subprocess.Popen(["pidof", "chemsh.x"], stdout=subprocess.PIPE)
        pid_out, pid_err = process.communicate()
        pid = pid_out.decode().strip().replace(" ", ",")
        while os.path.exists("rc_dlfind.log"):
            with open("rc_dlfind.log") as f:
                if "Terminated" in f.read():
                    print("RC Terminated by User")
                    exit()
            process = subprocess.Popen(["pidof", "chemsh.x"], stdout=subprocess.PIPE)
            pid_out, pid_err = process.communicate()
            new_pid = pid_out.decode().strip().replace(" ", ",")
            if pid != new_pid:
                break
        print("RC Terminated normally")
        print("RC Terminated. Now Running Define")
        
        
        define_input = '''a coord
        *
        no
        b all def2-SVP
        *
        eht
        y
        {}
        n
        u {}
        *
        n
        scf
        iter
        900

        dft
        on
        func b3-lyp

        *'''.format(charge, unp)
        
        define_output = subprocess.check_output("define", input=define_input, text=True)
        print(define_output)
        print("Executing RC Optimization")

        omit = subprocess.check_output(["pidof", "chemsh.x"])
        string = ",".join(omit.decode().split())

        subprocess.run(["tcsh", "-c", f"setenv PARNODES {nodes};nohup chemsh RC_dlfind.chm >& RC_dlfind.log &"])

        time.sleep(5)
        if not string:
            calc = subprocess.check_output(["pidof", "chemsh.x"])
        else:
            calc = subprocess.check_output(["pidof", "-o", string, "chemsh.x"])
            
        time.sleep(5)
        while psutil.pid_exists(int(calc)):
            time.sleep(1)

        if "Energy evaluation failed" in open("RC_dlfind.log").read():
            print("DSCF Failed. Now changing SCF iterlimit and Restarting")
            with open("control", "r+") as f:
                content = f.read()
                f.seek(0)
                f.write(re.sub(f"{scfiterlimit}\s+100", f"{scfiterlimit}      900", content))
                f.truncate()

            omit = subprocess.check_output(["pidof", "chemsh.x"])
            string = ",".join(omit.decode().split())
            subprocess.run(["tcsh", "-c", f"setenv PARNODES {nodes};nohup chemsh RC_dlfind.chm >& RC_dlfind.log &"])
            subprocess.run(f'echo "{job} {system} {frame} JOB SCF Error and Restarted" | mail -s "Job Restarted" simahjsr@gmail.com', shell=True)
            
            time.sleep(5)
            if not string:
                calc = subprocess.check_output(["pidof", "chemsh.x"])
            else:
                calc = subprocess.check_output(["pidof", "-o", string, "chemsh.x"])
                
            time.sleep(5)
            while psutil.pid_exists(int(calc)):
                time.sleep(1)
            subprocess.run(f'echo "Job Completed in {host} on `date` for {system} {jobname} at {job} " | mail -s "Job Completed {system}" {email}', shell=True)
        else:
            print("RC Completed")
            subprocess.run(f'echo "Job Completed in {host} on `date` for {system} {jobname} at {job} " | mail -s "Job Completed {system}" {email}', shell=True)


main()