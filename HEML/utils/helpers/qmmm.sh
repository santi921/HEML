#!/bin/bash
while getopts i:s:a:b:t:p: flag
do
    case "${flag}" in
	i) inp=${OPTARG};;
	s) step=${OPTARG};;
	a) A=${OPTARG};;
	b) B=${OPTARG};;
	t) transition=${OPTARG};;
	p) product=${OPTARG};;
esac
done

source ${inp}
user=$USER
host=$(hostname)

if [ "$step" = "0" ]; then

mkdir scratch
cd scratch

cp $parsefile .
file_name=$tleapinput

# Read the lines containing atom types
atom_type_lines=$(grep -oP '\{ "\K\w+(?=.*"sp3")' "$file_name")

# Initialize an associative array to store atom type values
declare -A atom_type_values

# Iterate over the atom type lines and extract the values
while read -r atom_type; do
    value=$(grep -oP "{ \"$atom_type\"  \"\K.*(?=\" \"sp3\")" "$file_name")
    atom_type_values["$atom_type"]="$value"
done <<< "$atom_type_lines"

​# Print the extracted values
for atom_type in "${!atom_type_values[@]}"; do
    echo "$atom_type: ${atom_type_values[$atom_type]}"
done
​
# Read the contents of the parse_amber.tcl file
parse_amber_file="./parse_amber.tcl"
content=$(cat "$parse_amber_file")
​
# Make the required substitutions
for atom_type in "${!atom_type_values[@]}"; do
    atom_value=${atom_type_values[$atom_type]}
    if [[ $atom_value == Fe ]]; then
        #sed -i "s/FN - F3 - F4 { return 26 }/FN - F3 - F4 - $atom_type { return 26 }/g" "$parse_amber_file"
        sed -i "s/\(.*\){ return 26 }/\1- $atom_type { return 26 }/" "$parse_amber_file"
    elif [[ $atom_value == N ]]; then
        sed -i "s/\(.*\){ return 7 }/\1- $atom_type { return 7 }/" "$parse_amber_file"
    elif [[ $atom_value == O ]]; then
        sed -i "s/\(.*\){ return 8 }/\1- $atom_type { return 8 }/" "$parse_amber_file"
    elif [[ $atom_value == S ]]; then
        sed -i "s/\(.*\){ return 16 }/\1- $atom_type { return 16 }/" "$parse_amber_file"
    elif [[ $atom_value == C ]]; then
        sed -i "s/\(.*\){ return 6 }/\1- $atom_type { return 6 }/" "$parse_amber_file"
    else
        echo "Atom type $atom_type not substituted in parse_amber.tcl file." || exit 1
    fi
done


#CREATING INPUT FILES FOR CPPTRAJ FOR PREPARING RC COMPLEX FILES

filename="$pdb"


# Save the original IFS value
OLDIFS=$IFS
# Set IFS to "_" and "-" to split the filename
IFS="_-."
# Use the read command to extract the desired part
read -ra parts <<< "$filename"

system="${parts[2]}"
run="${parts[3]}"
frame="${parts[4]}"


# Print the variables
echo "filename = $filename"
echo "system = $system"
echo "run = $run"
echo "frame = $frame"

cat > ${system}_${run}_${frame}.in << ENDOFFILE
parm ${parm}
trajin ${trajin} ${frame} ${frame}
trajout ${system}_${run}_${frame}.inpcrd restart
trajout ${system}_${run}_${frame}.pdb
run
exit
ENDOFFILE

#cat inpcrd_${frame}.in

cat > strip_${system}_${run}_${frame}.in <<ENDOFFILE
parm ${parm}
trajin ${system}_${run}_${frame}.inpcrd
reference ${system}_${run}_${frame}.inpcrd
strip :Na+,Cl-
strip !(:$numberofres<:10.0) outprefix stripped10
trajout stripped10.${system}_${run}_${frame}.inpcrd restart
run
exit
ENDOFFILE

#cat water_strip_${frame}.in

cat > rc_${system}_${run}_${frame}.in <<ENDOFFILE
parm stripped10.*.prmtop
trajin stripped10.${system}_${run}_${frame}.inpcrd
trajout rc.pdb
trajout rc.rst restart
run
exit
ENDOFFILE
cat rc_${system}_${run}_${frame}.in

#CREATING RC COMPLEX FILES
nohup cpptraj -i ${system}_${run}_${frame}.in > ${system}_${run}_${frame}.out &
process=$!
while ps -p ${process} > /dev/null;do sleep 1;done;
if [ "$(grep -c "Error" ${system}_${run}_${frame}.out)" -ge 1 ]; then
                echo "Cpptraj Error"
                exit
        else
                echo "Generated Frame PDB"
        fi
nohup cpptraj -i strip_${system}_${run}_${frame}.in > strip_${system}_${run}_${frame}.out &
process=$!
while ps -p ${process} > /dev/null;do sleep 1;done;
if [ "$(grep -c "Error" strip_${system}_${run}_${frame}.out)" -ge 1 ]; then
                echo "Cpptraj Error"
                exit
        else
                echo "Generated WaterStripped PDB and prmtop"
        fi
nohup cpptraj -i rc_${system}_${run}_${frame}.in > rc_${system}_${run}_${frame}.out &
process=$!
while ps -p ${process} > /dev/null;do sleep 1;done;
if [ "$(grep -c "Error" rc_${system}_${run}_${frame}.out)" -ge 1 ]; then
                echo "Cpptraj Error"
                exit
        else
                echo "Generated RC files for RC_OPT"
        fi
echo "Copying prmtop file: stripped10.${base}.prmtop to rc.prmtop"
cp stripped10.*.prmtop rc.prmtop
sed -i "9s/1/0/" rc.prmtop

#MAKING QMMM MODEL
echo "Creating QMMM Model"
cat > QM_MM_${system}_${run}_${frame}.tcl <<ENDOFFILE
mol load pdb rc.pdb
atomselect top "same residue as (within 8 of (resname $resname $substrate))"
atomselect0 num
atomselect0 writepdb MM_${system}_${run}_${frame}.pdb
set myfile [open mm_${system}_${run}_${frame}.txt w]
puts \$myfile [atomselect0 list]
close \$myfile
atomselect top "(resname $resname and not backbone and not type HA H) or (resname $substrate)"
atomselect1 num
atomselect1 writepdb QM_${system}_${run}_${frame}.pdb
atomselect1 writexyz QM_${system}_${run}_${frame}.xyz
set myfile1 [open qm_${system}_${run}_${frame}.txt w]
puts \$myfile1 [atomselect1 list]
close \$myfile1
atomselect top "resname $resname and type HA or resname FE1 or resname HM1 and name NA or resname CB1 and name C1 or resname $substrate and name C3"
set resid [atomselect1 get resid]
foreach elementid \$resid {dict set tmp \$elementid 1}
set id [dict keys \$tmp]
set resname [atomselect1 get resname]
foreach elementname \$resname {dict set tmp2 \$elementname 1}
set name [dict keys \$tmp2]
set myresidues [open qm_mm_${system}_${run}_${frame}.sh w]
puts \$myresidues "resid=(\$id)"
puts \$myresidues "resname=(\$name)"
puts \$myresidues "myresidues=()"
puts \$myresidues "n=\\\${#resname\\[@]}"
puts \$myresidues "for i in \\\$(seq 1 \\\$n);"
puts \$myresidues "do"
puts \$myresidues "myresidues+=(\\\${resname\\[i-1]}\\\${resid\\[i-1]})"
puts \$myresidues "done"
puts \$myresidues "cat > myresidues_${system}_${run}_${frame}.dat <<ENDOFFILE"
puts \$myresidues "set res \\[ pdb_to_res \\"rc.pdb\\"]"
puts \$myresidues "set myresidues  \\[ inlist function=combine residues= \\\\\\\$res sets= {\\\${myresidues\\[*]}} target=QM ]"
puts \$myresidues "ENDOFFILE"
close \$myresidues
exit
ENDOFFILE

echo "Using rc.pdb"
echo "Using Residues:${resname}"
echo "Using Substrate:${substrate}"

vmd -dispdev text -e QM_MM_${system}_${run}_${frame}.tcl

chmod +x qm_mm_${system}_${run}_${frame}.sh
./qm_mm_${system}_${run}_${frame}.sh

cat > addone.awk <<ENDOFFILE

BEGIN{
   RS = " "
}

{
a = \$1
++a
printf( "%d " , a )
}

ENDOFFILE
awk -f addone.awk mm_${system}_${run}_${frame}.txt > MM_${system}_${run}_${frame}.dat
sed -i '1s/^/set active {/' MM_${system}_${run}_${frame}.dat
echo "}" >> MM_${system}_${run}_${frame}.dat
awk -f addone-awk qm_${system}_${run}_${frame}.txt > QM_${system}_${run}_${frame}.dat
sed -i '1s/^/set qm_atoms {/' QM_${system}_${run}_${frame}.dat
echo "}" >> QM_${system}_${run}_${frame}.dat



#MAKING DIRECTORIES
if [[ ! -e ../../QMMM ]]; then
    mkdir ../../QMMM
elif [[ ! -e ../../QMMM ]]; then
    echo "QMMM already exists but is not a directory" 1>&2
fi

if [[ ! -e ../../QMMM/${system}_${run}_${frame} ]]; then
    mkdir ../../QMMM/${system}_${run}_${frame}
elif [[ ! -d ../../QMMM/${system}_${run}_${frame} ]]; then
    echo "QMMM already exists but is not a directory" 1>&2
fi

if [[ ! -e ../../QMMM/${system}_${run}_${frame}/1-rc-opt ]]; then
    mkdir ../../QMMM/${system}_${run}_${frame}/1-rc-opt
elif [[ ! -d ../../QMMM/${system}_${run}_${frame}/1-rc-opt ]]; then
    echo "1-rc-opt already exists but is not a directory" 1>&2
fi

#COPYING FILES TO THE RC_Opt DIRECTORIES
cp rc.pdb ../../QMMM/${system}_${run}_${frame}/1-rc-opt/.
cp rc.rst ../../QMMM/${system}_${run}_${frame}/1-rc-opt/.
cp rc.prmtop ../../QMMM/${system}_${run}_${frame}/1-rc-opt/.
cp QM_${system}_${run}_${frame}.dat ../../QMMM/${system}_${run}_${frame}/1-rc-opt/QM.dat
cp MM_${system}_${run}_${frame}.dat ../../QMMM/${system}_${run}_${frame}/1-rc-opt/MM.dat
cp myresidues_${system}_${run}_${frame}.dat ../../QMMM/${system}_${run}_${frame}/1-rc-opt/myresidues.dat
cp parse_amber.tcl ../../QMMM/${system}_${run}_${frame}/1-rc-opt/.
cp input.in ../../QMMM/${system}_${run}_${frame}/1-rc-opt/.
cd ../../QMMM/${system}_${run}_${frame}/1-rc-opt


cat > RC_dlfind.chm <<ENDOFFILE
# adenine - Amber example with polarisation turned off
# hybrid with electrostaic embedding
global sys_name_id
source parse_amber.tcl
source MM.dat
source QM.dat
source myresidues.dat
set sys_name_id rc
set prmtop rc.prmtop
set inpcrd rc.rst
load_amber_coords inpcrd=\$inpcrd prmtop=\$prmtop coords=rc.c
# # for the time being we have to calculate an energy to be able to call list_amber_atom_charges
energy energy=e coords=rc.c theory=dl_poly  : [ list \\
amber_prmtop_file=\$prmtop \\
scale14 = [ list [ expr 1 / 1.2 ] 0.5  ] \\
mxexcl=2000  \\
mxlist=40000 \\
cutoff=1000 \\
use_pairlist = no \\
save_dl_poly_files = yes \\
exact_srf=yes \\
list_option=none ]

set atom_charges [ list_amber_atom_charges ]

# optimize geometry with distance A-B fixed
dl-find coords=rc.c maxcycle=999 active_atoms= \$active residues= \$myresidues list_option=full result=\${sys_name_id}.opt.c \\
theory=hybrid : [ list \\
coupling= shift \\
qm_region= \$qm_atoms \\
atom_charges= \$atom_charges \\
qm_theory= turbomole : [list   \\
read_control= yes \\
scratchdir=ocean/projects/che160019p/$user/temp \\
hamiltonian= b3-lyp \\
scftype= uhf  ]  \\
mm_theory= dl_poly  : [ list \\
amber_prmtop_file= \$prmtop \\
exact_srf=yes \\
use_pairlist=no \\
mxlist=40000 \\
cutoff=1000 \\
mxexcl=2000  \\
debug_memory=no \\
scale14 = [ list [ expr 1 / 1.2 ] 0.5  ] \\
conn= rc.c \\
save_dl_poly_files = yes \\
list_option=none ]]

####
# save structure
read_pdb  file= \${sys_name_id}.pdb  coords=hybrid.dl_poly.coords
write_pdb file= \${sys_name_id}.opt.pdb coords= \${sys_name_id}.opt.c
write_xyz file= \${sys_name_id}.QMregion.opt.xyz coords=hybrid.\${qm_theory}.coords

  exit

ENDOFFILE

elif [ "$step" = "1" ]; then
job=$(pwd)
jobname="RC-Optimization"
nohup chemsh rc_dlfind.chm > rc_dlfind.log &
process=$!
while ps -p ${process} > /dev/null;do sleep 1;done;
if [ "$(grep -c "Terminated" rc_dlfind.log)" -ge 1 ]; then
		echo "RC Terminated by User"
		exit
	else
		echo "RC Terminated normally"
	fi
echo "RC Terminated.Now Running Define"

define <<EOF


a coord
*
no
b all def2-SVP
*
eht
y
$charge
n
u $unp
*
n
scf
iter
900

dft
on
func b3-lyp

*
EOF

echo "Executing RC Optimization"

omit=$(pidof chemsh.x)
string="${omit//${IFS:0:1}/,}"

tcsh -c "setenv PARNODES $nodes;nohup chemsh RC_dlfind.chm >& RC_dlfind.log &"

sleep 5
if [ -z "$string" ]
then
calc=$(pidof chemsh.x)
else
calc=$(pidof -o ${string} chemsh.x)
fi
sleep 5
while ps -p "${calc}" > /dev/null;do sleep 1;done;

if [ "$(grep -c "Energy evaluation failed" RC_dlfind.log)" -ge 1 ]; then
        echo "DSCF Failed. Now changing SCF iterlimit and Restarting"
        sed -i "s/$scfiterlimit      100/$scfiterlimit      900/" control
        omit=$(pidof chemsh.x)
	string="${omit//${IFS:0:1}/,}"
        tcsh -c "setenv PARNODES $nodes;nohup chemsh RC_dlfind.chm >& RC_dlfind.log &"
	echo "$job $system $frame JOB SCF Error and Restarted" | mail -s "Job Restarted" simahjsr@gmail.com
	sleep 5
	if [ -z "$string" ]
	then
	calc=$(pidof chemsh.x)
	else
	calc=$(pidof -o ${string} chemsh.x)
	fi
	sleep 5
	while ps -p "${calc}" > /dev/null;do sleep 1;done;
	echo "Job Completed in ${host} on `date` for ${system} ${jobname} at ${job} " | mail -s "Job Completed ${system}" simahjsr@gmail.com
else
        echo "RC Completed"
        echo "Job Completed in ${host} on `date` for ${system} ${jobname} at ${job} " | mail -s "Job Completed ${system}" simahjsr@gmail.com
fi
fi