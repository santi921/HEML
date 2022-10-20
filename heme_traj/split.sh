csplit -f ${1%.pdb}_ -b %04d.pdb --suppress-matched -z ${1} ./ENDMDL/ {*}
