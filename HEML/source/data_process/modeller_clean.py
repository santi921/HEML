# /home/santiagovargas/bin/modeller10.3/bin/modpy.sh modeller_clean.py

# TODO 

from modeller import *
from modeller.scripts import complete_pdb

from modeller import *
from modeller.automodel import *    # Load the AutoModel class

log.verbose()
env = Environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = AutoModel(env, alnfile = 'alignment.ali',
              knowns = '1qg8', sequence = '1qg8_fill')

a.starting_model= 1
a.ending_model  = 1
a.loop.starting_model = 1
a.loop.ending_model   = 2
a.loop.md_level       = refine.fast

a.make()
# TODO