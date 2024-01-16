# HEML

Installing

``pip install -e .``

In the current directory. This will give you the following scripts in your enviroment:
- ``get_charges.py`` - uses ChargeFW2 to compute the charges on a pdb with default paramters
- ``process_charges.py`` - constructs input files for CPET calculations. Added optionality for processing. 
- ``run_topology.py`` - runs cpet topology calculations on a folder
- ``run_box_calc.py`` - runs cpet box field calculations on a folder
- ``dist_and_compress_each_protein.py`` - specify code for analyzing MD trajectories
- ``run_sweep.py`` - runs a sweep of parameters for compression along with statistics to  pick the best clustering hyperparameters.
- ``dist_on_topo.py`` - runs a compression on the current folder


For Computing Fields/Boxs on an MD trajectory:

    1) Run``split_movies.py`` to separate your frames. If you need to remove solvent molecules, run ``remove_solvent_from_pdb.py`` subsequently.

    2) Setup an``./options.json`` file with the following entries:

* cpet_loc - currently CPET is a dependency and can be found at https://github.com/matthew-hennefarth/CPET. We are in the process of creating a pythonic variant that can be installed via pip
* chargefw2_loc - another dependency for computing electric fields.
* processed_pdb_folder - finalized location for pdb files.
* charges_folder - where raw charges are located.
* processed_charges_folder - where finalized charges are to be output.
* cpet_folder- location for output cpet calculations and input files.

3) You can now run ``get_charges.py`` and this will generate files in the pqr format in the destination (``charges_folder``).
4) Process the charges and generate the cpet input files via ``process_charges.py``
5) Run your box or topology calculation using ``run_topology.py`` or ``run_box_calc.py``
6) We have implemented helper methods to build your own compression workflows but you can use the ``dist_on_topo.py`` script to compress the current folder

If you use any of the scripts/utilities in this package it would be goated if you cited the following 
1) Todo: add HEML
2) Todo: add Proto
