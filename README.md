# HEML

Installing

``pip install -e .``

In the current directory. This will give you the following scripts in your enviroment:

- ``remove_solvent_from_pdb.py``
- ``split_movies.py``
- ``run_topology.py``
- ``run_box_calc.py``
- ``get_charges.py``
- ``process_charges.py``
- ``dist_and_compress_each_protein.py``
- ``sub_turbo.py``

For Computing Fields/Boxs on an MD trajectory:

    1) Run``split_movies.py`` to separate your frames. If you need to remove solvent molecules, run ``remove_solvent_from_pdb.py`` subsequently.

    2) Setup an``./options.json`` file with the following entries:

* cpet_loc - currently CPET is a dependency and can be found at https://github.com/matthew-hennefarth/CPET. We are in the process of creating a pythonic variant that can be installed via pip
* chargefw2_loc 
* processed_charges_folder
* cpet_folder
* charges_folder
* processed_pdb_folder

3) You can now run ``get_charges.py`` and this will generate files in the pqr format in the destination (``charges_folder``).
4) Process the charges and generate the cpet input files via ``process_charges.py``
5) Run your box or topology calculation using ``run_topology.py`` or ``run_box_calc.py``

If you use any of the scripts/utilities in this package it would be goated if you cited the following 
1) Todo: add HEML
2) Todo: add Proto
