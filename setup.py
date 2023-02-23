from setuptools import setup, find_packages
# to setup utils run the following
# pip install -e . 

setup(
    name = 'heml',
    version = '0.0.1',
    packages = find_packages(),
    scripts=['./HEML/utils/helpers/remove_solvent_from_pdb.py',
             './HEML/utils/helpers/split_movies.py',
             './HEML/source/data_process/run_topology.py',
             './HEML/source/data_process/run_box_calc.py',
             './HEML/source/data_process/get_charges.py',
             './HEML/source/data_process/process_charges.py',
             './HEML/source/data_process/dist_and_compress_each_protein.py',
             './HEML/source/data_process/sub_turbo.py'
             ],
)