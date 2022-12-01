import json
from HEML.utils.cpet import run_topology_calcs
from HEML.utils.data import get_options

def main():
    options = get_options("./options_local.json")
    cpet_folder = options["cpet_folder"]
    processed_charges_folder = options["processed_charges_folder"]
    run_topology_calcs(cpet_folder, processed_charges_folder)
    
main()