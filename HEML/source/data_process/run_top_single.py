import json
from HEML.utils.cpet import run_topology_calcs

def main():
    with open("./options.json") as f:
        options = json.load(f)
    cpet_folder = options["cpet_folder"]
    processed_charges_folder = options["processed_charges_folder"]
    run_topology_calcs(cpet_folder, processed_charges_folder)
    
main()