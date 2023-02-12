from HEML.utils.cpet import run_box_calcs
from HEML.utils.data import get_options
import os 

def main():
    options = get_options("./options/options_local.json")
    cpet_folder = options["cpet_folder"]
    processed_charges_folder = options["processed_charges_folder"]
    print(os.listdir(cpet_folder))
    print(os.listdir(processed_charges_folder))
    run_box_calcs(cpet_folder, processed_charges_folder)
main()