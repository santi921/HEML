import re
import numpy as np
from HEML.utils.data import get_fe_positions, get_ligand_info


def test_process_charges():
    charge_file = "./test_data/test.pqr"

    fe_dict = get_fe_positions(charge_file)
    assert fe_dict["id"] != None
    ligand_dict = get_ligand_info(charge_file, fe_dict["xyz"])
    ligand_identifier = ligand_dict["best_crit"].split(":")

    zero_active = True
    zero_everything_charged = True
    ligands_to_zero_radius = 5.0
    carbene_tf = True
    zero_radius = True

    ligand_zero_count = 0
    active_zero_count = 0
    carbene_zero_count = 0
    zero_charged_count = 0

    openfile = open(charge_file)
    readfile = openfile.readlines()

    for j in readfile:
        line_split = re.split(r"(\s+)", j)
        # strip whitespace from list
        line_split = [x.strip() for x in line_split]
        # remove empty strings from list
        line_split = list(filter(None, line_split))

        cond = (
            line_split[8] == ligand_identifier[0]
            and line_split[10] == ligand_identifier[1]
        )

        if zero_active and ("HETATM" in line_split[0] or cond):
            active_zero_count += 1

        elif carbene_tf and line_split[3] in ["CB1"]:
            carbene_zero_count += 1

        elif zero_everything_charged and line_split[3] in [
            "ASP",
            "GLU",
            "LYS",
            "ARG",
            "HIS",
        ]:
            zero_charged_count += 1

        elif zero_radius:
            xyz_str_x = j[30:38]
            xyz_str_y = j[38:46]
            xyz_str_z = j[46:54]

            distance = np.sqrt(
                (float(xyz_str_x) - float(fe_dict["xyz"][0])) ** 2
                + (float(xyz_str_y) - float(fe_dict["xyz"][1])) ** 2
                + (float(xyz_str_z) - float(fe_dict["xyz"][2])) ** 2
            )

            if distance < ligands_to_zero_radius:
                ligand_zero_count += 1
    assert ligand_zero_count == 17, "ligand_zero_count is wrong"
    assert active_zero_count == 292, "active_zero_count is wrong"
    assert carbene_zero_count == 0, "carbene_zero_count is wrong"
    assert zero_charged_count == 9578, "zero_charged_count is wrong"
