

import os 
def find_files_that_start_with_0_and_delete(folder):
    files = os.listdir(folder)
    for i in files:
        if(i[0] == "0"):
            os.remove(folder + i)
            print("removed {}".format(i))


def find_files_that_dont_start_with_pro_and_delete(folder):
    files = os.listdir(folder)
    for i in files:
        if(i[0:3] != "pro"):
            os.remove(folder + i)
            print("removed {}".format(i))
