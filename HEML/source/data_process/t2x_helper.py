import os 
root = './'
directories=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
for dir in directories:
    os.chdir(dir)
    os.chdir("./embedding/o")
    xyz_file = [f for f in os.listdir(os.getcwd()) if f.endswith('.xyz')][0]
    os.system("x2t {} > coord".format(xyz_file))
    # get the file in directory ending with .xyz
    os.chdir("../oh")
    xyz_file = [f for f in os.listdir(os.getcwd()) if f.endswith('.xyz')][0]
    os.system("x2t {} > coord".format(xyz_file))
    os.chdir("../normal")
    xyz_file = [f for f in os.listdir(os.getcwd()) if f.endswith('.xyz')][0]
    os.system("x2t {} > coord".format(xyz_file))
    os.chdir("../../../")