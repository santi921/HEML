import os 
root = './'
directories=[d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
for dir in directories:
    os.chdir(dir)
    os.chdir("./embedding/o")
    os.system("t2x > opt.xyz")
    os.system("t2x -c > final.xyz")
    os.chdir("../oh")
    os.system("t2x > opt.xyz")
    os.system("t2x -c > final.xyz")
    os.chdir("../normal")
    os.system("t2x > opt.xyz")
    os.system("t2x -c > final.xyz")
    os.chdir("../../../")