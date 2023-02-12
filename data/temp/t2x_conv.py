import os 
root = './'
for dir in os.listdir(root):
    os.system("./t2x {}/embedding/o/coord > {}/embedding/o/pos.xyz".format(dir, dir))
    os.system("./t2x {}/embedding/oh/coord > {}/embedding/oh/pos.xyz".format(dir, dir))
    os.system("./t2x {}/embedding/normal/coord > {}/embedding/normal/pos.xyz".format(dir, dir))
