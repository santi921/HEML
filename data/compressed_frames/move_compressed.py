#get all directory names 
import os 

directories = os.listdir()
target = "../frames_processed_2/"
# scan all *pdb files in another directory and they contain a folder name copy to that directory 
for directory in directories:
    if os.path.isdir(directory):
        files = os.listdir(target)
        for file in files:
            if file.endswith('.pdb'):
                os.system('cp %s/%s %s' % (target, file, directory))
