from glob import glob
from PIL import Image

files = glob("*svg")
for file in files:
    # convert file to png
    im = Image.open("./"+file)
    im.save(file.split(".")[0] + ".png")
    
