# find all the files in the current directory ending in *svg and convert them to png
from glob import glob
from PIL import Image


# downsample the image to 1/4 of the original size
def downsample_image(image, factor):
    width, height = image.size
    return image.resize((width // factor, height // factor), Image.ANTIALIAS)

for file in glob("*png"):
    im = Image.open(file)
    im = downsample_image(im, 4)
    im.save(file.split(".")[0] + "down.png", "PNG")
