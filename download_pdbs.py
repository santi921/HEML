import subprocess, os, sys
import unicodedata
import urllib.request

<<<<<<< HEAD
import pandas as pd 

=======
>>>>>>> c9b36e8852adf4014b6fe50e318edf3c55848e8c
def remove_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")

#-------------------
<<<<<<< HEAD
#with open("heme_list.txt") as f: 
#    lines = f.readlines()
#heme_names = [remove_control_characters(line.split(" ")[0]) for line in lines]
heme_names = pd.read_csv("protein_data.csv")['name']

for i in heme_names:
    print(i)
    try:
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', './pdbs/'+i+'1.pdb')
    except:
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', './pdbs/'+i+'.pdb')
=======
with open("heme_list.txt") as f: 
    lines = f.readlines()
heme_names = [remove_control_characters(line.split(" ")[0]) for line in lines]

for i in heme_names:
    print(i)
    #try:
    #    urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', './pdbs/'+i+'1.pdb')
    #except:
    urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', './pdbs/'+i+'.pdb')
>>>>>>> c9b36e8852adf4014b6fe50e318edf3c55848e8c
