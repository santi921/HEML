import subprocess, os, sys
import unicodedata
import urllib.request

import pandas as pd 

def remove_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")

#-------------------
#with open("heme_list.txt") as f: 
#    lines = f.readlines()
#heme_names = [remove_control_characters(line.split(" ")[0]) for line in lines]
heme_names = pd.read_csv("../../../data/protein_data.csv")['name']

for i in heme_names:
    print(i)
    try:
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', '../../data/pdbs/'+i+'.pdb1')
    except:
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', '../../data/pdbs/'+i+'.pdb')
