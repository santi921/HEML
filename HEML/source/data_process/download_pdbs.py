import subprocess, os, sys
import unicodedata
import urllib.request

import pandas as pd 

def remove_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")

def main():
    heme_names = pd.read_csv("../../../data/protein_data.csv")['name']

    for i in heme_names:
        print(i)
        try:    
            urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb1.gz', '../../../data/pdbs/'+i+'.pdb1')
            # uncompress gz file to pdb file
            subprocess.call(['gunzip', '../../../data/pdbs/'+i+'.pdb1'])
            print("bioassembly")
        except:
            print("NOT bioassembly" * 3)
            urllib.request.urlretrieve('https://files.rcsb.org/download/' + i + '.pdb', '../../../data/pdbs/'+i+'.pdb')


main()