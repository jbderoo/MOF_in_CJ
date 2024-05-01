# -*- coding: utf-8 -*-
"""
Created on Wed May  1 00:04:47 2024

@author: jderoo
"""

# 221_CJ_centered was create by `cat 221_fixed_CJ_centered_rechained_part1.pdb 221_fixed_CJ_centered_rechained_part2.pdb >  221_CJ_centered.pdb`
# they are split into part 1 and part 2 because it takes more than 26 monomers
# to create the pore. Thus chain overwrite happens, and pymol cannot correctly render. However, we still need
# all the coordinates in a single pdb object to do the analysis; thats what this script is for.

centered = '221_CJ_centered.pdb'

with open(centered) as file, open('221_CJ_only_chainA_centered.pdb', 'w') as writefile:
    for line in file:
        if line[:4] == 'ATOM':
            newline = line[:21] + 'A' + line[22:]
        else:
            newline = line
        
        writefile.write(newline)
