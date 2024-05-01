# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 23:19:24 2024

@author: jderoo
"""

import numpy as np
from scipy.spatial import distance
import sys
from pdb_work_v5 import ezPDB

            
mof_file        = 'trimmed_cubtc.pdb'
pore_file       = '221_CJ_only_chainA_centered.pdb'
cutoff_distance = 2.2
output_file     = 'output_CuBTC_pore_filling.pdb'

mof  = ezPDB(mof_file)
pore = ezPDB(pore_file)

mof_xyz  = mof.xyz
pore_xyz = pore.xyz

print('starting pore_mof analysis')
dists = distance.cdist(mof_xyz, pore_xyz)
pore_mof_tf = dists < 2.5


with open(mof_file) as file, open('tmp.pdb', 'w') as writefile:
    for i, line in enumerate(file):
        if not any(pore_mof_tf[i]):
            writefile.write(line)

print('starting mof_mof analysis')
trimmed_mof = ezPDB('tmp.pdb')
tmof_xyz    = trimmed_mof.xyz
dists       = distance.cdist(tmof_xyz, tmof_xyz)
mof_mof_tf  = (dists < 2.2) & (dists > 0.2)


cu_count = 0
with open('tmp.pdb') as file, open(output_file, 'w') as writefile:
    for i, line in enumerate(file):
        if np.sum(mof_mof_tf[i]) >= 2:
            writefile.write(line)
            if "Cu" in line:
                cu_count += 1


print(f"Filtered PDB written to real_{output_file}")
print(f'There are {cu_count} Cu atoms in the file')

