# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 23:01:42 2024

@author: jderoo
"""

import numpy as np

def read_pdb(file_path):
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append((x, y, z))
    return atoms

def center_atoms(atoms):
    # Calculate the centroid of the points
    centroid = np.mean(atoms, axis=0)
    # Center the atoms at the origin
    centered_atoms = atoms - centroid
    return centered_atoms

def write_pdb(original_file_path, new_file_path, centered_atoms):
    atom_idx = 0
    with open(original_file_path, 'r') as original, open(new_file_path, 'w') as new_file:
        for line in original:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x, y, z = centered_atoms[atom_idx]
                atom_idx += 1
                # Format the new coordinates into the line
                new_line = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
                new_file.write(new_line)
            else:
                new_file.write(line)


# center the large MOF structure

file_path = 'cubtc_775_no_dups.pdb'
new_file_path = 'cubtc_775_no_dups_centered.pdb'

atoms = read_pdb(file_path)
centered_atoms = center_atoms(np.array(atoms))
write_pdb(file_path, new_file_path, centered_atoms)
