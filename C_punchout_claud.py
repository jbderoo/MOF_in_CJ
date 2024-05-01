import numpy as np
from pdb_work_v5 import ezPDB

file_path = 'cubtc_775_no_dups_centered.pdb'
p = ezPDB(file_path)
cu_xyz = p.xyz

# Center the coordinates at (0, 0, 0)
center = np.mean(cu_xyz, axis=0)
cu_xyz_centered = cu_xyz - center


# Create boolean masks for each dimension
mask_x = (cu_xyz_centered[:, 0] > -200) & (cu_xyz_centered[:, 0] <= 200)
mask_y = (cu_xyz_centered[:, 1] > -200) & (cu_xyz_centered[:, 1] <= 200)
mask_z = (cu_xyz_centered[:, 2] > -25) & (cu_xyz_centered[:, 2] <= 25)
radius_mask = (cu_xyz_centered[:, 0]**2 + cu_xyz_centered[:, 1]**2)**0.5 <= 70

# Combine the masks to get the final mask
mask = mask_x & mask_y & mask_z & radius_mask

# Open the input and output files
with open(file_path, 'r') as f, open('trimmed_cubtc.pdb', 'w') as f1:
    # Write the start line
    start_line = f.readline()
    f1.write(start_line)

    # Iterate over the remaining lines
    for i, line in enumerate(f):
        # Check if the line corresponds to a coordinate within the desired bounds
        if mask[i]:
            x, y, z = cu_xyz_centered[i]
            f1.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            #f1.write(line)

    # Write the end line (if present)
    try:
        end_line = f.readline()
        f1.write(end_line)
    except:
        pass