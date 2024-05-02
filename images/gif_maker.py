# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:11:19 2024

@author: jderoo
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import imageio.v2 as imageio
import matplotlib.gridspec as gridspec


datafile = 'CuBTC_CJ_raw_data.npz'
if not os.path.isfile(datafile):
    import fabio
    from xrd_tools import xrd_extract 
    
    img_path = '2D-02_1_0001.img'
    img = fabio.open(img_path)
    data = img.data
    
    x_center, y_center, x_pixel, y_pixel, wl, distance = xrd_extract(img, 'synch', verbose=True)
    
    np.savez(datafile, data=data, x_center=x_center, y_center=y_center, 
         x_pixel=x_pixel, y_pixel=y_pixel, wl=wl, distance=distance)

    print(f"Data and extraction results saved to {datafile}")
    
else:
    loaded_data = np.load(datafile)
    x_center = loaded_data['x_center']
    y_center = loaded_data['y_center']
    x_pixel  = loaded_data['x_pixel']
    y_pixel  = loaded_data['y_pixel']
    wl       = loaded_data['wl']
    distance = loaded_data['distance']
    data     = loaded_data['data']
    
    print(f'loaded data from {datafile}')
    
def create_circular_mask(h, w, dp, center=None, radius=None):
    if center is None:
        center = [int(w/2), int(h/2)]
    if radius is None:
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y - center[1])**2)

    mask = (dist_from_center <= radius + dp) & (dist_from_center >= radius - dp)
    return mask

def calculate_frame_repeats(change_rate, threshold=0.5, max_repeats=10, min_repeats=1):
    """Calculate the number of repeats based on the rate of change."""
    if abs(change_rate) > threshold:
        return min_repeats  # Rapid change, fewer repeats
    else:
        return max_repeats  # Slow change, more repeats

# Load your data
with open('perfect_CuBTC.npy', 'rb') as f:
    x = np.load(f)
    y = np.load(f)

# just orientation stuff to get it all on the same field of view
y = y / max(y) * 60 
y += 65    


lower_lim = 100
upper_lim = 575
n_pts = 100  # n points used to integrate
dp = 10   # thickness of the circle
bnds = 575

wdata = data[y_center-bnds:y_center+bnds, x_center-bnds:x_center+bnds]
h, w = wdata.shape[:2]
radii = np.ceil(np.linspace(lower_lim, upper_lim, n_pts)).astype(int)

all_theta_r = []
all_int_r   = []

# Prepare to save frames
filename = 'azimuthal_integration_V5.gif'
writer = imageio.get_writer(filename, fps=5)  # Adjust fps to speed up or slow down the animation
last_value = None

# manual points added to slow gif down around the peaks

ra = list(radii)

extras = [251, 253, 256, 258, 261, 263, 266, 270, 271, 275, 276, 376, 377, 381, 382, 386, 390, 391, 395, 396, 461, 463, 466, 467, 468, 471, 473, 476, 477, 478, 482, 484]
for e in extras:
    ra.append(e)
ra = sorted(ra)



# Generate frames
for i, r in enumerate(ra):
    if i % 5 == 0:
        print(f'Completion: {i/len(ra)*100:1.2f}%')
    
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), dpi=600)
    
    # Set up the figure and GridSpec
    fig = plt.figure(figsize=(13, 5), dpi=600)  # Adjust total figsize to fit both plots
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 8])  # width_ratios controls the relative widths
    
    # Create subplots
    ax1 = fig.add_subplot(gs[0])  # First subplot
    ax2 = fig.add_subplot(gs[1])  # Second subplot, wider
    
    # Image with circle on the left
    ax1.imshow(wdata * 500, cmap="gray_r")
    circle = Circle((bnds, bnds), r, color='r', fill=False, linewidth=1.5)
    ax1.add_patch(circle)
    ax1.set_title('CuBTC@CJ scXRD')
    ax1.axis('off')  # Hide axes for a cleaner look
    

    # Plot on the right
    mask = create_circular_mask(h, w, dp, center=[bnds, bnds], radius=r)
    masked_data = wdata[mask]
    radial_averages = np.mean(masked_data)
    radial_std = np.std(masked_data)
    theta_r = 0.154 / wl * np.rad2deg(np.arctan(np.array(r) * x_pixel / distance))
    
    all_theta_r.append(theta_r)
    all_int_r.append(radial_averages)
    
    
    ax2.plot(x, y, color='black', label='expected pXRD', linewidth=0.75)
    ax2.plot(all_theta_r, all_int_r, color = 'r', label='computed pXRD', linewidth=2)
    ax2.axis([2, 15, 60, 140])
    ax2.legend(loc='upper center')
    ax2.set_xlabel('2\u03B8 (degrees)')
    ax2.set_ylabel('Count')
    ax2.set_title('CuBTC extracted pXRD')
    ax2.figure.tight_layout()

    # Save frame
    plt.savefig('temp_plot.png', dpi=600)
    plt.close(fig)  # Close the plot to free memory
    image = imageio.imread('temp_plot.png')
    writer.append_data(image)

    last_value = radial_averages  # Update last value for next iteration

# Finish the GIF
writer.close()