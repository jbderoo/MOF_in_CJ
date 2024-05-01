# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 17:42:05 2023

@author: jderoo
"""

import numpy as np
import matplotlib.pyplot as plt
import fabio
from xrd_tools import xrd_extract, find_max, theta_integration, pixel_to_theta, dspace_convert
import sys
import glob
import os

#path = 'day1/snap/'
path = 'results/'
all_images = sorted(glob.glob(path + '*.img'))
    
core_string = '1D-01_snap_0009.img'

filtered_images = []
for img_path in all_images:
    base_name = os.path.basename(img_path)
    if core_string not in base_name:
        filtered_images.append(img_path)

all_images = filtered_images

true_centers = [
    [1447, 1414],
    [1447, 1414],
    ]

plotter_num = 1
colors = ['tab:blue', 'tab:orange']
labels = [' ']
for i in range(0, len(all_images), 2):
    
    img1 = fabio.open(all_images[i])    
    img2 = fabio.open(all_images[i+1])
    
    #data1 = img1.data
    #data2 = img2.data
    
    all_data = [img1, img2]
    
    plt.figure(plotter_num, dpi=600)
    #plt.title(f'picture set {plotter_num}')

    for j, d in enumerate(all_data):
        #if plotter_num == 1:
        #    continue
        
        data1 = d.data
        x_center, y_center, x_pixel, y_pixel, wl, distance = xrd_extract(d, 'synch', True)
        x_og = x_center
        y_og = y_center
        n = 575
        wdata = data1[y_og-n:y_og+n, x_og-n:x_og+n]
        
        pixel_range = [0, max([x_center, y_center])] # yields full integration
        pixel_range = [100, 700] # big picture with meaningful data
        pixel_range = [200, 550] # most meaningful data
        dp = 12 # the thickness of the pixel's reading: e.g. at pixel 100, it reads 
                # pixel intensities at 95-105 and averages them to 100
                # affects smoothness/clarity of plot
         
        n  = 100 # number of azimuthal slices taken between the pixel range, affects 
                 # affects resolution of plot (and time to run)
        #x_center = 1432
        #y_center = 1416    
        if path == 'day2/snap/':
            x_center = 1448
            y_center = 1408
        else:
            x_center = 1432
            y_center = 1416
            
        if path == 'day1/snap/' and plotter_num == 1:
            x_center = 1448
            y_center = 1409
            
            
        pixels, pixel_average, pixel_std = theta_integration(data1, pixel_range, dp, n, x_center, y_center, True)
        r = pixels
        it = pixel_average
        #fig, ax = plt.subplots(dpi=600)
        #ax.plot(pixels, pixel_average, 'g')
        
        #plt.figure(dpi=600)
    
    
        # bound/extract the meaningful data together
        #wdata = data1[y_og-n:y_og+n, x_og-n:x_og+n]
        #plt.imshow(wdata, cmap="gray")
        
        maxpxrd = 100
        base    = 1
        converted_Braggs_2theta, rescaled_intensity_to_count, int_std = pixel_to_theta(d, 'synch', r, it, maxpxrd, base, std=pixel_std)
    
        rescaled_intensity_to_count += 10
    
        #plt.plot(converted_Braggs_2theta, rescaled_intensity_to_count, color=colors[i], label=d.filename.split('\\')[1].split('.')[0])
        plt.plot(converted_Braggs_2theta, rescaled_intensity_to_count, color=colors[j], label='          ')
        plt.plot(converted_Braggs_2theta, rescaled_intensity_to_count-int_std, color=colors[j], alpha=0.27)
        plt.plot(converted_Braggs_2theta, rescaled_intensity_to_count+int_std, color=colors[j], alpha=0.27)
        lower_std = rescaled_intensity_to_count - int_std
        upper_std = rescaled_intensity_to_count + int_std
        
        # Shade the area between the upper and lower standard deviation lines
        plt.fill_between(converted_Braggs_2theta, lower_std, upper_std, color=colors[j], alpha=0.2)

    plotter_num += 1
    plt.legend(loc=0)
    plt.savefig('../../Spring_2024/cds_notes/rotation_analysis.svg',format='svg', bbox_inches='tight' )
    if plotter_num == 4:
        sys.exit()
    
    
    
    #sys.exit()
    
