#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 14:47:20 2022

@author: jderoo
"""

import numpy as np
import matplotlib.pyplot as plt
import fabio
from xrd_tools import xrd_extract, find_max, theta_integration, pixel_to_theta, dspace_convert
import sys
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
from cycler import cycler


font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
colors = ['#18737a', '#a82e36','#68bd8f','#604380', '#19274f','#c6bf55','#576244' ]
plt.rcParams.update({'font.size': 14, 'font.weight':'bold','font.family':'normal'  }   )
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})

plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})
plt.rcParams.update({'axes.prop_cycle':cycler(color=colors)})


plt.rcParams.update({'xtick.major.width'   : 2.8 })
plt.rcParams.update({'xtick.labelsize'   : 10 })
plt.rcParams.update({'xtick.major.size'   : 6 })


plt.rcParams.update({'ytick.major.width'   : 2.8 })
plt.rcParams.update({'ytick.labelsize'   : 10})
plt.rcParams.update({'ytick.major.size'   : 6 })


plt.rcParams.update({'axes.linewidth':2.8})
plt.rcParams.update({'axes.labelpad':8})
plt.rcParams.update({'axes.titlepad':10})
plt.rcParams.update({'figure.dpi':600})


dp = 5
n = 300 # DO NOT XHANGE
store = np.zeros(shape=[12,n])


for i in range(11, 12):
    insert = str(i).zfill(2)
# ---- open img and convert to data matrix -----
    img_path = f'1E-01_snap_00{i}.img' # simplified to best image in the group
    
 

    img = fabio.open(img_path)
    data = img.data
    
    
    # ---- extract meaningful info from image metadata -----
    x_center, y_center, x_pixel, y_pixel, wl, distance = xrd_extract(img, 'synch', True)
    #distance = 40
    pixel_range = [215, 725] # big picture with meaningful data
    pixels, pixel_average = theta_integration(data, pixel_range, dp, n, x_center, y_center)

    plt.figure(dpi=600)
    plt.plot(pixels, pixel_average, 'g')
    plt.xlabel('pixels')
    plt.ylabel('pixel intensity')
    plt.title('raw picture/data output')
    
    plt.figure(2, dpi=600)
    plt.imshow(data*1000, cmap="gray_r")
    plt.tight_layout()
    plt.savefig('scXRD_pattern_UiO-67.png', dpi=600)
    #sys.exit()
    
    maxpxrd = 20
    base = 2
    r = pixels
    it = pixel_average
    converted_Braggs_2theta, rescaled_intensity_to_count = pixel_to_theta(img, 'synch', r, np.array(it), maxpxrd, base)
    #plt.figure(1)
    #plt.plot(converted_Braggs_2theta, rescaled_intensity_to_count)
    #store[i-1,:] = rescaled_intensity_to_count
   # plt.axis([])

    #break


    #converted_Braggs_2theta, rescaled_intensity_to_count = pixel_to_theta(img, 'rigaku', r, it, maxpxrd, base)
#y = store.sum(axis=0)
y = rescaled_intensity_to_count + 20



plt.figure(3, dpi=600)
plt.plot(converted_Braggs_2theta, y)
plt.xlabel('2theta')
plt.ylabel('count')
plt.title('extract pXRD of UiO-67@CJ from scXRD')
#plt.axis([4, 15, 18, 30])

dspace = dspace_convert(img, r, 'synch')


plt.figure(4, dpi=600)
plt.plot(dspace, y)
plt.xlabel('dspace')
plt.ylabel('count')
plt.title('extract dspace intensities of UiO-67@CJ from scXRD')
plt.gca().invert_xaxis()

compare = np.zeros([len(dspace), 2])
compare[:,0] = converted_Braggs_2theta
compare[:,1] = y


peak1 = [26, 46]
peak1 = np.arange(peak1[0], peak1[1]+1, 1)

peak2 = [57, 70]
peak2 = np.arange(peak2[0], peak2[1]+1, 1)

peak3 = [135, 150]
peak3 = np.arange(peak3[0], peak3[1]+1, 1)

peak4 = [198, 212]
peak4 = np.arange(peak4[0], peak4[1]+1, 1)

peak5 = [183, 196]
peak5 = np.arange(peak5[0], peak5[1]+1, 1)

plt.figure(3, dpi=600)
plt.plot(converted_Braggs_2theta, y, 'b', label='base error')
plt.plot(converted_Braggs_2theta[peak1], y[peak1], 'r--', label='ignored in error')
plt.plot(converted_Braggs_2theta[peak2], y[peak2], 'r--', label='ignored in error')
plt.plot(converted_Braggs_2theta[peak3], y[peak3], 'r--', label='ignored in error')
plt.plot(converted_Braggs_2theta[peak4], y[peak4], 'r--', label='ignored in error')
plt.plot(converted_Braggs_2theta[peak5], y[peak5], 'r--', label='ignored in error')
plt.xlabel('2theta')
plt.ylabel('count')
plt.title('pXRD UiO-67@CJ reigions ignored for error')



all_peaks = np.concatenate( (peak1, peak2, peak3, peak4, peak5) )
I = np.arange(0, n, 1)

a = list(set(I) - set(all_peaks))

curve_x = converted_Braggs_2theta[a]
curve_y = y[a]

plt.figure(5, dpi=600)
plt.plot(curve_x, curve_y)
plt.xlabel('2theta')
plt.ylabel('count')
plt.title('corresponding error')

def minme(s, x, y):
    A = s[0]
    B = s[1]
    diff = A/x + B - y
    return diff

x0 = [1, 15]
lst_sqs = least_squares(minme, x0, args = (curve_x, curve_y))
A = lst_sqs.x[0]
B = lst_sqs.x[1]

x_fit = np.linspace(4.5, 14, 300)
y_fit = A/x_fit + B

f = interp1d(curve_x, curve_y)
y_fit2 = f(x_fit)
'''
plt.figure(dpi=600)
plt.plot(curve_x, curve_y, 'b')
#plt.plot(x_fit, y_fit, 'r--')
plt.plot(x_fit, y_fit2, 'r--')
'''


proper_y = []

for pt in range(len(converted_Braggs_2theta)):
    my_y = y[pt]
    
    if my_y < 18.5: break
    proper_y.append( my_y - float(f( converted_Braggs_2theta[pt] )) + 1 )


UiOx = np.load('UiO-67_theta.npy')
UiOy = np.load('UiO-67_count.npy')

y = np.array(proper_y)
y = y/max(y)*100

plt.figure(dpi=600)
plt.plot(UiOx, UiOy+2, 'k', label='simulated')
plt.plot(converted_Braggs_2theta-.1, y-20,'c', label='observed', linewidth=2)
#plt.xlabel('2theta')
#plt.ylabel('intensity')
#plt.title('extract pXRD of UiO-67@CJ from scXRD without noise')
plt.axis([4, 12, 0, 105])
plt.legend(loc=0)
plt.tight_layout()
plt.savefig('publication_UiO-67@CJ.png',dpi=600)


