# -*- coding: utf-8 -*-
"""
Created on Thu May  6 20:17:04 2021

@author: jderoo

This "package" uses python 3.6 or later - there are f-strings in here!
"""

import numpy as np

def xrd_extract(image, machine, verbose):
    '''


    Parameters
    ----------
    image : an image file that has been opened by fabio, --> img = fabio.open(image_path)

    machine : string
        either rigaku or synchtrotron/synch, tells the function how you collected your data
    verbose : True/False
        Returns print statements detailing the extraction

    Returns
    -------
    x_center: the x coordinate of the center of the image. Returned in pixels
        
    y_center: the y coordinate of the center of the image. Returned in pixels
        
    x_pixel: the conversion of pixels to mm in x direction
        
    y_pixel: the conversion of pixels to mm in y direction
        
    wl: the wavelength in nm used by the instrument
        
    distance: how far the detector was from the sample at the time of shooting


    '''
     
    if machine == 'synch' or machine == 'synchtrotron':
        params   = list(image.header['RDI1_SPATIAL_DISTORTION_INFO'].split()[:])
        x_center = int(float(params[0]))
        y_center = int(float(params[1]))
        x_pixel  = float(params[2])
        y_pixel  = float(params[3])
        wl       = float(image.header['WAVELENGTH'])/10
        distance = float(image.header['DISTANCE'])
        filename = image.header['FILENAME'].split('/')[-1]
        
    elif machine == 'rigaku':  
        params   = list(image.header['PILT_SPATIAL_DISTORTION_INFO'].split()[:])
        x_center = int(float(params[0]))
        y_center = int(float(params[1]))
        x_pixel  = float(params[2])
        y_pixel  = float(params[3])
        wl       = float(image.header['SOURCE_WAVELENGTH'].split()[1])/10
        distance = float(image.header['PILT_GONIO_VALUES'].split()[-1])
        filename = image.header['FILENAME']
        
    else:
        return print(f'you chose {machine}, this function only accepts "rigaku" or "synchtrotron" as options')
    if verbose == True:
        print('------------')
        print(f'image {filename} has:')
        print(f'center coords {x_center}, {y_center}')
        print(f'pixel size {x_pixel} x {y_pixel} mm')
        print(f'has wavelength {wl} nm and distance {distance} mm')
        print('------------')
    return x_center, y_center, x_pixel, y_pixel, wl, distance 

def find_max(pt, data, y):
    """
    Finds the local maxima between two bounds values

    Parameters
    ----------
    pt : tuple
        [upperbound, lowerbound]
    data : list
        the independent list of your data. upperbound and lowerbound should exist in this list
    y : list
        the dependent list of your data. the maxima that is found will be on this ilst

    Returns
    -------
    xmax: the x coordinate of the maxima
    
    ymax : the y coordinate of the maxima
       

    """
    lb = pt[0]
    ub = pt[1]
    lbI = min(range(len(data)), key=lambda i: abs(data[i]-lb))
    ubI = min(range(len(data)), key=lambda i: abs(data[i]-ub))
    ymax = max(y[lbI:ubI])
    maxI = min(range(len(y)), key=lambda i: abs(y[i]-ymax))
    xmax = data[maxI]
    return xmax, ymax


def create_circular_mask(h, w, dr, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = (radius-dr <= dist_from_center ) & (dist_from_center <= radius)

    return mask


def theta_integration(data, pixel_range, dp, n, x_og, y_og, std=False):
    '''
    performs theta integration on the data
    data        = img.data extract from fabio
    pixel_range = what range of pixels would you like to integrate from. Looking at the 
                  image, this is starting in the center and working outward. E.g. for a data 
                  size of 500x500, an appropriate initial pixel range is [0, 250]
    dp          = the thickness of the ring that moves outward during the integration. The higher,
                  the smoother the peaks, the broader the bases. Play with until pretty graphs.              
    n           = number of steps the ring outward takes. Usually 100-500 is good.
    x_og, y_og  = center of the image. Can be obtained from "xrd_extract" easily.
    
    '''
    r = np.ceil(np.linspace(pixel_range[0], pixel_range[1], n)).astype(int)
    dr = dp
    radial_averages = np.zeros(n)
    radial_std      = np.zeros(n)
    
    for i in range(len(r)):
        h, w = data.shape[:2]
       # mask = create_circular_mask(h, w, dr, center = [origin_x-400, origin_y-407], radius = r[i])
        mask = create_circular_mask(h, w, dr, center = [x_og, y_og], radius = r[i])
        masked_img = data.copy()
        radial_averages[i] = np.mean(masked_img[mask])
        radial_std[i]      = np.std(masked_img[mask])
        

        if radial_averages[i] > 3000:
            print(masked_img[mask])
        #circle2 = plt.Circle((origin_x-400,  origin_y-407),r[i], color='r', fill=False)
        #circle2 = plt.Circle((x_og,  y_og),r[i], color='r', fill=False)
    if not std:
        return r, radial_averages
    else:
        return r, radial_averages, radial_std
    
def pixel_to_theta(img, machine, r, it, maxpxrd, base, std=False):
    '''
    

    Parameters
    ----------
    img : image extracted from fabio, e.g. img = fabio.open('/path/to/image')

    machine : either "rigaku" or "synchtrotron"/"synch". string

    r : list
        the pixel independent vector that is made in theta_integration
    it : list
        the pixel intensity vector (dependent vector) that is made in theta_integration
    maxpxrd : scalar 
        the maximum value of the known fingerprint/PXRD being compared to.
    base : scalar
        the average value that the plot is shifted to for overlay purposes

    Returns
    -------
    theta_r : list
        Bragg's 2theta 
    theta_it : list
        somewhat arbitrary rescaling of pixel intensity --> count

    '''
    
    x_center, y_center, x_pixel, y_pixel, wl, distance = xrd_extract(img, machine, False)
    #distance = 27
    theta_r     = 0.154/wl*np.rad2deg(np.arctan(np.array(r)*x_pixel/distance)) # 
     # 
    theta_it    = it/max(it)*maxpxrd
    diff        = (np.mean(theta_it) - base)
    theta_it    = theta_it - diff 
    
    if type(std) == bool:
        return theta_r, theta_it
    
    else:
        theta_std   = 0.154/wl*np.rad2deg(np.arctan(np.array(std)*x_pixel/distance))
        return theta_r, theta_it, theta_std
        
def dspace_convert(img, r, machine):
    '''
    wl given in nm
    theta2 is 2theta angle
    return d-space in A
    '''
    
    #dspace = wl*10/2/np.sin(theta2/(2) * np.pi/180)
    x_center, y_center, x_pixel, y_pixel, wl, distance = xrd_extract(img, machine, False)
    #distance = 27
    theta = np.rad2deg(np.arctan(r*x_pixel/distance)/2)
    dspace = (wl*10) / (2*np.sin(theta * np.pi/180))
    
    return dspace
    
    