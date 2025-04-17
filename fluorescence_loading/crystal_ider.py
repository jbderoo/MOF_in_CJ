# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 22:16:06 2021

@author: jderoo
"""

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

mydir = 'TopLeft_HCl++_Right_HCl+_Bottom_Control001'
mystr = 'TopLeft_HCl++_Right_HCl+_Bottom_Control001_T01_Z5.tif'
num_images = 41

def read_tiff(path):
    """
    path - Path to the multipage-tiff file
    """
    img = Image.open(path)
    images = []
    for i in range(img.n_frames):
        img.seek(i)
        images.append(np.array(img))
    return np.array(images) 

plt.rcParams.update({'figure.max_open_warning'   : num_images+1})



# top left
square1y = [100, 200]
square1x = [60, 175]

# bottom left
square2y = [350, 475]
square2x = [100, 250]

# bottom right
square3y = [275, 375]
square3x = [350, 450]

# background
square4y = [0, 200]
square4x = [300, 500]

cryst1 = []
cryst2 = []
cryst3 = []
bgrnd  = []

i=1

file = open(mydir + '/' + mystr,'r')
im = Image.open(mydir + '/' + mystr,'r')
imarray = np.array(im)
plt.figure()
plt.imshow(imarray, cmap='gray')
plt.colorbar()
plt.title(f'image {i}')

imarray[square1y[0]:square1y[1], square1x[0]:square1x[1]] = np.max(imarray)
imarray[square2y[0]:square2y[1], square2x[0]:square2x[1]] = np.max(imarray)
imarray[square3y[0]:square3y[1], square3x[0]:square3x[1]] = np.max(imarray)
imarray[square4y[0]:square4y[1], square4x[0]:square4x[1]] = np.max(imarray)



plt.figure()
plt.imshow(imarray, cmap='gray')
plt.colorbar()
plt.title(f'image {i}')
plt.text(67, 155, 'crystal 1')
plt.text(120, 420, 'crystal 2')
plt.text(350, 328, 'crystal 3')
plt.text(350, 100, 'background')


for i in range(1, num_images+1):
    mystr = f'TopLeft_HCl++_Right_HCl+_Bottom_Control001_T{i:02}_Z5.tif'
    

    path = mydir + '/' + mystr
    images = read_tiff(path)
    imarray = images[1]

    cryst1.append(np.mean(imarray[square1y[0]:square1y[1], square1x[0]:square1x[1]]))
    cryst2.append(np.mean(imarray[square2y[0]:square2y[1], square2x[0]:square2x[1]]))
    cryst3.append(np.mean(imarray[square3y[0]:square3y[1], square3x[0]:square3x[1]]))
    bgrnd.append(np.mean(imarray[square4y[0]:square4y[1], square4x[0]:square4x[1]]))

init = [cryst1[0],cryst2[0],cryst3[0], bgrnd[0]]
for i in range(len(cryst1)):
    cryst1[i] = cryst1[i]-init[0]
    cryst2[i] = cryst2[i]-init[1]
    cryst3[i] = cryst3[i]-init[2]
    bgrnd[i]  = bgrnd[i] -init[3]

t = np.arange(0,len(cryst1)/2, 0.5)
plt.figure(dpi=600)
plt.plot(t, cryst1, label = 'crystal 1') 
plt.plot(t, cryst2, label = 'crystal 2') 
plt.plot(t, cryst3, label = 'crystal 3')
plt.plot(t, bgrnd, label = 'background') 
plt.legend()   
plt.xlabel('min')
plt.ylabel('d fluoresence')
plt.title('change in fluoresence over time')
plt.xticks(ticks=np.arange(0, 31, 5))
plt.tight_layout()
    
