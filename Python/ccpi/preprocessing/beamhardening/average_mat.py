#
# Code to average over number of tif images of Crown samples
# and apply flat and dark corrections.
# Use I/I0  = (image-dark) / (flat-dark)
# image, flat and dark are averages over the set of tiffs in the
# corresponding directory.
# Input required is a set of directories containing one or more tiff
# images each. This list of directories is provided in a file, one
# directory per line. The first must be the dark images, the second the
# flat images and subsequent directories are of the calibaration materials.
# e.g.
#    python average_mat.py matdirs.txt output.raw
#
# Where "matdirs.txt" might contain:
#   //Rcu029/MXIF/data/temp/Sara/ForRon/SecondBatch/dark
#   //Rcu029/MXIF/data/temp/Sara/ForRon/SecondBatch/flat
#   //Rcu029/MXIF/data/temp/Sara/ForRon/SecondBatch/Material1
#   //Rcu029/MXIF/data/temp/Sara/ForRon/SecondBatch/Material2
#   //Rcu029/MXIF/data/temp/Sara/ForRon/SecondBatch/Material3
#
# Within each of the above directories there must be one or more files
# with the .tif extension.
# Each image must be of the same resolution (e.g. 2000x2000) and with the
# same exposure time as each other.
# Present practice is to take 11 tiff images per material sample, but this
# program should read any number of tiffs in each directory.
# The final output is the averaged and normalised image of each material.
# The I/I0 values are floating point values in the range 0 to 1.0.
# To save space, without significant loss of accuracy, values are
# multiplied by 65535 and then converted from float32 to uint16.
#
# This format must be labelled in the <name>.data file used for fitting
# as: uint16_65535
#
# An option exists to save only a subsection of each image by adding
# four integers to the end of the argument list:
#
#  python average_mat.py output.raw xmin ymin xmax ymax
#
# e.g.
#  python average_mat.py output.raw 500 500 1500 1500
#  
# which will only write the corrected image for the central part
# of the image between (500,500) and (1500,1500) in the case of
# 2000x2000 pixel images.
#
# This makes the image data much smaller, though at the expense of preventing
# a correction fit that varies with the line number (as may be needed if
# the X-ray take-off angle is important). However, with the X-tek set up
# using the crown to image samples, the source may be moved closer to the
# object for imaging the samples. Varying this distance will prevent
# correction on take-off angle for the sample. This is because the take-off
# angle for a line in the material image is not necessarily the same as that
# in the case of the actual sample.
#
# If a subset of the full image is used, e.g. (500,500) to (1500,1500),
# it is important to set the reduced image size in the .data file used
# to define the images (1000 by 1000 in this case).
#
# Note that raw format is used for the output, so it is important that
# the same architecture is used in both this preprocessing step as in the
# fitting process. Most machines use Intel processors so this is unlikely
# to be a problem.
#
from __future__ import print_function
import sys
import numpy as np
import os
from os import walk
from PIL import Image

def getaverage(dir):
    """ Given a directory of tiff files, assumed to be images of
    uint16 and same size, sum all and return average as float32 """
    print(dir)
    for (dirpath,dirname,filenames) in walk(dir):
        count = 0
        for f in filenames:
            if f.endswith(".tif"):
                count = count+1
                with open(dir+'/'+f,'rb') as fptr:
                    im = Image.open(fptr)
                    imarray = np.array(im)
                    print("mean/max/min: ",np.mean(imarray),np.max(imarray),np.min(imarray))
                    if count==1:
                        farray = imarray.astype('float32')
                        new = False
                    else:
                        farray = farray+imarray
                    del imarray
        farray = farray/count
        print("ave mean/max/min: ",np.mean(farray),np.max(farray),np.min(farray))
        return farray

print("Average and dark and flat field correct the Crown images")
print("Output in raw uint16 with normalisation to 65535")
if len(sys.argv)>2:
    indirfile = sys.argv[1]
    outfile = sys.argv[2]
    xmin = -1
    if len(sys.argv)==7 :
        xmin=int(sys.argv[3])
        ymin=int(sys.argv[4])
        xmax=int(sys.argv[5])
        ymax=int(sys.argv[6])
        print("output limited to: ",xmin,ymin,xmax,ymax)
else:
    print("syntax: python ",sys.argv[0]," <fileOfDirs> <out.raw> [xmin ymin xmax ymax]")
    sys.exit(1)

dirs=[]
with open(indirfile) as f:
    dirlist = f.readlines()

dirlist = [x.strip() for x in dirlist]
for dr in dirlist:
    if not os.path.isdir(dr):
        print("Error: Directory ",dr," not found")
        sys.exit(1)

means = []
for dir in dirlist:
    means.append(getaverage(dir))

dark = means[0]
flat = means[1]
divider = flat-dark
if len(divider[divider==0.])>0 :
    print("Warning: zeros in flat-dark: ",len(divider[divider==0.]) )
    divider[divider==0.] = 1000

fout = open(outfile,'wb')
for i in range(len(means)-2):
    corrected = (means[i+2]-dark)/divider
    corrected[corrected>1.0] = 1.0
    corrected = corrected*65535
    if xmin==-1:
        ui16 = corrected.astype('uint16')
    else:
        ui16 = corrected[xmin:xmax,ymin:ymax].astype('uint16')
    print("Image ",i," mean = ",np.mean(ui16)," max = ",np.max(ui16)," min = ",np.min(ui16))
    fout.write(ui16.tobytes())
fout.close()
