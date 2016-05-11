#
# Program to read in the polynomial fit(s) from carousel analysis and apply them to
# a set of data. The Data can be float values of I0/I or int16 values of I0/I with I0 defined
# in some way. e.g. single I0 value or matrix of I0 values
# Corrected values will be written in same format as input.
#
from __future__ import print_function
import os.path
import sys, getopt
import numpy as np
import numpy.polynomial.polynomial as nppoly
import pdb


def readPolyCoeff(fileParam):
    # read the "param.log" type data of polynomial coeffs
    # retunrs numpy array polyC
    if not os.path.isfile(fileParam):
        return 0
    else:
        fl = open(fileParam,"rb")
        polyfits = np.load(fl)
        #print("polyfits="+str(polyfits))
        return polyfits

def checkArgs(argstr):
    # set default lines and rows for QMUL data
    lines=800
    rows=600
    whiteLevel = 0
    polyf = "polyfit.npz"
    bhtFile = ""
    syntax = "applyTrans.py [-r rows -l lines -p poly.npz -w whiteLevel -x file.bht] [-d] [file1.ext] [filen.ext]"
    try:
        opts, args = getopt.getopt(argstr,"hr:l:dp:w:x:",["rows=","lines=","debug","polyfit","whiteLevel","xtekbht"])
    except getopt.GetoptError:
        print(syntax)
        sys.exit(2)
    debug = False
    for opt, arg in opts:
        if opt == "-h":
            print(syntax)
            sys.exit()
        elif opt in ("-r", "--rows"):
            rows = int(arg)
        elif opt in ("-l", "--lines"):
            lines = int(arg)
        elif opt in ("-p", "--polyfit"):
            polyf = arg
        elif opt in ("-w", "--whiteLevel"):
            whiteLevel = int(arg)
        elif opt in ("-x", "--xtekbht"):
            bhtFile = arg
        elif opt in ("-d", "--debug"):
            debug = True
    return rows,lines,polyf,whiteLevel,bhtFile,debug,args

def processRaw(infile,outfile,polyC,rows,lines,debug):
    #
    # Read a binary raw file containing float32 data of the ln(I0/I) image data
    # and apply the correction given by the polynomial(s) in polyC.
    # Results output written to outfile in same raw format.
    #
    infi = open(infile,"rb")
    outfi = open(outfile,"wb")
    count = 0
    while True:
        img=np.fromfile(infi,np.float32,rows*lines)
        count=count+1
        if len(img)<1:
            break
        print("read block len="+str(len(img)))
        try:
            img = img.reshape(lines,rows)
        except:
            print("Error reshaping block "+str(count))
            sys.exit(1)
        if len(polyC[:,0])==lines:
            print("Using "+str(lines)+" lines for correction")
        else:
            print("Using 1 line for correction")
        for l in range(lines):
            ln = 0
            if len(polyC[:,0])==1:
                ln = 0
            elif len(polyC[:,0])==lines:
                ln = l
            img[l,:] = correct(img[l,:],polyC[ln,::-1])
        img.tofile(outfi)
        
    infi.close()
    outfi.close()

def processTif(infile,outpre,polyC,rows,lines,whiteLevel,debug):
    #
    # Read a tif file containing uint16 data of "I" image data
    # and apply the correction given by the polynomial(s) in polyC.
    # White level "whiteLevel" must be provided.
    # Have to convert to float format for correction and retur to uint16
    # afterwards. A correction table in the reconstruction process would
    # be more efficient.
    # Results output written to outfile in same tif format.
    #
    try:
        #from PIL import Image
        import tifffile as tf
    except:
        #print("Import PIL failed")
        print("Import tifffile failed")
        sys.err(1)
    #im = Image.open(infile)
    imArr = tf.imread(infile)
    if len(imArr[0,:]) != lines or len(imArr[:,0]) != rows:
        print("* Error: tiff image size is ",len(imArr[0,:])," by ",len(imArr[:,0]))
        print("Must set lines and rows the same")
        sys.exit(1)
    #imArr = np.array(im)
    if whiteLevel==0:
        whiteLevel = imArr.max()
        if debug:
           print("Using whiteLevel=",whiteLevel)
    #
    imLn = np.log(float(whiteLevel)/imArr)

    for l in range(lines):
        ln = 0
        if len(polyC[:,0])==1:
            ln = 0
        elif len(polyC[:,0])==lines:
            ln = l
        imLn[l,:] = correct(imLn[l,:],polyC[ln,::-1])
    imCor = whiteLevel/np.exp(imLn)
    imCor16 = np.array(imCor,dtype='uint16')
    if debug:
        print("imCor16: ",imCor16[0:2,0:2])
    #imOut = Image.fromarray(imCor16)
    #print("infile=",infile) #," outpre=",outpre)
    outfile = outpre+infile
    #imOut.save(outfile)
    tf.imsave(outfile,imCor16)

def correct(lineData,polyData):
    return nppoly.polyval(lineData,polyData)

def genbht(bhtFile,whiteLevel,polyC):
    # Generate the bht file from the first polynomial data. Only one correction
    # curve can be used in a bht file, so line wise correction is not possible.
    # Any values above the whiteLevel are set as 1.0, i.e. no attenuation.
    # Use first correction polynomial at present.
    bhtf = open(bhtFile,"w")
    fwl = float(whiteLevel)
    values = 2**16
    oatt = np.ones(values)
    for l in range(values):
        if l>=whiteLevel:
            #bhtf.write("1.0\n")
            oatt[l] = 1.0
        elif l>0:
            oatt[l] = fwl/l
        else:
            oatt[l] = fwl/0.5
    oatt = np.log(oatt)
    oatt = correct(oatt,polyC[0,::-1])
    oatt = np.exp(oatt)
    for l in range(values):
        bhtf.write(str(oatt[l])+"\n")
    bhtf.close()
    print("bht data written to "+bhtFile)

#

print("Apply beam hardening corrections to image file OR generate bht file")
rows,lines,polyf,whiteLevel,bhtFile,debug,args = checkArgs(sys.argv[1:])

if debug:
    pdb.set_trace()

print("\nReading correction polynomials:")
polyCoeffs = readPolyCoeff(polyf)
if type(polyCoeffs)==type(0):
    print("failed to read polyfit file")
    sys.exit()
# polynomial is usually fitted as 8th order, but as the leading term is forced to zero, it is
# actually 9th order.
print("Number of polynomials: ",str(len(polyCoeffs[:,0]))," order: ",str(len(polyCoeffs[0,:])-1) )
nargs = len(args)

if bhtFile <> "":
    print("generating .bht file")
    if nargs>0:
        print("Ignoring file arguments")
    genbht(bhtFile,whiteLevel,polyCoeffs)
    sys.exit(0)

print("Number of input files:", nargs, " rows=",rows," lines=",lines)
if nargs==0:
    sys.exit("No files to process")
if not os.path.isfile(args[0]):
    sys.exit("File not found: "+args[0])
print("First file:", args[0])
filename,fileext = os.path.splitext(args[0])
print("Extension = ",fileext)
if fileext==".raw" and nargs==1:
    print("Processing raw file "+args[0])
    if whiteLevel>0:
        print("* Ignoring whiteLevel for .raw file")
    outfile = "bhc_"+args[0]
    processRaw(args[0],outfile,polyCoeffs,rows,lines,debug)
elif fileext==".tif":
    print("Processing .tif file(s)")
    #if whiteLevel<1:
    #    print("* Must set whiteLevel to process tif files")
    #    sys.exit(1)
    if whiteLevel>65536:
        print("WhiteLevel too high for 16 bit data")
        sys.exit(1)
    for tfile in args:
        processTif(tfile,"bhc_",polyCoeffs,rows,lines,whiteLevel,debug)
else:
    print("extension not recognised")
