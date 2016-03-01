#
# simple program to read E,S(E) and calculate mean and variance of energy with
# various filers. Input parameters are energy and take-off angle, the filter material
# name density and thickness and the thickness of th CsI detector.
# e.g
# pythom meanrms.py  50 1 Sn 7.3 0.4 0.3
# for 50KeV beam at 1deg take off through Sn filter of 7.3 density and .4 width (cgs)
# with a .3 thickness of CsI. At the moment an additional 1mm of Al is used as a pre
# filter to prove that beam softening is possible under some conditions.
# This program must be run in the "test" directory of the bhc code to ensure it
# has access to the spectra data for W and the xcom data for attenuation.
#
import os
import sys
import numpy as np

sums=0.
sume=0.
sume2=0.
# default if nothing specified on cmd line, but other parts will fail.
filename="spectra/W/070220.spc"
if len(sys.argv)>2:
    if int(sys.argv[2])>9:
        if int(sys.argv[1])<100:
            filename="spectra/W/0"+str(sys.argv[1])+str(sys.argv[2])+"0.spc"
        else:
            filename="spectra/W/"+str(sys.argv[1])+str(sys.argv[2])+"0.spc"
    else:
        if int(sys.argv[1])<100:
            filename="spectra/W/0"+str(sys.argv[1])+"0"+str(sys.argv[2])+"0.spc"
        else:
            filename="spectra/W/"+str(sys.argv[1])+"0"+str(sys.argv[2])+"0.spc"

print "filename=",filename
if not os.path.isfile(filename):
    print "file not found: ", filename, " in specData"
else:
    with open(filename, 'r') as fl:
        en, amp = np.loadtxt(fl, unpack=True)
if len(sys.argv)>3:
    filename="xcom/"+sys.argv[3]+".txt"
    density=float(sys.argv[4])
    width=float(sys.argv[5])
    if not os.path.isfile(filename):
        print "file not found: ", filename, " in xcom"
    else:
        with open(filename, 'r') as fl:
            en2, mu = np.loadtxt(fl, unpack=True)
        slen=len(en)
        mut=mu[0:slen]*density*width

if True:
    # add Al filter
    filename="xcom/Al.txt"
    density=2.64
    width=2.0
    if not os.path.isfile(filename):
        print "file not found: ", filename, " in xcom"
    else:
        with open(filename, 'r') as fl:
            ena, mua = np.loadtxt(fl, unpack=True)
        slen=len(en)
        muat=mua[0:slen]*density*width
    print("applying AL pre filter at "+str(width)+"cm")
    amp=amp*np.exp(-muat)

filename="xcom/CsI.txt"
if not os.path.isfile(filename):
    print "file not found: ", filename, " in xcom"
else:
    with open(filename, 'r') as fl:
        en3, mucsi = np.loadtxt(fl, unpack=True)
    slen=len(en)
    density=4.51
    if len(sys.argv)>5:
        width=float(sys.argv[6])
    else:
        width=0.02
    print "width=",str(width)
    mucsit=mucsi[0:slen]*density*width
sumAmp = sum(amp)
sumE = sum(en*amp)
sumE2 = sum(en*en*amp)
meanE = sumE/sumAmp
var2 = sumE2/sumAmp-meanE*meanE
print "file = ",filename
print "mean E = ",str(meanE)
print "var = ",str(var2)
print "RMS = ",str(np.sqrt(var2))

ampF=amp*np.exp(-mut)
sumAmpF = sum(ampF)
sumEF = sum(en*ampF)
sumEF2 = sum(en*en*ampF)
meanEF = sumEF/sumAmpF
varF2 = sumEF2/sumAmpF-meanEF*meanEF
print "mean E F = ",str(meanEF)
print "var F = ",str(varF2)
print "RMS F = ",str(np.sqrt(varF2))
print "ratio of Amps = ",str(sumAmp/sumAmpF)

ampR = ampF*en*(1-np.exp(-mucsit))
sumAmpR = sum(ampR)
sumER = sum(en*ampR)
sumER2 = sum(en*en*ampR)
meanER= sumER/sumAmpR
varR2 = sumER2/sumAmpR-meanER*meanER
print "mean E R = ",str(meanER)
print "var F = ",str(varR2)
print "RMS F = ",str(np.sqrt(varR2))
print "ratio of Amps = ",str(sumAmp/sumAmpR)
