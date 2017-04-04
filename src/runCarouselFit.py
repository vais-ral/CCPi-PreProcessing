""" Control file to run simple carousel analysis and return
    attenuation correction function. Based on IDL code by Anthony
    Evershed.

    This is the command line interface which is executed by:
       python ../src/runCarouselFit.py

    within the "test" directory. The location is important since the
    code assumes that spectral data is available in the "spectra"
    subdirectory which is currently in test.
    Most of the fitting operations are implemented in the file carouselUtils.py.
    This code should run on both Python 2 and Python 3.
"""
from __future__ import print_function
import sys
import os.path
import logging
import timeit
import time
import numpy as np
import itertools as it
from numpy.polynomial.polynomial import polyval
# import pyplot, checking if NOX11 env var set; if so use Agg plotting (no X11)
try:
    import matplotlib
    if 'NOX11' in os.environ:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Error: cant find matplotlib")
import carouselUtils as cu
import pdb

markerErr = it.cycle((',', '+', '.', 'o', '*'))

def loadAll(string):
    """ read both the carousel definition file and the data file with the
        calibration data """
    global carouselData, carouselCal, xSpec
    if len(string)<3:
        print("syntax: load <cardef> <carrun>")
        return

    if debug:
        pdb.set_trace()
    carouselData = cu.carousel(string[1])
    if not carouselData.isValid():
        print("** failed to load carousel data")
        return
    carouselCal = cu.carouselCalibrationData(string[2], carouselData)
    if not carouselCal.isValid():
        print("** failed to load calibration data")
        return
    xSpec = carouselCal.spec
    # set guess for spectra peak to half the maximum X-ray voltage
    startX[4] = carouselCal.voltage/2.
    if carouselCal.spec.getS() is None and vary[4]==-1:
        print("Using 'vary spectra 0' as pre-computed spectra not found")
        vary[4:] = 0
    if not xSpec.isValid():
        sys.exit("** failed to load spectrum")


def showImg(string):
    """ plot the n calibration images on one plot; user must kill
        window to continue"""
    if carouselCal == None:
        print("must load data first")
        return
    try:
        width = carouselCal.width
    except:
        print("** must read calibration data first")
        return
    plt.figure(FIG_IMG)
    carouselCal.plotCalData(False, width)
    plt.show(block=False)


def showSpec(string):
    """ plot spectra of source along with filtered spectra and response spectra.
        Note that these just use input data, not fitted data.
    """
    if carouselCal == None:
        print("must load data first")
        return
    # if no SpeKCalc data available, need to allow fitted spectra
    if carouselCal.spec.getS() is None and vary[4]==-1:
            print("** No spectrum is defined. Must use e.g. 'vary spectra 0'")
            return
    if xSpec.valid:
        plt.figure(FIG_SPEC)
        line = 0
        if len(string) > 1 :
            try:
                line = int(string[1])
            except ValueError:
                print("Unrecognised line number")
                return
            if line<0 or line>=fit.nlines:
                print("line must be >=0 and < ",fit.nlines)
                return
        if carouselCal.filterCount>0:
            n = len(xSpec.getE())
            try:
                # this is PATCHED UP to get the Gaussian spectra if this option
                # has been selected and to push this into the calculation pipeline
                # HOWEVER the above plot has not been fixed up, neither has contant
                # values for others??? Or is that done by calcWidths?
                tw,dw,fw,ec,spectra = fit.calcWidths(res,fit.nlines,xSpec.getE())
                varf = fit.varFilter
                if isinstance(spectra,np.ndarray):
                    norm = np.sum(spectra)
                    attSpec = spectra/norm
                else:
                    norm = np.sum(xSpec.getS())
                    attSpec = xSpec.getS()/norm
                plt.plot(xSpec.getE(), attSpec,label='raw')
            except:
                line = 0
                print("No fit, using preset values")
                tw= [0.]
                dw= [0.1]
                varf = 1 # maybe?
                fw= [carouselCal.filterWidth[varf]]
                ec = xSpec.getE()
            for i in range(carouselCal.filterCount):
                # allow for unphysical -ve filter thicknesses
                if varf==i:
                    expo = -carouselCal.filterAtten[i].getMu()[0:n]*fw[line]
                else:
                    expo = -carouselCal.filterAtten[i].getMu()[0:n]*carouselCal.filterWidth[i]
                expo[expo>600] = 600
                attSpec = attSpec*np.exp(expo)
                #attSpec = attSpec*np.exp(-carouselCal.filterAtten[i].getMu()[0:n]*carouselCal.filterWidth[i])
            #
            norma = np.sum(attSpec)
            attSpec = attSpec/norma
            plt.plot(xSpec.getE(),attSpec,label='filtered')
            meanE = np.sum(attSpec*xSpec.getE())
            dev2 = np.sum(attSpec*xSpec.getE()*xSpec.getE()) - meanE*meanE
            print("For filtered spectrum:")
            print("mean E =",meanE," std dev = ",np.sqrt(dev2)," total atten ratio = ",1.0/norma)
            #
            dE = xSpec.getE()[1]-xSpec.getE()[0]
            nmean = int(meanE/dE)
            # This expression is not valid if using fitted spectra
            #print(" atten ratio at mean energy = ",xSpec.getS()[nmean]/(attSpec[nmean]*norma*norm))
            #
            expo = -carouselCal.targetAtten.getMu()[0:n]*tw[line]
            expo[expo>600] = 600
            attSpec = attSpec*np.exp(expo)
            detWid = dw[line]
            attDet = detWid*carouselCal.detectorAtten.getMu()[:len(attSpec)]
            resSpec = attSpec*ec*(1.-np.exp(-attDet))
            resSpec = resSpec/np.sum(resSpec)
            plt.plot(xSpec.getE(),resSpec,label='response')
            meanE = np.sum(resSpec*xSpec.getE())
            dev2 = np.sum(resSpec*xSpec.getE()*xSpec.getE()) - meanE*meanE
            print("For response spectrum:")
            print("mean E =",meanE," std dev = ",np.sqrt(dev2)," total atten ratio = ",1.0/norma)
            #
            if carouselCal.whiteLevel>0:
                print("Calculation used whiteLevel = ",carouselCal.whiteLevel)

        #else:
        #    print("string=",string)
        plt.legend()
        plt.xlabel('KeV')
        plt.ylabel('S(E) (normalised)')
        plt.draw()
        plt.show(block=False)
    else:
        print("must load data first")

def quitCarousel(string):
    """ exit the command level """
    sys.exit(0)

def showAtt(string):
    """ 1D plots of attenuation of sample n"""
    if carouselCal == None:
        print("must load data first")
        return
    defline=400
    plt.figure(FIG_ATT1D)
    if len(string)>1:
        try:
            samp = int(string[1])
        except:
            print("failed to read sample value")
            return
        if len(string)>2:
            try:
                defline = int(string[2])
            except:
                print("failed to read line number")
        if samp > -1 and samp < carouselCal.samples and defline > -1 and defline < carouselCal.lines:
            z = carouselCal.getImage(samp)
            plt.plot(z[defline,:])
            plt.xlabel("Column number at line "+str(defline))
            plt.ylabel("attenuation")
            plt.draw()
            plt.show(block=False)
        else:
            print("sample number out of range")
    else:
        for i in range(carouselCal.samples):
            z = carouselCal.getImage(i)
            plt.plot(z[defline,:])
        plt.xlabel("Column number at line "+str(defline))
        plt.ylabel("attenuation")
        plt.draw()
        plt.show(block=False)
    return

def showCor(string):
    """ plot the fitted correction curve from the polynomial data """
    if not 'xtab' in globals():
        print("No correction data available; run fitatt first")
        return
    if len(string)==1:
        linevals = [ 0, int((len(xtab)-1)/2), len(xtab)-1 ]
    else:
        linevals = []
        for i in range(len(string)-1):
            try:
                lineno=int(string[i+1])
            except:
                print("Bad integer value: ",string[i+1])
                return
            if lineno>=0 and lineno<=len(xtab):
                linevals.append(lineno)
            else:
                print("skipping out of range value: ",lineno)
    # check if we actual have more than one line fitted
    if len( polyfit[:,0] ) > 1 :
        print("lines=",linevals)
    else:
        print("Only 1 line fitted to data; using line 0")
        linevals = [0]
    #print(ytab[::30])
    #for line in linevals:
    #    xvals = polyval( xtab[0,:], polyfit[line,::-1])
    #    print("pv=",polyfit[line,::-1])
    #    print(xvals[::30])
    #    print(" ")
    mymax=np.max(xtab[:,-1])
    xval = np.linspace(0.0,mymax,num=300)
    plt.figure(FIG_COR)
    for line in linevals:
        yval = polyval( xval,polyfit[line,::-1])
        plt.plot(xtab[line,:],ytab,xval,yval,':')
        plt.xlabel('Observed log attenuation ratio')
        plt.ylabel('Effective single energy log attenuation ratio')
    # add the x=y line for comparison with circles for fit points
    nsamp = len(carouselData.mask)-1
    xarr = np.zeros(nsamp)
    count = 0
    for i in range(nsamp):
        if not carouselData.mask[i]:
            xarr[count] = carouselCal.getAvAtten(linevals[0],i)
            count = count+1
            #print("mask ",i,"t",carouselCal.getAvAtten(linevals[0],i))
    plt.plot([0.0,mymax],[0.0,mymax],'r--')
    plt.plot(xarr,xarr,'ro')
    plt.draw()
    plt.show(block=False)

def setFilter(string):
    """ List filters defined or change settings of an existing filter. Can not at
        present add a new filter """
    try:
        if carouselCal.isValid():
            #print("have carouselCal - str: ",string)
            if len(string) == 3:
                #print("try and set ",string[1:])
                try:
                    mat = string[1]
                    val = float(string[2])
                except:
                    print("** Bad values")
                    return
                for i in range(carouselCal.filterCount):
                    if carouselCal.filterMaterial[i] == mat:
                        carouselCal.filterWidth[i] = val
                        print("set ",mat," filter width to ",val)
                        return
                print("filter type not found")
            else:
                print("Filters set:")
                print("filter, width:")
                for i in range(carouselCal.filterCount):
                    print(carouselCal.filterMaterial[i], ",", carouselCal.filterWidth[i])
    except:
        print("no carousel data file loaded")


def fitAtt(string):
    """ Check necessary data has been set then fit model to carousel data.
        Finally generate curves for attenuation over each line using the
        correction material (corMat/corEn) and then fit a polynomial to this for
        correction purposes.
        """
    global res,xtab,ytab,fit,polyfit,xpolyfit
    # lstep is the line step; e.g. 1 for every line, 2 for every other line in fitting
    lstep = 1
    if carouselData == None or carouselCal == None:
        print("must load data first")
        return
    if not carouselData.valid or not carouselCal.valid:
        print("data not correctly loaded")
        return
    if corMat.name=="":
        print(" ** Must define corrrection material and energy using 'setcormat'")
        return
    # if no SpeKCalc data available, need to allow fitted spectra
    if carouselCal.spec.getS() is None and vary[4]==-1:
            print("** No spectrum is defined. Must use e.g. 'vary spectra 0'")
            return
    defMat = "Cu"
    fit = cu.fitData(carouselData, carouselCal, defMat)
    fit.verbose = debug
    if debug:
        np.seterr(over='warn',invalid='warn')
    else:
        np.seterr(over='ignore',invalid='ignore')

    if np.max(vary)<0:
        print("** Error: no parameters to fit, check setvary")
        return
    x = np.zeros(7+np.sum(vary))
    if len(string) == 2 or len(string) == 3:
        print("Fitting variables: ",np.sum(vary)+len(vary))
        # The fit function consists of 3 polynomial expressions in the
        # the line number, plus a possible polynomial in the energy.
        # Initial values for the zero order terms are
        # given here, the higher terms (if any) are set to zero.
        # Updated to allow any of the variables to be excluded from the fit (-1)
        # In this case the initial value, in startX, should be used, which is passed
        # to fit.
        offset = vary[0]
        if vary[0]>-1:
            x[offset] = startX[0]
        offset = offset+1+vary[1]
        if vary[1]>-1:
            x[offset] = startX[1]
        offset = offset+1+vary[2]
        if vary[2]>-1:
            x[offset] = startX[2]
        offset = offset+2+vary[3]+vary[4]
        if vary[4]>-1:
            x[offset] = startX[4]
        offset = offset+1+vary[5]
        if vary[5]>-1:
            x[offset] = startX[5]
        offset = offset+1+vary[6]
        if vary[6]>-1:
            x[offset] = startX[6]
        fit.defaults = startX
        try:
            nlines = int(string[1])
            if len(string) == 3:
                lstep = int(string[2])
        except:
            print("Wrong arguments: fitatt nlines [linestep]")
            return
    else:
        print("wrong number of args: need fitatt nlines [linestep]")
        print("where nlines=number of lines to fit and lstep is step between")
        print("lines, default 1")
        return
    if nlines < 1 or nlines > carouselCal.lines:
        print("fit lines out of range, must be 1 to ",carouselCal.lines)
        return
    fit.vary_target = vary[0]
    fit.vary_detector = vary[1]
    fit.vary_filter = vary[2]
    fit.vary_energy = vary[3]
    fit.vary_epk = vary[4]
    fit.vary_ewidlow = vary[5]
    fit.vary_ewidhigh = vary[6]

    fit.solver = solverChoice

    t0 = timeit.default_timer()

    try:
        res,cov,infodict,mesg,ier = fit.dofit(nlines,lstep,x)
    except Exception as experr:
        print("** Fit failed due to exception: ",experr)
        return
    
    tim = timeit.default_timer()-t0

    print("time=",tim)
    print("dofit returned: ")
    print(" best fit values = ",res)
    if ier>0 and ier<5:
        print(" ier = ",ier)
    else:
        print("** Fit failed: ier=",ier)
        print("   message=",mesg)
    print(" iterations = ",infodict["nfev"])
    # measure error
    samples = carouselCal.samples
    ofile = open('fit.log','w')
    ofile.write('time={0:12.6f}\n'.format(tim))
    ofile.write('dofit returned: ')
    ofile.write(' best fit values = \n')
    ofile.write(str(res)+'\n')
    ofile.write(' ier = {0:5d}\n'.format(ier))
    ofile.write(' iterations = {0:5d}\n'.format(infodict["nfev"]))
    ofile.write(' cov = ')
    ofile.write(str(cov)+'\n')
    ofile.write(' mesg = '+mesg+'\n')
    rfile = open('param.log','w')
    rfile.write('fit input: lines={0:5d}\n'.format(nlines))
    rfile.write('guess: ')
    rfile.write(str(x))
    rfile.write('\n')
    rfile.write('solution: ')
    rfile.write(str(res))
    rfile.write('\n')
    sumtot=0.
    summax=0.

    # calculate the attenuation(corEn) vs attenuation(observed) and return
    # polynomial fit to these curves for each line.
    attLnWid = 14.0
    attPts = 300
    xtab,ytab,polyfit,xpolyfit = fit.linesPolyFit(res,corMat,corEn,attPts,attLnWid)
    # find average and max error for each line
    rfile.write('polyfits '+str(polyfit.shape[0])+" "+str(polyfit.shape[1])+"\n")
    lsumsq = []
    avatt = np.zeros((2,samples))
    for line in range(nlines):
        sumsq = 0.
        for sample in range(samples):
            if carouselData.mask[sample]:
                continue
            sumsq += (fit.atten[line,sample] - carouselCal.getAvAtten(line,sample) ) ** 2
            if line==0:
                avatt[0,sample] = carouselCal.getAvAtten(line,sample)
                avatt[1,sample] = fit.atten[line,sample]
        ofile.write(' {0:5d}  {1:12.6f}\n'.format(line,sumsq))
        lsumsq.append(sumsq)
        sumtot += sumsq
        if sumsq>summax:
            summax = sumsq
        if debug:
            print("Line: ",line," ave error:",sumsq)
        # save poly data to param.log
        if len(polyfit[:,0])>line:
            rfile.write('{0:5d} '.format(line)+str(polyfit[line,:])+'\n')
    # write data for xtek
    rfile.write('xpoly:\n')
    for line in range(len(xpolyfit[:,0])):
        rfile.write('{0:5d} '.format(line)+str(xpolyfit[line,:])+'\n')
    # write data in binary file
    bfile = open("polyfit.npz","wb")
    np.save(bfile,polyfit)
    bfile.close()
    #
    print("average error: ",sumtot/nlines)
    print("max error: ",summax)
    rfile.write("average error: {0:12.6f}\nmax error: {1:12.6f}\n".format(sumtot/nlines,summax))
    try:
        plt.figure(FIG_ERR)
        plt.plot(lsumsq)
        plt.xlabel('line number')
        plt.ylabel('mean sq error')
        plt.draw()
        plt.show(block=False)
        #
        plt.figure(FIG_ATTCOMP)
        nsamp = len(avatt[0,:])
        xnum = np.array(range(nsamp))+1
        plt.subplot(211)
        plt.plot(xnum,avatt[0,:],'bx')
        mark = next(markerErr)
        plt.plot(xnum,avatt[1,:],marker=mark)
        plt.xlabel('sample number')
        plt.ylabel('log(I0/I)')
        plt.subplot(212)
        plt.plot(xnum,avatt[0,:]-avatt[1,:],marker=mark)
        plt.ylabel('err log(I0/I)')
        plt.draw()
        plt.show(block=False)
    except:
        print("Plotting failed")
    ofile.close()
    rfile.close()

def initGuess(words):
    """ Set initial values to use for the variables of the target absortion width, detector
    width and filter width """
    try:
        startX[0] = float(words[1])
        startX[1] = float(words[2])
        startX[2] = float(words[3])
        if len(words)==8:
            startX[3] = float(words[4])
            startX[4] = float(words[5])
            startX[5] = float(words[6])
            startX[6] = float(words[7])
        elif len(words)>4:
            print("Some values ignored!")
    except:
        print("initguess requires 3 or 7 floating point values: dt, ln(dd) and df. dt is the target "+
              "absortion width, ln(dd) is the log of detector width, and df the filter width.")
        print("If 7 values used additional values are: energy coeff, spectral peak, low/high width")
        print("Current guess = ",startX[0:7])


def setWidth(words):
    """ set the half width of area along row to be averaged"""
    if len(words)>1:
        try:
            width = float(words[1])
            carouselCal.width = width
            carouselCal.setWidthAve(width)
        except:
            print("load carousel data before setting width")
            return
    else:
        try:
            print("width= ",carouselCal.width, " (No. of pixels about centre of line to average)")
        except:
            print("width not set until carousel data loaded")


def showCalConf(string):
    """ prints out some calibration data"""
    try:
        filterCount = carouselCal.filterCount
    except NameError:
        print("** must read calibration data first")
        return
    print("filter count = ", filterCount)
    print("filter  material width    density   fittable")
    for i in range(filterCount):
        print('{0:4d}     {1:7s} {2:7f} {3:7f}'.format(i,
               carouselCal.filterMaterial[i],
               carouselCal.filterWidth[i], carouselCal.filterDensity[i]))
    print("Detector:")
    print('         {0:7s} {1:7f} {2:7f}'.format(
               carouselCal.detectorMaterial,
               carouselCal.detectorWidth, carouselCal.detectorDensity))
    print("Source filter:")
    print('         {0:7s} {1:7f} {2:7f}'.format(carouselCal.targetMat,
               0.0, carouselCal.targetDensity))
    print("Voltage=",  carouselCal.voltage, " angle=", carouselCal.angle)

def setVary(strlst):
    """ define polynomial order of parameters to fit """
    if len(strlst)==1:
        print("Control order of polynomial used for fitting across lines")
        print(" - 3 widths are fitted for target, detector and filter")
        print(" - Experimental is to set energy dependence away from linear")
        print(" - If using energy, suggest only 0 order and check results")
        print(" - using -1 implies do not vary parameter, ONLY FOR ENERGY")
        print("current settings:")
        print("target: ",vary[0])
        print("detector: ",vary[1])
        print("filter: ",vary[2])
        print("energy dependence: ",vary[3])
        print("spectra: ",vary[4])
        return
    if len(strlst)==3:
        try:
            npo = int(strlst[2])
        except:
            print("** failed to read npo")
            npo = 0
        if npo<-1 or npo>3:
            print("** Order must be range -1 to 3")
            return
        if strlst[1]=="target":
            vary[0] = npo
        elif strlst[1]=="detector":
            vary[1] = npo
        elif strlst[1]=="filter":
            vary[2] = npo
        elif strlst[1]=="energy":
            vary[3] = npo
        elif strlst[1]=="spectra":
            if npo>0:
                print("Error: spectra can only take values -1 and 0 at present")
                return
            vary[4] = npo
            vary[5] = npo
            vary[6] = npo
        else:
            print("Not recognised: ",strlst[1])
    else:
        print("syntax: vary [target|detector|filter|energy n]")

def debugToggle(cmd):
    """ toggle debug """
    global debug
    debug = not debug
    print("debug set to ",debug)


def helpCar(cmd, string):
    """ just print list of commands"""
    print("Carousel calibration utility")
    print(" ")
    print("cmds:")
    for i in cmd:
        print("  ", i)
    print(" ")
    print("To execute script file use: read <filename>")
    print(" ")
    print("The required input is a set of images of test pieces imaged at the exact same")
    print("Xray settings as the sample. These may be in a carousel or a crown device. The")
    print("details of each test item (material/density/thickness) must be provided in a")
    print("file in './carouselData'. Using these images a fit is made to the effective")
    print("energy response function of the Xray source/detector. Using this fit a correction")
    print("curve can be determined to map from observed attenuation to true attenuation of the")
    print("dominant Xray material of the sample. This will be inaccurate on samples with muliple")
    print("material types in them. The output is a file of polynomial fits giving the corrections")
    print("which can be used in the python program 'applyTransForm.py'")
    print(" ")

def setCorMat(words):
    """ Input the name of the material that will be the target of the attenuation correction
        e.g. Calcium hydroxyapatite which might be defined in a file cahypa.txt, precalculated
        using the program xcom. Without arguments, list current setting, if any."""
    global corMat,corEn
    if len(words)>2:
        name = words[1]
        try:
            corEn = float(words[2])
        except:
            print("failed to read energy value")
            return
    else:
        if corMat.name != "":
            print("corMat is set to ",corMat.name," and energy ",corEn,"KeV")
        else:
            print("corMat is not set")
        return
    name = words[1]
    print("reading correction material definition from file: ",name)
    try:
        corMat = cu.materialAtt(name,1.0)
    except:
        print("error reading material type")

def mask(words):
    """ show or set a mask array which is used to define if some of the sample data
        is not to be used in the fit. e.g. "mask 7" will mean sample 7 is omitted
        from the fit. "mask off" will restore all samples to the fit.
    """
    if carouselData== None:
        print("must load data first")
        return
    if len(words) == 1:
        print("Mask = ",carouselData.mask)
    elif words[1] == "off":
        carouselData.mask.fill(False)
    else:
        try:
            for i in range(len(words)-1):
                num = int(words[i+1])
                if num > 0:
                    carouselData.mask[num-1] = True
                elif num < 0:
                    carouselData.mask[-num-1] = False
                else:
                    print("Warning: label 0 ignored")
        except:
            print("Error: bad value in list")

def transform(words):
    """ TEST code to try and fix up some data; should get right data
        in first place, Assume we have "I" data and want log(I0/I) values, where I0 is
        provided by the user.
    """
    if len(words)<2:
        print("transform command is a test command. Was used to map image intensity")
        print(" data I to log(I0/I), in the case where only I is provided. Now redundant")
        print(" as uint16 data is assumed to include I0 image at start and transform applied.")
        return
    I0 = float(words[1])
    nsamp = len(carouselData.mask)-1
    for i in range(nsamp):
        z = carouselCal.getImage(i)
        zerocount = sum(sum(z<=0))
        if zerocount>0:
            print("zeros for sample ",i," are ",zerocount)
        z[z<=0] = 1e-4
        z = np.log(z)
        z = np.log(I0) - z
        carouselCal.getImage(i)[:,:] = z

def setOptions(words):
    """ Set options controlling the fit process. Currently just select between the new and old
        versions of the least squares solver.
    """
    global solverChoice
    if len(words)<2:
        print("Options:")
        print(" solver = ",solverChoice)
    elif words[1]=="solver=new":
        solverChoice="new"
    elif words[1]=="solver=old":
        solverChoice="old"
    else:
        print("Option not recognised")
       

def checkVersions():
    """ Check version of matplotlib and exit if too old
    """
    import matplotlib as mpl
    def versiontuple(v):
        """ split . separated string for comparison
        """
        return tuple(map(int, (v.split("."))))
    if versiontuple(mpl.__version__) < versiontuple('1.1'):
        print("Matplotlib version too old, need at least 1.1; have ",mpl.__version__)
        sys.exit(1)
        


# Set of commands that are implemented and the corresponding function names.
# Additional commands and functions can be added here.
cmd_switch = { "load":loadAll,
               "showimg":showImg,
               "showspec":showSpec,
               "showconf":showCalConf,
               "showatt":showAtt,
               "showcor":showCor,
               "setfilter":setFilter,
               "setwidth":setWidth,
               "fitatt":fitAtt,
               "initguess":initGuess,
               "vary":setVary,
               "setcormat":setCorMat,
               "debug":debugToggle,
               "mask":mask,
               "quit":quitCarousel,
               "help":helpCar,
               "transform":transform,
               "setoptions":setOptions,
               }

# set figures to use for different plots
FIG_COR = "Correction"
FIG_ATT1D = "CarouselLines"
FIG_SPEC = "Spectra"
FIG_IMG = "CarouselImages"
FIG_ERR = "ErrorInFit"
FIG_ATTCOMP = "ObservedVsFittedAtten"
# simple command line loop to allow loading of data and run
# of fitting code.

carouselData = None
carouselCal = None
# initial guesses unless user changes them
startX = np.array([0.01,-6.0,0.01,0.,30.,0.05,0.05])
checkVersions()

if __name__ == "__main__":
    logging.basicConfig(filename='runCarouselFit.log',level=logging.INFO)
    localtime = time.asctime( time.localtime(time.time()) )
    logging.info('started at '+localtime)

    nargs = len(sys.argv)
    debug = False
    carouselCal = None

    print(" ")
    print(" *** Carousel data fitting program ***")
    print(" ")
    print(" Code based on algorithms developed at Queen Mary University")
    print(" of London by Graham Davis, Anthony Evershed et al")
    print(" This implementation by Ron Fowler at STFC")
    print(" contact ronald.fowler@stfc.ac.uk")
    print(" Funded by the CCPi project")
    print(" ")
    # set the polynomial order of the fitting variables. Variables are
    # function of line number.
    vary = np.zeros(7,dtype='int')
    vary[3:] = -1
    # set the default solver type to be "old"
    solverChoice = "old"
    # set an object for the material to which attenuation is to be corrected to;
    # this is null until the user provides one
    corMat = cu.materialAtt("",1.0)
    # command loop
    filein = False
    while True:
        try:
            if filein:
                cmd = infile.readline().strip()
                if cmd=='':
                    filein = False
                    infile.close()
                    continue
                print(" bhc: ",cmd)
            else:
                #if sys.version_info.major == 2:
                if sys.version_info[0] == 2:
                    cmd = raw_input("bhc: ").strip()
                else:
                    cmd = input("bhc: ").strip()
        except EOFError as ex:
            logging.info('bhc: EOF')
            print("EOF")
            sys.exit(0)
        logging.info('bhc: '+cmd)
        words = cmd.split(" ")
        try:
            if words[0] == "help":
                cmd_switch[words[0]](cmd_switch, words)
            elif words[0] == "read":
                if len(words)>1 and os.path.isfile(words[1]):
                    rfile = words[1]
                    infile = open(rfile,'r')
                    filein = True
                else:
                    print("Error - syntax: read file")
                    continue
            elif words[0] == "#":
                continue
            else:
                cmd_switch[words[0]](words)
        except SystemExit as ex:
            sys.exit(0)
        except KeyError as ex:
            if not ( words[0] == "" ):
                print("** command not found: ", words[0])
        except:
            print("** internal error: ", sys.exc_info()[0])
            raise

