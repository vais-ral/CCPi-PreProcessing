""" Control file to run simple carousel analysis and return
    attenuation correction function. Based on IDL code by Anthony
    Evershed.
"""
import sys
import logging
import timeit
import numpy as np
from numpy.polynomial.polynomial import polyval
from carouselUtils import *
import pdb

def getAtt(dw,dcu,dal,dti,dcsi,me):
    """ test function to evaluate attenuation """
    cuatt = np.zeros(me)
    alatt = np.zeros(me)
    watt = np.zeros(me)
    csiatt = np.zeros(me)
    tiatt = np.zeros(me)
    cuatt[0:me-1] = cuAtten.getMu()[0:me-1]
    tiatt[0:me-1] = tiAtten.getMu()[0:me-1]
    alatt[0:me-1] = alAtten.getMu()[0:me-1]
    watt[0:me-1] = wAtten.getMu()[0:me-1]
    csiatt[0:me-1] = csiAtten.getMu()[0:me-1]
    at_se = se*np.exp(-cuatt*dcu-watt*dw-alatt*dal)*xe*(1-np.exp(-csiatt*dcsi))
    return at_se


def loadAll(string):
    """ read both the carousel definition file and the data file with the
        calibration data """
    global carouselData, carouselCal, xSpec
    if len(string)<3:
        print "syntax: load <cardef> <carrun>"
        return
    
    carouselData = carousel(string[1])
    if not carouselData.isValid():
        print "** failed to load carousel data"
        return
    carouselCal = carouselCalibrationData(string[2], carouselData)
    if not carouselCal.isValid():
        print "** failed to load calibration data"
        return
    xSpec = specData(carouselCal.targetMat, carouselCal.voltage, carouselCal.angle)
    if not xSpec.isValid():
        sys.exit("** failed to load spectrum")

 
def showImg(string):
    """ plot the n calibration images on one plot; user must kill
        window to continue"""
    try:
        width = carouselCal.width
    except NameError:
        print "** must read calibration data first"
        return
    carouselCal.plotCalData(True, width)


def showSpec(string):
    """ plot spectra """
    if xSpec.valid:
        yval = np.zeros(xSpec.getS().size)
        norm = np.sum(xSpec.getS())
        yval = xSpec.getS()/norm
        plt.plot(xSpec.getE(), yval,label='raw')
        if carouselCal.filterCount>0:
            n = len(xSpec.getS())
            attSpec = np.copy(yval)
            for i in range(carouselCal.filterCount):
                attSpec = attSpec*np.exp(-carouselCal.filterAtten[i].getMu()[0:n]*carouselCal.filterWidth[i])
                #print "attn0-1=",carouselCal.filterAtten[i].getMu()[0:9]
                #print "width=",carouselCal.filterWidth[i]
            norm = np.sum(attSpec)
            attSpec = attSpec/norm
            plt.plot(xSpec.getE(),attSpec,label='filtered')
            meanE = np.sum(attSpec*xSpec.getE())
            dev2 = np.sum(attSpec*xSpec.getE()*xSpec.getE()) - meanE*meanE
            print "mean E =",meanE," std dev = ",np.sqrt(dev2)
        if len(string) > 1 and string[1] == 'log':
            plt.yscale('log')
        #else:
        #    print "string=",string
        plt.legend()
        plt.xlabel('KeV')
        plt.ylabel('S(E) (normalised)')
        plt.draw()
        plt.show()
    else:
        print "must load data first"
        
def quitCarousel(string):
    """ exit the command level """
    sys.exit("quit called")

def showAtt(string):
    """ 1D plots of attenuation of sample n"""
    if len(string)>1:
        samp = int(string[1])
        print "samp: ",samp
        defline=400
        if samp > -1 and samp < carouselCal.samples:
            z = carouselCal.getImage(samp)
            plt.plot(z[defline,:])
            plt.draw()
            plt.show()
            return
    print "syntax: showatt n"
    print "where n is the carousel number range 0:",carouselCal.samples-1

def showCor(string):
    """ plot the fitted correction curve from the polynomial data """
    if not 'xtab' in globals():
        print "No correction data available; run fitatt first"
        return
    if len(string)==1:
        linevals = [ 0, int((len(xtab)-1)/2), len(xtab)-1 ]
    else:
        linevals = []
        for i in range(len(string)-1):
            try:
                lineno=int(string[i+1])
            except:
                print "Bad integer value: ",string[i+1]
                return
            if lineno>=0 and lineno<=len(xtab):
                linevals.append(lineno)
            else:
                print "skipping out of range value: ",lineno
    print "lines=",linevals
    #print ytab[::30]
    #for line in linevals:
    #    xvals = polyval( xtab[0,:], polyfit[line,::-1])
    #    print "pv=",polyfit[line,::-1]
    #    print xvals[::30]
    #    print " "
    mymax=np.max(xtab[:,-1])
    xval = np.linspace(0.0,mymax,num=300)
    for line in linevals:
        yval = polyval( xval,polyfit[line,::-1])
        plt.plot(xtab[line,:],ytab,xval,yval,':')
    # add the x=y line for comparison
    plt.plot([0.0,mymax],[0.0,mymax],'r--')
    plt.draw()
    plt.show()

def setFilters(string):
    try:
        if carouselCal.isValid():
            #print "have carouselCal - str: ",string
            if len(string) == 3:
                #print "try and set ",string[1:]
                try:
                    mat = string[1]
                    val = float(string[2])
                except:
                    print "** Bad values"
                    return
                for i in range(carouselCal.filterCount):
                    if carouselCal.filterMaterial[i] == mat:
                        carouselCal.filterWidth[i] = val
                        print "set ",mat," filter width to ",val
                        return
            else:
                print "filter, width:"
                for i in range(carouselCal.filterCount):
                    print carouselCal.filterMaterial[i], ",", carouselCal.filterWidth[i]
    except:
        print "no carousel data file loaded"

def calcAtt(string):
    pass

def fitAtt(string):
    """ Check necessary data has been set then fit model to carousel data.
        Finally generate curves for attenuation over each line using the
        correction material (corMat/corEn) and then fit a polynomial to this for
        correction purposes.
        """
    global carouselData, carouselCal, xSpec, debug, vary
    global res,xtab,ytab,polyfit
    defMat = "Cu"
    fit = fitData(carouselData, carouselCal, defMat)
    fit.verbose = debug
    if debug:
        np.seterr(over='warn',invalid='warn')
    else:
        np.seterr(over='ignore',invalid='ignore')

    if len(string) == 5:
        print "Fitting variables: ",np.sum(vary)+4
        # The fit function consists of 3 polynomial expressions in the
        # the line number, plus a possible polynomial in the energy.
        # Initial values for the zero order terms are
        # given here, the higher terms (if any) are set to zero.
        x = np.zeros(4+np.sum(vary))
        offset = vary[0]
        x[offset] = float(string[2])
        offset = offset+1+vary[1]
        x[offset] = float(string[3])
        offset = offset+1+vary[2]
        x[offset] = float(string[4])
        nlines = int(string[1])
    else:
        print "wrong number of args: need fitatt nlines x1 x2 x3"
        print "where nlines=number of lines to fit and x1/2/3 are initial values"
        print "for the unknowns: target width, ln(detector width), filter width"
        return
    fit.vary_target = vary[0]
    fit.vary_detector = vary[1]
    fit.vary_filter = vary[2]
    fit.vary_energy = vary[3]
    t0 = timeit.default_timer()
    res,cov,infodict,mesg,ier = fit.dofit(nlines,x)
    tim = timeit.default_timer()-t0
    print "time=",tim
    print "dofit returned: "
    print " best fit values = ",res
    print " ier = ",ier
    print " iterations = ",infodict["nfev"]
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
    #if debug:
    #    pdb.set_trace()
    # calcualte the attenuation(corEn) vs attenuation(observed) and return
    # polynomial fit to these curves for each line.
    xtab,ytab,polyfit = fit.linesPolyFit(res,corMat,corEn,300,12.0)
    # find average and max error for each line
    rfile.write('polyfits:')
    for line in range(nlines):
        sumsq = 0.
        for sample in range(samples):
            sumsq += (fit.atten[line,sample] - carouselCal.getAvAtten(line,sample) ) ** 2
        ofile.write(' {0:5d}  {1:12.6f}\n'.format(line,sumsq))
        sumtot += sumsq
        if sumsq>summax:
            summax = sumsq
        if debug:
            print "Line: ",line," ave error:",sumsq
        rfile.write(str(polyfit[line,:])+'\n')
    print "average error: ",sumtot/nlines
    print "max error: ",summax
    rfile.write("average error: {0:12.6f}\nmax error: {1:12.6f}\n".format(sumtot/nlines,summax))
    print "writing polynomial correction files"
    #(tableRes,polyFit) = fit.linesPolyFit(x,corMat,corEn,300,12.0)
    #print "result="
    #for i in range(len(ytab)):
    #    if len(xtab) == 2:
    #        print xtab[0,i],xtab[1,i],ytab[i]
    #    else:
    #        print xtab[0,i],ytab[i]
    print "polyfit="
    print polyfit
    rfile.write('polyfits:')
    
    ofile.close()
    rfile.close()


def setWidth(words):
    """ set the half width of area along row to be averaged"""
    if len(words)>1:
        try:
            width = float(words[1])
        except:
            print "value not recognised"
            return
        carouselCal.width = width
    else:
        print "width= ",carouselCal.width


def showCalConf(string):
    """ prints out some calibration data"""
    try:
        filterCount = carouselCal.filterCount
    except NameError:
        print "** must read calibration data first"
        return
    print "filter count = ", filterCount
    print "filter  material width    density   fittable"
    for i in range(filterCount):
        print  '{0:4d}     {1:7s} {2:7f} {3:7f}'.format(i,
               carouselCal.filterMaterial[i],
               carouselCal.filterWidth[i], carouselCal.filterDensity[i])
    print "Detector:"
    print '         {0:7s} {1:7f} {2:7f}'.format(
               carouselCal.detectorMaterial,
               carouselCal.detectorWidth, carouselCal.detectorDensity)
    print "Source filter:"
    print '         {0:7s} {1:7f} {2:7f}'.format(carouselCal.targetMat,
               0.0, carouselCal.targetDensity)
    print "Voltage=",  carouselCal.voltage, " angle=", carouselCal.angle

def setVary(strlst):
    """ define polynomial order of parameters to fit """
    global vary
    if len(strlst)==1:
        print "settings:"
        print "target: ",vary[0]
        print "detector: ",vary[1]
        print "filter: ",vary[2]
        print "energy dependence: ",vary[3]
        return
    if len(strlst)==3:
        try:
            np = int(strlst[2])
        except:
            print "** failed to read np"
            np = 0
        if strlst[1]=="target":
            vary[0] = np
        elif strlst[1]=="detector":
            vary[1] = np
        elif strlst[1]=="filter":
            vary[2] = np
        elif strlst[1]=="energy":
            vary[3] = np
        else:
            print "Not recognised: ",strlst[1]
    else:
        print "syntax: vary [target|detector|filter|energy n]"

def debugToggle(cmd):
    """ toggle debug """
    global debug
    debug = not debug
    

def helpCar(cmd, string):
    """ just print list of commands"""
    print "Carousel calibration utility"
    print " "
    print "cmds:"
    for i in cmd:
        print "  ", i
    print " "

def setCorMat(words):
    """ Inout the name of the material that will be the target of the attenuation correction
        e.g. Calcium hydroxyapatite which might be defined in a file cahypa.txt, precalculated
        using the program xcom."""
    global corMat,corEn
    if len(words)>2:
        name = words[1]
        corEn = float(words[2])
    else:
        if corMat.isValid:
            print "corMat is set"
        else:
            print "corMat is not set"
        return
    name = words[1]
    print "reading correction material definition from file: ",name
    try:
        corMat = materialAtt(name,1.0)
    except:
        print "error reading material type"
    

# Set of commands that are implemented and the corresponding function names.
# Additional commands and functions can be added here.
cmd_switch = { "load":loadAll,
               "showimg":showImg,
               "showspec":showSpec,
               "showconf":showCalConf,
               "showatt":showAtt,
               "showcor":showCor,
               "setfilters":setFilters,
               "setwidth":setWidth,
               "calcatt":calcAtt,
               "fitatt":fitAtt,
               "vary":setVary,
               "setcormat":setCorMat,
               "debug":debugToggle,
               "quit":quitCarousel,
               "help":helpCar,
               }


# simple command line loop to allow loading of data and run
# of fitting code.
if __name__ == "__main__":
    logging.basicConfig(level = logging.WARNING)
    nargs = len(sys.argv)
    debug = False
#    if nargs > 1 and sys.argv[1] == "-debug":
#        debug = True
#        print "set debug True"
#    else:
#        debug = False
    # test reading attenuation data for standard materials
    print " "
    print " *** Carousel test program ***"
    print " "
    # set the polynomial order of the fitting variables. Variables are
    # function of line number.
    vary = np.zeros(4,dtype='int')
    vary[3] = -1
    # set an object for the material to which attenuation is to be corrected to; this is null until the user provides one
    corMat = materialAtt("",1.0)
    # command loop
    while True:
        try:
            cmd = raw_input("cmd: ").strip()
        except EOFError as ex:
            sys.exit("EOF")
        words = cmd.split(" ")
        try:
            if words[0] == "help":
                cmd_switch[words[0]](cmd_switch, words)
            else:
                cmd_switch[words[0]](words)
        except SystemExit as ex:
            sys.exit("quit")
        except KeyError as ex:
            if not ( words[0] == "" ):
                print "** command not found: ", words[0]
        except:
            print "** internal error: ", sys.exc_info()[0]
            raise
        
