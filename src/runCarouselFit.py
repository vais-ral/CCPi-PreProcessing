""" Control file to run simple carousel analysis and return
    attenuation correction function. Based on IDL code by Anthony
    Evershed.
"""
import sys
import logging
import timeit
import numpy as np
from carouselUtils import *

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
    print "where n is the carousel data number, max=",carouselCal.samples

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
    """ meant to do repeated fitting"""
    global carouselData, carouselCal, xSpec, debug, vary
    global res
    defMat = "Cu"
    fit = fitData(carouselData, carouselCal, defMat)
    fit.verbose = debug
    if debug:
        np.seterr(over='warn',invalid='warn')
    else:
        np.seterr(over='ignore',invalid='ignore')

    if len(string) == 5:
        print "Fitting variables: ",np.sum(vary)+3
        # The fit function consists of 3 polynomial expressions in the
        # the line number. Initial values for the zero order terms are
        # given here, the higher terms (if any) are set to zero.
        x = np.zeros(3+np.sum(vary))
        offset = vary[0]
        x[offset] = float(string[2])
        offset = offset+1+vary[1]
        x[offset] = float(string[3])
        offset = offset+1+vary[2]
        x[offset] = float(string[4])
        nlines = int(string[1])
#    elif len(string) == 7:
#        print "Fitting 5 variables"
#        x = np.zeros(5)
#        x[0] = float(string[2])
#        x[1] = float(string[3])
#        x[2] = float(string[4])
#        x[3] = float(string[5])
#        x[5] = float(string[6])
#        nlines = int(string[1])
#        print "args=",nlines,x
    else:
        print "wrong number of args: need fitatt nlines x1 x2 x3"
        return
    fit.vary_target = vary[0]
    fit.vary_detector = vary[1]
    fit.vary_filter = vary[2]
    if debug:
        tw,dw,fw = fit.calcWidths(x,nlines)
        print "initial guesses:"
        print "tw=",tw
        print "dw=",dw
        print "fw=",fw
    t0 = timeit.default_timer()
    res,cov,infodict,mesg,ier = fit.dofit(nlines,x)
    print "time=",timeit.default_timer()-t0
    print "dofit returned: "
    print " best fit values = ",res
    print " ier = ",ier
    print " iterations = ",infodict["nfev"]
    if debug:
        print " cov = ",cov
        print " mesg = ",mesg
    # measure error
    samples = carouselCal.samples
    ofile = open('fit.log','w')
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
    for line in range(nlines):
        sumsq = 0.
        for sample in range(samples):
            sumsq += (fit.atten[sample,line] - carouselCal.getAvAtten(line,sample) ) ** 2
        ofile.write(' {0:5d}  {1:12.6f}\n'.format(line,sumsq))
        sumtot += sumsq
        if sumsq>summax:
            summax = sumsq
        if debug:
            print "Line: ",line," ave error:",sumsq
    print "average error: ",sumtot/nlines
    print "max error: ",summax
    rfile.write("average error: {0:12.6f}\nmax error: {1:12.6f}\n".format(sumtot/nlines,summax))
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
    print "Voltage=",  carouselCal.voltage, " angel=", carouselCal.angle

def setVary(strlst):
    """ define polynomial order of parameters to fit """
    global vary
    if len(strlst)==1:
        print "settings:"
        print "target: ",vary[0]
        print "detector: ",vary[1]
        print "filter: ",vary[2]
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
        else:
            print "Not recognised: ",strlst[1]
    else:
        print "syntax: vary [target|detector|filter n]"

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

# Set of commands that are implemented and the corresponding function names.
# Additional commands and functions can be added here.
cmd_switch = { "load":loadAll,
               "showimg":showImg,
               "showspec":showSpec,
               "showconf":showCalConf,
               "showatt":showAtt,
               "setfilters":setFilters,
               "setwidth":setWidth,
               "calcatt":calcAtt,
               "fitatt":fitAtt,
               "vary":setVary,
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
    vary = np.zeros(3,dtype='int')
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
        
