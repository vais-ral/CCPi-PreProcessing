""" Utilities for beam hardening correction method. Classes are provided
to load the various data types that are needed by the algorithm.
classes:
    specData: object to load and store X-ray spectrum for specified case
    carousel: object to load and store description of test carousel
    carouselCalibrationData: object to load and store the calibration data
    fitData: object with methods and data related to the fitting process
    
"""
import sys
import os
import logging
import pdb
try:
    import numpy as np
except ImportError:
    sys.exit("Error: cant find numpy")
try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Error: cant find matplotlib")
try:
    from scipy.optimize import minimize
except ImportError:
    sys.exit("Error: cant find scipy.minimize")


class specData:
    """
    Spectral data for given condition

    Data referes to the initial values of target, voltage and angle.
    This object contains the spectral data for the given
    target material, applied voltage and take off angle.
    Currently this is only available for Tungsten in the range
    1-150KeV from data obtained from Spekcalc at fixed set of values.
    Interpolation may be used for intermeidate points.
    The resolution is as set in the data files.
    Note that only one spectrum is loaded as we do not vary
    the voltage in this fitting procedure.

    Attributes
    ----------
    target : string
        X-ray tube target material, chemical symbol e.g. W
    voltage : float
        applied voltage to the X-ray tube in KeV
    angle : float
        take off angle in degrees for the X-ray beam. This is kept constant

    Methods
    -------
    getE()
        Get array of energies over which the spectrum is defined (KeV)
    getS()
        Get array of spectral intensities corresponding to energies of getE()
    getEnergyRes()
        Return the step size between spectral points, assumed constant.
    getMaxEnergy()
        Return the applied voltage (KeV)
    isValid()
        Flag to indicate if object correctly set up.
    """

    def __init__(self, target, voltage, angle):
        """ set up spectrum from file for given state """
        self.target = target
        self.voltage = voltage
        self.angle = angle
        # initialise these for pylint
        self.en = None
        self.amp = None
        self.valid = False
        self.readData(target, voltage, angle)

    def readData(self, target, voltage, angle):
        """ Read the spectra for this voltage/angle pair from file.
            If file not found, return spec=0 """
        filename = "./spectra/%s/%03d%03d.spc" % (target, voltage, angle*10)
        if not os.path.isfile(filename):
            print "file not found: ", filename, " in specData"
            self.valid = False
        else:
            with open(filename, 'r') as fl:
                self.en, self.amp = np.loadtxt(fl, unpack=True)
            self.valid = True

    def getE(self):
        """ return energy array; should remove """
        return self.en

    def getS(self):
        """ return intensity of spectra at energy points"""
        return self.amp

    def getEnergyRes(self):
        """ return the (assumed constant) energy resolution"""
        return self.en[1]-self.en[0]

    def getMaxEnergy(self):
        """ voltage is max energy"""
        return self.voltage

    def isValid(self):
        """ has object been set up OK"""
        return self.valid

class materialAtt:
    """Attenuation for a material expressed as formula

    For each material of interest will create this object to represent attenuation
    as a function of energy. Material is idenitified by chemical symbol and data is
    read from file produced by xcom.f.
    """

    def __init__(self, formula, density):
        self.name = formula
        # for pylint
        self.mu = None
        self.energy = None
        self.valid = False
        if len(formula) != 0:
            self.readFile(formula, density)

    def readFile(self, formula, density):
        """ load the mu data for this material from a simple ascii file in ./xcom
             mu is mulitplied by density """
        filename = "./xcom/%s.txt" % (formula)
        if os.path.isfile(filename):
            with open(filename, 'r') as fl:
                self.energy, self.mu = np.loadtxt(fl, unpack=True)
                self.mu = self.mu*density
            self.valid = True
        else:
            print "failed to find attenuation file: ", filename
            self.valid = False

    def getE(self):
        """ return array of energy data"""
        return self.energy

    def getMu(self):
        """ attenuation array"""
        return self.mu

    def getEnergyRes(self):
        """ energy resolution, assumed constant"""
        return self.energy[1]-self.energy[0]

    def getMaxEnergy(self):
        """ highest energy value recorded"""
        return self.energy[-1]

    def getMuByE(self,energyVal):
        """ return attenuation for the given energy """
        if energyVal<0:
            print "error: bad value in getMuByE"
            return -1.
        for i in range(len(self.energy)):
            if energyVal >= self.energy[i]:
                if energyVal == self.energy[i]:
                    return(self.mu[i])
                else:
                    frac = (energyVal-self.energy[i-1])/(self.energy[i]-self.energy[i-1])
                    muVal = self.mu[i-1]+(self.mu[i]-self.mu[i-1])*frac
        print "error: failed to match energy in getMuByE"
        return -1.

    def isValid(self):
        """ if object ok"""
        return self.valid

class carousel:
    """ class for data describing the test carousel """

    def __init__(self, defFile):
        self.defFile = defFile
        self.info = None
        self.density = None
        self.numSamples = None
        self.sampWidth = None
        self.materialTypes = None
        self.mask = None
        self.valid = False
        self.__readFile(defFile)

    def __readFile(self, defFile):
        """ read a simple ascii file describing carousel from file
        """
        if os.path.isfile(defFile):
            with open(defFile, 'r') as fl:
                self.info = fl.readline()
                self.numSamples = int(fl.readline())
                self.materialTypes = []
                mline = fl.readline()
                self.materialTypes = mline.split(',')
                self.density = np.fromstring(fl.readline(), dtype = float, sep=',')
                self.sampWidth = np.fromstring(fl.readline(), dtype = float,
                                   sep=',')
                self.filterAtt = {}
                # assume last sample is labeled "Nothing"
                for i in range(self.numSamples - 1):
                    try:
                        self.filterAtt[i] = materialAtt(self.materialTypes[i],self.density[i])
                    except:
                        print "** failed to set carousel attenuation for ",self.materialTypes[i]
                self.mask = np.zeros((self.numSamples),dtype=bool)
            self.valid = True
        else:
            print "failed to find carousel file: ", defFile
            self.valid = False

    def getSamples(self):
        """ elements in the carousel """
        return self.numSamples

    def isValid(self):
        """ is object ok"""
        return self.valid

class carouselCalibrationData:
    """ This class reads the calibration data from a carousel calibration file.

    The calibration file contains information about the X-ray voltage, take-off angle,
    target material and image resolution.
    It also gives the name & possibly the format of the image data file for the carousel.
    """
    def __init__(self, calFile, carouselInfo):
        self.calFile = calFile
        self.samples = carouselInfo.getSamples()-1 # note that last sample "Nothing" has no image
        self.voltage = None
        self.targeMat = None
        self.rows = 0
        self.lines = 0
        self.angle = 0
        self.filterDensity = 0
        self.valid = False
        self.__readCalFile(calFile)
        if self.valid:
            self.__readImageFile(self.imageFile)
            self.__setAverages()
            self.width = 100
            self.__cacheAveSet = False
            self.__cacheAve = np.zeros(shape=(self.samples, self.lines))
            logging.debug('initialised cacheAve')

    def __readCalFile(self, calFile):
        """ read calibration data file, and from that the actual image data"""
        if os.path.isfile(calFile):
            with open(calFile, 'r') as fl:
                try:
                    pos = "start"
                    self.info = self.__readLineStrip(fl)
                    pos = "Voltage"
                    self.voltage = float(self.__readLineStrip(fl))
                    self.angle = float(self.__readLineStrip(fl))
                    self.targetMat = self.__readLineStrip(fl).rstrip()
                    self.targetDensity = float(self.__readLineStrip(fl).rstrip())
                    if self.targetMat != "W":
                        print "Warning: only W (Tungsten) target supported at present: not '", self.targetMat, "'"
                    try:
                        self.targetAtten = materialAtt(self.targetMat,self.targetDensity)
                    except:
                        print "** failed to set target attenuation for ",self.targetMat
                    try:
                        self.spec = specData(self.targetMat,self.voltage,self.angle)
                    except:
                        print "** failed to load spectra for ",self.targetMat,self.voltage,self.angle
                    pos = "rows"
                    self.rows = int(self.__readLineStrip(fl))
                    self.lines = int(self.__readLineStrip(fl))
                    self.imageFile = self.__readLineStrip(fl)
                    self.imageFileFormat = self.__readLineStrip(fl)
                    pos = "filterCount"
                    self.filterCount = int(self.__readLineStrip(fl))
                    self.filterMaterial = {}
                    self.filterWidth = {}
                    self.filterDensity = {}
                    self.filterAtten = {}
                    pos = "filtermatcnt"
                    for i in range(self.filterCount):
                        self.filterMaterial[i] = self.__readLineStrip(fl)
                        self.filterWidth[i] = float(self.__readLineStrip(fl))
                        self.filterDensity[i] = float(self.__readLineStrip(fl))
                        try:
                            self.filterAtten[i] = materialAtt(self.filterMaterial[i],self.filterDensity[i])
                        except:
                            print "** failed to set attenuation for ",self.filterMaterial[i]

                    pos = "detector"
                    self.detectorMaterial = self.__readLineStrip(fl)
                    self.detectorWidth = float(self.__readLineStrip(fl))
                    self.detectorDensity = float(self.__readLineStrip(fl))
                    try:
                        self.detectorAtten = materialAtt(self.detectorMaterial,self.detectorDensity)
                    except:
                        print "** failed to set detector attenuation for ",self.detectorMaterial
                    self.valid = True
                except (ValueError, IOError):
                    print "Read Calibration file failed beyond: ", pos
                    self.valid = False
        else:
            print "failed to find calibration file: ", calFile
            self.valid = False

    def __readLineStrip(self, fl):
        strng = fl.readline()
        if strng.find("#") == -1:
            return strng.rstrip()
        else:
            return strng[:strng.find("#")].rstrip()

    def __readImageFile(self, imageFile):
        if os.path.isfile(imageFile):
            with open(imageFile, 'rb') as fl:
                self.image = np.fromfile(fl, dtype = self.imageFileFormat,
                                       count = self.rows*self.lines*self.samples).reshape(self.samples, self.lines, self.rows)
        else:
            print "Image file not found!: ", imageFile

    def printImageStats(self, carInf):
        """ print out some data for each frame in the set of images"""
        for i in range(self.samples):
            nancount = np.count_nonzero(np.isnan(self.image[i,:,:]))
            if nancount>0:
                print "*** img ", i, " contains: ", nancount, " NaNs (", nancount*100./(self.rows*self.lines), "%)"
                maskedimage = np.ma.array(self.image[i,:,:],mask = np.isnan(self.image[i,:,:]))
                print "    average(masked)= ", np.ma.average(maskedimage), "  max= ", np.ma.max(maskedimage)
            else:
                print "img ", i, " ", carInf.filterMaterial[i]," ", carInf.filterWidth[i], " average= ", \
                      np.average(self.image[i,:,:])

    def __setAverages(self):
        """ pre compute mean centre of each image row
        """
        self.__centre = np.zeros(self.samples*self.lines).reshape(self.samples,
                               self.lines)
        for samp in range(self.samples):
            for li in range(self.lines):
                self.__centre[samp, li] = self.__getCentrePos(li, samp)

    def __getCentrePos(self, line, sample):
        """ do a weighted average of position by signal to find the centre of
            the signal, in pixels"""
        vals = self.image[sample,line,:]
        meanpt = np.sum(np.arange(len(vals))*vals)/np.sum(vals)
        return meanpt

    def getCentrePos(self, line, sample):
        """ get centre point"""
        return int(self.__centre[sample, line])

    def getLines(self):
        """ return the number of lines of the data"""
        return self.lines

    def getRows(self):
        """ return the number of rows of the data"""
        return self.rows

    def getImage(self, imgNum):
        """ access a given image by number"""
        if imgNum >= 0 and imgNum <self.samples:
            return self.image[imgNum,:,:]
        else:
            return 0

    def getAvAtten(self, line, sample):
        """ average over a named range of rows """
        if self.__cacheAveSet:
            return self.__cacheAve[sample, line]
        logging.debug('calc cache values of Ave')
        self.__cacheAveSet =True
        for s in range(self.samples):
            for l in range(self.lines):
                logging.debug('s= %d l=%d c=%f', s, l, self.__centre[s, l])
                rowStart = int(self.__centre[s, l]-self.width)
                rowEnd = int(rowStart+2*self.width)
                self.__cacheAve[s, l] = np.average(self.image[s, l, rowStart:rowEnd])
        return self.__cacheAve[sample, line]

    def setWidthAve(self, width):
        """ set the (half) width to be used when calculating the average
            attenuation along a row
            This is measured either size of the centre point. All points should
            lie on the sample image.
        """
        self.width = width
        self.__cacheAveSet = False

    def isValid(self):
        """ is object ok"""
        return self.valid

    def plotCalData(self, showplt, markWidth):
        """ plot calibration images in one window using matplotlib. Only
            return when window closed. Mark image with width of "averaging",
            markWidth """
        nsample = self.samples
        if nsample<5:
            pltgrid=2
        elif nsample<10:
            pltgrid=3
        else:
            pltgrid=4
            nsample = min(nsample, 15)
        for isam in range(nsample):
            # deep copy image so can modify
            z = np.copy(self.getImage(isam))
            zzmax = 0
            # mark width of each line used for average
            for li in range(self.lines):
                cen = self.getCentrePos(li, isam)
                if markWidth>0:
                    z[li, cen-markWidth] = 0
                    z[li, cen+markWidth] = 0
            maskedimage = np.ma.array(z, mask = np.isnan(z))
            zmax = np.max(maskedimage)
            zzmax = max(zmax, zzmax)
            plt.subplot(pltgrid, pltgrid, isam+1, axisbg='y')
            plt.imshow(maskedimage, cmap='RdBu', vmin = 0, vmax = zzmax,
                       aspect='auto')
            plt.draw()
        if showplt:
            plt.show()

    def plotCalProf(self):
        """ will be used to plot profile along a row; not yet complete"""
        nsample = self.samples
        if self.samples < 5:
            pltgrid = 2
        elif self.samples < 10:
            pltgrid = 3
        else:
            pltgrid = 4
            nsample = max(self.samples)
        for isam in range(nsample):
            plt.subplot(pltgrid, pltgrid, isam+1, axisbg = 'y')
            


class fitData:
    """ This object contain fit related data and functions
    """

    def __init__(self, carInfo, carCal, defMat):
        if not ( carInfo.isValid and carCal.isValid ):
            self.isValid = False
            return
        self.carInfo = carInfo
        self.carCal = carCal
        #
        # set default fitting parameters; values are order of polynomials
        # in line number. Hence "0" is a constant, not dependent on line number
        # default is for 3rd polys for detector and target width, with global
        # width for Cu filters.
        self.vary_target=1  #3
        self.vary_detector=1  #3
        self.vary_filter=0
        self.vary_energy=-1
        # since there may be several filters, define which one should vary
        self.vary_filter_name="Cu"
        self.verbose = False
        # define the fitting parameters and their global/local status
        self.varFilter = -1
        # check if we have filter of material defMat
        for i in carCal.filterMaterial:
            #print carCal.filterMaterial[i]
            if defMat == carCal.filterMaterial[i]:
                self.varFilter = i
                self.defFilterMat = defMat
                # this value not used
                #self.defFilter = 'global'
        # if defMat not found, use first material in list, if any
        if self.varFilter==-1:
            print "default filter: ", defMat," not found"
            if len(carCal.filterMaterial)>0:
                self.defFilterMat=carCal.filterMaterial[0]
                self.varFilter = 0
                print "using material = ",self.defFilterMat," for fitting"
            else:
                print "no filter material present, hence cannot fit this"
            # self.isValid = False
        self.atten = np.zeros([self.carCal.lines,self.carInfo.getSamples()])
        #if len([m for m in carCal.filterMaterial if defMat in m]) > 0 :
        #    self.defMat=defMat
        #    self.defFilter='global'
        #else:
        #    self.defMat=""
        #    self.defFilter="none"
        #    self.isValid=False
        #    return
        #self.isValid=True

    def calc_atten(self,widths,atten,energy,steps,detec_wid,detec_atten,spec):
        """ find the attenuation from set of filters

            widths - array of filter widths
            atten  - 2D array, energy steps by filter count
            steps  - scalar count of energy steps
            detec_wid - scalar width of detector
            detec_atten - array of detector attenuation values
            spec   - array of spectral intensity
        """
        nfilter = len(widths)
        at = np.ones(len(spec))
        for i in range(nfilter):
            at = at*widths(i)*atten(i)
        at = np.exp(-at)*spec*(1.-np.exp(-detec_wid*detec_atten))
        atsum = np.sum(at)
        return atsum
            
    def calcWidths(self,x0,nlines,xe):
        """ simple function to return the 3 widths for the target,
            the detector and the filter from the set of fitting
            variables. Each is assumed to be a polnomial of some
            order in the line number. Global varables are zero
            order polynomials """
        lines=np.array(range(nlines),dtype="double")
        nt=self.vary_target+1
        nd=self.vary_detector+1
        nf=self.vary_filter+1
        ne=self.vary_energy+1
        if len(x0)<nt+nd+nf+ne:
            print "** calcWidthd called with too few values in x0"
            return
        # Polynomial expressions: highest order term is first in the array.
        twidth=np.polyval(x0[:nt],lines)
        #twidth=np.polyval(x0[nt-1:0:-1],lines)
        dwidth=np.polyval(x0[nt:nt+nd],lines)
        #dwidth=np.polyval(x0[nt+nd-1:nt:-1],lines)
        dwidth = np.exp( dwidth ) # force >0 by working in log space
        fwidth=np.polyval(x0[nt+nd:nt+nd+nf],lines)
        if ne>0:
            ecoeffs = xe + xe*xe*np.polyval(x0[nt+nd+nf:],xe)
        else:
            ecoeffs = xe
        return twidth,dwidth,fwidth,ecoeffs

    def dofit(self,nlines,xin):
        """ perform fit """
        try:
            from scipy.optimize import minimize
            from scipy.optimize import leastsq
        except:
            print "** cannot find scipy.minmize - check python has scipy"
            return
        if self.verbose:
            pdb.set_trace()
        x = xin
        self.nlines = nlines
        #self.objFun(x)
        #res = minimize(self.objFun, x, method='Nelder-Mead')
        res = leastsq(self.objFunSq, x, full_output = True)
        #for field, val in  leastsq(self.objFunSq, x, full_output = True).items():
        #    print "field=",field," value=",value
        print "atten=",self.atten[0,:]
        expt = np.zeros(self.carCal.samples+1)
        for i in range(self.carCal.samples):
             expt[i] = self.carCal.getAvAtten(0,i)
        print "expt=",expt
        #plt.plot(self.atten[0,:])
        #plt.plot(expt)
        #plt.draw()
        #plt.show()
        return res

    def objFunSq(self,x):
        """ The function to minimize; returns the squared error for every point on each
            selected line  """
        # Get the 3 widths: target(e.g. W), detector(e.g. CsI), global filter(e.g. Cu)
        # target and detector widths depend on line number, filter is a global value
        # for flexiblity all 3 are dimesioned by nlines
        xe = self.carCal.spec.getE()
        tw,dw,fw,ec = self.calcWidths(x,self.nlines,xe)
        nsamples = self.carInfo.numSamples - 1 # ignore null sample
        ans = np.zeros(nsamples*self.nlines)
        tarAtt = self.carCal.targetAtten
        se = self.carCal.spec.getS()
        #
        for line in range(self.nlines):
            # compute filter attenuation and hence i0 as signal level with no sample for this line
            attDet = dw[line]*self.carCal.detectorAtten.getMu()[:len(xe)]
            attSum = np.zeros(len(xe))
            for filt in range(self.carCal.filterCount):
                if filt == self.varFilter:
                    fwid = fw[line]
                else:
                    fwid = self.carCal.filterWidth[filt]
                attSum = attSum + fwid*self.carCal.filterAtten[filt].getMu()[:len(xe)]
            # this is the key integral done as a simple sum. Can ignore width of each value
            # as constant energy steps, so cancels in I/I0
            at_se = se*np.exp(-attSum-tarAtt.getMu()[:len(xe)]*tw[line])*ec*(1-np.exp(-attDet))
            # remove nan's - why are nan's present? exp overflow gives inf, multiply by 0 gives nan
            # in most cases nans are OK to ignore.
            at_se_finite = at_se[np.logical_not(np.isnan(at_se))]
            i0 = np.sum(at_se_finite)
            #
            for sample in range(nsamples):
                # skip masked samples
                if self.carInfo.mask[sample]:
                    continue
                widSam = self.carInfo.sampWidth[sample]
                #if sample < nsamples-1:
                attSam = widSam*self.carInfo.filterAtt[sample].getMu()[:len(xe)]
                #else:
                #    attSam = np.ones(len(xe))
                #attDet = dw[line]*self.carCal.detectorAtten.getMu()[:len(xe)] #csiAtten.getMu()[0:me-1]
                at_se_sample = at_se * np.exp(-attSam)
                #
                # remove nan's - see above
                at_se_sample_finite = at_se_sample[np.logical_not(np.isnan(at_se_sample))]
                i_sample = np.sum(at_se_sample_finite)
                if i0==0. and self.verbose:
                    print "warn: i0 zero at ",line
                    i0=1.
                if i_sample == 0 and self.verbose:
                    print "i_sample=0"
                if i_sample < 0. and self.verbose:
                    print "i_sample<0",at_se_sample[:8]
                #sumSq = sumSq + ( (i_sample/i0) - self.carCal.getAvAtten(line,sample) ) ** 2
                ans[line*nsamples+sample] = ( np.log(i0/i_sample) - self.carCal.getAvAtten(line,sample) ) ** 2
                self.atten[line,sample] = np.log(i0/i_sample)
        if self.verbose:
            print "tw,dw,fw,sumSq: ",tw[0],dw[0],fw[0],np.sum(ans)
        #
        # return vector of squared errors: length=samples*lines
        return ans

    def linesPolyFit(self,soln,corMat,corEn,npoints,attrange):
        """ Function to calculate the attenuation over "npoints" for attenuation up to "attrange" for each line.
            Uses the fitted parameters for attenuation in "soln". Having calculated apparent attenuation for the correction
            material "corMat", map the observed attenuation to the true attenuation at the corEn energy. Then fit a
            polynomial to the data and save the coefficients.
            """
        # Get the 3 widths: target(e.g. W), detector(e.g. CsI), global filter(e.g. Cu)
        # target and detector widths depend on line number, filter is a global value
        # for flexiblity all 3 are dimesioned by nlines. Also the energy parameter, if fitted.
        xe = self.carCal.spec.getE()
        tw,dw,fw,ec = self.calcWidths(soln,self.nlines,xe)
        #nsamples = self.carInfo.numSamples
        # allocate space to store all calculated points and polynomials fitted to them
        attout = np.zeros(shape=(self.nlines,npoints+1))
        odpoly = 8
        fitdata = np.zeros(odpoly*self.nlines)
        #
        # find the actual attenuation of the correction material at the correction energy
        corrAtt = corMat.getMuByE(corEn)
        
        tarAtt = self.carCal.targetAtten
        se = self.carCal.spec.getS()

        # generate points to evaluate attenuation at.
        mulist = np.arange(npoints+1,dtype='float')*attrange/(npoints*corMat.getMuByE(corEn))
        # attin is the observed attenuation; for each line want to find corresponding attout
        # the "true" attenuation at energy corEn
        attin = mulist*corMat.getMuByE(corEn)
        polyfit = np.zeros(shape=(self.nlines,odpoly+2))
        #
        # for each line generate npoints values of attenuation from fit data
        #
        for line in range(self.nlines):
            # compute filter attenuation and hence i0 as signal level with no sample for this line
            attDet = dw[line]*self.carCal.detectorAtten.getMu()[:len(xe)]
            attSum = np.zeros(len(xe))
            for filt in range(self.carCal.filterCount):
                if filt == self.varFilter:
                    fwid = fw[line]
                else:
                    fwid = self.carCal.filterWidth[filt]
                attSum = attSum + fwid*self.carCal.filterAtten[filt].getMu()[:len(xe)]
            # this is the key integral done as a simple sum. Can ignore width of each value
            # as constant energy steps, so cancels in I/I0
            at_se = se*np.exp(-attSum-tarAtt.getMu()[:len(xe)]*tw[line])*ec*(1-np.exp(-attDet))
            # remove nan's - why are nan's present? exp overflow gives inf, multiply by 0 gives nan
            # in most cases nans are OK to ignore.
            at_se_finite = at_se[np.logical_not(np.isnan(at_se))]
            i0 = np.sum(at_se_finite)
            #
            count = 0
            # loop over required attenuation values
            for muwid in mulist:
                attSam = muwid*corMat.getMu()[:len(xe)]
                at_se_sample = at_se * np.exp(-attSam)
                #
                # remove nan's - see above
                #  should not be needed here if fit is OK
                at_se_sample_finite = at_se_sample[np.logical_not(np.isnan(at_se_sample))]
                i_sample = np.sum(at_se_sample_finite)
                if i0==0. and self.verbose:
                    print "warn: i0 zero at ",line
                    i0=1.
                if i_sample == 0 and self.verbose:
                    print "i_sample=0"
                if i_sample < 0. and self.verbose:
                    print "i_sample<0",at_se_sample[:8]
                #sumSq = sumSq + ( (i_sample/i0) - self.carCal.getAvAtten(line,sample) ) ** 2
                attout[line,count] = np.log(i0/i_sample)
                count = count+1
            #
            polyfit[line,0:odpoly+1] = np.polyfit(attout[line,1:],attin[1:]/attout[line,1:],odpoly)
        #
        # following carousel.pro, attout is the apparent attenuation or the x-axis of our correction
        # graph. the y-axis should be the actual attenuation at monochromatic energy corEn for the
        # correction material. Since the density of the correction material is unknown, which is the
        # main point of this code, width and density are not used here. For fitting of a polynomial
        # through (0,0) can divide y values by x values, ignoring origin (first point in this case).
        #
        
        return attout,attin,polyfit
