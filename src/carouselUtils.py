""" Utilities for beam hardening correction method. Classes are provided
to load the various data types that are needed by the algorithm.
classes:

    specData  object to load and store X-ray spectrum for specified case
    carousel  object to load and store description of test carousel
    carouselCalibrationData  object to load and store the calibration data
    fitData  object with methods and data related to the fitting process

"""
# just for reading 16 bit images and conversion to log(I0/I)
from __future__ import division
from __future__ import print_function

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


class specData(object):
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
   
    target : string
        X-ray tube target material, chemical symbol e.g. W
    voltage : float
        applied voltage to the X-ray tube in KeV
    angle : float
        take off angle in degrees for the X-ray beam. This is kept constant

    Methods
    
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
            print("File not found: ", filename, " in spectra")
            print("Only fitting with function spectra possible")
            # just set energy array 0 to voltage, in 0.5KeV steps
            estep=0.5
            self.en = np.arange(voltage/estep+estep*0.5)*estep
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

class materialAtt(object):
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
            print("failed to find attenuation file: ", filename)
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
            print("error: bad value in getMuByE")
            return -1.
        for i in range(len(self.energy)):
            if energyVal <= self.energy[i]:
                if energyVal == self.energy[i]:
                    return self.mu[i]
                else:
                    frac = (energyVal-self.energy[i-1])/(self.energy[i]-self.energy[i-1])
                    muVal = self.mu[i-1]+(self.mu[i]-self.mu[i-1])*frac
                    return muVal
        print("error: failed to match energy in getMuByE")
        return -1.

    def isValid(self):
        """ if object ok"""
        return self.valid

class carousel(object):
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
        self.filterAtt = {}
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
                        print("** failed to set carousel attenuation for ",self.materialTypes[i])
                self.mask = np.zeros((self.numSamples),dtype=bool)
            self.valid = True
        else:
            print("failed to find carousel file: ", defFile)
            self.valid = False

    def getSamples(self):
        """ elements in the carousel """
        return self.numSamples

    def isValid(self):
        """ is object ok"""
        return self.valid

class carouselCalibrationData(object):
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

        self.filterMaterial = {}
        self.filterWidth = {}
        self.filterDensity = {}
        self.filterAtten = {}

        self.detectorMaterial = ''
        self.detectorWidth = 0.
        self.detectorDensity = 0.
        self.detectorAtten = {}

        self.targetMat = ''
        self.targetDensity = 0.
        self.targetAtten = {}

        self.info = ''
        self.imageFileFormat = ''
        self.imageFile = ''
        self.image = {}
        self.spec = {}
        self.filterCount = 0
        self.__centre = None

        self.__readCalFile(calFile)
        if self.valid:
            try:
                self.whiteLevel = 0
                self.__readImageFile(self.imageFile)
                self.__setAverages()
                self.width = 100
                self.__cacheAveSet = False
                self.__cacheAve = np.zeros(shape=(self.samples, self.lines))
                logging.debug('initialised cacheAve')
            except:
                self.valid = False
                logging.debug('reading data failed')

    def __readCalFile(self, calFile):
        """ read calibration data file, and from that the actual image data"""
        if os.path.isfile(calFile):
            with open(calFile, 'r') as fl:
                try:
                    self.info = self.__readLineStrip(fl)
                    self.voltage = float(self.__readLineStrip(fl))
                    self.angle = float(self.__readLineStrip(fl))
                    self.targetMat = self.__readLineStrip(fl).rstrip()
                    self.targetDensity = float(self.__readLineStrip(fl).rstrip())
                    if self.targetMat != "W":
                        print("Warning: only W (Tungsten) target supported at present: not '", self.targetMat, "'")
                    try:
                        self.targetAtten = materialAtt(self.targetMat,self.targetDensity)
                    except:
                        print("** failed to set target attenuation for ",self.targetMat)
                    try:
                        self.spec = specData(self.targetMat,self.voltage,self.angle)
                    except:
                        print("** failed to load spectra for ",self.targetMat,self.voltage,self.angle)
                    self.rows = int(self.__readLineStrip(fl))
                    self.lines = int(self.__readLineStrip(fl))
                    self.imageFile = self.__readLineStrip(fl)
                    self.imageFileFormat = self.__readLineStrip(fl)
                    self.filterCount = int(self.__readLineStrip(fl))
                    for i in range(self.filterCount):
                        self.filterMaterial[i] = self.__readLineStrip(fl)
                        self.filterWidth[i] = float(self.__readLineStrip(fl))
                        self.filterDensity[i] = float(self.__readLineStrip(fl))
                        try:
                            self.filterAtten[i] = materialAtt(self.filterMaterial[i],self.filterDensity[i])
                        except:
                            print("** failed to set attenuation for ",self.filterMaterial[i])

                    self.detectorMaterial = self.__readLineStrip(fl)
                    self.detectorWidth = float(self.__readLineStrip(fl))
                    self.detectorDensity = float(self.__readLineStrip(fl))
                    try:
                        self.detectorAtten = materialAtt(self.detectorMaterial,self.detectorDensity)
                    except:
                        print("** failed to set detector attenuation for ",self.detectorMaterial)
                    self.valid = True
                except (ValueError, IOError):
                    print("Read Calibration file failed")
                    self.valid = False
        else:
            print("failed to find calibration file: ", calFile)
            self.valid = False

    def __readLineStrip(self, fl):
        """ read line: strip newline and anything beyond #, if present
        """
        strng = fl.readline()
        if strng.find("#") == -1:
            return strng.rstrip()
        else:
            return strng[:strng.find("#")].rstrip()

    def __readImageFile(self, imageFile):
        """ read the image data based on specificed format """
        if os.path.isfile(imageFile):
            if self.imageFileFormat == "uint16":
                # if raw int16 data, assume first image is flat field and normalise rst of image by this and do log
                nimages = self.samples+1
                with open(imageFile, 'rb') as fl:
                    tmpimage = np.fromfile(fl, dtype = self.imageFileFormat,
                                           count = self.rows*self.lines*nimages)
                    tmpimage = tmpimage.reshape(nimages, self.lines, self.rows)
                self.image = np.zeros(self.rows*self.lines*self.samples,
                                      dtype=float).reshape(self.samples, self.lines, self.rows)
                # note - imported division from _future_ to avoid int div

                # assume that the first image is the white level, I0, a constant
                # over the image. To impose this assumption we take the average
                # value over the first (flat field, shading corrected) image.
                # This is only set for uint16; it is not known for float32
                whiteLev = np.average(tmpimage[0,:,:])
                self.whiteLevel = whiteLev
                for i in range(nimages-1):
                    # catch zero data:
                    if np.min(tmpimage[i+1,:,:])==0:
                        imt = tmpimage[i+1,:,:]
                        imt[imt==0] = np.mean(imt)
                        logging.warning('zero data replaced in image %d',i+1)
                    self.image[i,:,:] = np.log( whiteLev / tmpimage[i+1,:,:] )

            elif self.imageFileFormat == "float32":
                # assume float data already transformed by I0 and log
                with open(imageFile, 'rb') as fl:
                    try:
                        self.image = np.fromfile(fl, dtype = self.imageFileFormat,
                                 count = self.rows*self.lines*self.samples).reshape(self.samples, self.lines, self.rows)
                    except:
                        print("Failed in reading/converting image file: check data")
            else:
                print("** error: image format name ",self.imageFileFormat," not recognised")

        else:
            print("Image file not found!: ", imageFile)

    def printImageStats(self, carInf):
        """ print out some data for each frame in the set of images"""
        for i in range(self.samples):
            nancount = np.count_nonzero(np.isnan(self.image[i,:,:]))
            if nancount>0:
                print("*** img ", i, " contains: ", nancount, " NaNs (", nancount*100./(self.rows*self.lines), "%)")
                maskedimage = np.ma.array(self.image[i,:,:],mask = np.isnan(self.image[i,:,:]))
                print("    average(masked)= ", np.ma.average(maskedimage), "  max= ", np.ma.max(maskedimage))
            else:
                print("img ", i, " ", carInf.filterMaterial[i]," ", carInf.filterWidth[i], " average= ",
                      np.average(self.image[i,:,:]))

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
        if meanpt>0.75*np.size(vals) or meanpt<0.25*np.size(vals):
            meanpt = np.size(vals)*0.5
            logging.info('reset getCentrePos to %f for sample %d',meanpt,sample)
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
                    z[li, int(round(cen-markWidth))] = 0
                    z[li, int(round(cen+markWidth))] = 0
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



class fitData(object):
    """ This object contains fit related data and functions
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
        # width for Cu filters. Also Gaussian for spectra, if no pre-defined values.
        #
        self.vary_target = 1
        self.vary_detector = 1
        self.vary_filter = 0
        self.vary_energy = -1
        self.vary_epk = -1
        self.vary_ewidlow = -1
        self.vary_ewidhigh = -1
        # make space for default values to use if vary==-1
        self.defaults = np.zeros(4)
        #
        self.nlines = 0
        # since there may be several filters, define which one should vary
        self.vary_filter_name="Cu"
        self.verbose = False
        # define the fitting parameters and their global/local status
        self.varFilter = -1
        # check if we have filter of material defMat
        for i in carCal.filterMaterial:
            if defMat == carCal.filterMaterial[i]:
                self.varFilter = i
                self.defFilterMat = defMat
        # if defMat not found, use first material in list, if any
        if self.varFilter==-1:
            print("default filter: ", defMat," not found")
            if len(carCal.filterMaterial)>0:
                self.defFilterMat=carCal.filterMaterial[0]
                self.varFilter = 0
                print("using material = ",self.defFilterMat," for fitting")
            else:
                print("no filter material present, hence cannot fit this")
        self.atten = np.zeros([self.carCal.lines,self.carInfo.getSamples()])
        self.objFnCalls = 0
        self.nlines = 0
        self.lineStep = 1
        self.linestep = 1
        self.bounds = False
        self.boundsValues = {}
        self.solver = "old"

    def calcWidths(self,x0,nlines,xe):
        """ Function to return the 3 widths for the target,
            the detector and the filter from the set of fitting
            variables. Each is assumed to be a polynomial of some
            order in the line number. Global varables are zero
            order polynomials.
            Also returns the energy array which can be a polynomial:
            E+aE**2+... ; this should be constrained >=0, not done at
            present."""
        lines=np.array(range(nlines),dtype="double")
        nt = self.vary_target+1
        nd = self.vary_detector+1
        nf = self.vary_filter+1
        ne = self.vary_energy+1
        ns = self.vary_epk+self.vary_ewidlow+self.vary_ewidhigh+3
        if len(x0)<nt+nd+nf+ne+ns:
            print("** calcWidthd called with too few values in x0")
            sys.exit(1)
        # Polynomial expressions: highest order term is first in the array.
        if nt>0:
            twidth = np.polyval(x0[:nt],lines)
        else:
            twidth = np.polyval(self.defaults[0:1],lines)
        #twidth=np.polyval(x0[nt-1:0:-1],lines)
        if nd>0:
            dwidth = np.polyval(x0[nt:nt+nd],lines)
        else:
            dwidth = np.polyval(self.defaults[1:2],lines)
        #dwidth=np.polyval(x0[nt+nd-1:nt:-1],lines)
        dwidth = np.exp( dwidth ) # force >0 by working in log space
        if nf>0:
            fwidth = np.polyval(x0[nt+nd:nt+nd+nf],lines)
        else:
            fwidth = np.polyval(self.defaults[2:3],lines)
        if ne>0:
            # This term should be constrained as >=0 for all xe but is not at present.
            # -Ve values will give errors in output stage.
            ecoeffs = xe + xe*xe*(np.polyval(x0[nt+nd+nf:nt+nd+nf+ne],xe))
        else:
            ecoeffs = xe
        if ns>0:
            i0 = nt+nd+nf+ne
            indarr = xe>x0[i0:i0+1]
            spectra = xe-x0[i0:i0+1]
            spectra[indarr] = spectra[indarr]*x0[i0+1:i0+2]
            spectra[~indarr] = spectra[~indarr]*x0[i0+2:i0+3]
            spectra = np.exp(-spectra**2)
            # mask out lowest 10% of spectra, as in spekCalc
            masklen = int(len(spectra)*0.1)
            spectra[0:masklen] = 0.
        else:
            spectra = 0.
        return twidth,dwidth,fwidth,ecoeffs,spectra

    def dofit(self,nlines,lstep,xin):
        """ perform fit """
        got=0
        try:
            # from scipy.optimize import minimize
            from scipy.optimize import leastsq
            got=1
            from scipy.optimize import least_squares
        except:
            # if found old lib try and continue
            if got==0:
                print("** cannot find scipy leastsq or least_squares - check python has scipy")
                return

        if self.verbose:
            pdb.set_trace()
        x = xin
        self.nlines = nlines
        self.lineStep = lstep

        # use either old or new solver interface from scipy for least squares
        if self.solver=="old":
            res = leastsq(self.objFunSq, x, full_output = True)
        else:
            if self.bounds:
                resobj = least_squares(self.objFunSq, x, verbose = 1, bounds=self.boundsValues)
            else:
                resobj = least_squares(self.objFunSq, x, verbose = 1, method='lm')
            infodict = {"nfev":resobj.nfev}
            cov = [0]
            res = (resobj.x,cov,infodict,resobj.message,resobj.status)

        print("Line 0 atten=",self.atten[0,:])
        expt = np.zeros(self.carCal.samples+1)
        for i in range(self.carCal.samples):
            expt[i] = self.carCal.getAvAtten(0,i)
        print("Line 0 expt=",expt)
        # Do final calulation on all lines:
        self.lineStep = 1
        if self.solver=="old":
            self.objFunSq(res[0])
        else:
            self.objFunSq(resobj.x)
        #
        return res

    def objFunSq(self,x):
        """ The function to minimize; returns the squared error for every point on each
            selected line  """
        # Get the 3 widths: target(e.g. W), detector(e.g. CsI), global filter(e.g. Cu)
        # target and detector widths depend on line number, filter is a global value
        # for flexiblity all 3 are dimesioned by nlines
        #
        # mask out low en spectral points which may get undue weight if -ve filter widths occur
        minpt = int(0.1*len(self.carCal.spec.getE()))
        #
        xe = self.carCal.spec.getE()
        tw,dw,fw,ec,spectra = self.calcWidths(x,self.nlines,xe)
        nsamples = self.carInfo.numSamples - 1 # ignore null sample
        ans = np.zeros(nsamples*self.nlines)
        tarAtt = self.carCal.targetAtten
        if isinstance(spectra,np.ndarray):
            se = spectra
        else:
            se = self.carCal.spec.getS()
        #
        for line in range(0,self.nlines,self.lineStep):
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
            # drop low energy terms below 10%
            at_se_t = at_se[minpt:]
            # remove nan's - why are nan's present? exp overflow gives inf, multiply by 0 gives nan
            # in most cases nans are OK to ignore.
            at_se_finite = at_se_t[np.logical_not(np.isnan(at_se_t))]
            i0 = np.sum(at_se_finite)
            #
            for sample in range(nsamples):
                # skip masked samples
                if self.carInfo.mask[sample]:
                    continue
                widSam = self.carInfo.sampWidth[sample]
                #
                attSam = widSam*self.carInfo.filterAtt[sample].getMu()[:len(xe)]
                #
                # drop low energy terms < 10%
                at_se_sample = (at_se * np.exp(-attSam))[minpt:]
                #
                # remove nan's - see above
                at_se_sample_finite = at_se_sample[np.logical_not(np.isnan(at_se_sample))]
                i_sample = np.sum(at_se_sample_finite)
                if i0==0. and self.verbose:
                    print("warn: i0 zero at ",line)
                    i0=1.
                if i_sample == 0 and self.verbose:
                    print("i_sample=0")
                if i_sample < 0. and self.verbose:
                    print("i_sample<0",at_se_sample[:8])
                #sumSq = sumSq + ( (i_sample/i0) - self.carCal.getAvAtten(line,sample) ) ** 2
                ans[line*nsamples+sample] = ( np.log(i0/i_sample) - self.carCal.getAvAtten(line,sample) ) ** 2
                self.atten[line,sample] = np.log(i0/i_sample)
        if self.verbose:
            print("tw,dw,fw,sumSq: ",tw[0],dw[0],fw[0],np.sum(ans))
            # for debugging we print the normalised detector response for the current parameters
            plotFreq=10
            if self.objFnCalls%plotFreq == 0:
                plt.figure('NormRespone')
                plt.xlabel('energy')
                plt.ylabel('Normalised response')
                # map nans to zero for plotting
                at_se[np.isnan(at_se)] = 0.
                plt.plot(xe,at_se/i0)
                plt.draw()
                plt.show(block=False)
        self.objFnCalls=self.objFnCalls+1

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
        tw,dw,fw,ec,spectra = self.calcWidths(soln,self.nlines,xe)
        #nsamples = self.carInfo.numSamples
        # allocate space to store all calculated points and polynomials fitted to them
        attout = np.zeros(shape=(self.nlines,npoints+1))
        # set order of polynomial fits to use
        # for xtek a 4th order polynomial can be used by the reconstruction - note
        # that as constant term is forced to zero, 3 gives 4th order
        odpoly = 8
        xtekodpoly = 3
        # determine if the solution varies with line number; if not only one fit required
        vary_line = (self.vary_target>0 or self.vary_detector>0 or self.vary_filter>0)
        #
        # find the actual attenuation of the correction material at the correction energy
        corrAtt = corMat.getMuByE(corEn)

        tarAtt = self.carCal.targetAtten
        #
        if isinstance(spectra,np.ndarray):
            se = spectra
        else:
            se = self.carCal.spec.getS()
        #

        # generate points to evaluate attenuation at.
        mulist = np.arange(npoints+1,dtype='float')*attrange/(npoints*corrAtt)
        # attin is the observed attenuation; for each line want to find corresponding attout
        # the "true" attenuation at energy corEn
        attin = mulist*corrAtt
        #
        # for each line generate npoints values of attenuation from fit data
        #
        if vary_line:
            nlines = self.nlines
        else:
            nlines = 1
        #
        polyfit = np.zeros(shape=(nlines,odpoly+2))
        xpolyfit = np.zeros(shape=(nlines,xtekodpoly+2))
        #
        for line in range(nlines):
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
                    print("warn: i0 zero at ",line)
                    i0=1.
                if i_sample == 0 and self.verbose:
                    print("i_sample=0")
                if i_sample < 0. and self.verbose:
                    print("i_sample<0",at_se_sample[:8])
                #sumSq = sumSq + ( (i_sample/i0) - self.carCal.getAvAtten(line,sample) ) ** 2
                attout[line,count] = np.log(i0/i_sample)
                count = count+1
            #
            try:
                polyfit[line,0:odpoly+1] = np.polyfit(attout[line,1:],attin[1:]/attout[line,1:],odpoly)
                xpolyfit[line,0:xtekodpoly+1] = np.polyfit(attout[line,1:],attin[1:]/attout[line,1:],xtekodpoly)
            except:
                print("*** Polynomial fit of result failed")
        #
        # following carousel.pro, attout is the apparent attenuation or the x-axis of our correction
        # graph. the y-axis should be the actual attenuation at monochromatic energy corEn for the
        # correction material. Since the density of the correction material is unknown, which is the
        # main point of this code, width and density are not used here. For fitting of a polynomial
        # through (0,0) can divide y values by x values, ignoring origin (first point in this case).
        #

        return attout,attin,polyfit,xpolyfit
