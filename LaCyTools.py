#! /usr/bin/env python
from datetime import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from scipy.interpolate import interp1d
import scipy.optimize
from scipy.optimize import curve_fit
from Tkinter import *
import base64
import collections
import glob
import itertools
import linecache
import math
import matplotlib.pyplot as plt
import matplotlib
import numpy
import os
import struct
import tkFileDialog
import tkMessageBox
import ttk
import zlib
import tables
# Dev Imports
#import timeit
#import inspect

tables.parameters.MAX_NUMEXPR_THREADS = None
tables.parameters.MAX_BLOSC_THREADS = None

# File Parameters
EXTENSION = ".mzXML"
EXTRACTION = "aligned"
OUTPUT = "Summary.txt"

# Alignment Parameters
ALIGNMENT_TIME_WINDOW = 10      # The +/- time window that the program is allowed to look for the feature for alignment (EIC time axis)
ALIGNMENT_MASS_WINDOW = 0.1     # The +/- m/z window (not charge state corrected) that is used to detect the feature used for alignment. Afterwards a spline fit is used to detect the measured time
ALIGNMENT_BACKGROUND_MULTIPLIER = 2 # The multiplier of the timewindow used for background determination
ALIGNMENT_S_N_CUTOFF = 9        # The minimum S/N value of a feature to be used for alignment
ALIGNMENT_MIN_PEAK = 5          # The minimum number of features used for alignment

# Calibration Parameters
SUM_SPECTRUM_RESOLUTION = 100   # Number of data points per 1 whole m/z unit
CALIB_MASS_WINDOW = 0.3         # This +/- mass window used to detect the accurate mass of a calibra
CALIB_S_N_CUTOFF = 9            # The minimum S/N value of a feature to be used for calibration
CALIB_MIN_PEAK = 3              # Minimum number of calibrants

# PARAMETERS
MASS_MODIFIERS = []             # The mass modifiers refer to changes to the analyte.
                                # Charge carrier should NOT go here (the program assumes that your analytes are protonated)

# Extraction Parameters
EXTRACTION_TYPE = 2             # 1 = Max, 0 = Total and 2 = Area
MASS_WINDOW = 0.2               # The +/- m/z window used around each feature for extraction
TIME_WINDOW = 8                 # The +/- time window that will be used around a cluster, to create the sum spectrum
MIN_CHARGE = 2                  # The minimum charge state that the program will integrate for all features (unless overwritten in the composition file)
MAX_CHARGE = 3                  # The maximum charge state that the program will integrate for all features (unless overwritten in the composition file)
#MIN_CONTRIBUTION = 0.01        # Minimum contribution to isotopic distrubition to be included (NOT BEING USED ATM)
MIN_TOTAL = 0.99                # Desired contribution of extracted isotopes of total isotopic pattern
BACKGROUND_WINDOW = 10          # Total m/z window (+ and -) to search for background

# The maximum distance between distinct isotopic masses to be 'pooled'
EPSILON = 0.5                   # DO NOT TOUCH THIS UNLESS YOU KNOW WTF YOU ARE DOING! Read below if you truly want to know the meaning:
                                # This value represents the maximum distance (in Da) for which the element specific isotopic mass defect will be combined

# Isotopic Mass Differences
C = [('13C',0.0107,1.00335)]
H = [('2H',0.00012,1.00628)]
N = [('15N',0.00364,0.99703)]
O18 = [('18O',0.00205,2.00425)]
O17 = [('17O',0.00038,1.00422)]
S33 = [('33S',0.0076,0.99939)]
S34 = [('34S',0.0429,1.9958)]
S36 = [('36S',0.0002,3.99501)]

# Building block properties
BLOCKS = {  #######################
			# Structural Features #
			#######################
				###################
				# Monosaccharides #
				###################
					'F':{'mass':146.05790879894,
						'carbons':6,
						'hydrogens':10,
						'nitrogens':0,
						'oxygens':4,
						'sulfurs':0},
					'H':{'mass':162.0528234185,
						'carbons':6,
						'hydrogens':10,
						'nitrogens':0,
						'oxygens':5,
						'sulfurs':0},
					'N':{'mass':203.07937251951,
						'carbons':8,
						'hydrogens':13,
						'nitrogens':1,
						'oxygens':5,
						'sulfurs':0},
					'S':{'mass':291.09541650647,
						'carbons':11,
						'hydrogens':17,
						'nitrogens':1,
						'oxygens':8,
						'sulfurs':0},
					'L':{'mass':273.08485182277,
						'carbons':11,
						'hydrogens':15,
						'nitrogens':1,
						'oxygens':7,
						'sulfurs':0},
					'M':{'mass':305.11106657061,
						'carbons':12,
						'hydrogens':19,
						'nitrogens':1,
						'oxygens':8,
						'sulfurs':0},
					'E':{'mass':319.12671663475,
						'carbons':13,
						'hydrogens':21,
						'nitrogens':1,
						'oxygens':8,
						'sulfurs':0},
				#########################
				# Mouse Monosaccharides #
				#########################
					'G':{'mass':307.0903311261,
						'carbons':11,
						'hydrogens':17,
						'nitrogens':1,
						'oxygens':9,
						'sulfurs':0},
					'Gl':{'mass':289.0797664424,
						'carbons':11,
						'hydrogens':15,
						'nitrogens':1,
						'oxygens':8,
						'sulfurs':0},
					'Ge':{'mass':335.1216312544,
						'carbons':13,
						'hydrogens':21,
						'nitrogens':1,
						'oxygens':9,
						'sulfurs':0},
				#######################
				# Sugar Modifications #
				#######################
					'P':{'mass':79.96633088875,
						'carbons':0,
						'hydrogens':1,
						'nitrogens':0,
						'oxygens':3,
						'sulfurs':0},
					'Su':{'mass':79.95681485868,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':3,
						'sulfurs':1},
					'Ac':{'mass':42.0105646837,
						'carbons':2,
						'hydrogens':2,
						'nitrogens':0,
						'oxygens':1,
						'sulfurs':0},
				##############################
				# Reducing End Modifications #
				##############################
					'aa':{'mass':139.06332853255,
						'carbons':7,
						'hydrogens':9,
						'nitrogens':1,
						'oxygens':2,
						'sulfurs':0},
					'ab':{'mass':138.07931294986,
						'carbons':7,
						'hydrogens':10,
						'nitrogens':2,
						'oxygens':1,
						'sulfurs':0},
					'free':{'mass':18.0105646837,
						'carbons':0,
						'hydrogens':2,
						'nitrogens':0,
						'oxygens':1,
						'sulfurs':0},
			###################
			# Charge Carriers #
			###################
					#################
					# Positive Mode #
					#################
					'sodium':{'mass':22.9897692809,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'potassium':{'mass':38.96370668,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'proton':{'mass':1.007276466812,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					#################
					# Negative Mode #
					#################
					'protonLoss':{'mass':-1.007276466812,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'electron':{'mass':0.00054857990946,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
			############
			# Elements #
			############
					'_H':{'mass':1.007825,
						'carbons':0,
						'hydrogens':1,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'_C':{'mass':12.000000,
						'carbons':1,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'_N':{'mass':14.003074,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':1,
						'oxygens':0,
						'sulfurs':0},
					'_O':{'mass':15.994915,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':1,
						'sulfurs':0},
					'_S':{'mass':31.972071,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':1},
					'_P':{'mass':30.973761,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':1},
					'_F':{'mass':18.998403,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'_Na':{'mass':22.989770,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
					'_K':{'mass':38.963708,
						'carbons':0,
						'hydrogens':0,
						'nitrogens':0,
						'oxygens':0,
						'sulfurs':0},
			#################
			# Custom Blocks #
			#################
					####################
					# Immunoglobulin G #
					####################
					'IgGI':{'mass':1188.5047,	# Get exacter mass
						'carbons':50,
						'hydrogens':72,
						'nitrogens':14,
						'oxygens':20,
						'sulfurs':0},
					'IgGIV':{'mass':1172.5098,	# Get exacter mass
						'carbons':50,
						'hydrogens':72,
						'nitrogens':14,
						'oxygens':19,
						'sulfurs':0},
					'IgGII':{'mass':1156.5149,	# Get exacter mass
						'carbons':50,
						'hydrogens':72,
						'nitrogens':14,
						'oxygens':18,
						'sulfurs':0},
					####################
					# Immunoglobulin A #
					####################
					'Q':{'mass':4135.882086,
						'carbons':177,
						'hydrogens':270,
						'nitrogens':50,
						'oxygens':59,
						'sulfurs':3},
					'R':{'mass':2962.590442,
						'carbons':128,
						'hydrogens':219,
						'nitrogens':37,
						'oxygens':41,
						'sulfurs':1},
					'T':{'mass':2346.1348023,
						'carbons':101,
						'hydrogens':163,
						'nitrogens':27,
						'oxygens':33,
						'sulfurs':2},
					'U':{'mass':2183.0709257,
						'carbons':92,
						'hydrogens':154,
						'nitrogens':26,
						'oxygens':33,
						'sulfurs':2}}
UNITS = BLOCKS.keys()

###################
# DATA STRUCTURES #
###################
class Analyte():
    def __init__(self):
        self.composition = None
        self.mass = None
        self.massWindow = None
        self.time = None
        self.timeWindow = None
        self.minCharge = None
        self.maxCharge = None
        self.isotopes = None

class Isotope():
    def __init__(self):
        self.isotope = None
        self.charge = None
        self.mass = None
        self.obsInt = None
        self.obsMax = None
        self.expInt = None
        self.qc = None
        self.background = None
        self.backgroundPoint = None
        self.noise = None

###############################
# Start of actual application #
###############################
class App():
    def __init__(self,master):
        # VARIABLES
        self.master = master
        self.inputFile = ""
        self.inputFileIdx = 0
        self.refFile = ""
        self.alFile = ""
        self.calFile = IntVar()
        self.ptFile = None
        self.rmMZXML = IntVar()
        self.batchFolder = ""
        self.batchProcessing = 0
        self.batchWindow = 0
        self.dataWindow = 0
        self.outputWindow = 0
        self.analyteIntensity = IntVar()
        self.analyteRelIntensity = IntVar()
        self.analyteBackground = IntVar()
        self.analyteNoise = IntVar()
        self.analytePerCharge = IntVar()
        self.analyteBckSub = IntVar()
        self.alignmentQC = IntVar()
        self.ppmQC = IntVar()
        self.qualityControl = IntVar()
        self.SN = IntVar()
        self.log = True
        self.noise = "RMS"
        #self.noise = "MM"
        self.batch = False
        self.fig = matplotlib.figure.Figure(figsize=(8, 6))

        # The LacyTools Logo (created by ....)
        if os.path.isfile('./UI.png'):
            image = plt.imread('./UI.png')
            plt.axis('off')
            plt.tight_layout()
            im = plt.imshow(image)
        # The Canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master = master)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, root)
        self.canvas.get_tk_widget().pack(fill=BOTH,expand=YES)
        self.canvas.draw()

        # FRAME
        frame = Frame(master)
        master.title("LaCyTools")

        # MENU
        menu = Menu(root)
        root.config(menu = menu)

        filemenu = Menu(menu,tearoff=0)
        menu.add_cascade(label="File", menu=filemenu)
        filemenu.add_command(label="Open Input File", command = self.openFile)

        extractmenu = Menu(menu,tearoff=0)
        menu.add_cascade(label="Extraction", menu=extractmenu)
        extractmenu.add_command(label="Open ref file", command = self.openRefFile)
        extractmenu.add_command(label="Extract", command = self.extractData)

        menu.add_command(label="Batch Process", command = lambda: self.batchPopup(self))

        menu.add_command(label="Data Storage", command = lambda: self.dataPopup(self))

    def feature_reader(self,file):
        """ reads the contents of the 'features.txt' file and stores
        the relevant values in a list.
        """
        features = []
        with open(file,'r') as fr:
            for line in fr:
                if line[0][0].isdigit():
                    line = line.rstrip().split()
                    features.append(map(float,line))
                else:
                    continue
        return features

    def fitFunc(self, x,a,b,c):
        penalty = 0
        if b > 2.:
            penalty = abs(b-1.)*10000
        if b < 0.:
            penalty = abs(2.-b)*10000
        return a*x**b + c + penalty

    def calcQuadratic(self,data):
        """ This function fits the specified function in 'fitFunc'
        to the data, using the curve_fit package from scipy.optimize.

        INPUT: A list of (m/z,int) tuples
        OUTPUT: The parameters for the fitted function
        """
        expected = []
        observed = []
        for i in data:
            expected.append(i[0])
            observed.append(i[1])
        try:
            z = curve_fit(self.fitFunc, observed, expected)#,maxfev=10000)
            name = self.inputFile.split(".")[0]
            name = os.path.join(self.batchFolder,name)
            #############
            # Plot Code #
            #############
            minX = min(expected)-0.1*min(expected)
            maxX = max(expected)+0.1*max(expected)
            newX = numpy.linspace(minX,maxX,2500*(maxX-minX))
            linY = newX
            yNew = self.fitFunc(newX,*z[0])
            minY = self.fitFunc(minX,*z[0])
            maxY = self.fitFunc(maxX,*z[0])
            fig =  plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            plt.scatter(expected,observed,c='b',label='Raw',alpha=0.5)
            observedCalibrated = []
            for index, j in enumerate(observed):
                observedCalibrated.append(self.fitFunc(j,*z[0]))
            plt.scatter(expected,observedCalibrated,c='r',label='Calibrated',marker='s',alpha=0.5)
            numbers = ["%.2f" % number for number in z[0]]
            if float(numbers[2]) > 0.0:
                plt.plot(newX,yNew,label='Fit, Function: '+str(numbers[0])+"x"+"$^{"+str(numbers[1])+"}$+"+str(numbers[2]),c='b')
            else:
                plt.plot(newX,yNew,label='Fit, Function: '+str(numbers[0])+"x"+"$^{"+str(numbers[1])+"}$"+str(numbers[2]),c='b')
            plt.plot(newX,linY,label='Target',c='r',linestyle='--')
            plt.legend(loc='best')
            plt.xlabel("Expected rt (s.)")
            plt.ylabel("Observed rt (s.)")
            plt.xlim(minX,maxX)
            plt.ylim(minY,maxY)
            fig.savefig(name,dpi=800)
            plt.close()
            ###############
            # end of plot #
            ###############
        except RuntimeError:
            z = None
        return z

    def dataPopup(self,master):
        if master.dataWindow == 1:
            return
        master.dataWindow = 1
        self.folder = StringVar()
        self.ptFileName = StringVar()
        def close(self):
            master.dataWindow = 0
            top.destroy()
        def batchButton():
            master.openBatchFolder()
            self.folder.set(master.batchFolder)
        top = self.top = Toplevel()
        top.protocol( "WM_DELETE_WINDOW", lambda: close(self))
        self.batchDir = Button(top, text = "Batch Directory", width = 25, command = lambda: batchButton())
        self.batchDir.grid(row = 0, column = 0, sticky = W)
        self.batch = Label(top, textvariable = self.folder, width = 25)
        self.batch.grid(row = 0, column = 1)
        self.remove = Checkbutton(top, text = "Remove mzXML files", variable = master.rmMZXML, onvalue = 1, offvalue = 0)
        self.remove.grid(row = 1, column = 0, sticky = W)
        self.convertButton = Button(top, text = "Batch Convert to pyTables", width = 25, command = lambda: master.batchConvert(master))
        self.convertButton.grid(row = 2, column = 0,columnspan = 2)

    def batchConvert(self,master):
        """ TODO: COMMENT THIS FUNCTION PLEASE.
        This function does x, using Y

        INPUT: stuff
        OUTPUT: stuff
        """
        import time
        start_time = time.time()

        filenames = glob.glob(str(self.batchFolder)+"/*" + EXTENSION)

        print "Converting..."

        filename = filenames[0]
        array = []
        self.inputFile = filename
        self.readData(array,None)

        nscans = len(array)
        size = 0
        for rt, spectrum in array:
            size = max(size, len(spectrum))

        SCAN_SIZE = int(size * 1.1)

        try:
            rawfile = tables.open_file(os.path.join(self.batchFolder, "pytables.h5"), "w", filters=tables.Filters(complevel=4, complib="blosc:lz4"))
        except tables.HDF5ExtError:
            print "Error creating pyTables file"
            raise

        class Scan(tables.IsDescription):
            sample = tables.Int64Col(pos=0)
            scan = tables.Int64Col(pos=1)
            rt = tables.Float64Col(pos=2)
            art = tables.Float64Col(pos=3)      # aligned retention time
            idx = tables.Int64Col(pos=4)
            size = tables.Int64Col(pos=5)

        rawfile.create_vlarray('/', 'filenames', atom=tables.VLUnicodeAtom(), expectedrows=len(filenames))
        rawfile.create_table('/', 'scans', description=Scan, expectedrows=len(filenames)*nscans)
        rawfile.create_earray('/', 'mzs', atom=tables.Float64Atom((SCAN_SIZE,)), shape=(0,), chunkshape=(1,))
        rawfile.create_earray('/', 'Is', atom=tables.Int64Atom((SCAN_SIZE,)), shape=(0,), chunkshape=(1,))

        row = rawfile.root.scans.row

        idx = 0

        # main loop
        for count, filename in enumerate(filenames):
            self.inputFile = filename
            if count >= 1:
                array = []
                self.readData(array,None)

            mzs = numpy.zeros((len(array), SCAN_SIZE), numpy.float64)
            Is = numpy.zeros((len(array), SCAN_SIZE), numpy.int64)

            # loop over spectra
            for scan, spectrum in enumerate(array):
                rt, spectrum = spectrum
                size = min(len(spectrum), SCAN_SIZE)

                spectrum = numpy.array(spectrum).T
                spectrum[0, 1:] = numpy.diff(spectrum[0])

                mzs[scan, :size], Is[scan, :size] = spectrum[:, :size]

                row['sample'] = count
                row['scan'] = scan
                row['rt'] = rt
                row['idx'] = idx
                row['size'] = size

                row.append()

                idx += 1

            rawfile.root.mzs.append(mzs)
            rawfile.root.Is.append(Is)

            rawfile.root.filenames.append(filename)

            if self.rmMZXML.get() == 1:
                try:
                    os.remove(filename)
                except:
                    raise

        rawfile.close()
        print "Finished converting."
        end_time = time.time()
        print "Batch convertion lasted for", (end_time - start_time) / 60., "minutes, or", (end_time - start_time) / len(filenames), "seconds per sample."
        tkMessageBox.showinfo("Status Message","Batch Convert finished on "+str(datetime.now()))

    def batchProcess(self,master):
        import time
        start = time.time()
        self.batch = True
        # Check if reference or alignment file was selected
        if self.refFile == "" and self.alFile == "" and self.calFile == "":
            tkMessageBox.showinfo("File Error","No reference or alignment file selected")
        # Check for pytables file
        if os.path.isfile(os.path.join(self.batchFolder,"pytables.h5")):
            ptFileName = os.path.join(self.batchFolder,"pytables.h5")
            if self.ptFile is None:
                self.ptFile = tables.open_file(ptFileName, mode='a')
            filenames = self.ptFile.root.filenames[:]
            self.readData = self.readPTData
            self.transform_mzXML = self.alignRTs
            filenames2idx = dict([(filename, idx) for idx, filename in enumerate(filenames)])
            print 'Found "pytables.h5" in batch folder.'
        else:
            filenames = glob.glob(os.path.join(str(self.batchFolder),"*"+EXTENSION))
            filenames2idx = dict([(filename, idx) for idx, filename in enumerate(filenames)])
        # ALIGNMENT
        if self.alFile != "":
            features = []
            features = self.feature_reader(self.alFile)
            features = sorted(features, key = lambda tup: tup[1])
            # reset aligned rts to 0
            if self.ptFile is not None and self.ptFile.isopen:
                for scan in self.ptFile.root.scans:
                    scan['art'] = 0
                    scan.update()
                self.ptFile.flush()
            for file in filenames:
                array = []
                timePairs = []
                self.inputFile = file
                self.inputFileIdx = filenames2idx[file]
                readTimes = self.matchFeatureTimes(features)
                self.readData(array,readTimes)
                for i in features:
                    peakTime = 0
                    peakIntensity = 0
                    dataPoints = []
                    leftPoints = []
                    rightPoints = []
                    signalPoints = []
                    for j in array:
                        if j[0] > i[1] - 2*ALIGNMENT_TIME_WINDOW and j[0] < i[1] - ALIGNMENT_TIME_WINDOW:
                            dataPoints.append(self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW))
                            leftPoints.append(self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW))
                        if j[0] > i[1] + ALIGNMENT_TIME_WINDOW and j[0] < i[1] + 2*ALIGNMENT_TIME_WINDOW:
                            dataPoints.append(self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW))
                            rightPoints.append(self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW))
                        if j[0] < i[1] + ALIGNMENT_TIME_WINDOW and j[0] > i[1] - ALIGNMENT_TIME_WINDOW:
                            signalPoints.append(self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW))
                            if self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW) > peakIntensity:
                                peakIntensity = self.feature_finder(j[1],i[0]-ALIGNMENT_MASS_WINDOW,i[0]+ALIGNMENT_MASS_WINDOW)
                                peakTime = j[0]
                    ############################################################################################
                    # Awesome cool method (with milk and cookies!) to determine background and noise in an EIC #
                    ############################################################################################
                    sortedData = sorted(dataPoints)
                    startSize = int(0.25 * float(len(sortedData)))
                    currSize = startSize
                    currAverage = numpy.average(sortedData[0:currSize])
                    if self.noise == "MM":
                        currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                    elif self.noise == "RMS":
                        currNoise = numpy.std(sortedData[0:currSize])
                    for k in range(0,len(sortedData)-(startSize+1)):
                        if sortedData[currSize+1] < currAverage + 3 * currNoise:
                            currSize += 1
                            currAverage =  numpy.average(sortedData[0:currSize])
                            if self.noise == "MM":
                                currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                            elif self.noise == "RMS":
                                currNoise = numpy.std(sortedData[0:currSize])
                        else:
                            break
                    background = currAverage
                    noise = currNoise
                    ######################
                    # End of awesomeness #
                    ######################
                    # Plot Code #
                    #############
                    """plotPoints = leftPoints + signalPoints + rightPoints
                    fig =  plt.figure()
                    ax = fig.add_subplot(111)
                    plt.plot(sorted(plotPoints))
                    #plt.plot(plotPoints)
                    plt.axhline(y=currAverage, color ='k')
                    plt.axhline(y=currAverage + 3* currNoise, color = 'r')
                    #plt.axvline(x=len(plotPoints)/3, color = 'r')
                    #plt.axvline(x=(len(plotPoints)/3)*2, color = 'r')
                    plt.show()"""
                    ###############
                    # end of plot #
                    ###############
                    if (peakIntensity - background)/noise > ALIGNMENT_S_N_CUTOFF:
                        timePairs.append((i[1],peakTime))
                    else:
                        print "Feature: "+str(i)+" was not above ALIGNMENT_S_N_CUTOFF: "+str(ALIGNMENT_S_N_CUTOFF)
                        if self.log == True:
                            with open('LaCyTools.log', 'a') as flog:
                                flog.write(str(datetime.now())+"\tFeature: "+str(i)+" was not above alignment S/N cutoff: "+str(ALIGNMENT_S_N_CUTOFF)+" in file: "+str(file)+"\n")
                # Make sure that enough features are used for alignment
                if len(timePairs) >= ALIGNMENT_MIN_PEAK:
                    warnUser = False
                    alignFunction = self.calcQuadratic(timePairs)
                    #except TypeError:
                    if alignFunction == None:
                        if self.log == True:
                            with open('LaCyTools.log', 'a') as flog:
                                flog.write(str(datetime.now())+"\tFile: "+str(file)+" could not be aligned. Curve_fit exceeded maximum number of iterations\n")
                        outFile = os.path.split(file)[-1]
                        outFile = "unaligned_"+outFile
                        outFile = os.path.join(self.batchFolder,outFile)
                        open(outFile,'w').close()
                        continue
                    # Create alignment output file
                    alignmentOutput = self.inputFile.split(".")[0]
                    alignmentOutput = alignmentOutput + ".alignment"
                    with open(alignmentOutput,'w') as falign:
                        lsq = 0
                        falign.write("Peak\tExpected RT\tOriginal RT\tAligned RT\n")
                        for index,timePair in enumerate(timePairs):
                            falign.write(str(features[index][0])+"\t"+str(timePair[0])+"\t"+str(timePair[1])+"\t"+str(self.fitFunc(float(timePair[1]),*alignFunction[0]))+"\n")
                            lsq += float(features[index][0]) - self.fitFunc(float(timePair[1]),*alignFunction[0])
                    self.transform_mzXML(file,alignFunction[0])
                else:
                    print "File not aligned due to lack of features"
                    outFile = os.path.split(file)[-1]
                    outFile = "unaligned_"+outFile
                    outFile = os.path.join(self.batchFolder,outFile)
                    open(outFile,'w').close()
        # (CALIBRATION AND) EXTRACTION
        if self.refFile != "":
            if self.analyteIntensity.get() == 0 and self.analyteRelIntensity.get() == 0 and self.analyteBackground.get() == 0 and self.analyteNoise.get() == 0 and self.alignmentQC.get() == 0 and self.qualityControl.get() == 0 and self.ppmQC.get() == 0 and self.SN.get() == 0:
                tkMessageBox.showinfo("Output Error","No outputs selected")
            self.initCompositionMasses(self.refFile)
            ref = []
            self.refParser(ref)
            times = []
            for i in ref:
                times.append(i[4])
            chunks = collections.OrderedDict()
            for i in times:
                if i not in chunks.keys():
                    chunks['%s' % i] = []
            for i in ref:
                chunks['%s' % i[4]].append(i)
            if os.path.isfile(os.path.join(self.batchFolder,"pytables.h5")) == False:
                filenames = glob.glob(os.path.join(str(self.batchFolder),EXTRACTION+"*"+EXTENSION))
                filenames2idx = dict([(filename, idx) for idx, filename in enumerate(filenames)])
            for file in filenames:
                results = []
                self.inputFile = file
                self.inputFileIdx = filenames2idx[file]
                array = []
                readTimes = self.matchAnalyteTimes(ref)
                self.readData(array, readTimes)
                for index,i in enumerate(chunks.keys()):
                    spectrum = self.sumSpectrum(i,array)
                    calibrants = []
                    # Calibrate the sum spectrum
                    if self.calFile.get() == 1:
                        for j in ref:
                            if j[6] == "True" and int(round(float(j[4]))) == int(round(float(i))):
                                charge = j[0].split("_")[1]
                                calibrants.append((float(j[1]),int(charge)))
                        measuredMaxima = self.getLocalMaxima(calibrants,spectrum)
                        presentCalibrants = self.getObservedCalibrants(measuredMaxima,calibrants)
                        measuredMaximaMZ = []
                        # Strip the m/z values from the maxima
                        for j in measuredMaxima:
                            measuredMaximaMZ.append(j[0])
                        # Perform 2d degree polynomial fit
                        if len(measuredMaximaMZ) >= CALIB_MIN_PEAK:
                            z = numpy.polyfit(measuredMaximaMZ,presentCalibrants,2) # This should be the correct one
                        else:
                            print "Unable to calibrate the sum spectrum at "+str(i)+" seconds"
                            outFile = os.path.split(str(self.inputFile))[1]
                            outFile = outFile.split(".")[0]
                            outFile = "Uncalibrated_sumSpectrum_"+str(i)+"_"+str(outFile)+".xy"
                            outFile = os.path.join(str(self.batchFolder),outFile)
                            with open(outFile,'w') as fw:
                                fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                            continue
                        #z = numpy.polyfit(presentCalibrants,measuredMaximaMZ,2)
                        f = numpy.poly1d(z)
                        calOut = str(file.split(".")[0])+"_"+str(index)+".calibration"
                        with open(calOut,'w') as fw2:
                            for index,j in enumerate(measuredMaximaMZ):
                                fw2.write("accurate mass: "+str(presentCalibrants[index])+" measured at "+str(j) +" being calibrated to: "+str(f(j))+"\n")
                        #y = f(measuredMaximaMZ)
                        outputBatch = []
                        mzList = []
                        intList = []
                        for j in spectrum:
                            mzList.append(float(j[0]))
                            intList.append(int(j[1]))
                        # Transform python list into numpy array
                        mzArray = numpy.array(mzList)
                        newArray = f(mzArray)
                        newSpectrum = []
                        for index,j in enumerate(newArray):
                            newSpectrum.append((j,intList[index]))
                        spectrum = newSpectrum
                        outFile = os.path.split(str(self.inputFile))[1]
                        outFile = outFile.split(".")[0]
                        outFile = "Calibrated_sumSpectrum_"+str(i)+"_"+str(outFile)+".xy"
                        outFile = os.path.join(str(self.batchFolder),outFile)
                        with open(outFile,'w') as fw:
                            fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                    else:
                        outFile = os.path.split(str(self.inputFile))[1]
                        outFile = outFile.split(".")[0]
                        outFile = "sumSpectrum_"+str(i)+"_"+str(outFile)+".xy"
                        outFile = os.path.join(str(self.batchFolder),outFile)
                        with open(outFile,'w') as fw:
                            fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                    self.extractData(chunks[i],spectrum,results)
                self.writeResults(results,file)
                master.batchProcess = 0
            self.combineResults()
        if self.ptFile is not None:
            self.ptFile.close()
        end = time.time()
        print "Batch process lasted for", (end - start) / 60., "minutes"
        tkMessageBox.showinfo("Status Message","Batch Process finished on "+str(datetime.now()))

    def writeCalibration(self,function,array):
        """ TODO
        """
        endian = "!"
        started = False
        with open(self.inputFile,'r') as fr:
            name = os.path.split(str(self.inputFile))[-1]
            name = name.split(".")[0]
            name = "calibrated_"+name+".mzXML" # TODO: Make the extension dynamic
            with open(name,'w') as fw:
                counter = 0
                mzList = []
                intList = []
                values = []
                for line in fr:
                    if 'zlib' in line:
                        fw.write(line)
                        compression = True
                    elif 'byteOrder' in line:
                        fw.write(line)
                        byteOrder = line.split("byteOrder")[1]
                        byteOrder = byteOrder.split("\"")[1]
                        endian = "!"
                        if byteOrder == 'little':
                            endian = '<'
                        elif byteOrder == 'big':
                            endian = '>'
                    elif 'precision' in line:
                        fw.write(line)
                        precision = line.split("precision")[1]
                        precision = precision.split("\"")[1]
                        if int(precision) == 64:
                            precision = 'd'
                        else:
                            precision = 'f'
                    elif 'contentType="m/z-int">' in line:
                        mzList = []
                        intList = []
                        values = []
                        for i in array[counter][1]:
                            mzList.append(i[0])
                            intList.append(i[1])
                        mzArray = numpy.array(mzList)
                        newArray = function(mzArray)
                        for index, i in enumerate(newArray):
                            values.append(i)
                            values.append(intList[index])
                        format = str(endian)+str(len(values))+precision
                        data = struct.pack(format, *values)
                        if compression == True:
                            data = zlib.compress(data)
                        data = base64.b64encode(data)
                        fw.write('contentType="m/z-int">'+str(data)+'</peaks>\n')
                        counter += 1
                    else:
                        fw.write(line)

    def getObservedCalibrants(self,maxima,potentialCalibrants):
        """ This function compares the list of local maxima with the
        expected calibrants. The function will perceive the observed
        local maxima that is closest to a desired calibrant as being the
        m/z where the calibrant was observed in the spectrum. The
        function then appends the theoretical m/z value of a calibrants
        that were actually observed to a list (actualCalibrants) which
        is returned at the end of the function.

        INPUT 1: A list of floats containg the observed local maxima (of
        the spline fit within each inclusion range, assuming that they
        were above user specified S/N cut off).
        INPUT 2: A list of floats containing the theoretical m/z of all
        calibrants.
        OUTPUT: A list of floats containing the theoretical m/z of the
        calibrants which were near an oberved local maxima.
        """
        actualCalibrants = []
        for i in maxima:
            diff = 4.0
            closest = 0
            for j in potentialCalibrants:
                if abs(float(j[0])-float(i[0])) < diff:
                    diff = abs(float(j[0])-float(i[0]))
                    closest = float(j[0])
            actualCalibrants.append(closest)
        return actualCalibrants

    def getLocalMaxima(self,features,spectrum):
        """ TODO
        """
        maxima = []
        for i in features:
            mass, charge = i
            window = CALIB_MASS_WINDOW / charge
            lowMz = self.binarySearch(spectrum,float(mass)-float(window),len(spectrum)-1,'left')
            highMz = self.binarySearch(spectrum,float(mass)+float(window),len(spectrum)-1,'right')
            maximum = (0,0)
            x_points = []
            y_points = []
            for j in spectrum[lowMz:highMz]:
                x_points.append(j[0])
                y_points.append(j[1])
            newX = numpy.linspace(x_points[0],x_points[-1],2500*(x_points[-1]-x_points[0]))
            f = interp1d(x_points,y_points, kind='cubic')
            ySPLINE = f(newX)
            # Plot Code (for testing purposes)
            """fig =  plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(x_points, y_points, 'b*')
            #plt.plot(newX,newY, 'b--')
            plt.plot(newX,ySPLINE,'r--')
            #plt.legend(['Raw Data','Guassian (All Points)','Cubic Spline'], loc='best')
            plt.legend(['Raw Data','Cubic Spline'], loc='best')
            plt.show()"""
            for index, j in enumerate(ySPLINE):
                if j > maximum[1]:
                    maximum = (newX[index],j)
            # Check if maxima above S/N cut-off
            values = self.getBackground(spectrum, maximum[0], 1, MASS_WINDOW)
            background,noise = values[0], values[2]
            if maximum[1] > background + CALIB_S_N_CUTOFF * noise:
                maxima.append(maximum)
        return maxima

    def readCalibrationFeatures(self):
        """ This function reads the calibration file and returns the
        features in a tuple containg the lower time and upper time values
        followed by a list of the m/z coordinates

        INPUT: None
        OUTPUT: Tuple containing (lowTime, highTime, [m/z coordinatse])
        """
        with open(self.calFile,'r') as fr:
            firstLine = fr.readline()
            lowTime, highTime = firstLine.strip().split("\t")
            mz = []
            for line in fr:
                mz.append(float(line.strip()))
        return (float(lowTime), float(highTime), mz)

    def sumSpectrum(self,time,array):
        """ This function creates a summed spectrum and returns the
        resulting spectrum back to the calling function.
        """
        # This is returning None's now
        lowTime = self.binarySearch(array,float(time)-TIME_WINDOW,len(array)-1,'left')
        highTime = self.binarySearch(array,float(time)+TIME_WINDOW,len(array)-1,'right')
        LOW_MZ = 25000.0
        HIGH_MZ = 0.0
        for i in array[lowTime:highTime]:
            if i[1][0][0] < LOW_MZ:
                LOW_MZ = i[1][0][0]
            if i[1][-1][0] > HIGH_MZ:
                HIGH_MZ = i[1][-1][0]
        # This should be dynamically determined
        arraySize = (float(HIGH_MZ) - float(LOW_MZ)) * float(SUM_SPECTRUM_RESOLUTION)
        combinedSpectra = numpy.zeros(shape=(arraySize+2,2))
        bins = []
        for index, i in enumerate(combinedSpectra):
            i[0] = float(LOW_MZ) + index*(float(1)/float(SUM_SPECTRUM_RESOLUTION))
            bins.append(float(LOW_MZ) + index*(float(1)/float(SUM_SPECTRUM_RESOLUTION)))
        fullSet = []
        mz = []
        start = datetime.now()
        for i in array[lowTime:highTime]:
            for j in i[1]:
                fullSet.append(j)
                mz.append(j[0])
        fullSet.sort(key = lambda tup: tup[0])
        mz.sort()
        mzArray = numpy.asarray(mz)
        binsArray = numpy.asarray(bins)
        test = numpy.searchsorted(binsArray,mzArray)
        for index, i in enumerate(fullSet):
            try:
                combinedSpectra[test[index]][1] += i[1]
            except:
                # We ignore the data points at m/z edge, if they are important
                # then the user should do a proper measurement.
                pass
        return combinedSpectra

    def findNearest(self,array,value):
        if value >= array[0][0] and value <= array[-1][0]:
            diff = 1
            # First Pass
            a = 0
            b = len(array)
            while a < b:
                mid = (a+b)//2
                if array[mid][0] > value:
                    b = mid
                else:
                    a = mid+1
            if array[a][0] - value < diff:
                diff = array[a][0] - value
                index = a
            # Second Pass
            a = 0
            b = len(array)
            while a < b:
                mid = (a+b)//2
                if array[mid][0] < value:
                    a=mid+1
                else:
                    b=mid
            if array[a][0] - value < diff:
                diff = array[a][0] - value
                index = a
            return a

    def transform_mzXML(self,file,polynomial):
        """Reads the mzXML file and transforms the reported retention
        time by the specified polynomial function.
        """
        with open(file,'r') as fr:
            outFile = os.path.split(file)[-1]
            outFile = "aligned_"+outFile
            outFile = os.path.join(self.batchFolder,outFile)
            print "Writing output file: "+outFile
            with open(outFile,'w') as fw:
                for line in fr:
                    if 'retentionTime' in line:
                        time = line.strip()
                        time = time.split("\"")
                        for index,i in enumerate(time):
                            if 'retentionTime' in i:
                                time=time[index+1]
                                break
                        if time[0] == 'P':
                            time = time[2:-1]
                        # The below line is only to make it work with mzMine
                        if  self.fitFunc(float(time),*polynomial) < 0:
                            newTime = str(0)
                        else:
                            newTime = str(self.fitFunc(float(time),*polynomial))
                        line = line.replace(time,newTime)
                        fw.write(line)
                    else:
                        fw.write(line)

    def alignRTs(self,file,polynomial):
        """Reads the mzXML file and transforms the reported retention
        time by the specified polynomial function.
        """
        print "Aligning file", self.inputFile
        i = self.inputFileIdx
        for row in self.ptFile.root.scans.where("sample == i"):
            time = row['rt']
            # The below line is only to make it work with mzMine
            if  self.fitFunc(time,*polynomial) > 0:
                row['art'] = self.fitFunc(time,*polynomial)
            row.update()
        self.ptFile.flush()

    def feature_finder(self,data,lowMass,highMass):
        low, high = 0, len(data)
        while low != high:
            mid = (low + high) // 2
            mz = data[mid][0]
            if mz < lowMass:
                low = mid + 1
            else:
                high = mid
        start, high = low, len(data)
        while low != high:
            mid = (low + high) // 2
            mz = data[mid][0]
            if mz > highMass:
                high = mid
            else:
                low = mid + 1
        end = high
        return max(I for mz, I in data[start:end])

    def combineResults(self):
        """ This function reads all the raw files and creates the summary
        output file.
        """
        total = []
        for file in glob.glob(os.path.join(str(self.batchFolder),"*.raw")):
            compositions = []
            trigger = 0
            results = []
            with open(file,'r') as fr:
                name = str(file)
                name = os.path.split(str(name))[-1]
                current = None
                for _ in xrange(2):
                    next(fr)
                for line in fr:
                    if not line:
                        break
                    line = line.strip().split("\t")
                    if current:
                        if current.composition == line[0] and current.time == line[5]:
                            foo = Isotope()
                            foo.isotope = line[2]
                            foo.mass = float(line[3])
                            foo.measMass = float(line[4])
                            foo.charge = line[1]
                            foo.obsInt = float(line[9])
                            foo.obsMax = float(line[13])
                            foo.expInt = float(line[6])
                            foo.background = float(line[10])
                            foo.backgroundPoint = float(line[11])
                            foo.noise = float(line[12])
                            current.isotopes.append(foo)
                        else:
                            results.append(current)
                            current = Analyte()
                            current.composition = line[0]
                            current.time = line[5]
                            current.timeWindow = line[8]
                            current.massWindow = line[7]
                            current.isotopes = []
                            foo = Isotope()
                            foo.isotope = line[2]
                            foo.mass = float(line[3])
                            foo.measMass = float(line[4])
                            foo.charge = line[1]
                            foo.obsInt = float(line[9])
                            foo.obsMax = float(line[13])
                            foo.expInt = float(line[6])
                            foo.background = float(line[10])
                            foo.backgroundPoint = float(line[11])
                            foo.noise = float(line[12])
                            current.isotopes.append(foo)
                    else:
                        current = Analyte()
                        current.composition = line[0]
                        current.time = line[5]
                        current.timeWindow = line[8]
                        current.massWindow = line[7]
                        current.isotopes = []
                        foo = Isotope()
                        foo.isotope = line[2]
                        foo.mass = float(line[3])
                        foo.measMass = float(line[4])
                        foo.charge = line[1]
                        foo.obsInt = float(line[9])
                        foo.obsMax = float(line[13])
                        foo.expInt = float(line[6])
                        foo.background = float(line[10])
                        foo.backgroundPoint = float(line[11])
                        foo.noise = float(line[12])
                        current.isotopes.append(foo)
            results.append(current)
            total.append((name,results))
        for file in glob.glob(str(self.batchFolder)+"/unaligned*"+EXTENSION):
            name = str(file)
            name = os.path.split(str(name))[-1]
            total.append((name,[]))
        total.sort()

        # Test chunk to see if class conversion worked
        """for i in total:
            for j in i[1]:
                print j.composition
                if j.composition == "IgGI1H3N4F1":
                    for k in j.isotopes:
                        print k.isotope, k.charge, k.obsInt"""

        #################################
        # Generate the summaryFile name #
        #################################
        utc_datetime = datetime.utcnow()
        s = utc_datetime.strftime("%Y-%m-%d-%H%MZ")
        filename = s +"_"+OUTPUT
        summaryFile = os.path.join(self.batchFolder,filename)

        ####################################################################
        # Get list of analytes, required for correct alignment of clusters #
        ####################################################################
        compositions = []
        with open(self.refFile,'r') as fr:
            for line in fr:
                if line[0] == "#":
                    continue
                parts=line.rstrip('\n').split('\t')
                compositions.append((parts[0],parts[1]))

        #############################
        # Start writing the results #
        #############################
        with open(summaryFile,'w') as fw:
            ##############
            # Parameters #
            ##############
            fw.write("Parameter Settings\n")
            if self.alFile != "":
                fw.write("Alignment Parameters\n")
                fw.write("ALIGNMENT_TIME_WINDOW\t"+str(ALIGNMENT_TIME_WINDOW)+"\n")
                fw.write("ALIGNMENT_MASS_WINDOW\t"+str(ALIGNMENT_MASS_WINDOW)+"\n")
                fw.write("ALIGNMENT_S_N_CUTOFF\t"+str(ALIGNMENT_S_N_CUTOFF)+"\n")
                fw.write("ALIGNMENT_MIN_PEAK\t"+str(ALIGNMENT_MIN_PEAK)+"\n")
            if self.calFile.get() == 1:
                fw.write("Calibration Parameters\n")
                fw.write("CALIB_MASS_WINDOW\t"+str(CALIB_MASS_WINDOW)+"\n")
                fw.write("CALIB_S_N_CUTOFF\t"+str(CALIB_S_N_CUTOFF)+"\n")
                fw.write("CALIB_MIN_PEAK\t"+str(CALIB_MIN_PEAK)+"\n")
            if self.refFile != "":
                fw.write("Extraction Parameters\n")
                fw.write("SUM_SPECTRUM_RESOLUTION\t"+str(SUM_SPECTRUM_RESOLUTION)+"\n")
                fw.write("MASS_WINDOW\t"+str(MASS_WINDOW)+"\n")
                fw.write("TIME_WINDOW\t"+str(TIME_WINDOW)+"\n")
                fw.write("MIN_CHARGE\t"+str(MIN_CHARGE)+"\n")
                fw.write("MAX_CHARGE\t"+str(MAX_CHARGE)+"\n")
                fw.write("MIN_TOTAL\t"+str(MIN_TOTAL)+"\n")
                fw.write("BACKGROUND_WINDOW\t"+str(BACKGROUND_WINDOW)+"\n\n")

            ##############################
            # Analyte Absolute Intensity #
            ##############################
            if self.analyteIntensity.get() == 1 and self.analyteBckSub.get() == 0:
                ##########################
                # Combined charge states #
                ##########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Abs Int")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += l.obsInt
                                except AttributeError:
                                    pass
                            if sumInt > 0:
                                fw.write("\t"+str(sumInt))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("Abs Int ("+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.obsInt
                                    except AttributeError:
                                        pass
                                if sumInt > 0:
                                    fw.write("\t"+str(sumInt))
                                else:
                                    fw.write("\t")
                            fw.write("\n")
                        fw.write("\n")

            ######################################################
            # Analyte Absolute Intensity (Background subtracted) #
            ######################################################
            if self.analyteIntensity.get() == 1 and self.analyteBckSub.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Abs Int (Bck Sub)")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += max(0, l.obsInt - l.background)
                                except AttributeError:
                                    pass
                            if sumInt > 0:
                                fw.write("\t"+str(sumInt))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("Abs Int (Bck Sub, "+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += max(0, m.obsInt - m.background)
                                    except AttributeError:
                                        pass
                                if sumInt > 0:
                                    fw.write("\t"+str(sumInt))
                                else:
                                    fw.write("\t")
                            fw.write("\n")
                        fw.write("\n")

            ##############################
            # Analyte Relative Intensity #
            ##############################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 0:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Rel Int")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        totalIntensity = 1
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            totalIntensity += l.obsInt
                                except AttributeError:
                                    pass
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += l.obsInt
                                except AttributeError:
                                    pass
                            if sumInt > 0:
                                fw.write("\t"+str(float(sumInt)/float(totalIntensity)))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("Rel Int ("+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            totalIntensity = 1
                            for k in compositions:
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    totalIntensity += m.obsInt
                                    except AttributeError:
                                        pass
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.obsInt
                                    except AttributeError:
                                        pass
                                if sumInt > 0:
                                    fw.write("\t"+str(float(sumInt)/float(totalIntensity)))
                                else:
                                    fw.write("\t")
                            fw.write("\n")
                        fw.write("\n")

            ##############################################
            # Relative Intensity (Background Subtracted) #
            ##############################################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Rel Int (Bck Sub)")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        totalIntensity = 1
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            totalIntensity += max(0, l.obsInt - l.background)
                                except AttributeError:
                                    pass
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += max(0, l.obsInt - l.background)
                                except AttributeError:
                                    pass
                            if sumInt > 0:
                                fw.write("\t"+str(float(sumInt)/float(totalIntensity)))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("Rel Int (Bck Sub, "+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            totalIntensity = 1
                            for k in compositions:
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    totalIntensity += max(0, m.obsInt - m.background)
                                    except AttributeError:
                                        pass
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += max(0, m.obsInt - m.background)
                                    except AttributeError:
                                        pass
                                if sumInt > 0:
                                    fw.write("\t"+str(float(sumInt)/float(totalIntensity)))
                                else:
                                    fw.write("\t")
                            fw.write("\n")
                        fw.write("\n")

            ################################
            # Analyte Background Intensity #
            ################################
            if self.analyteBackground.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("BCK")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += l.background
                                except AttributeError:
                                    pass
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("BCK ("+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.background
                                    except AttributeError:
                                        pass
                                fw.write("\t"+str(sumInt))
                            fw.write("\n")
                        fw.write("\n")


            #######################
            # Analyte Noise Value #
            #######################
            if self.analyteNoise.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Noise")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            sumInt += l.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        charge = 0
                        for j in total:
                            if len(j[1]) == len(compositions):
                                for k in j[1]:
                                    if k.composition == i[0] and float(k.time) == float(i[1]):
                                        for l in k.isotopes:
                                            if l.charge != charge:
                                                charge = l.charge
                                                if masses == "":
                                                    masses ="["+str(l.mass)+"]"
                                                else:
                                                    masses += " ["+str(l.mass)+"]"
                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for i in total:
                        fw.write(str(i[0]))
                        for j in compositions:
                            sumInt = 0
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            sumInt += l.noise
                                except AttributeError:
                                    pass
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                    fw.write("\n")

                ####################
                # Per charge state #
                ####################
                if self.analytePerCharge.get() == 1:
                    minCharge = sys.maxint
                    maxCharge = 0
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    if k.composition == j[0] and float(k.time) == float(j[1]):
                                        for l in k.isotopes:
                                            if int(l.charge) < minCharge:
                                                minCharge = int(l.charge)
                                            elif int(l.charge) > maxCharge:
                                                maxCharge = int(l.charge)
                                except AttributeError:
                                    pass
                    for i in xrange(minCharge,maxCharge+1):
                        # This is a time intensive function
                        # Header
                        fw.write("Noise ("+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas (not being correct?)
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.expInt
                                    break
                            fw.write("\t"+str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in total:
                                if len(k[1]) == len(compositions):
                                    for l in k[1]:
                                        if l.composition == j[0] and float(l.time) == float(j[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    masses ="["+str(m.mass)+"]"
                                                    break
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        for j in total:
                            fw.write(str(j[0]))
                            for k in compositions:
                                sumInt = 0
                                for l in j[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == float(k[1]):
                                            for m in l.isotopes:
                                                if int(m.charge) == i:
                                                    sumInt += m.noise
                                    except AttributeError:
                                        pass
                                fw.write("\t"+str(sumInt))
                            fw.write("\n")
                        fw.write("\n")

            #######################
            # Alignment Residuals #
            #######################
            if self.alignmentQC.get() == 1:
                # Get results
                totalResults = []
                for file in glob.glob(os.path.join(str(self.batchFolder),"*.alignment")):
                    resultBuffer = []
                    with open (file,'r') as fr:
                        for line in fr:
                            line = line.strip().split()
                            resultBuffer.append(line)
                    totalResults.append((file,resultBuffer))
                # Header
                header = []
                for i in totalResults:
                    if len(i[1]) > len(header):
                        header = i[1][:]
                fw.write("Alignment Residuals")
                for i in header[1:]:
                    fw.write("\t"+str(i[0]))
                fw.write("\tRMS\n")
                # Actual Data
                for i in totalResults:
                    RMS = 0
                    fw.write(str(i[0]))
                    for j in header[1:]:
                        flag = 0
                        for k in i[1]:
                            if j[0] == k[0]:
                                fw.write("\t"+str(float(k[3])-float(k[1])))
                                RMS += (float(k[3])-float(k[1]))**2
                                flag = 1
                        if flag == 0:
                            fw.write("\t")
                    fw.write("\t"+str(math.sqrt(RMS))+"\n")
                fw.write("\n")

            ###############################
            # Analyte Mass Error (in PPM) #
            ###############################
            print self.ppmQC.get()
            if self.ppmQC.get() == 1:
                minCharge = sys.maxint
                maxCharge = 0
                for i in total:
                    for j in compositions:
                        for k in i[1]:
                            try:
                                if k.composition == j[0] and float(k.time) == float(j[1]):
                                    for l in k.isotopes:
                                        if int(l.charge) < minCharge:
                                            minCharge = int(l.charge)
                                        elif int(l.charge) > maxCharge:
                                            maxCharge = int(l.charge)
                            except AttributeError:
                                pass
                for i in xrange(minCharge,maxCharge+1):
                    # This is a time intensive function
                    # Header
                    fw.write("PPM Error ("+str(i)+"+)")
                    for j in compositions:
                        fw.write("\t"+str(j[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                sumInt += m.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                masses ="["+str(m.mass)+"]"
                                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual Data
                    for j in total:
                        fw.write(str(j[0]))
                        for k in compositions:
                            relContribution = 0.0
                            targetMass = 0.0
                            actualMass = 0.0
                            for l in j[1]:
                                try:
                                    if l.composition == k[0] and float(l.time) == float(k[1]):
                                        for m in l.isotopes:
                                            if m.expInt > relContribution and int(m.charge) == i:
                                                relContribution = m.expInt
                                                targetMass = m.mass
                                                actualMass = m.measMass
                                except AttributeError:
                                    pass
                            try:
                                ppm = ((actualMass - targetMass) / targetMass) * 1000000
                                fw.write("\t"+str(ppm))
                            except ZeroDivisionError:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

            ###############
            # Isotopic QC #
            ###############
            if self.qualityControl.get() == 1:
                minCharge = sys.maxint
                maxCharge = 0
                for i in total:
                    for j in compositions:
                        for k in i[1]:
                            try:
                                if k.composition == j[0] and float(k.time) == float(j[1]):
                                    for l in k.isotopes:
                                        if int(l.charge) < minCharge:
                                            minCharge = int(l.charge)
                                        elif int(l.charge) > maxCharge:
                                            maxCharge = int(l.charge)
                            except AttributeError:
                                pass
                for i in xrange(minCharge,maxCharge+1):
                    # This is a time intensive function
                    # Header
                    fw.write("QC ("+str(i)+"+)")
                    for j in compositions:
                        fw.write("\t"+str(j[0]))
                    fw.write("\n")
                    # List of theoretical areas (not being correct?)
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                sumInt += m.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                masses ="["+str(m.mass)+"]"
                                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for j in total:
                        fw.write(str(j[0]))
                        for k in compositions:
                            sumInt = 0
                            qc = 0
                            for l in j[1]:
                                try:
                                    if l.composition == k[0] and float(l.time) == float(k[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                sumInt += max(float(m.obsInt) - float(m.background),0)
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                try:
                                                    maxIntensityBackCorrected = max(float(m.obsInt) - float(m.background),0)
                                                    qc += abs((maxIntensityBackCorrected  / float(sumInt)) - m.expInt)
                                                except ZeroDivisionError:
                                                    pass
                                except AttributeError:
                                    pass
                            if qc > 0:
                                fw.write("\t"+str(qc))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

            ##############################
            # Signal to Background ratio #
            ##############################
            if self.SN.get() == 1:
                minCharge = sys.maxint
                maxCharge = 0
                for i in total:
                    for j in compositions:
                        for k in i[1]:
                            try:
                                if k.composition == j[0] and float(k.time) == float(j[1]):
                                    for l in k.isotopes:
                                        if int(l.charge) < minCharge:
                                            minCharge = int(l.charge)
                                        elif int(l.charge) > maxCharge:
                                            maxCharge = int(l.charge)
                            except AttributeError:
                                pass
                for i in xrange(minCharge,maxCharge+1):
                    # This is a time intensive function
                    # Header
                    fw.write("S/N ("+str(i)+"+)")
                    for j in compositions:
                        fw.write("\t"+str(j[0]))
                    fw.write("\n")
                    # List of theoretical areas (not being correct?)
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                sumInt += m.expInt
                                break
                        fw.write("\t"+str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in total:
                            if len(k[1]) == len(compositions):
                                for l in k[1]:
                                    if l.composition == j[0] and float(l.time) == float(j[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                masses ="["+str(m.mass)+"]"
                                                break
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for j in total:
                        fw.write(str(j[0]))
                        for k in compositions:
                            expInt = 0
                            SN = 0
                            for l in j[1]:
                                try:
                                    if l.composition == k[0] and float(l.time) == float(k[1]):
                                        for m in l.isotopes:
                                            if m.expInt > expInt and int(m.charge) == i:
                                                try:
                                                    SN = (m.obsMax - m.backgroundPoint) / m.noise
                                                except ZeroDivisionError:
                                                    pass
                                                expInt = m.expInt
                                except AttributeError:
                                    pass
                            if SN > 0:
                                fw.write("\t"+str(SN))
                            else:
                                fw.write("\t")
                        fw.write("\n")
                    fw.write("\n")

    def writeResults(self,results,file):
        """ This function writes the resultes per file away to a raw
        file.
        """
        outFile = os.path.split(file)[-1]
        outFile = outFile.split(".")[0]
        outFile = outFile+".raw"
        outFile = os.path.join(self.batchFolder,outFile)
        with open(outFile,'w') as fw:
            fw.write(str(file)+"\n")
            fw.write("Composition\tCharge\tIsotope\tExact Mass\tAccurate Mass\tTime\tTheor Area\tMass Window\tTime Window\tArea\tBackground Area\tBackground Point\tNoise\tMax Intensity\n")
            for index,i in enumerate(results):
                composition,charge,isotope = i[4][0].split("_")
                fw.write(str(composition)+"\t"+str(charge)+"\t"+str(isotope)+"\t"+str(i[4][1])+"\t"+str(i[2])+"\t"+str(i[4][4])+"\t"+str(i[4][2])+"\t"+str(i[4][3])+"\t"+str(i[4][5])+"\t"+str(i[0])+"\t"+str(i[1][1])+"\t"+str(i[1][0])+"\t"+str(i[1][2])+"\t"+str(i[3])+"\n")

    def batchPopup(self,master):
        """ TODO

        INPUT: None
        OUTPUT: None
        """
        if master.batchWindow == 1:
            return
        master.batchWindow = 1
        self.al = StringVar()
        self.ref = StringVar()
        self.folder = StringVar()
        def alButton():
            master.openAlFile()
            self.al.set(master.alFile)
        def refButton():
            master.openRefFile()
            self.ref.set(master.refFile)
        def batchButton():
            master.openBatchFolder()
            self.folder.set(master.batchFolder)
        def close(self):
            master.batchWindow = 0
            top.destroy()
        top = self.top = Toplevel()
        top.protocol( "WM_DELETE_WINDOW", lambda: close(self))
        self.aligns = Button(top, text = "Alignment File", widt = 25, command = lambda: alButton())
        self.aligns.grid(row = 2, column = 0, sticky = W)
        self.alLabel = Label(top, textvariable = self.al, width = 25)
        self.alLabel.grid(row = 2, column = 1)
        self.calibrate = Checkbutton(top, text = "Calibration", variable = master.calFile, onvalue = 1, offvalue = 0)
        self.calibrate.grid(row = 3, column = 0, sticky = W)
        self.compos = Button(top, text = "Reference File", width = 25, command = lambda: refButton())
        self.compos.grid(row = 4, column = 0, sticky = W)
        self.com = Label(top, textvariable = self.ref, width = 25)
        self.com.grid(row = 4, column = 1)
        self.batchDir = Button(top, text = "Batch Directory", width = 25, command = lambda: batchButton())
        self.batchDir.grid(row = 5, column = 0, sticky = W)
        self.batch = Label(top, textvariable = self.folder, width = 25)
        self.batch.grid(row = 5, column = 1)
        self.output = Button(top, text = "Output Format", width = 25, command = lambda: master.outputPopup(master))
        self.output.grid(row = 6, column = 0,columnspan = 2)
        self.run = Button(top, text = "Run Batch Process", width = 25, command = lambda: master.batchProcess(master))
        self.run.grid(row = 7, column = 0, columnspan = 2)
        #top.lift()
        # Couple the attributes to button presses
        top.attributes("-topmost", True)

    def outputPopup(self,master):
        """ This function creates a pop up box to specify what output
        should be shown in the final summary. The default value for all
        variables is off (0) and by ticking a box it is set to on (1).

        INPUT: None
        OUTPUT: None
        """
        if master.outputWindow == 1:
            return
        master.outputWindow = 1
        def select_all(self):
            master.analyteIntensity.set(1)
            master.analyteRelIntensity.set(1)
            master.analyteBackground.set(1)
            master.analyteNoise.set(1)
            master.analytePerCharge.set(1)
            master.analyteBckSub.set(1)
            master.alignmentQC.set(1)
            master.qualityControl.set(1)
            master.ppmQC.set(1)
            master.SN.set(1)
        def select_none(self):
            master.analyteIntensity.set(0)
            master.analyteRelIntensity.set(0)
            master.analyteBackground.set(0)
            master.analyteNoise.set(0)
            master.analytePerCharge.set(0)
            master.analyteBckSub.set(0)
            master.alignmentQC.set(0)
            master.qualityControl.set(0)
            master.ppmQC.set(0)
            master.SN.set(0)
        def close(self):
            master.outputWindow = 0
            top.destroy()
        top = self.top = Toplevel()
        top.protocol( "WM_DELETE_WINDOW", lambda: close(self))
        self.all = Button(top, text = "Select All", command = lambda: select_all(self))
        self.all.grid(row = 0, column = 0, sticky = W)
        self.none = Button(top, text = "Select None", command = lambda: select_none(self))
        self.none.grid(row = 0, column = 1, sticky = E)
        self.text1 = Label(top, text = "Base Outputs", font="bold")
        self.text1.grid(row = 1, column = 0, sticky = W)
        self.text2 = Label(top, text = "Output Modifiers", font="bold")
        self.text2.grid(row = 1, column = 1, sticky = W)
        # Analyte Intensity (*,#)
        self.ai = Checkbutton(top, text = u"Analyte Intensity\u00B9\u00B7\u00B2", variable = master.analyteIntensity, onvalue = 1, offvalue = 0)
        self.ai.grid(row = 2, column = 0, sticky = W)
        self.ri = Checkbutton(top, text = u"Relative Intensity\u00B9\u00B7\u00B2", variable = master.analyteRelIntensity, onvalue = 1, offvalue = 0)
        self.ri.grid(row = 3, column = 0, sticky = W)
        self.back = Checkbutton(top, text = u"Analyte Background\u00B9", variable = master.analyteBackground, onvalue = 1, offvalue = 0)
        self.back.grid(row = 4, column = 0, sticky = W)
        self.analNoise = Checkbutton(top, text = u"Analyte Noise\u00B9", variable = master.analyteNoise, onvalue = 1, offvalue = 0)
        self.analNoise.grid(row = 5, column = 0, sticky = W)
        self.chargeState = Checkbutton(top, text = u"\u00B9Areas per Charge State", variable = master.analytePerCharge, onvalue = 1, offvalue = 0)
        self.chargeState.grid(row = 2, column = 1, sticky = W)
        self.bckSub = Checkbutton(top, text = u"\u00B2Background subtracted Areas", variable = master.analyteBckSub, onvalue = 1, offvalue = 0)
        self.bckSub.grid(row = 3, column = 1, sticky = W)
        self.align = Checkbutton(top, text="Alignment Residuals", variable=master.alignmentQC, onvalue=1, offvalue=0)
        self.align.grid(row = 6, column=0, sticky=W)
        self.qc = Checkbutton(top, text = "QC", variable = master.qualityControl, onvalue = 1, offvalue = 0)
        self.qc.grid(row = 7, column = 0, sticky = W)
        self.ppm = Checkbutton(top, text = "PPM Error", variable = master.ppmQC, onvalue = 1, offvalue = 0)
        self.ppm.grid(row = 8, column = 0, sticky = W)
        self.snratio = Checkbutton(top, text = "Signal to Noise", variable = master.SN, onvalue = 1, offvalue = 0)
        self.snratio.grid(row = 9, column = 0, sticky = W)
        self.button = Button(top,text='Ok',command = lambda: close(self))
        self.button.grid(row = 10, column = 0, columnspan = 2)
        top.lift()

    def binarySearch(self, array, target, high, direction):
        """Returns element number directly to the left in array 'array'
        of specified element 'target', assuming 'array[x][0]' is sorted,
        if direction is set as 'left'.

        The return value a is such that all elements in array[:a] have
        element < target, and all e in array[a:] have element >= target.

        Returns element number directly to the right in array 'array'
        of specified element 'target', assuming 'array[x][0]' is sorted,
        if direction is set as 'right'

        The return value a is such that all elements in array[:a] have
        element <= target, and all e in array[a:] have element > target.

        former left"""
        if target >= array[0][0] and target <= array[high][0]:
            a = 0
            b = high
            while a < b:
                mid=(a+b)//2
                if direction == 'left':
                    if array[mid][0] < target:
                        a=mid+1
                    else:
                        b=mid
                if direction == 'right':
                    if array[mid][0] > target:
                        b = mid
                    else:
                        a = mid+1
            return a

    def openFile(self):
        """ This function opens a Tkinter filedialog, asking the user
        to select a file. The chosen file is then read (by the readData
        function) and the read data is used to plot the selected spectrum
        on the screen (by the plotData function).

        INPUT: None
        OUTPUT: None
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'inputFile',file_path)

    def openCalFile(self):
        """ TODO
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'calFile',file_path)

    def openBatchFolder(self):
        """ This function opens a Tkinter filedialog, asking the user
        to select a directory. The chosen directory is then set to the
        self.batchFolder variable.

        INPUT: None
        OUTPUT: None
        """
        folder_path = tkFileDialog.askdirectory()
        if not folder_path:
            pass
        else:
            setattr(self,'batchFolder',folder_path)

    def openRefFile(self):
        """ PLACE HOLDER.
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'refFile',file_path)

    def openAlFile(self):
        """ PLACE HOLDER.
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'alFile',file_path)

    def processBlock(self, block, array, readTimes):
        """ This function processes a data block as taken from the input
        file.
        """
        """if "scan num" in block:
            scan = block.split("scan num")[1]
            scan = scan.split("\"")[1]
        """

        if "retentionTime" in block:
            rt = block.split("retentionTime")[1]
            rt = rt.split("\"")[1]
            if rt[0] == 'P':
                rt = rt[2:-1]

        """if "peaksCount" in block:
            peaks = block.split("peaksCount")[1]
            peaks = peaks.split("\"")[1]
        """

        # FIX not to catch zlib in encoded data
        if '"zlib"' in block:
            compression = True
        # FIX for implicit no compression
        else:
            compression = False

        if "byteOrder" in block:
            byteOrder = block.split("byteOrder")[1]
            byteOrder = byteOrder.split("\"")[1]

        if "precision" in block:
            precision = block.split("precision")[1]
            precision = precision.split("\"")[1]

        # FIX pairOrder is Bruker format bending
        if "contentType" in block or "pairOrder" in block:
            peaks = block.split('"m/z-int">')[1]
            peaks = peaks.split("</peaks>")[0]

        if peaks:
            if readTimes:
                flag = 0
                for i in readTimes:
                    if float(rt) >= i[0] and float(rt) <= i[1]:
                        self.mzXMLDecoder(rt, peaks, precision, compression, byteOrder, array)
                        flag = 1
                if flag == 0:
                    array.append((float(rt),None))
            else:
                self.mzXMLDecoder(rt, peaks, precision, compression, byteOrder, array)


    ######################################################
    # START OF FUNCTIONS RELATED TO PARSING ANALYTE FILE #
    ######################################################
    def getChanceNetwork(self,(mass,carbons,hydrogens,nitrogens,oxygens17,oxygens18,sulfurs33,sulfurs34,sulfurs36)):
        """ This function calculates the total chance network based on
        all the individual distributions. The function multiplies all
        the chances to get a single chance for a single option.

        INPUT: A list containing the Analyte m/z followed by several
        other lists (1 for each isotopic state).
        OUTPUT: A list of float tuples (isotopic m/z, isotopic chance)
        """
        totals = []
        for x in itertools.product(carbons,hydrogens,nitrogens,oxygens17,oxygens18,sulfurs33,sulfurs34,sulfurs36):
            i, j, k, l, m, n, o, p = x
            totals.append((mass+i[0]+j[0]+k[0]+l[0]+m[0]+n[0]+o[0]+p[0],
                           i[1]*j[1]*k[1]*l[1]*m[1]*n[1]*o[1]*p[1]))
        return totals

    def mergeChances(self,totals):
        """ This function merges all the isotopic chances based on the
        specified resolution of the machine.

        INPUT: A list of float tuples (isotopic m/z, isotopic chance)
        OUTPUT: A sorted list of float tuples (isotopic m/z, isotopic
        chance).
        """
        results = []
        newdata = {d: True for d in totals}
        for k, v in totals:
            if not newdata[(k,v)]: continue
            newdata[(k,v)] = False
            # use each piece of data only once
            keys,values = [k*v],[v]
            for kk, vv in [d for d in totals if newdata[d]]:
                if abs(k-kk) < EPSILON:
                    keys.append(kk*vv)
                    values.append(vv)
                    newdata[(kk,vv)] = False
            results.append((sum(keys)/sum(values),sum(values)))
        return results

    def calcDistribution(self, element, number):
        """ This function calculates the fraction of the total intensity
        that is present in each isotope of the given element based on
        a binomial distribution. The function takes the name of the
        element and the number of atoms of said element as an input and
        returns a list of (m/z,fraction) tuples. The number of isotopes
        that is returned is dependant on the distribution, once fractions
        fall below 0.01 the function stops.

        INPUT1: A string containing the code for the element (ie 33S)
        INPUT2: An integer listing the number of atoms
        OUTPUT: A list of float tuples (isotope m/z, isotope fraction).
        """
        fractions = []
        for i in element:
            j = 0
            while j <= number:
                nCk = math.factorial(number) / (math.factorial(j) * math.factorial(number - j))
                f = nCk * i[1]**j * (1 - i[1])**(number-j)
                fractions.append((i[2]*j,f))
                j+= 1
                if f < 0.01:
                    break
        return fractions

    def parseAnalyte(self,Analyte):
        """ This function splits the Analyte input string into a parts
        and calculates the total number of each element of interest per
        Analyte. The function will then attach further elements based on
        the user specified mass modifiers before calling the isotopic
        distribution function. The function finally returns a list
        containing the analyte mass and distribution lists for each
        isotopic state.

        INPUT: A string containing the Analyte (ie 'H4N4')
        OUTPUT: A list containing the Analyte m/z followed by several
        other lists (1 for each isotopic state).
        """
        results = []
        mass = 0
        numCarbons = 0
        numHydrogens = 0
        numNitrogens = 0
        numOxygens = 0
        numSulfurs = 0
        totalElements = 0
        units = ["".join(x) for _,x in itertools.groupby(Analyte,key=str.isdigit)]
        # Calculate the bass composition values
        for index,j in enumerate(units):
            for k in UNITS:
                    if j == k:
                        mass += float(BLOCKS[k]['mass']) * float(units[index+1])
                        numCarbons += int(BLOCKS[k]['carbons']) * int(units[index+1])
                        numHydrogens += int(BLOCKS[k]['hydrogens']) * int(units[index+1])
                        numNitrogens += int(BLOCKS[k]['nitrogens']) * int(units[index+1])
                        numOxygens += int(BLOCKS[k]['oxygens']) * int(units[index+1])
                        numSulfurs += int(BLOCKS[k]['sulfurs']) * int(units[index+1])
        # Attach the mass modifier values
        for j in MASS_MODIFIERS:
            mass += float(BLOCKS[j]['mass'])
            numCarbons += float(BLOCKS[j]['carbons'])
            numHydrogens += int(BLOCKS[j]['hydrogens'])
            numNitrogens += int(BLOCKS[j]['nitrogens'])
            numOxygens += int(BLOCKS[j]['oxygens'])
            numSulfurs += int(BLOCKS[j]['sulfurs'])
        # Calculate the distribution for the given value
        carbons = self.calcDistribution(C,numCarbons)
        hydrogens = self.calcDistribution(H,numHydrogens)
        nitrogens = self.calcDistribution(N,numNitrogens)
        oxygens17 = self.calcDistribution(O17,numOxygens)
        oxygens18 = self.calcDistribution(O18,numOxygens)
        sulfurs33 = self.calcDistribution(S33,numSulfurs)
        sulfurs34 = self.calcDistribution(S34,numSulfurs)
        sulfurs36 = self.calcDistribution(S36,numSulfurs)
        return ((mass,carbons,hydrogens,nitrogens,oxygens17,oxygens18,sulfurs33,sulfurs34,sulfurs36))

    def initCompositionMasses(self, file):
        """ This function reads the composition file. Calculates the
        masses for the compositions read from the composition file.
        The function then calculates the mass and fraction of total
        ions that should be theoretically present in. The final output
        is a modified reference list containing each analyte's structure
        and window followed by a list of isotope m/z and isotopic
        fraction.

        INPUT: A string containing the path of the composition file
        OUTPUT: None
        """
        lines = []
        with open(file,'r') as fr:
            for line in fr:
                line = line.rstrip()
                lines.append(line)
        # Chop composition into sub units and get exact mass & carbon count
        analyteFile = os.path.join(self.batchFolder,"analytes.ref")
        #if os.path.exists(analyteFile) == True:
        #   print "USING EXISTING REFERENCE FILE"
        #   return
        print "PRE-PROCESSING REFERENCE FILE"
        print "THIS MAY TAKE A WHILE"
        with open(analyteFile,'w') as fw:
            fw.write("# Peak\tm/z\tRel Area\twindow\trt\ttime window\Calibration\n")
            for i in lines:
                if i[0] == "#":
                    continue
                i = i.split("\t")
                # Initialize variables
                massWindow = MASS_WINDOW
                timeWindow = TIME_WINDOW
                minCharge = MIN_CHARGE
                maxCharge = MAX_CHARGE
                calibration = False
                # Check optional variables
                if len(i) >= 2:
                    time = float(i[1])
                if len(i) > 2:
                    if i[2]:
                        massWindow = float(i[2])
                if len(i) > 3:
                    if i[3]:
                        timeWindow = int(i[3])
                if len(i) > 4:
                    if i[4]:
                        minCharge = int(i[4])
                if len(i) > 5:
                    if i[5]:
                        maxCharge = int(i[5])
                if len(i) > 6:
                    if i[6]:
                        calibration = True
                # End of variable check
                values =  self.parseAnalyte(i[0])
                totals = self.getChanceNetwork(values)
                results = self.mergeChances(totals)
                results.sort(key=lambda x: x[0])
                # Make the range inclusive
                for j in range(minCharge,maxCharge+1):
                    contribution = 0.0
                    maxIsotope = max(results,key=lambda tup:tup[1])[1]
                    for index,k in enumerate(results):
                        contribution += float(k[1])
                        if calibration == True and k[1] == maxIsotope:
                            fw.write(str(i[0])+"_"+str(j)+"_"+str(index)+"\t"+str((k[0]+j*C[0][2])/j)+"\t"+str(k[1])+"\t"+str(massWindow)+"\t"+str(time)+"\t"+str(timeWindow)+"\tTrue\n")
                        else:
                            fw.write(str(i[0])+"_"+str(j)+"_"+str(index)+"\t"+str((k[0]+j*C[0][2])/j)+"\t"+str(k[1])+"\t"+str(massWindow)+"\t"+str(time)+"\t"+str(timeWindow)+"\tFalse\n")
                        #if k[1] < MIN_CONTRIBUTION:
                        #   break
                        if contribution > MIN_TOTAL:
                            break
        print "PRE-PROCESSING COMPLETE"

    ####################################################
    # END OF FUNCTIONS RELATED TO PARSING ANALYTE FILE #
    ####################################################

    def extractData(self,ref,array,results):
        """ PLACE HOLDER.
        """
        if self.refFile == "":
            tkMessageBox.showinfo("Error Message","No reference file selected")
            return
        if self.inputFile == "":
            tkMessageBox.showinfo("Error Message","No input file selected")
            return
        for i in ref:
            intensity = 0
            x_points = []
            y_points = []
            maximum = (0,0)
            maxInt = 0
            # Not pretty but it fixes the charge not being taken into account
            # with extraction and background determination
            charge = int(i[0].split("_")[-2])
            if '#' in i[0]:
                pass
            else:
                lowMz = self.binarySearch(array,float(i[1])-float(i[3]),len(array)-1,'left')
                highMz = self.binarySearch(array,float(i[1])+float(i[3]),len(array)-1,'right')
                if int(i[0].split("_")[-1]) == 0:
                    background = self.getBackground(array, float(i[1]), charge, float(i[3]))
                if lowMz and highMz:
                    try:
                        range(lowMz,highMz)
                    except TypeError:
                        print "\nReference: "+str(i[0])+" has incorrect m/z parameters"
                        raw_input("Press ENTER to exit")
                        sys.exit()
                    for k in range(lowMz, highMz):
                        # Get maximum point for S/N calculation
                        if int(array[k][1]) > maxInt:
                            maxInt = int(array[k][1])
                        if EXTRACTION_TYPE == 1:
                            if int(array[k][1]) > intensity:
                                intensity = int(array[k][1])
                        elif EXTRACTION_TYPE == 0:
                            intensity += int(array[k][1])
                        elif EXTRACTION_TYPE == 2:
                            intensity += array[k][1] * ((array[highMz][0] - array[lowMz][0]) / (highMz - lowMz))
                        # We need these to get the local maxima
                        x_points.append(array[k][0])
                        y_points.append(array[k][1])
                    ######################################################################
                    # Only spend time on doing this if we actually wanted the PPM Errors #
                    # This is not being used yet, but should!                            #
                    ######################################################################
                    if self.ppmQC.get() == 1:
                        try:
                            newX = numpy.linspace(x_points[0],x_points[-1],2500*(x_points[-1]-x_points[0]))
                            f = interp1d(x_points,y_points, kind='cubic')
                            ySPLINE = f(newX)
                            for index, j in enumerate(ySPLINE):
                                if j > maximum[1]:
                                    maximum = (newX[index],j)
                        except ValueError:
                            data = zip(x_points,y_points)
                            print "Guassian Curve Fit failed for analyte: "+str(i[0])+", reverting to non fitted local maximum"
                            for j in data:
                                if j[1] > maximum[1]:
                                    maximum = (j[0],j[1])
                        except:
                            print "Analyte: "+str(i[0])+" is being troublesome, kill it"
                else:
                    intensity = 0
                    background = (0,0,0)
                    maximum = (0,0)
            # Check if maxima above S/N cut-off
            results.append((intensity,background,maximum[0],maxInt,i))
        if self.batch == True:
            return results
        else:
            for i in results:
                print i

    def getBackground(self, array, target, charge, width):
        """ TODO
        """
        backgroundPoint = 1000000000000000000000000000000000000000000000000000      # Ridiculous start value
        totals = []
        for i in numpy.arange(-BACKGROUND_WINDOW,BACKGROUND_WINDOW,1.0/charge):
            windowAreas = []
            windowIntensities = []
            begin = self.binarySearch(array,(float(target)-i*C[0][2])-float(width),len(array)-1,'left')
            end = self.binarySearch(array,(float(target)-i*C[0][2])+float(width),len(array)-1,'right')
            if begin == None or end == None:
                print "Specified m/z value of " +str((float(target)-i*C[0][2])-float(width)) + " or " + str((float(target)-i*C[0][2])+float(width))+ " outside of spectra range"
                raw_input("Press enter to exit")
                sys.exit()
            for j in array[begin:end]:
                windowAreas.append(j[1] * ((array[end][0] - array[begin][0]) / (end - begin)))
                windowIntensities.append(j[1])
            totals.append((windowAreas,windowIntensities))
        # Find the set of 5 consecutive windows with lowest average intensity
        for i in range(0,(2*BACKGROUND_WINDOW)-4):
            mix = totals[i][1]+totals[i+1][1]+totals[i+2][1]+totals[i+3][1]+totals[i+4][1]
            avgBackground = numpy.average([sum(totals[i][0]),sum(totals[i+1][0]),sum(totals[i+2][0]),sum(totals[i+3][0]),sum(totals[i+4][0])])
            avg = numpy.average(mix)
            if avg < backgroundPoint:
                if self.noise == "RMS":
                    noise = numpy.std(mix)
                elif self.noise == "MM":
                    minNoise = 10000000000000000000000
                    maxNoise = 0
                    for k in mix:
                        if k > maxNoise:
                            maxNoise = k
                        if k < minNoise:
                            minNoise = k
                    noise = maxNoise - minNoise
                backgroundPoint = avg
                backgroundArea = avgBackground
        return (backgroundPoint,backgroundArea,noise)

    def matchFeatureTimes(self, features):
        """ This function takes a list of features/times and combines
        them into a singe list, useful for reading only relevant
        scans later in the program.

        INPUT: A list of (m/z,rt) tuples
        OUTPUT: A list of (rt,rt) tuples
        """
        wanted = []
        features = sorted(features, key=lambda x:x[1])
        current = (float(features[0][1])-ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW, float(features[0][1])+ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW)
        for i in features:
            if float(i[1])-ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW >= current[0] and float(i[1])-ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW < current[1]:
                if float(i[1])+ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW > current[1]:
                    current = (current[0],float(i[1])+ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW)
            else:
                wanted.append(current)
                current = (float(i[1])-ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW, float(i[1])+ALIGNMENT_BACKGROUND_MULTIPLIER*ALIGNMENT_TIME_WINDOW)
        wanted.append(current)
        return wanted

    def matchAnalyteTimes(self, ref):
        """ This function takes a list of references and creates a list
        of time tuples, that is needed to read only relevant scans later
        in the program.

        INPUT: A list of references (name, mz, int, window and so forth)
        OUTPUT: A list of (rt,rt) tuples
        """
        wanted = []
        ref = sorted(ref,key=lambda x:x[4])
        current = (float(ref[0][4])-float(ref[0][5]), float(ref[0][4])+float(ref[0][5]))
        for i in ref:
            if float(i[4])-float(i[5]) >= current[0] and float(i[4])-float(i[5]) < current[1]:
                if float(i[4])+float(i[5]) > current[1]:
                    current = (current[0],float(i[4])+float(i[5]))
            else:
                wanted.append(current)
                current= (float(i[4])-float(i[5]), float(i[4])+float(i[5]))
        wanted.append(current)
        return wanted

    def readData(self, array, readTimes):
        """ This function reads mzXML files and has the scans decoded on
        a per scan basis. The scans are identified by getting the line
        number of the beginning and ending tag for a scan.

        INPUT: file handle
        OUTPUT: TBA
        """
        header = True
        started = False
        block = ""
        with open(self.inputFile,'r') as fr:
            print "Processing "+str(self.inputFile)
            for number, line in enumerate(fr):
                if '</dataProcessing>' in line:
                    header = False
                if '<scan num="' in line and header == False:
                    started = True
                if started == True:
                    block +=line
                if '</scan>' in line and header == False and started == True:
                    self.processBlock(block, array, readTimes)
                    started = False
                    block = ""
            #print "Finished processing "+str(self.inputFile)

    def readPTData(self, array, readTimes):
        print "Processing "+str(self.inputFile)
        i = self.inputFileIdx

        for row in self.ptFile.root.scans.where("sample == i"):
            rt, art, idx, size = row[2:]
            if art != 0:
                rt = art
            if readTimes:
                for s, e in readTimes:
                    if rt >= s and rt <= e:
                        break
                else:
                    array.append((rt, None))
                    continue
            mzs = self.ptFile.root.mzs[idx][:size]
            Is = self.ptFile.root.Is[idx][:size]
            array.append((rt, numpy.vstack((numpy.cumsum(mzs), Is)).T))

    def refParser(self, ref):
        """Reads the reference file and fills the list 'ref' with names
        and parameters inherent to the chosen analytes to be integrated.
        """
        with open(os.path.join(self.batchFolder,"analytes.ref"),'r') as fr:
            for line in fr:
                if line[0] == "#":
                    continue
                parts=line.rstrip('\n').split('\t')
                ref.append(parts)

    def mzXMLDecoder(self, rt, peaks, precision, compression, byteOrder, array):
        """ This function parses the encoded string from an mzXML file.
        The decoded data is finally added to the
        """
        endian = ">"
        if byteOrder == 'little':
            endian = '<'
        elif byteOrder == 'big':
            endian = '>'

        # get precision
        if int(precision) == 64:
            precision = '8'
        else:
            precision = '4'

        # decode data
        data = base64.b64decode(peaks)

        # decompression
        if compression == True:
            data = zlib.decompress(data)

        data = numpy.frombuffer(data, dtype=endian + 'f' + precision)

        # format
        new = numpy.vstack((data[::2], data[1::2])).T

        # list notation
        array.append((float(rt),new))


# Call the main app
root = Tk()
app = App(root)
root.mainloop()

