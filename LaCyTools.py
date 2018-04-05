#! /usr/bin/env python
#
# Copyright 2014-2016 Bas C. Jansen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# You should have received a coyp of the Apache 2.0 license along
# with this program; if not, see 
# http://www.apache.org/licenses/LICENSE-2.0

from datetime import datetime
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.optimize
from scipy.optimize import curve_fit
from Tkinter import *
from PIL import Image, ImageTk
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
EXTENSION = ".mzXML"            # File types that will be used by MassyTools
EXTRACTION = "aligned"          # Pre-fix required for files to be extracted
OUTPUT = "Summary.txt"          # Name of the output file
SETTINGS_FILE = "Settings.txt"  # Name of the settings file

# Alignment Parameters
ALIGNMENT_TIME_WINDOW = 10      # The +/- time window that the program is allowed to look for the feature for alignment (EIC time axis)
ALIGNMENT_MASS_WINDOW = 0.1     # The +/- m/z window (not charge state corrected) that is used to detect the feature used for alignment. Afterwards a spline fit is used to detect the measured time
ALIGNMENT_BACKGROUND_MULTIPLIER = 2 # The multiplier of the timewindow used for background determination
ALIGNMENT_S_N_CUTOFF = 9        # The minimum S/N value of a feature to be used for alignment
ALIGNMENT_MIN_PEAK = 5          # The minimum number of features used for alignment

# Calibration Parameters
SUM_SPECTRUM_RESOLUTION = 100   # Number of data points per 1 whole m/z unit
CALIB_MASS_WINDOW = 0.5         # This +/- mass window (in Dalton) used to detect the accurate mass of a calibra
CALIB_S_N_CUTOFF = 9            # The minimum S/N value of a feature to be used for calibration
CALIB_MIN_PEAK = 3              # Minimum number of calibrants

# PARAMETERS
MASS_MODIFIERS = []             # The mass modifiers refer to changes to the analyte
CHARGE_CARRIER = ['proton']     # The charge carrier that is used for ionization

# Extraction Parameters
EXTRACTION_TYPE = 2             # 1 = Max, 0 = Total and 2 = Area
MASS_WINDOW = 0.2               # The +/- m/z window used around each feature for extraction
TIME_WINDOW = 8                 # The +/- time window that will be used around a cluster, to create the sum spectrum
MIN_CHARGE = 2                  # The minimum charge state that the program will integrate for all features (unless overwritten in the composition file)
MAX_CHARGE = 3                  # The maximum charge state that the program will integrate for all features (unless overwritten in the composition file)
#MIN_CONTRIBUTION = 0.01        # Minimum contribution to isotopic distrubition to be included (NOT BEING USED ATM)
MIN_TOTAL = 0.95                # Desired contribution of extracted isotopes of total isotopic pattern
BACKGROUND_WINDOW = 10          # Total m/z window (+ and -) to search for background
S_N_CUTOFF = 9                  # Minimum signal to noise value of an analyte to be included in the percentage QC

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
                        'available_for_charge_carrier':0,
                        'carbons':6,
                        'hydrogens':10,
                        'nitrogens':0,
                        'oxygens':4,
                        'sulfurs':0},
                    'H':{'mass':162.0528234185,
                        'available_for_charge_carrier':0,
                        'carbons':6,
                        'hydrogens':10,
                        'nitrogens':0,
                        'oxygens':5,
                        'sulfurs':0},
                    'N':{'mass':203.07937251951,
                        'available_for_charge_carrier':0,
                        'carbons':8,
                        'hydrogens':13,
                        'nitrogens':1,
                        'oxygens':5,
                        'sulfurs':0},
                    'S':{'mass':291.09541650647,
                        'available_for_charge_carrier':0,
                        'carbons':11,
                        'hydrogens':17,
                        'nitrogens':1,
                        'oxygens':8,
                        'sulfurs':0},
                    'L':{'mass':273.08485182277,
                        'available_for_charge_carrier':0,
                        'carbons':11,
                        'hydrogens':15,
                        'nitrogens':1,
                        'oxygens':7,
                        'sulfurs':0},
                    'M':{'mass':305.11106657061,
                        'available_for_charge_carrier':0,
                        'carbons':12,
                        'hydrogens':19,
                        'nitrogens':1,
                        'oxygens':8,
                        'sulfurs':0},
                    'E':{'mass':319.12671663475,
                        'available_for_charge_carrier':0,
                        'carbons':13,
                        'hydrogens':21,
                        'nitrogens':1,
                        'oxygens':8,
                        'sulfurs':0},
                #########################
                # Mouse Monosaccharides #
                #########################
                    'G':{'mass':307.0903311261,
                        'available_for_charge_carrier':0,
                        'carbons':11,
                        'hydrogens':17,
                        'nitrogens':1,
                        'oxygens':9,
                        'sulfurs':0},
                    'Gl':{'mass':289.0797664424,
                        'available_for_charge_carrier':0,
                        'carbons':11,
                        'hydrogens':15,
                        'nitrogens':1,
                        'oxygens':8,
                        'sulfurs':0},
                    'Ge':{'mass':335.1216312544,
                        'available_for_charge_carrier':0,
                        'carbons':13,
                        'hydrogens':21,
                        'nitrogens':1,
                        'oxygens':9,
                        'sulfurs':0},
                #######################
                # Sugar Modifications #
                #######################
                    'P':{'mass':79.96633088875,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':1,
                        'nitrogens':0,
                        'oxygens':3,
                        'sulfurs':0},
                    'Su':{'mass':79.95681485868,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':3,
                        'sulfurs':1},
                    'Ac':{'mass':42.0105646837,
                        'available_for_charge_carrier':0,
                        'carbons':2,
                        'hydrogens':2,
                        'nitrogens':0,
                        'oxygens':1,
                        'sulfurs':0},
                ##############################
                # Reducing End Modifications #
                ##############################
                    'aa':{'mass':139.06332853255,
                        'available_for_charge_carrier':0,
                        'carbons':7,
                        'hydrogens':9,
                        'nitrogens':1,
                        'oxygens':2,
                        'sulfurs':0},
                    'ab':{'mass':138.07931294986,
                        'available_for_charge_carrier':0,
                        'carbons':7,
                        'hydrogens':10,
                        'nitrogens':2,
                        'oxygens':1,
                        'sulfurs':0},
                    'free':{'mass':18.0105646837,
                        'available_for_charge_carrier':0,
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
                        'available_for_charge_carrier':1,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    'potassium':{'mass':38.96370668,
                        'available_for_charge_carrier':1,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    'proton':{'mass':1.007276466812,
                        'available_for_charge_carrier':1,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    #################
                    # Negative Mode #
                    #################
                    'protonLoss':{'mass':-1.007276466812,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    'electron':{'mass':0.00054857990946,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
            ############
            # Elements #
            ############
                    '_H':{'mass':1.007825,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':1,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    '_C':{'mass':12.000000,
                        'available_for_charge_carrier':0,
                        'carbons':1,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    '_N':{'mass':14.003074,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':1,
                        'oxygens':0,
                        'sulfurs':0},
                    '_O':{'mass':15.994915,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':1,
                        'sulfurs':0},
                    '_S':{'mass':31.972071,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':1},
                    '_P':{'mass':30.973761,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':1},
                    '_F':{'mass':18.998403,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    '_Na':{'mass':22.989770,
                        'available_for_charge_carrier':0,
                        'carbons':0,
                        'hydrogens':0,
                        'nitrogens':0,
                        'oxygens':0,
                        'sulfurs':0},
                    '_K':{'mass':38.963708,
                        'available_for_charge_carrier':0,
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
                    'IgGI':{'mass':1188.5047307674,
                        'available_for_charge_carrier':0,
                        'carbons':50,
                        'hydrogens':72,
                        'nitrogens':14,
                        'oxygens':20,
                        'sulfurs':0},
                    'IgGIV':{'mass':1172.5098161478,
                        'available_for_charge_carrier':0,
                        'carbons':50,
                        'hydrogens':72,
                        'nitrogens':14,
                        'oxygens':19,
                        'sulfurs':0},
                    'IgGII':{'mass':1156.5149015283,
                        'available_for_charge_carrier':0,
                        'carbons':50,
                        'hydrogens':72,
                        'nitrogens':14,
                        'oxygens':18,
                        'sulfurs':0},
                    ##########################    
                    # Mouse Immunoglobulin G #
                    ##########################
                    'MIgGI':{'mass':1156.514906,
                        'available_for_charge_carrier':0,
                        'carbons':50,
                        'hydrogens':72,
                        'nitrogens':14,
                        'oxygens':18,
                        'sulfurs':0},
                    'MIgGII':{'mass':996.451243,
                        'available_for_charge_carrier':0,
                        'carbons':41,
                        'hydrogens':64,
                        'nitrogens':12,
                        'oxygens':17,
                        'sulfurs':0},
                    'MIgGIII':{'mass':1114.504341,
                        'available_for_charge_carrier':0,
                        'carbons':48,
                        'hydrogens':70,
                        'nitrogens':14,
                        'oxygens':17,
                        'sulfurs':0},
                    ####################
                    # Immunoglobulin A #
                    ####################
                    'Q':{'mass':4135.882086,
                        'available_for_charge_carrier':0,
                        'carbons':177,
                        'hydrogens':270,
                        'nitrogens':50,
                        'oxygens':59,
                        'sulfurs':3},
                    'R':{'mass':2962.590442,
                        'available_for_charge_carrier':0,
                        'carbons':128,
                        'hydrogens':219,
                        'nitrogens':37,
                        'oxygens':41,
                        'sulfurs':1},
                    'T':{'mass':2346.1348023,
                        'available_for_charge_carrier':0,
                        'carbons':101,
                        'hydrogens':163,
                        'nitrogens':27,
                        'oxygens':33,
                        'sulfurs':2},
                    'U':{'mass':2183.0709257,
                        'available_for_charge_carrier':0,
                        'carbons':92,
                        'hydrogens':154,
                        'nitrogens':26,
                        'oxygens':33,
                        'sulfurs':2},
                    ###############
                    # Haptoglobin #
                    ###############
                    'HNLT':{'mass':1736.8516,
                        'available_for_charge_carrier':0,
                        'carbons':50,
                        'hydrogens':72,
                        'nitrogens':14,
                        'oxygens':20,
                        'sulfurs':0},
                    'HNLTMC':{'mass':2678.3851,
                        'available_for_charge_carrier':0,
                        'carbons':118,
                        'hydrogens':191,
                        'nitrogens':33,
                        'oxygens':36,
                        'sulfurs':1},
                    'HNHS':{'mass':972.4665,
                        'available_for_charge_carrier':0,
                        'carbons':43,
                        'hydrogens':64,
                        'nitrogens':12,
                        'oxygens':14,
                        'sulfurs':0},
                    'HNAT':{'mass':503.2704,
                        'available_for_charge_carrier':0,
                        'carbons':20,
                        'hydrogens':37,
                        'nitrogens':7,
                        'oxygens':8,
                        'sulfurs':0},
                    'HNHSNAT':{'mass':1457.7263,
                        'available_for_charge_carrier':0,
                        'carbons':63,
                        'hydrogens':99,
                        'nitrogens':19,
                        'oxygens':21,
                        'sulfurs':0},
                    'HNYS':{'mass':1269.6354,
                        'available_for_charge_carrier':0,
                        'carbons':57,
                        'hydrogens':87,
                        'nitrogens':15,
                        'oxygens':18,
                        'sulfurs':0}}
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

################################################################################################
# Tooltip code - Taken from http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml #
################################################################################################
class ToolTip(object):
    
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 27
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        try:
            # For Mac OS
            tw.tk.call("::tk::unsupported::MacWindowStyle",
                       "style", tw._w,
                       "help", "noActivates")
        except TclError:
            pass
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      wraplength=500, font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def createToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
    
###############################
# Start of actual application #
###############################
class App():
    def __init__(self,master):
        # VARIABLES
        self.master = master
        self.version = "1.0.1"
        self.build = "20180405b"
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
        self.normalizeCluster = IntVar()
        self.alignmentQC = IntVar()
        self.ppmQC = IntVar()
        self.qualityControl = IntVar()
        self.spectraQualityControl = IntVar()
        self.SN = IntVar()
        self.log = True
        # Background can be determined in two ways
        # Options are 'MIN', 'MEDIAN' and 'NOBAN'
        self.background = "MIN"
        # Nose can be determined in multiple ways
        # Options are 'RMS' and 'MM'
        self.noise = "RMS"
        self.fig = matplotlib.figure.Figure(figsize=(12, 6))

        # Attempt to retrieve previously saved settings from settingsfile
        if os.path.isfile('./'+str(SETTINGS_FILE)):
            self.getSettings()

        # The LacyTools Logo (Placeholder figure)
        if os.path.isfile('./UI/LaCyTools.png'):
            background_image = self.fig.add_subplot(111)
            image = matplotlib.image.imread('./ui/LaCyTools.png')
            background_image.axis('off')
            self.fig.set_tight_layout(True)
            background_image.imshow(image)
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
        
        menu.add_command(label="Settings", command = lambda: self.settingsPopup(self))

    def settingsPopup(self,master):
        """ This function creates a window in which the user can change
        all the parameters that are for normal use. Certain advanced
        settings such as the extraction type and noise determination
        method remain hidden from the user through this window.
        """

        def close(self):
            """ This function closes the settings popup and applies
            all the entered values to the parameters.
            """
            global ALIGNMENT_TIME_WINDOW
            global ALIGNMENT_MASS_WINDOW
            global ALIGNMENT_S_N_CUTOFF
            global ALIGNMENT_MIN_PEAK
            global CALIB_MASS_WINDOW
            global CALIB_S_N_CUTOFF
            global CALIB_MIN_PEAK
            global SUM_SPECTRUM_RESOLUTION
            global MASS_WINDOW
            global TIME_WINDOW
            global MIN_CHARGE
            global MAX_CHARGE
            global CHARGE_CARRIER
            global MIN_TOTAL
            global BACKGROUND_WINDOW
            global S_N_CUTOFF
            ALIGNMENT_TIME_WINDOW = float(self.alignTimeWindow.get())
            ALIGNMENT_MASS_WINDOW = float(self.alignMassWindow.get())
            ALIGNMENT_S_N_CUTOFF = int(self.alignSn.get())
            ALIGNMENT_MIN_PEAK = int(self.alignMin.get())
            CALIB_MASS_WINDOW = float(self.calibMassWindow.get())
            CALIB_S_N_CUTOFF = int(self.calibSn.get())
            CALIB_MIN_PEAK = int(self.calibMin.get())
            SUM_SPECTRUM_RESOLUTION = int(self.sumSpec.get())
            MASS_WINDOW = float(self.extracMassWindow.get())
            TIME_WINDOW = float(self.extracTimeWindow.get())
            MIN_CHARGE = int(self.extracMinCharge.get())
            MAX_CHARGE = int(self.extracMaxCharge.get())
            CHARGE_CARRIER = []
            for i in UNITS:
                if str(i) == master.chargeCarrierVar.get() and BLOCKS[i]['available_for_charge_carrier'] == 1:
                    CHARGE_CARRIER.append(i)
            MIN_TOTAL = float(self.extracMinTotal.get())
            BACKGROUND_WINDOW = int(self.extracBack.get())
            S_N_CUTOFF = int(self.extracSnCutoff.get())
            master.measurementWindow = 0
            top.destroy()
            
        def save(self):
            """ This function saves all changed settings to the 
            settings file.
            """
            global CHARGE_CARRIER
            CHARGE_CARRIER = []            
            for i in UNITS:
                if str(i) == master.chargeCarrierVar.get() and BLOCKS[i]['available_for_charge_carrier'] == 1:
                    CHARGE_CARRIER.append(i)
            with open(SETTINGS_FILE,'w') as fw:
                fw.write("ALIGNMENT_TIME_WINDOW\t"+str(float(self.alignTimeWindow.get()))+"\n")
                fw.write("ALIGNMENT_MASS_WINDOW\t"+str(float(self.alignMassWindow.get()))+"\n")
                fw.write("ALIGNMENT_S_N_CUTOFF\t"+str(int(self.alignSn.get()))+"\n")
                fw.write("ALIGNMENT_MIN_PEAK\t"+str(int(self.alignMin.get()))+"\n")
                fw.write("CALIB_MASS_WINDOW\t"+str(float(self.calibMassWindow.get()))+"\n")
                fw.write("CALIB_S_N_CUTOFF\t"+str(int(self.calibSn.get()))+"\n")
                fw.write("CALIB_MIN_PEAK\t"+str(int(self.calibMin.get()))+"\n")
                fw.write("SUM_SPECTRUM_RESOLUTION\t"+str(int(self.sumSpec.get()))+"\n")
                fw.write("MASS_WINDOW\t"+str(float(self.extracMassWindow.get()))+"\n")
                fw.write("TIME_WINDOW\t"+str(float(self.extracTimeWindow.get()))+"\n")
                fw.write("MIN_CHARGE\t"+str(int(self.extracMinCharge.get()))+"\n")
                fw.write("MAX_CHARGE\t"+str(int(self.extracMaxCharge.get()))+"\n")
                fw.write("CHARGE_CARRIER\t"+str(CHARGE_CARRIER[0])+"\n")
                fw.write("MIN_TOTAL\t"+str(float(self.extracMinTotal.get()))+"\n")
                fw.write("BACKGROUND_TOTAL\t"+str(int(self.extracBack.get()))+"\n")
                fw.write("S_N_CUTOFF\t"+str(int(self.extracSnCutoff.get()))+"\n")
        
        master.measurementWindow = 1
        top = self.top = Toplevel()
        self.chargeCarrierVar = StringVar()
        self.chargeCarrierVar.set(CHARGE_CARRIER[0])
        options = []
        top.protocol( "WM_DELETE_WINDOW", lambda: close(self))
        self.alignmentLabel = Label(top, text="Alignment parameters", font="bold")
        self.alignmentLabel.grid(row=0, columnspan=2, sticky=W)
        self.alignTimeWindowLabel = Label(top, text="Alignment time window")
        self.alignTimeWindowLabel.grid(row=1, column=0, sticky=W)
        self.alignTimeWindow = Entry(top)
        self.alignTimeWindow.insert(0, ALIGNMENT_TIME_WINDOW)
        self.alignTimeWindow.grid(row=1, column=1, sticky=W)
        self.alignMassWindowLabel = Label(top, text="Alignment m/z window")
        self.alignMassWindowLabel.grid(row=2, column=0, sticky=W)
        self.alignMassWindow = Entry(top)
        self.alignMassWindow.insert(0, ALIGNMENT_MASS_WINDOW)
        self.alignMassWindow.grid(row=2, column=1, sticky=W)
        self.alignSnLabel = Label(top, text="Minimal S/N for alignment")
        self.alignSnLabel.grid(row=3, column=0, sticky=W)
        self.alignSn = Entry(top)
        self.alignSn.insert(0, ALIGNMENT_S_N_CUTOFF)
        self.alignSn.grid(row=3, column=1, sticky=W)
        self.alignMinLabel = Label(top, text="Minimal features for alignment")
        self.alignMinLabel.grid(row=4, column=0, sticky=W)
        self.alignMin = Entry(top)
        self.alignMin.insert(0, ALIGNMENT_MIN_PEAK)
        self.alignMin.grid(row=4, column=1, sticky=W)
        self.calibrationLabel = Label(top, text="Calibration parameters", font="bold")
        self.calibrationLabel.grid(row=5, columnspan=2, sticky=W)
        self.calibMassWindowLabel = Label(top, text="Calibration mass window")
        self.calibMassWindowLabel.grid(row=6, column=0, sticky=W)
        self.calibMassWindow = Entry(top)
        self.calibMassWindow.insert(0, CALIB_MASS_WINDOW)
        self.calibMassWindow.grid(row=6, column=1, sticky=W)
        self.calibSnLabel = Label(top, text="Minimal S/N for calibration")
        self.calibSnLabel.grid(row=7, column=0, sticky=W)
        self.calibSn = Entry(top)
        self.calibSn.insert(0, CALIB_S_N_CUTOFF)
        self.calibSn.grid(row=7, column=1, sticky=W)
        self.calibMinLabel = Label(top, text="Minimal number of calibrants")
        self.calibMinLabel.grid(row=8, column=0, sticky=W)
        self.calibMin = Entry(top)
        self.calibMin.insert(0, CALIB_MIN_PEAK)
        self.calibMin.grid(row=8, column=1, sticky=W)
        self.extractionLabel = Label(top, text="Extraction parameters", font="bold")
        self.extractionLabel.grid(row=9, columnspan=2, sticky=W)
        self.sumSpecLabel = Label(top, text="Data points per 1 m/z")
        self.sumSpecLabel.grid(row=10, column=0, sticky=W)
        self.sumSpec = Entry(top)
        self.sumSpec.insert(0, SUM_SPECTRUM_RESOLUTION)
        self.sumSpec.grid(row=10, column=1, sticky=W)
        self.extracMassWindowLabel = Label(top, text="Extraction m/z window")
        self.extracMassWindowLabel.grid(row=12, column=0, sticky=W)
        self.extracMassWindow = Entry(top)
        self.extracMassWindow.insert(0, MASS_WINDOW)
        self.extracMassWindow.grid(row=12, column=1, sticky=W)
        self.extracTimeWindowLabel = Label(top, text="Extraction time window")
        self.extracTimeWindowLabel.grid(row=13, column=0, sticky=W)
        self.extracTimeWindow = Entry(top)
        self.extracTimeWindow.insert(0, TIME_WINDOW)
        self.extracTimeWindow.grid(row=13, column=1, sticky=W)
        self.extracMinChargeLabel = Label(top, text="Minimum charge state")
        self.extracMinChargeLabel.grid(row=14, column=0, sticky=W)
        self.extracMinCharge = Entry(top)
        self.extracMinCharge.insert(0, MIN_CHARGE)
        self.extracMinCharge.grid(row=14, column=1, sticky=W)
        self.extracMaxChargeLabel = Label(top, text="Maximum charge state")
        self.extracMaxChargeLabel.grid(row=15, column=0, sticky=W)
        self.extracMaxCharge = Entry(top)
        self.extracMaxCharge.insert(0, MAX_CHARGE)
        self.extracMaxCharge.grid(row=15, column=1, sticky=W)
        for i in UNITS:
            if BLOCKS[i]['available_for_charge_carrier'] == 1:
                options.append(i)
        self.chargeCarrierLabel = Label(top, text="Charge carrier")
        self.chargeCarrierLabel.grid(row=16, column=0, sticky=W)
        self.chargeCarrier = OptionMenu(top, self.chargeCarrierVar, *options)
        self.chargeCarrier.grid(row=16, column=1, sticky=W)
        self.extracMinTotalLabel = Label(top, text="Minimum isotopic fraction")
        self.extracMinTotalLabel.grid(row=17, column=0, sticky=W)
        self.extracMinTotal = Entry(top)
        self.extracMinTotal.insert(0, MIN_TOTAL)
        self.extracMinTotal.grid(row=17, column=1, sticky=W)
        self.extracBackLabel = Label(top, text="Background detection window")
        self.extracBackLabel.grid(row=18, column=0, sticky=W)
        self.extracBack = Entry(top)
        self.extracBack.insert(0, BACKGROUND_WINDOW)
        self.extracBack.grid(row=18, column=1, sticky=W)
        self.extracSnCutoffLabel = Label(top, text="Spectra QC S/N cutoff")
        self.extracSnCutoffLabel.grid(row=19, column=0, sticky=W)
        self.extracSnCutoff = Entry(top)
        self.extracSnCutoff.insert(0, S_N_CUTOFF)
        self.extracSnCutoff.grid(row=19,column=1, sticky=W)
        self.ok = Button(top,text = 'Ok', command = lambda: close(self))
        self.ok.grid(row = 20, column = 0, sticky = W)
        self.save = Button(top, text = 'Save', command = lambda: save(self))
        self.save.grid(row = 20, column = 1, sticky = E)
        # Tooltips
        createToolTip(self.alignTimeWindowLabel,"The time window in seconds around the specified time of an alignment feature that LaCyTools is allowed to look for the maximum intensity of each feature.")
        createToolTip(self.alignMassWindowLabel,"The m/z window in Thompson around the specified exact m/z of an alignment feature, that LaCyTools will use to find the maximum of each feature.")
        createToolTip(self.alignSnLabel,"The minimum S/N of an alignment feature to be included in the alignment.")
        createToolTip(self.alignMinLabel,"The minimum number of features that have a S/N higher than the minimum S/N for alignment to occur.")
        createToolTip(self.calibMassWindowLabel,"The mass window in Dalton around the specified exact m/z of a calibrant, that LaCyTools uses to determine the uncalibrated accurate mass. This value will be charge state corrected, i.e. for a triple charged analyte the used window will be the value specified here divided by 3.")
        createToolTip(self.calibSnLabel,"The minimum S/N of a calibrant to be included in the calibration.")
        createToolTip(self.calibMinLabel,"The minimum number of calibrants that have a S/N higher than the minimum S/N for calibration to occur.")
        createToolTip(self.sumSpecLabel,"The number of bins per m/z that will be used in the sum spectrum. A value of 100 means that each data point in the sum spectrum is spaced at 0.01 m/z.")
        createToolTip(self.extracMassWindowLabel,"The m/z window in Thompson around the specified exact m/z of a feature that LaCyTools will use for quantitation. For example, a value of 0.1 results in LaCyTools quantifying 999.9 to 1000.1 for a feature with an m/z value of 1000.")
        createToolTip(self.extracTimeWindowLabel,"The rt window in seconds around the specified elution time of each cluster that contains features for quantitation. For example, a value of 10 will result in LaCyTools creating a sum spectrum from 90 s. to 110 s. for a cluster eluting at 100s.")
        createToolTip(self.extracMinChargeLabel,"The minimum charge state that LaCyTools will attempt to use in calibration and quantitation for all features listed in the analyte reference file.")
        createToolTip(self.extracMaxChargeLabel,"The maximum charge state that LaCyTools will attempt to use in calibration and quantitation for all features listed in the analyte reference file.")
        createToolTip(self.chargeCarrierLabel,"The charge carrier that is applied to all specified analytes for quantitation.")
        createToolTip(self.extracMinTotalLabel,"The minimum fraction of the theoretical isotopic pattern that LaCyTools will use for quantitation. For example, a value of 0.95 means that LaCyTools will quantify isotopes until the sum of the quantified isotopes exceeds 0.95 of the total theoretcal isotopic pattern.")
        createToolTip(self.extracBackLabel,"The mass window in Dalton that LaCyTools is allowed to look for the local background and noise for each analyte. For example, a value of 10 means that LaCyTools will look from 990 m/z to 1010 m/z for an analyte with an m/z of 1000.")
        createToolTip(self.extracSnCutoffLabel,"The minimum S/N of an analyte to be included in the spectral QC. Specifically, for the output that lists what fraction of the total quantified analytes passed the here specified S/N value.")

    def getSettings(self):
        """ This function reads the settings file as specified in the
        program, applying them to the program.
        """
        with open(SETTINGS_FILE,'r') as fr:
            for line in fr:
                line = line.rstrip('\n')
                chunks = line.split()
                if chunks[0] == "ALIGNMENT_TIME_WINDOW":
                    global ALIGNMENT_TIME_WINDOW
                    ALIGNMENT_TIME_WINDOW = float(chunks[1])
                if chunks[0] == "ALIGNMENT_MASS_WINDOW":
                    global ALIGNMENT_MASS_WINDOW
                    ALIGNMENT_MASS_WINDOW = float(chunks[1])
                if chunks[0] == "ALIGNMENT_S_N_CUTOFF":
                    global ALIGNMENT_S_N_CUTOFF
                    ALIGNMENT_S_N_CUTOFF = int(chunks[1])
                if chunks[0] == "ALIGNMENT_MIN_PEAK":
                    global ALIGNMENT_MIN_PEAK
                    ALIGNMENT_MIN_PEAK = int(chunks[1])
                if chunks[0] == "CALIB_MASS_WINDOW":
                    global CALIB_MASS_WINDOW
                    CALIB_MASS_WINDOW = float(chunks[1])
                if chunks[0] == "CALIB_S_N_CUTOFF":
                    global CALIB_S_N_CUTOFF
                    CALIB_S_N_CUTOFF = int(chunks[1])
                if chunks[0] == "CALIB_MIN_PEAK":
                    global CALIB_MIN_PEAK
                    CALIB_MIN_PEAK = int(chunks[1])
                if chunks[0] == "SUM_SPECTRUM_RESOLUTION":
                    global SUM_SPECTRUM_RESOLUTION
                    SUM_SPECTRUM_RESOLUTION = int(chunks[1])
                if chunks[0] == "MASS_WINDOW":
                    global MASS_WINDOW
                    MASS_WINDOW = float(chunks[1])
                if chunks[0] == "TIME_WINDOW":
                    global TIME_WINDOW
                    TIME_WINDOW = float(chunks[1])
                if chunks[0] == "MIN_CHARGE":
                    global MIN_CHARGE
                    MIN_CHARGE = int(chunks[1])
                if chunks[0] == "MAX_CHARGE":
                    global MAX_CHARGE
                    MAX_CHARGE = int(chunks[1])
                if chunks[0] == "MIN_TOTAL":
                    global MIN_TOTAL
                    MIN_TOTAL = float(chunks[1])
                if chunks[0] == "BACKGROUND_TOTAL":
                    global BACKGROUND_TOTAL
                    BACKGROUND_TOTAL = int(chunks[1])
                if chunks[0] == "S_N_CUTOFF":
                    global S_N_CUTOFF
                    S_N_CUTOFF = int(chunks[1])
                        
    def feature_reader(self,file):
        """ This reads the contents of the alignmen features file and 
        stores the relevant values in a list.
        
        INPUT: A filename
        OUTPUT: A list of m/z,retention lists (elements are type float)
        """
        features = []
        with open(file,'r') as fr:
            for line in fr:
                try:
                    if line[0][0].isdigit():
                        line = line.rstrip().split()
                        features.append(map(float,line))
                except IndexError:
                    print "Incorrect line observed in: "+str(file)
                    if self.log == True:
                        with open('LaCyTools.log', 'a') as flog:
                            flog.write(str(datetime.now())+ "\tIncorrect line observed in: "+str(analyteFile)+"\n")
                except:
                    print "Unexpected Error: ", sys.exc_info()[0]
        return features

    def fitFunc(self, x,a,b,c):
        penalty = 0
        if b > 2.:
            penalty = abs(b-1.)*10000
        if b < 0.:
            penalty = abs(2.-b)*10000
        return a*x**b + c + penalty
        
    def fitFuncLin(self, x,a,b):
        return a*x + b

    def calcQuadratic(self,data,func):
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
            if func == "PowerLaw":
                z = curve_fit(self.fitFunc, observed, expected)#,maxfev=10000)
            elif func == "Linear":
                z = curve_fit(self.fitFuncLin,observed,expected)
            name = self.inputFile.split(".")[0]
            name = os.path.join(self.batchFolder,name)
            #############
            # Plot Code #
            #############
            minX = min(expected)-0.1*min(expected)
            maxX = max(expected)+0.1*max(expected)
            newX = numpy.linspace(minX,maxX,2500*(maxX-minX))
            linY = newX
            if func == "PowerLaw":
                yNew = self.fitFunc(newX,*z[0])
                minY = self.fitFunc(minX,*z[0])
                maxY = self.fitFunc(maxX,*z[0])
            elif func == "Linear":
                yNew = self.fitFuncLin(newX,*z[0])
                minY = self.fitFuncLin(minX,*z[0])
                maxY = self.fitFuncLin(maxX,*z[0])
            fig =  plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            plt.scatter(expected,observed,c='b',label='Raw',alpha=0.5)
            observedCalibrated = []
            for index, j in enumerate(observed):
                if func == "PowerLaw":
                    observedCalibrated.append(self.fitFunc(j,*z[0]))
                elif func == "Linear":  
                    observedCalibrated.append(self.fitFuncLin(j,*z[0])) 
            plt.scatter(expected,observedCalibrated,c='r',label='Calibrated',marker='s',alpha=0.5)
            numbers = ["%.2f" % number for number in z[0]]
            if func == "PowerLaw":
                if float(numbers[2]) > 0.0:
                    plt.plot(newX,yNew,label="Fit, Function: "+str(numbers[0])+"x"+"$^{"+str(numbers[1])+"}$+"+str(numbers[2]),c='b')
                else:
                    plt.plot(newX,yNew,label="Fit, Function: "+str(numbers[0])+"x"+"$^{"+str(numbers[1])+"}$"+str(numbers[2]),c='b')
            elif func == "Linear":
                if float(numbers[1]) > 0.0:
                    plt.plot(newX,yNew, label="Fit, Function: "+str(numbers[0])+"x+"+str(numbers[1]),c='b')
                else:
                    plt.plot(newX,yNew, label="Fit, Function: "+str(numbers[0])+"x"+str(numbers[1]),c='b')
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
        """ This function creates a popup window belonging to the HD5
        data format. The window has a button where the user has to select
        the location of his mzXML files, a checkbox indicating if the
        mzXML files can be deleted afterwards and lastly a run button.
        
        INPUT: None
        OUTPUT: None
        """
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
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tFinished converting\n")
        end_time = time.time()
        print "Batch conversion lasted for", (end_time - start_time) / 60., "minutes, or", (end_time - start_time) / len(filenames), "seconds per sample."
        tkMessageBox.showinfo("Status Message","Batch Convert finished on "+str(datetime.now()))

    def batchProcess(self,master):
        """ This is the main controller function for batch processing.
        First, the function checks if any reference or alignment file 
        was selected, producing a message box if this is not the case.
        Afterwards, the function checks whether or not it has to read
        from HD5 files or other accepted file formats. Subsequently,
        it performs alignment if an alignment file is selected, followed
        by quantitation (and calibration) if a reference list is
        selected. Finally, it will combine all the individual results
        into a summary file before cleaning up.
        
        INPUT: None
        OUTPUT: A summary file
        """
        import time
        start = time.time()
        # Safety feature (prevents batchProcess from being started multiple times)
        if self.batchProcessing == 1:
            tkMessageBox.showinfo("Error Message", "Batch Process already running")
            return
        self.batchProcessing = 1
        #####################
        # PROGRESS BAR CODE #
        #####################
        self.alPerc = StringVar()
        self.extPerc = StringVar()
        self.alPerc.set("0%")
        self.extPerc.set("0%")
        # barWindow = Tk()
        barWindow = self.top = Toplevel()
        barWindow.title("Progress Bar")
        al = Label(barWindow, text="Alignment", padx=25)
        al.grid(row=0, column=0, sticky="W")
        ft = ttk.Frame(barWindow)
        ft.grid(row=1, columnspan=2, sticky="")
        perc1 = Label(barWindow, textvariable=self.alPerc)
        perc1.grid(row=0, column=1, padx=25)
        progressbar = ttk.Progressbar(ft, length=100, mode='determinate')
        progressbar.grid(row=1, columnspan=2, sticky="")
        ext = Label(barWindow, text="Quantitation", padx=25)
        ext.grid(row=2, column=0, sticky="W")
        ft2 = ttk.Frame(barWindow)
        ft2.grid(row=3, columnspan=2, sticky="")
        perc2 = Label(barWindow, textvariable=self.extPerc)
        perc2.grid(row=2, column=1, padx=25)
        progressbar2 = ttk.Progressbar(ft2, length=100, mode='determinate')
        progressbar2.grid(row=3, columnspan=2, sticky="")
        ###################
        # END OF BAR CODE #
        ###################
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
            for index,file in enumerate(filenames):
                self.alPerc.set(str(int( (float(index) / float(len(filenames) ) ) *100))+"%")
                progressbar["value"] = int( (float(index) / float(len(filenames) ) ) *100)
                progressbar.update()
                array = []
                timePairs = []
                self.inputFile = file
                self.inputFileIdx = filenames2idx[file]
                readTimes = self.matchFeatureTimes(features)
                self.readData(array,readTimes)
                strippedFeatures = []
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
                    directionFlag = 0
                    for k in range(0,len(sortedData)-(startSize+1)):
                        if sortedData[currSize+1] < currAverage + 3 * currNoise:
                            directionFlag == 1
                            currSize += 1
                            currAverage =  numpy.average(sortedData[0:currSize])
                            if self.noise == "MM":
                                currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                            elif self.noise == "RMS":
                                currNoise = numpy.std(sortedData[0:currSize])
                        else:
                            if sortedData[currSize-1] > currAverage + 3 * currNoise and directionFlag == 0:
                                currSize -= 1
                                currAverage = numpy.average(sortedData[0:currSize])
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
                    if peakIntensity > background + ALIGNMENT_S_N_CUTOFF * noise:
                        timePairs.append((i[1],peakTime))
                        strippedFeatures.append(i)
                    else:
                        if self.log == True:
                            with open('LaCyTools.log', 'a') as flog:
                                flog.write(str(datetime.now())+"\tFeature: "+str(i)+" was not above alignment S/N cutoff: "+str(ALIGNMENT_S_N_CUTOFF)+" in file: "+str(file)+"\n")
                # Make sure that enough features are used for alignment
                if len(timePairs) >= ALIGNMENT_MIN_PEAK:
                    warnUser = False
                    # Attempt advanced alignment (PowerLaw)
                    alignFunction = self.calcQuadratic(timePairs,"PowerLaw")
                    # Fall back to basic alignment (Linear)
                    if alignFunction == None:
                        if self.log == True:
                            with open('LaCyTools.log', 'a') as flog:
                                flog.write(str(datetime.now())+ "\tAdvanced alignment failed on file: "+str(file)+", switching to basic alignment\n")
                        alignFunction = self.calcQuadratic(timePairs,"Linear")
                        if alignFunction == None:
                            if self.log == True:
                                with open('LaCyTools.log', 'a') as flog:
                                    flog.write(str(datetime.now())+"\tFile: "+str(file)+" could not be aligned. Both advanced and basic alignment fits failed\n")
                            outFile = os.path.split(file)[-1]
                            outFile = "unaligned_"+outFile
                            outFile = os.path.join(self.batchFolder,outFile)
                            open(outFile,'w').close()
                            continue
                    # Bind correct fit function to fit (PowerLaw or Linear)
                    if len(alignFunction[0]) == 3:
                        fit = self.fitFunc
                    elif len(alignFunction[0]) == 2:
                        fit = self.fitFuncLin                               
                    # Create alignment output file
                    alignmentOutput = self.inputFile.split(".")[0]
                    alignmentOutput = alignmentOutput + ".alignment"
                    with open(alignmentOutput,'w') as falign:
                        lsq = 0
                        falign.write("Peak\tExpected RT\tOriginal RT\tAligned RT\n")
                        for index,timePair in enumerate(timePairs):
                            falign.write(str(strippedFeatures[index][0])+"\t"+str(timePair[0])+"\t"+str(timePair[1])+"\t"+str(fit(float(timePair[1]),*alignFunction[0]))+"\n")
                            lsq += float(strippedFeatures[index][0]) - fit(float(timePair[1]),*alignFunction[0])
                    self.transform_mzXML(file,fit,alignFunction[0])
                else:
                    if self.log == True:
                        with open('LaCyTools.log', 'a') as flog:
                            flog.write(str(datetime.now())+ "\tFile not aligned due to lack of features\n")
                    outFile = os.path.split(file)[-1]
                    outFile = "unaligned_"+outFile
                    outFile = os.path.join(self.batchFolder,outFile)
                    open(outFile,'w').close()
        self.alPerc.set("100%")
        progressbar["value"] = 100
        # (CALIBRATION AND) EXTRACTION
        if self.refFile != "":
            if self.analyteIntensity.get() == 0 and self.analyteRelIntensity.get() == 0 and self.analyteBackground.get() == 0 and self.analyteNoise.get() == 0 and self.alignmentQC.get() == 0 and self.qualityControl.get() == 0 and self.ppmQC.get() == 0 and self.SN.get() == 0 and self.spectraQualityControl.get() == 0:
                tkMessageBox.showinfo("Output Error","No outputs selected")
            self.initCompositionMasses(self.refFile)
            ref = []
            self.refParser(ref)
            times = []
            for i in ref:
                times.append((i[4],i[5]))
            chunks = collections.OrderedDict()
            for i in times:
                if i not in chunks.keys():
                    chunks['%s' % '-'.join(i)] = []
            for i in ref:
                chunks['%s' % '-'.join((i[4],i[5]))].append(i)
            if os.path.isfile(os.path.join(self.batchFolder,"pytables.h5")) == False:
                filenames = glob.glob(os.path.join(str(self.batchFolder),EXTRACTION+"*"+EXTENSION))
                filenames2idx = dict([(filename, idx) for idx, filename in enumerate(filenames)])
            for index,file in enumerate(filenames):
                self.extPerc.set(str(int( (float(index) / float(len(filenames) ) ) *100))+"%")
                progressbar2["value"] = int( (float(index) / float(len(filenames) ) ) *100)
                progressbar2.update()
                results = []
                self.inputFile = file
                self.inputFileIdx = filenames2idx[file]
                array = []
                readTimes = self.matchAnalyteTimes(ref)
                self.readData(array, readTimes)
                for index,i in enumerate(chunks.keys()):
                    spectrum = self.sumSpectrum(i,array)
                    # Dirty hack to now get rid of the time window again
                    rt = tuple(i.split('-'))[0]
                    calibrants = []
                    # Calibrate the sum spectrum
                    if self.calFile.get() == 1:
                        for j in ref:
                            if j[6] == "True" and int(round(float(j[4]))) == int(round(float(rt))):
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
                            if self.log == True:
                                with open('LaCyTools.log', 'a') as flog:
                                    flog.write(str(datetime.now())+ "\tUnable to calibrate the sum spectrum at "+str(i)+" seconds\n")
                            # Adjust filename
                            (old, new) = os.path.split(self.inputFile)
                            old = os.path.abspath(old)
                            new = os.path.splitext(new)[0]
                            new = "Uncalibrated_sumSpectrum_"+str(i)+"_"+str(new)+".xy"
                            new = os.path.join(old,new)
                            outFile = "\\\\?\\"+new
                            # Write
                            with open(outFile,'w') as fw:
                                fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                            continue
                        f = numpy.poly1d(z)
                        calOut = str(file.split(".")[0])+"_"+str(index)+".calibration"
                        with open(calOut,'w') as fw2:
                            for index,j in enumerate(measuredMaximaMZ):
                                fw2.write("accurate mass: "+str(presentCalibrants[index])+" measured at "+str(j) +" being calibrated to: "+str(f(j))+"\n")
                        mzList = []
                        intList = []
                        for j in spectrum:
                            mzList.append(float(j[0]))
                            intList.append(int(j[1]))
                        # Additional batches still work here
                        # Transform python list into numpy array
                        mzArray = numpy.array(mzList)
                        newArray = f(mzArray)
                        newSpectrum = []
                        for index,j in enumerate(newArray):
                            newSpectrum.append((j,intList[index]))
                        spectrum = newSpectrum
                        # Adjust filename
                        (old, new) = os.path.split(self.inputFile)
                        old = os.path.abspath(old)
                        new = os.path.splitext(new)[0]
                        new = "sumSpectrum_"+str(i)+"_"+str(new)+".xy"
                        new = os.path.join(old,new)
                        outFile = "\\\\?\\"+new
                        # Write
                        with open(outFile,'w') as fw:
                            fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                    else:
                        # Adjust filename
                        (old, new) = os.path.split(self.inputFile)
                        old = os.path.abspath(old)
                        new = os.path.splitext(new)[0]
                        new = "sumSpectrum_"+str(i)+"_"+str(new)+".xy"
                        new = os.path.join(old,new)
                        outFile = "\\\\?\\"+new
                        # Write
                        with open(outFile,'w') as fw:
                            fw.write("\n".join(str(j[0])+"\t"+str(j[1]) for j in spectrum))
                    self.extractData(chunks[i],spectrum,results)
                self.writeResults(results,file)
            # Wrap up stuff
            self.extPerc.set("100%")
            progressbar2["value"] = 100
            barWindow.destroy()
            self.combineResults()
        if self.ptFile is not None:
            self.ptFile.close()
        self.batchProcessing = 0
        end = time.time()
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tBatch process lasted for "+str((end - start) / 60.)+"minutes\n")
        tkMessageBox.showinfo("Status Message","Batch Process finished on "+str(datetime.now()))

    def writeCalibration(self,function,array):
        """ This function creates a calibrated mzXML file. However, the
        function is currently not being used and might be removed in the
        future.
        
        INPUT: Calibration function and the raw data in an array
        OUTPUT: A calibrated mzXML file
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
        """ This function takes a list of potential calibrants and will
        identify the m/z value that shows the maximum intensity. The
        function will determine the accurate mass from a interpolated 
        univariate spline that is fitted through the data points, 
        yielding improved post calibration mass accuracy.
        
        INPUT: A spectrum and a list of features (mass,charge)
        OUTPUT: A containing (accurate mass, intensity) tuples for the
                calibrants that passed the user specified S/N cutoff.
        """
        maxima = []
        for i in features:
            mass, charge = i
            window = CALIB_MASS_WINDOW / charge
            lowMz = self.binarySearch(spectrum,float(mass)-float(window),len(spectrum)-1,'left')
            highMz = self.binarySearch(spectrum,float(mass)+float(window),len(spectrum)-1,'right')
            x_points = []
            y_points = []
            for j in spectrum[lowMz:highMz]:
                x_points.append(j[0])
                y_points.append(j[1])
            newX = numpy.linspace(x_points[0],x_points[-1],2500*(x_points[-1]-x_points[0]))
            maximum = (newX[int(len(newX)/2)],0)
            try:
                f = InterpolatedUnivariateSpline(x_points,y_points)
                ySPLINE = f(newX)
                for index, j in enumerate(ySPLINE):
                    if j > maximum[1]:
                        maximum = (newX[index],j)
            except ValueError:
                data = zip(x_points,y_points)
                if self.log == True:
                    with open('LaCyTools.log', 'a') as flog:
                        flog.write(str(datetime.now())+ "\tGuassian Curve Fit failed for analyte: "+str(i[0])+", reverting to non fitted local maximum\n")
                for j in data:
                    if j[1] > maximum[1]:
                        maximum = (j[0],j[1])
            except:
                print "Analyte: "+str(i[0])+" is being troublesome, kill it"
            # Plot Code (for testing purposes)
            """fig =  plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(x_points, y_points, 'b*')
            #plt.plot(newX,newY, 'b--')
            #plt.plot(newX,ySPLINE,'r--')
            plt.plot(newX,yINTER,'r--')
            plt.plot(newX,yUNIVAR,'g--')
            #plt.legend(['Raw Data','Guassian (All Points)','Cubic Spline'], loc='best')
            #plt.legend(['Raw Data','Cubic Spline'], loc='best')
            plt.legend(['Raw Data','Interp1d','Univariate Spline'], loc='best')
            plt.show()"""
            # Check if maxima above S/N cut-off
            values = self.getBackground(spectrum, maximum[0], charge, window)
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
        
        INPUT: The retention time-time window and an array containing
               the entire measurement
        OUTPUT: A sum spectrum in array form (m/z, intensity)
        """

        time = tuple(time.split('-'))
        # This is returning None's now
        lowTime = self.binarySearch(array,float(time[0])-float(time[1]),len(array)-1,'left')
        highTime = self.binarySearch(array,float(time[0])+float(time[1]),len(array)-1,'right')
        LOW_MZ = 25000.0
        HIGH_MZ = 0.0
        for i in array[lowTime:highTime]:
            if i[1][0][0] < LOW_MZ:
                LOW_MZ = i[1][0][0]
            if i[1][-1][0] > HIGH_MZ:
                HIGH_MZ = i[1][-1][0]
        # This should be dynamically determined
        arraySize = (float(HIGH_MZ) - float(LOW_MZ)) * float(SUM_SPECTRUM_RESOLUTION)
        combinedSpectra = numpy.zeros(shape=(int(arraySize+2),2))
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
        #from scipy.signal import savgol_filter
        #new = savgol_filter(combinedSpectra,21,3)
        #return new
        return combinedSpectra

    def findNearest(self,array,value):
        """ A depracated function, will most likely be removed in the
        near future.
        """
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

    def transform_mzXML(self,file,fit,alignFunction):
        """Reads the mzXML file and transforms the reported retention
        time by the specified polynomial function.
        
        INPUT: A filename, alignment function and fitting model
        OUTPUT: An aligned mzXML file
        """
        with open(file,'r') as fr:
            outFile = os.path.split(file)[-1]
            # Use * to indicate files that were aligned using the basic alignment
            if len(alignFunction) == 3:
                outFile = "aligned_"+outFile
            elif len(alignFunction) == 2:
                outFile = "alignedLin_"+outFile
            outFile = os.path.join(self.batchFolder,outFile)
            if self.log == True:
                with open('LaCyTools.log', 'a') as flog:
                    flog.write(str(datetime.now())+ "\tWriting output file: "+outFile+"\n")
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
                        if  fit(float(time),*alignFunction) < 0:
                            newTime = str(0)
                        else:
                            newTime = str(fit(float(time),*alignFunction))
                        line = line.replace(time,newTime)
                        fw.write(line)
                    else:
                        fw.write(line)

    def alignRTs(self,file,polynomial):
        """Reads the mzXML file and transforms the reported retention
        time by the specified polynomial function.
        
        INPUT: A filename and the alignment function
        OUTPUT: An aligned mzXML file
        """
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tAligning file: "+str(self.inputFile)+"\n")
        i = self.inputFileIdx
        for row in self.ptFile.root.scans.where("sample == i"):
            time = row['rt']
            # The below line is only to make it work with mzMine
            if  self.fitFunc(time,*polynomial) > 0:
                row['art'] = self.fitFunc(time,*polynomial)
            row.update()
        self.ptFile.flush()

    def feature_finder(self,data,lowMass,highMass):
        # Proper fix
        intensity = 0
        start = self.binarySearch(data,lowMass,len(data)-1,'left')
        end = self.binarySearch(data,highMass,len(data)-1,'right')
        for i in data[start:end]:
            if i[1] > intensity:
                intensity = i[1]
        return intensity


    def combineResults(self):
        """ This function reads all the raw files and creates the summary
        output file.
        
        INPUT: None
        OUTPUT: A summary file
        """
        total = []
        ref = []
        self.refParser(ref)
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
            fw.write("LaCyTools Version\t"+str(self.version)+"\n")
            fw.write("LaCyTools Build\t"+str(self.build)+"\n")
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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

            ####################################################
            # Analyte Relative Intensity (Total Normalization) #
            ####################################################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 0 and self.normalizeCluster.get() == 0:
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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

            ######################################################
            # Analyte Relative Intensity (Cluster Normalization) #
            ######################################################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 0 and self.normalizeCluster.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Rel Int (Cluster Norm)")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual Data
                    clusters = []
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    currentCluster = "-".join((k.time, k.timeWindow))
                                    if currentCluster not in clusters:
                                        clusters.append(currentCluster)
                                except AttributeError:
                                    continue
                    for i in total:
                        clusterValues = []
                        fw.write(str(i[0]))
                        for j in clusters:
                            clusterTime = float(j.split("-")[0])
                            clusterWindow = float(j.split("-")[1])
                            totalIntensity = 1
                            for k in compositions:
                                for l in i[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == clusterTime and float(l.timeWindow) == clusterWindow and float(k[1]) == float(l.time):
                                            for m in l.isotopes:
                                                totalIntensity += m.obsInt
                                    except AttributeError:
                                        pass
                            clusterValues.append((clusterTime, clusterWindow, totalIntensity))
                        for j in compositions:
                            flag = 0
                            sumInt = 0
                            for k in i[1]:
                                for l in clusterValues:
                                    try:
                                        if k.composition == j[0] and float(k.time) == l[0] and float(k.timeWindow) == l[1] and float(j[1]) == float(k.time):
                                            flag = 1
                                            for m in k.isotopes:                                             
                                                sumInt += m.obsInt
                                            if sumInt > 0:
                                                fw.write("\t"+str(float(sumInt)/float(l[2])))
                                            else:
                                                fw.write("\t")
                                    except AttributeError:
                                            pass
                            if flag == 0:
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
                        fw.write("Rel Int (Cluster Norm, "+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        clusters = []
                        for j in total:
                            for k in compositions:
                                for l in j[1]:
                                    try:
                                        currentCluster = "-".join((l.time, l.timeWindow))
                                        if currentCluster not in clusters:
                                            clusters.append(currentCluster)
                                    except AttributeError:
                                        continue
                        for j in total:
                            clusterValues = []
                            fw.write(str(j[0]))
                            for k in clusters:
                                clusterTime = float(k.split("-")[0])
                                clusterWindow = float(k.split("-")[1])
                                totalIntensity = 1
                                for l in compositions:
                                    for m in j[1]:
                                        try:
                                            if m.composition == l[0] and float(m.time) == clusterTime and float(m.timeWindow) == clusterWindow and float(l[1]) == float(m.time):
                                                for n in m.isotopes:
                                                    if int(n.charge) == i:
                                                        totalIntensity += n.obsInt
                                        except AttributeError:
                                            pass
                                clusterValues.append((clusterTime, clusterWindow, totalIntensity))   
                            for k in compositions:
                                flag = 0
                                sumInt = 0
                                for l in j[1]:
                                    for m in clusterValues:
                                        try:
                                            if l.composition == k[0] and float(l.time) == m[0] and float(l.timeWindow) == m[1] and float(k[1]) == float(l.time):
                                                flag = 1
                                                for n in l.isotopes:
                                                    if int(n.charge) == i:
                                                        sumInt += n.obsInt
                                                if sumInt > 0:
                                                    fw.write("\t"+str(float(sumInt)/float(m[2])))
                                                else:
                                                    fw.write("\t")
                                        except AttributeError:
                                                pass
                                if flag == 0:
                                    fw.write("\t")
                            fw.write("\n")
                        fw.write("\n")

            ##################################################################
            # Background Subtracted Relative Intensity (Total Normalization) #
            ##################################################################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 1 and self.normalizeCluster.get() == 0:
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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

            #####################################################################
            # Background Subtracted Relative Intensity (Cluster Normalization)  #
            #####################################################################
            if self.analyteRelIntensity.get() == 1 and self.analyteBckSub.get() == 1 and self.normalizeCluster.get() == 1:
                #########################
                # Combined charge state #
                #########################
                if self.analytePerCharge.get() == 0:
                    # Header
                    fw.write("Rel Int (Bck Sub, Cluster Norm)")
                    for i in compositions:
                        fw.write("\t"+str(i[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for i in compositions:
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual Data
                    clusters = []
                    for i in total:
                        for j in compositions:
                            for k in i[1]:
                                try:
                                    currentCluster = "-".join((k.time, k.timeWindow))
                                    if currentCluster not in clusters:
                                        clusters.append(currentCluster)
                                except AttributeError:
                                    continue
                    for i in total:
                        fw.write(str(i[0]))
                        clusterValues = []
                        for j in clusters:
                            clusterTime = float(j.split("-")[0])
                            clusterWindow = float(j.split("-")[1])
                            totalIntensity = 1
                            for k in compositions:
                                for l in i[1]:
                                    try:
                                        if l.composition == k[0] and float(l.time) == clusterTime and float(l.timeWindow) == clusterWindow and float(k[1]) == float(l.time):
                                            for m in l.isotopes:
                                                totalIntensity += max(0, m.obsInt - m.background)
                                    except AttributeError:
                                        pass
                            clusterValues.append((clusterTime, clusterWindow, totalIntensity))
                        for j in compositions:
                            flag = 0
                            sumInt = 0
                            for k in i[1]:
                                for l in clusterValues:
                                    try:
                                        if k.composition == j[0] and float(k.time) == l[0] and float(k.timeWindow) == l[1] and float(j[1]) == float(k.time):
                                            flag = 1
                                            for m in k.isotopes:                                             
                                                sumInt += max(0, m.obsInt - m.background)
                                            if sumInt > 0:
                                                fw.write("\t"+str(float(sumInt)/float(l[2])))
                                            else:
                                                fw.write("\t")
                                    except AttributeError:
                                            pass
                            if flag == 0:
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
                        fw.write("Rel Int (Bck Sub, Cluster Norm, "+str(i)+"+)")
                        for j in compositions:
                            fw.write("\t"+str(j[0]))
                        fw.write("\n")
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
                            fw.write("\t"+masses)
                        fw.write("\n")
                        # Actual data
                        clusters = []
                        for j in total:
                            for k in compositions:
                                for l in j[1]:
                                    try:
                                        currentCluster = "-".join((l.time, l.timeWindow))
                                        if currentCluster not in clusters:
                                            clusters.append(currentCluster)
                                    except AttributeError:
                                        continue
                        for j in total:
                            clusterValues = []
                            fw.write(str(j[0]))
                            for k in clusters:
                                clusterTime = float(k.split("-")[0])
                                clusterWindow = float(k.split("-")[1])
                                totalIntensity = 1
                                for l in compositions:
                                    for m in j[1]:
                                        try:
                                            if m.composition == l[0] and float(m.time) == clusterTime and float(m.timeWindow) == clusterWindow  and float(l[1]) == float(m.time):
                                                for n in m.isotopes:
                                                    if int(n.charge) == i:
                                                        totalIntensity += max(0, n.obsInt - n.background)
                                        except AttributeError:
                                            pass
                                clusterValues.append((clusterTime, clusterWindow, totalIntensity))
                            for k in compositions:
                                flag = 0
                                sumInt = 0
                                for l in j[1]:
                                    for m in clusterValues:
                                        try:
                                            if l.composition == k[0] and float(l.time) == m[0] and float(l.timeWindow) == m[1] and float(k[1]) == float(l.time):
                                                flag = 1
                                                for n in l.isotopes:
                                                    if int(n.charge) == i:
                                                        sumInt += max(0, n.obsInt - n.background)
                                                if sumInt > 0:
                                                    fw.write("\t"+str(float(sumInt)/float(m[2])))
                                                else:
                                                    fw.write("\t")
                                        except AttributeError:
                                                pass
                                if flag == 0:
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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
                        sumInt = 0.      
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]):
                                sumInt += float(j[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for i in compositions:
                        masses = ""
                        for j in ref:
                            glycan = j[0].split("_")[0]
                            charge = j[0].split("_")[1]
                            isotope = j[0].split("_")[2]
                            if str(glycan) == str(i[0]) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(j[1])+"]"
                                else:
                                    masses += " ["+str(j[1])+"]"
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
                        # List of theoretical areas
                        fw.write("Fraction")
                        for j in compositions:
                            sumInt = 0.      
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i):
                                    sumInt += float(k[2])
                            fw.write("\t")
                            if sumInt > 0.:
                                fw.write(str(sumInt))
                        fw.write("\n")
                        # List of monoisotopic masses
                        fw.write("Monoisotopic Mass")
                        for j in compositions:
                            masses = ""
                            for k in ref:
                                glycan = k[0].split("_")[0]
                                charge = k[0].split("_")[1]
                                isotope = k[0].split("_")[2]
                                if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                    if masses == "":
                                        masses ="["+str(k[1])+"]"
                                    else:
                                        masses += " ["+str(k[1])+"]"
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

            ##################################
            # Analyte Mass Accuracy (in PPM) #
            ##################################
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
                    fw.write("Mass Accuracy [ppm] ("+str(i)+"+)")
                    for j in compositions:
                        fw.write("\t"+str(j[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0.      
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i):
                                sumInt += float(k[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(k[1])+"]"
                                else:
                                    masses += " ["+str(k[1])+"]"
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
                    fw.write("IPQ ("+str(i)+"+)")
                    for j in compositions:
                        fw.write("\t"+str(j[0]))
                    fw.write("\n")
                    # List of theoretical areas
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0.      
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i):
                                sumInt += float(k[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(k[1])+"]"
                                else:
                                    masses += " ["+str(k[1])+"]"
                        fw.write("\t"+masses)
                    fw.write("\n")
                    # Actual data
                    for j in total:
                        fw.write(str(j[0]))
                        for k in compositions:
                            sumInt = 0
                            totalExpInt = 0
                            qc = 0
                            for l in j[1]:
                                try:
                                    if l.composition == k[0] and float(l.time) == float(k[1]):
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                sumInt += max(float(m.obsInt) - float(m.background),0)
                                                totalExpInt += float(m.expInt)
                                        for m in l.isotopes:
                                            if int(m.charge) == i:
                                                try:
                                                    maxIntensityBackCorrected = max(float(m.obsInt) - float(m.background),0)
                                                    qc += abs((maxIntensityBackCorrected  / float(sumInt)) - (m.expInt/totalExpInt))
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
                    # List of theoretical areas
                    fw.write("Fraction")
                    for j in compositions:
                        sumInt = 0.      
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i):
                                sumInt += float(k[2])
                        fw.write("\t")
                        if sumInt > 0.:
                            fw.write(str(sumInt))
                    fw.write("\n")
                    # List of monoisotopic masses
                    fw.write("Monoisotopic Mass")
                    for j in compositions:
                        masses = ""
                        for k in ref:
                            glycan = k[0].split("_")[0]
                            charge = k[0].split("_")[1]
                            isotope = k[0].split("_")[2]
                            if str(glycan) == str(j[0]) and int(charge) == int(i) and int(isotope) == 0:
                                if masses == "":
                                    masses ="["+str(k[1])+"]"
                                else:
                                    masses += " ["+str(k[1])+"]"
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
                    

            #########################################
            # Fraction of analytes above S/N cutoff #
            #########################################
            if self.spectraQualityControl.get() == 1:
                ref = []
                self.refParser(ref)
                times = []
                for i in ref:
                    times.append((i[4],i[5]))
                chunks = collections.OrderedDict()
                for i in times:
                    if i not in chunks.keys():
                        chunks['%s' % '-'.join(i)] = []
                for i in ref:
                    chunks['%s' % '-'.join((i[4],i[5]))].append(i)
                # Header
                fw.write("Fraction of analytes above S/N Cutoff")
                for index,i in enumerate(chunks.keys()):
                    fw.write("\t"+str(i))
                fw.write("\n")
                #for index,i in enumerate(chunks.keys()):    
                # Actual data
                for j in total:
                    fw.write(str(j[0]))
                    for index,i in enumerate(chunks.keys()):  
                        numberTotal = 0
                        numberPass = 0
                        for k in compositions:
                            expInt = 0
                            SN = 0
                            for l in j[1]:
                                try:
                                    if l.composition == k[0] and float(l.time) == float(i.split("-")[0]) and float(l.timeWindow) == float(i.split("-")[1]):
                                        numberTotal += 1
                                        for m in l.isotopes:
                                            if m.expInt > expInt: # and int(m.charge) == i:
                                                try:
                                                    SN = (m.obsMax - m.backgroundPoint) / m.noise
                                                except ZeroDivisionError:
                                                    pass
                                                expInt = m.expInt
                                        if SN > S_N_CUTOFF:
                                            numberPass += 1
                                except AttributeError:
                                    pass
                        if numberTotal > 0:
                            fw.write("\t"+str(float(numberPass)/float(numberTotal)))
                        else:
                            fw.write("\t")
                    fw.write("\n")
                fw.write("\n")

    def writeResults(self,results,file):
        """ This function writes the resultes per file away to a raw
        file.
        
        INPUT: A file name and a list of results
        OUTPUT: A raw file per measurement
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
        """ This function creates a pop up box in which all the parameters
        for a batch process can be set and visualized. The window can
        access and set the masters alFile, refFile and batchFolder.
        The window can also call the outputPopup function (to specify 
        the contents of final summary) and start the actual
        batchProcess function.

        INPUT: None
        OUTPUT: None
        """
        if master.batchWindow == 1:
            return
        master.batchWindow = 1
        self.al = StringVar()
        self.ref = StringVar()
        self.folder = StringVar()

        if master.alFile:
            self.al.set(master.alFile)

        if master.refFile:
            self.ref.set(master.refFile)

        if master.batchFolder:
            self.folder.set(master.batchFolder)

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

        def run():
            master.batchWindow = 0
            top.destroy()
            master.batchProcess(master)

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
        self.run = Button(top, text = "Run Batch Process", width = 25, command = lambda: run())
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
            master.normalizeCluster.set(1)
            master.alignmentQC.set(1)
            master.qualityControl.set(1)
            master.ppmQC.set(1)
            master.SN.set(1)
            master.spectraQualityControl.set(1)
        def select_none(self):
            master.analyteIntensity.set(0)
            master.analyteRelIntensity.set(0)
            master.analyteBackground.set(0)
            master.analyteNoise.set(0)
            master.analytePerCharge.set(0)
            master.analyteBckSub.set(0)
            master.normalizeCluster.set(0)
            master.alignmentQC.set(0)
            master.qualityControl.set(0)
            master.ppmQC.set(0)
            master.SN.set(0)
            master.spectraQualityControl.set(0)
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
        self.ri = Checkbutton(top, text = u"Relative Intensity\u00B9\u00B7\u00B2\u00B7\u00B3", variable = master.analyteRelIntensity, onvalue = 1, offvalue = 0)
        self.ri.grid(row = 3, column = 0, sticky = W)
        self.back = Checkbutton(top, text = u"Analyte Background\u00B9", variable = master.analyteBackground, onvalue = 1, offvalue = 0)
        self.back.grid(row = 4, column = 0, sticky = W)
        self.analNoise = Checkbutton(top, text = u"Analyte Noise\u00B9", variable = master.analyteNoise, onvalue = 1, offvalue = 0)
        self.analNoise.grid(row = 5, column = 0, sticky = W)
        self.chargeState = Checkbutton(top, text = u"\u00B9Intensities per Charge State", variable = master.analytePerCharge, onvalue = 1, offvalue = 0)
        self.chargeState.grid(row = 2, column = 1, sticky = W)
        self.bckSub = Checkbutton(top, text = u"\u00B2Background subtracted Intensities", variable = master.analyteBckSub, onvalue = 1, offvalue = 0)
        self.bckSub.grid(row = 3, column = 1, sticky = W)
        self.norClus = Checkbutton(top, text = u"\u00B3Normalization per cluster", variable = master.normalizeCluster, onvalue = 1, offvalue = 0)
        self.norClus.grid(row = 4, column = 1, sticky = W)
        self.align = Checkbutton(top, text="Alignment Residuals", variable=master.alignmentQC, onvalue=1, offvalue=0)
        self.align.grid(row = 6, column=0, sticky=W)
        self.qc = Checkbutton(top, text = "Isotopic Pattern Quality", variable = master.qualityControl, onvalue = 1, offvalue = 0)
        self.qc.grid(row = 7, column = 0, sticky = W)
        self.ppm = Checkbutton(top, text = "Mass Accuracy [ppm]", variable = master.ppmQC, onvalue = 1, offvalue = 0)
        self.ppm.grid(row = 8, column = 0, sticky = W)
        self.snratio = Checkbutton(top, text = "Signal-to-Noise Ratio", variable = master.SN, onvalue = 1, offvalue = 0)
        self.snratio.grid(row = 9, column = 0, sticky = W)
        self.specQC = Checkbutton(top, text="Spectral QC", variable = master.spectraQualityControl, onvalue=1, offvalue=0)
        self.specQC.grid(row = 10, column = 0, sticky = W)
        self.button = Button(top,text='Ok',command = lambda: close(self))
        self.button.grid(row = 11, column = 0, columnspan = 2)
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
        """ This function opens a Tkinter filedialog, asking the user
        to select a file. The chosen file is then set to the
        self.calFile variable.
        
        INPUT: None
        OUTPUT: None
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
        """ This function opens a Tkinter filedialog, asking the user
        to select a file. The chosen file is then set to the
        self.refFile variable.
        
        INPUT: None
        OUTPUT: None
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'refFile',file_path)

    def openAlFile(self):
        """ This function opens a Tkinter filedialog, asking the user
        to select a file. The chosen file is then set to the
        self.alFile variable.
        
        INPUT: None
        OUTPUT: None
        """
        file_path = tkFileDialog.askopenfilename()
        if not file_path:
            pass
        else:
            setattr(self,'alFile',file_path)

    def processBlock(self, block, array, readTimes):
        """ This function processes a data block as taken from the input
        file.
        
        INPUT: A data block from the mzXML file
        OUTPUT: None
        """
        #if "scan num" in block:
        #    scan = block.split("scan num")[1]
        #    scan = scan.split("\"")[1]

        if "retentionTime" in block:
            rt = block.split("retentionTime")[1]
            rt = rt.split("\"")[1]
            if rt[0] == 'P':
                rt = rt[2:-1]

        #if "peaksCount" in block:
        #    peaks = block.split("peaksCount")[1]

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
        # Calculate the bass composition values\
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
        #   TODO: MAKE THIS CHECK EXISTING ANALYTE FILE IF IT MATCHES SUPPLIED REF FILE
        #   print "USING EXISTING REFERENCE FILE"
        #   return
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tPRE-PROCESSING REFERENCE FILE\n")
        with open(analyteFile,'w') as fw:
            fw.write("# Peak\tm/z\tRel Area\twindow\trt\ttime window\tCalibration\n")
            for i in lines:
                try:
                    if i[0] == "#":
                        continue
                except IndexError:
                    if self.log == True:
                        with open('LaCyTools.log', 'a') as flog:
                            flog.write(str(datetime.now())+ "\tIncorrect line observed in: "+str(analyteFile)+"\n")
                except:
                    if self.log == True:
                        with open('LaCyTools.log', 'a') as flog:
                            flog.write(str(datetime.now())+ "\tUnexpected Error: "+str(sys.exc_info()[0])+"\n")
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
                            fw.write(str(i[0])+"_"+str(j)+"_"+str(index)+"\t"+str((k[0]+j*BLOCKS[CHARGE_CARRIER[0]]['mass'])/j)+"\t"+str(k[1])+"\t"+str(massWindow)+"\t"+str(time)+"\t"+str(timeWindow)+"\tTrue\n")
                        else:
                            fw.write(str(i[0])+"_"+str(j)+"_"+str(index)+"\t"+str((k[0]+j*BLOCKS[CHARGE_CARRIER[0]]['mass'])/j)+"\t"+str(k[1])+"\t"+str(massWindow)+"\t"+str(time)+"\t"+str(timeWindow)+"\tFalse\n")
                        #if k[1] < MIN_CONTRIBUTION:
                        #   break
                        if contribution > MIN_TOTAL:
                            break
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tPRE-PROCESSING COMPLETE\n")

    ####################################################
    # END OF FUNCTIONS RELATED TO PARSING ANALYTE FILE #
    ####################################################

    def extractData(self,ref,array,results):
        """ This is the controller function for quantitation of data.
        
        INPUT: A ref file and input file
        OUTPUT: A list of results consisting of (area, 'background',
                accurate mass, maximum intensity point and the current
                analyte)
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
                            f = InterpolatedUnivariateSpline(x_points,y_points)
                            ySPLINE = f(newX)
                            for index, j in enumerate(ySPLINE):
                                if j > maximum[1]:
                                    maximum = (newX[index],j)
                        except ValueError:
                            data = zip(x_points,y_points)
                            if self.log == True:
                                with open('LaCyTools.log', 'a') as flog:
                                    flog.write(str(datetime.now())+ "\tGuassian Curve Fit failed for analyte: "+str(i[0])+", reverting to non fitted local maximum\n")
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
        if self.batchProcessing == 1:
            return results
        else:
            for i in results:
                print i

    def getBackground(self, array, target, charge, width):
        """ This functin will determine the background and noise for a
        given analyte.
        
        INPUT: The spectrum in array form, the exact m/z of the analyte,
               the charge of the analyte and the m/z window
        OUTPUT: A list of (the average background, the background area
                and the noise)
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
        if self.background == "MIN":
            for i in range(0,(2*BACKGROUND_WINDOW)-4):
                mix = totals[i][1]+totals[i+1][1]+totals[i+2][1]+totals[i+3][1]+totals[i+4][1]
                avgBackground = numpy.average([sum(totals[i][0]),sum(totals[i+1][0]),sum(totals[i+2][0]),sum(totals[i+3][0]),sum(totals[i+4][0])])
                dev = numpy.std(mix)
                avg = numpy.average(mix)
                if avg < backgroundPoint:
                    backgroundPoint = avg
                    backgroundArea = avgBackground
                    if self.noise == "RMS":
                        noise = dev
                    elif self.noise == "MM":
                        noise = max(mix) - min(mix)
        # Find the set of 5 consecutive windows with median average intensity
        elif self.background == "MEDIAN":
            values = []
            for i in range(0, (2*OUTER_BCK_BORDER)-4):
                mix = totals[i][1]+totals[i+1][1]+totals[i+2][1]+totals[i+3][1]+totals[i+4][1]
                avgBackground = numpy.average([sum(totals[i][0]), sum(totals[i+1][0]), sum(totals[i+2][0]), sum(totals[i+3][0]), sum(totals[i+4][0])])
                dev = numpy.std(mix)
                avg = numpy.average(mix)
                if self.noise == "RMS":
                    noise = dev
                elif self.noise == "MM":
                    noise = max(mix) - min(mix)
                values.append((avg, avgBackground, noise))
            sortedValues = sorted(values, key=lambda x: x[0])
            a, b, c = zip(*sortedValues)
            backgroundPoint = a[len(a)//2]
            backgroundArea = b[len(b)//2]
            noise = c[len(c)//2]
        # NOBAN METHOD
        elif self.background == "NOBAN":
            dataPoints = []
            for i in range(0, (2*OUTER_BCK_BORDER)):
                dataPoints.extend(totals[i][1])
            sortedData = sorted(dataPoints)
            startSize = int(0.25 * float(len(sortedData)))
            currSize = startSize
            currAverage = numpy.average(sortedData[0:currSize])
            if self.noise == "MM":
                currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
            elif self.noise == "RMS":
                currNoise = numpy.std(sortedData[0:currSize])
            directionFlag = 0
            for k in range(0,len(sortedData)-(startSize+1)):
                if sortedData[currSize+1] < currAverage + 3 * currNoise:
                    directionFlag == 1
                    currSize += 1
                    currAverage =  numpy.average(sortedData[0:currSize])
                    if self.noise == "MM":
                        currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                    elif self.noise == "RMS":
                        currNoise = numpy.std(sortedData[0:currSize])
                else:
                    if sortedData[currSize-1] > currAverage + 3 * currNoise and directionFlag == 0:
                        currSize -= 1
                        currAverage = numpy.average(sortedData[0:currSize])
                        if self.noise == "MM":
                            currNoise = max(sortedData[0:currSize]) - min(sortedData[0:currSize])
                        elif self.noise == "RMS":
                            currNoise = numpy.std(sortedData[0:currSize])
                    else:
                        break
            # Get Area
            # Get length of window
            windowLength = 0
            for i in range(0, (2*OUTER_BCK_BORDER)):
                if len(totals[i][1]) > windowLength:
                    windowLength = len(totals[i][1])
            # Get spacing in window
            begin = self.search_right(data, lowEdge, len(data))
            end = self.search_left(data, highEdge, len(data))
            for j in data[begin:end]:
                spacing =  (data[end][0] - data[begin][0]) / (end - begin)
            currArea = windowLength * (currAverage * spacing)
            # Assign values to generic names
            backgroundPoint = currAverage
            backgroundArea = currArea
            noise = currNoise
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
        for i in ref:
            if (float(i[4])-float(i[5])) not in wanted:
                wanted.append((float(i[4])-float(i[5]),float(i[4])+float(i[5])))
        return list(self.merge_ranges(wanted))

    def merge_ranges(self,ranges):
        """ 
        Merge overlapping and adjacent ranges and yield the merged ranges
        in order. The argument must be an iterable of pairs (start, stop).

        >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
        [(-1, 7)]
        >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
        [(1, 2), (3, 4), (5, 6)]
        >>> list(merge_ranges([]))
        []
        
        Source = http://stackoverflow.com/questions/24130745/convert-generator-object-to-list-for-debugging
        """
        ranges = iter(sorted(ranges))
        current_start, current_stop = next(ranges)
        for start, stop in ranges:
            if start > current_stop:
                # Gap between segments: output current segment and start a new one.
                yield current_start, current_stop
                current_start, current_stop = start, stop
            else:
                # Segments adjacent or overlapping: merge.
                current_stop = max(current_stop, stop)
        yield current_start, current_stop
 
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
            if self.log == True:
                with open('LaCyTools.log', 'a') as flog:
                    flog.write(str(datetime.now())+ "\tProcessing "+str(self.inputFile)+"\n")
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
        """ TODO by Genadij Razdorov
        """
        if self.log == True:
            with open('LaCyTools.log', 'a') as flog:
                flog.write(str(datetime.now())+ "\tProcessing "+str(self.inputFile)+"\n")
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
        
        INPUT: An empty ref list
        OUTPUT: A filled ref list
        """
        with open(os.path.join(self.batchFolder,"analytes.ref"),'r') as fr:
            for line in fr:
                if line[0] == "#":
                    continue
                parts=line.rstrip('\n').split('\t')
                ref.append(parts)

    def mzXMLDecoder(self, rt, peaks, precision, compression, byteOrder, array):
        """ This function parses the encoded string from an mzXML file.
        The decoded data is finally added to the data array containing
        the entire measurement.
        
        INPUT: An encoded string and the data array containing all of 
                the measuremen that has been processed up to this point
        OUTPUT: A data array containing all of the measuremen that has
                been processed up to this point
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

