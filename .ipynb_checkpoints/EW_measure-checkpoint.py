"""
Uses PHEW on target exposures

returns: 
outputs PHEW figures for each target exposure
saves PHEW perc to .csv file for each target exposure
"""

from __future__ import print_function, division, absolute_import

import os

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.io.ascii as at
import astropy.io.fits as fits
import astropy.units as u
from scipy.interpolate import interp1d
from cycler import cycler

import sys
sys.path.insert(0, '/Users/amandaash/Desktop/Research/scripts/')

import spectra_functions
import pyspeckit

from PHEW import EqW
import pyspeckit as p

import csv
import glob as glob

target_list = []
with open('/Users/amandaash/Desktop/Research/data/Target_names.csv', 'r') as csvFile:
    target_data = csv.reader(csvFile)
    for row in target_data:
        target_list.append(row[1])
csvFile.close()

targets = target_list[2:]
days=[111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,204,205,206,207,208,209]

for star in targets:
    for day in days:
        flist_single = '/Users/amandaash/Desktop/Research/data/CSCU_reductions/20180{0}/finals/trim.{1}.fits'.format(day,star)
        if os.path.exists(flist_single) is True:
            flist = glob.glob('/Users/amandaash/Desktop/Research/data/CSCU_reductions/20180{0}/finals/trim.{1}.fits'.format(day,star))
            for n in sorted(flist):
                day = n.split('/')[7]
                target_info1 = n.split('/')[9]
                target_name = n.split('.')[1]
                #exposure = n.split('.')[2]
                eperc = EqW.measure_equivalent_width(str(n),
                                                     6550,6576,6560,6566,
                                                     1000,xunit="Angstrom",to_plot=True,    filebase="/Users/amandaash/Desktop/Research/plots/figures/{0}_{1}.pdf".format(target_name,day))
                
        else:
            pass
    for day in days:
        flist_multi = '/Users/amandaash/Desktop/Research/data/CSCU_reductions/20180{0}/finals/trim.{1}.1.fits'.format(day,star)
        if os.path.exists(flist_multi) is True:
            flist = glob.glob('/Users/amandaash/Desktop/Research/data/CSCU_reductions/20180{0}/finals/trim.{1}.?.fits'.format(day,star))
            for n in sorted(flist):
                day = n.split('/')[7]
                target_info1 = n.split('/')[9]
                target_name = n.split('.')[1]
                exposure = n.split('.')[2]
                eperc = EqW.measure_equivalent_width(str(n),
                                                     6550,6576,6560,6566,
                                                     1000,xunit="Angstrom",to_plot=True,
                                                    filebase="/Users/amandaash/Desktop/Research/plots/figures/{0}_{1}_{2}.pdf".format(target_name,day, exposure))
            else:
                continue
        