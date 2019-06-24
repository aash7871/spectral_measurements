target_list_MDM = ['08403953 + 1849', 'A 575', 'AD 2642', 'AD 4269', 'HSHJ272', 'HSHJ385', 'JC143', 'JS230', 'JS244', 'JS267', 'JS281', 'JS283', 'JS301', 'JS315', 'JS317', 'JS329', 'JS349', 'JS352', 'JS391', 'JS414', 'JS441', 'JS452','JS455', 'JS457', 'JS468', 'JS473', 'JS488', 'JS536', 'JS547', 'JS561', 'JS566', 'JS706', 'JS726', 'KW563', 'KW569']

import matplotlib.pyplot as plt, matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm
from astropy import log
import csv
import glob
from astropy.io import fits
import sys
sys.path.insert(0, '/Users/amandaash/Desktop/Research/scripts/')

import spectra_functions
from PHEW import EqW
import pyspeckit as p

def obs_table(target):
    
    
    if target == 'A 575':
        target_name = 'A575'
    
    elif target == 'AD 4269':
        target_name = 'AD4269'
    
    elif target == 'AD 2642':
        target_name = 'AD2642'
    
    else:
        target_name = target
        
    target_files_single = glob.glob('/Users/amandaash/Desktop/Research/data/CSCU_reductions/*/finals/trim.{0}.fits'.format(target_name))
    target_files_multi = glob.glob('/Users/amandaash/Desktop/Research/data/CSCU_reductions/*/finals/trim.{0}.?.fits'.format(target_name))
    target_files = sorted(target_files_single + target_files_multi)
    
    for filename in target_files:
        
        #first up is the spectral data
        local_date = filename.split('/')[7]
        
        eperc = EqW.measure_equivalent_width(filename,
                                     6550,6576,6560,6566,
                                     1000,xunit="Angstrom",to_plot=True,
                                     filebase="/Users/amandaash/Desktop/Research/plots/figures_6_21_19/{0}_{1}.pdf".format(target,local_date))
        perc16 = eperc[0]
        perc50 = eperc[1]
        perc84 = eperc[2]
  
    
    
        #alright now we're going to get the observational data by reading the .fits files
    
        obs_info = fits.open(filename)

        date = obs_info[0].header['DATE-OBS']
        time = obs_info[0].header['TIME-OBS']
        secz = obs_info[0].header["AIRMASS"]
        HA = obs_info[0].header["HA"]
        
    #last thing we're going to get is the EPIC ID, the RA, and the DEC from the target_epicID file
        
        epic_info = []
        with open('/Users/amandaash/Desktop/Research/observation_info/target_epicID.csv', 'r') as csvFile:
            epic_data = csv.reader(csvFile)
            for row in epic_data:
                epic_info.append(row)
        csvFile.close()

        epic_targets = []

        for n in np.arange(0,(len(epic_info[1:])),1):
            epic_targets.append(epic_info[1:][n][0])
    
        index = epic_targets.index(target)
        EPICID = epic_info[1:][index][1]
        RA = epic_info[1:][index][2]
        DEC = epic_info[1:][index][3]
        
        row = [target, EPICID, RA, DEC, local_date, date, time, HA, secz, perc16, perc50, perc84]
       
        
        with open('/Users/amandaash/Desktop/Research/observation_info/observation_table.csv', 'a') as csvFile:
            wr = csv.writer(csvFile, dialect = 'excel')
            wr.writerow(row)
        csvFile.close()
    
with open('/Users/amandaash/Desktop/Research/observation_info/observation_table.csv', 'w') as csvFile:
    row =['litname', 'EPIC ID', 'RA J2000', 'DEC J2000', 'local date obs.', 'UTC date obs.', 'UTC time obs.', 
          'HA', 'sec(z)', 'EW 16th', 'EW 50th', 'EW8 4th']
    wr = csv.writer(csvFile, dialect = 'excel')
    wr.writerow(row)
csvFile.close()

for t in target_list_MDM:

    obs_table(t)