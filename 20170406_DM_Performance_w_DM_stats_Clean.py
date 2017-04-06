# coding: utf-8

#------------------------------------------------------
# Script to read in text output from https://eew.geo.berkeley.edu:8443/eew-bk-prod1/dmreview/
# Written by Angela Chung < aichung@berkeley.edu >
# From webpage, output text file. Copy and paste text into document.
# Save text file into one of the following subdirectories: CA, LA, BAY_AREA, ALL_REGIONS'
# Things that need to be changed are directory of text file and region.
# User can choose to output results for DM, E2, ON, and VS by uncommenting respective sections
#------------------------------------------------------

#------------------------------------------------------
# Import modules:
#------------------------------------------------------
from pylab import *
import os
import numpy as np
import obspy
from obspy import *
from obspy.signal import *
from obspy.mseed import *
from obspy.core import *
import fnmatch
from obspy.signal.trigger import classicSTALTA, triggerOnset, recSTALTA
from obspy.signal.trigger import plotTrigger
from datetime import date, datetime
import time
import pylab as P
import matplotlib.colors as mcolors
#import matplotlib.mlab as mlab
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.cm as cm
import pandas as pd

import glob
from collections import defaultdict
#
#
#
close('all')


####### PUT DM REVIEW OUTPUT FILES IN DIRECTORY HERE:

os.chdir('/Users/Angie/Documents/BSL/RESEARCH_REPORTS/EVENT_PERFORMANCE_WRITEUPS/20170406_PERFORMANCE/')
cwd = os.getcwd()

####### CHOOSE REGION FROM ONE OF THE FOLLOWING: 'CA', 'LA', 'BAY_AREA', 'ALL_REGIONS'
d_dir = 'CA'

os.chdir('%s/%s' % (cwd, d_dir))
filelist = glob.glob('*.txt')

p = 0
for filename in filelist:
    print('*************************************')
    
    headerlist = ['anss_origin','anss_lat','anss_lon','anss_depth','anss_mag','DM_id','DM_ver','DM_origin',\
                  'DM_alert','DM_lat','DM_lon','DM_depth','DM_mag','E2_id','E2_ver','E2_origin','E2_alert',\
                  'E2_lat','E2_lon','E2_depth','E2_mag','ON_id','ON_ver','ON_origin','ON_alert','ON_lat',\
                  'ON_lon','ON_depth','ON_mag','VS_id','VS_ver','VS_origin','VS_alert','VS_lat','VS_lon',\
                  'VS_depth','VS_mag']
                  
    print(len(headerlist))
    
    title_text = pd.read_csv('%s/%s/%s' % (cwd, d_dir, filename), nrows=1)
    title_text = list(title_text)
    date_range = title_text[0][-82:-58]
    mag = title_text[0][-56:-53]
    mag = float(mag)


    title_text = title_text[0][0:39]


    
    DM_cat = pd.read_csv('%s/%s/%s' % (cwd, d_dir, filename), delim_whitespace=True, skiprows=3, names=headerlist)
    DM_cat = DM_cat[DM_cat.anss_origin != '-']
    DM_cat = DM_cat.replace(to_replace='-',value='0')
    
    
    DM_cat.anss_lat = DM_cat.anss_lat.astype(float)
    DM_cat.anss_lon = DM_cat.anss_lon.astype(float)
    DM_cat.anss_depth = DM_cat.anss_depth.astype(float)
    DM_cat.anss_mag = DM_cat.anss_mag.astype(float)
    DM_cat.DM_lat = DM_cat.DM_lat.astype(float)
    DM_cat.DM_lon = DM_cat.DM_lon.astype(float)
    DM_cat.DM_depth = DM_cat.DM_depth.astype(float)
    DM_cat.DM_mag = DM_cat.DM_mag.astype(float)
    DM_cat.E2_lat = DM_cat.E2_lat.astype(float)
    DM_cat.E2_lon = DM_cat.E2_lon.astype(float)
    DM_cat.E2_depth = DM_cat.E2_depth.astype(float)
    DM_cat.E2_mag = DM_cat.E2_mag.astype(float)
    DM_cat.ON_lat = DM_cat.ON_lat.astype(float)
    DM_cat.ON_lon = DM_cat.ON_lon.astype(float)
    DM_cat.ON_depth = DM_cat.ON_depth.astype(float)
    DM_cat.ON_mag = DM_cat.ON_mag.astype(float)
    DM_cat.VS_lat = DM_cat.VS_lat.astype(float)
    DM_cat.VS_lon = DM_cat.VS_lon.astype(float)
    DM_cat.VS_depth = DM_cat.VS_depth.astype(float)
    DM_cat.VS_mag = DM_cat.VS_mag.astype(float)
    
    DM_cat = DM_cat.reset_index()


    
    
    
    
    
    
        #print(filelist[p][22:-4])
    E2_min = filelist[p][2]
    ON_min = filelist[p][2]
    VS_min = filelist[p][2]
    DM_dist_err = []
    E2_dist_err = []
    ON_dist_err = []
    VS_dist_err = []
    DM_mag_err = []
    E2_mag_err = []
    ON_mag_err = []
    VS_mag_err = []
    DM_abs_mag_err = []
    E2_abs_mag_err = []
    ON_abs_mag_err = []
    VS_abs_mag_err = []
    DM_t_diff = []
    E2_t_diff = []
    ON_t_diff = []
    VS_t_diff = []
    DM_lons = []
    E2_lons = []
    ON_lons = []
    VS_lons = []
    DM_lats = []
    E2_lats = []
    ON_lats = []
    VS_lats = []
    E2_mag_array = []
    
    ALL_false = 0
    ALL_missed = 0
    DM_false  = 0
    DM_missed = 0
    E2_false  = 0
    E2_missed = 0
    ON_false  = 0
    ON_missed = 0
    VS_false  = 0
    VS_missed = 0
    E2_ON_false  = 0
    E2_ON_missed = 0
    E2_VS_false  = 0
    E2_VS_missed = 0
    ON_VS_false  = 0
    ON_VS_missed = 0
    event_count = 0
    
    ######## BELOW CODE CALCULATES NUMBER OF MATCHED, MISSED, FALSE EVENTS, AND STATISTICS FOR MATCHED EVENTS:
    
    for i in range(0,len(DM_cat)):
        
        # Filter out events in Northern California Only:
        #if (DM_cat['anss_lat'][i] > (0.8395*(DM_cat['anss_lon'][i])+136.29) or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['DM_lat'][i] > (0.8395*(DM_cat['DM_lon'][i])+136.29)))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['E2_lat'][i] > (0.8395*(DM_cat['E2_lon'][i])+136.29)))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['ON_lat'][i] > (0.8395*(DM_cat['ON_lon'][i])+136.29)))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['VS_lat'][i] > (0.8395*(DM_cat['VS_lon'][i])+136.29)))):
        #if ((DM_cat['anss_lat'][i] < (0.8395*(DM_cat['anss_lon'][i])+136.29) and (DM_cat['anss_lat'][i] >1))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['DM_lat'][i] < (0.8395*(DM_cat['DM_lon'][i])+136.29)) and (DM_cat['DM_lat'][i] >1))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['E2_lat'][i] < (0.8395*(DM_cat['E2_lon'][i])+136.29)) and (DM_cat['E2_lat'][i] >1))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['ON_lat'][i] < (0.8395*(DM_cat['ON_lon'][i])+136.29)) and (DM_cat['ON_lat'][i] >1))\
        #    or ((DM_cat['anss_lat'][i] ==0) and (DM_cat['VS_lat'][i] < (0.8395*(DM_cat['VS_lon'][i])+136.29)) and (DM_cat['VS_lat'][i] >1))):
            
           
            if (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                
                event_count = event_count + 1
                    
                ALL_missed = ALL_missed + 1
                DM_missed = DM_missed + 1
                E2_missed = E2_missed + 1
                ON_missed = ON_missed + 1
                VS_missed = VS_missed + 1
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                
                E2_missed = E2_missed + 1
                DM_missed = DM_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
                    
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                
                E2_missed = E2_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                    
                DM_missed = DM_missed + 1
                ON_missed = ON_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                    
                ON_missed = ON_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                
    
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
                    
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                DM_missed = DM_missed + 1
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
    
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                    
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
    
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                    
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                    
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                    
                DM_missed = DM_missed + 1
                E2_ON_missed = E2_ON_missed + 1
                E2_missed = E2_missed + 1
                ON_missed = ON_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                    
                E2_ON_missed = E2_ON_missed + 1
                E2_missed = E2_missed + 1
                ON_missed = ON_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
                    
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                DM_missed = DM_missed + 1
                E2_VS_missed = E2_VS_missed + 1
                E2_missed = E2_missed + 1
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                E2_VS_missed = E2_VS_missed + 1
                E2_missed = E2_missed + 1
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                    
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
            
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                DM_missed = DM_missed + 1
                ON_VS_missed = ON_VS_missed + 1
                ON_missed = ON_missed + 1
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                event_count = event_count + 1
                    
                ON_VS_missed = ON_VS_missed + 1
                ON_missed = ON_missed + 1
                VS_missed = VS_missed + 1
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                    
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                ALL_false = ALL_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                #print('ONSITE FALSE',DM_cat.ON_origin[i][0:22])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                ALL_false = ALL_false + 1
                DM_false = DM_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                #print('ONSITE FALSE',DM_cat.ON_origin[i][0:22])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                DM_false = DM_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                DM_false = DM_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                #print('ONSITE FALSE',DM_cat.ON_origin[i][0:22])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                DM_false = DM_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                E2_ON_false = E2_ON_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == True):
                    
                DM_false = DM_false + 1
                E2_ON_false = E2_ON_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                #print('ONSITE FALSE',DM_cat.ON_origin[i][0:22])
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                E2_VS_false = E2_VS_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                DM_false = DM_false + 1
                E2_VS_false = E2_VS_false + 1
                if (DM_cat.E2_mag[i] < E2_min):
                    E2_false = E2_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                ON_VS_false = ON_VS_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                DM_false = DM_false + 1
                ON_VS_false = ON_VS_false + 1
                if (DM_cat.ON_mag[i] < ON_min):
                    ON_false = ON_false + 1
                if (DM_cat.VS_mag[i] < VS_min):
                    VS_false = VS_false + 1
                
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == True) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                DM_missed = DM_missed + 1
                    
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
    
                
                if (DM_cat.E2_mag[i] >= E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                    
                if (DM_cat.ON_mag[i] >= ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] >= VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
        
                
            elif (fnmatch.fnmatch(DM_cat.anss_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.DM_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.E2_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.ON_origin[i][0:22], '0') == False) and \
                (fnmatch.fnmatch(DM_cat.VS_origin[i][0:22], '0') == False):
                    
                event_count = event_count + 1
                    
                t_origin = UTCDateTime(DM_cat.anss_origin[i][0:22])
                
                d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.DM_lat[i],DM_cat.DM_lon[i])
                DM_dist_err.append(abs(d[0]/1000))
                DM_mag_err.append((DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_abs_mag_err.append(abs(DM_cat.DM_mag[i] - DM_cat.anss_mag[i]))
                DM_t_alert = UTCDateTime(DM_cat.DM_alert[i][0:22])
                DM_t_diff.append(DM_t_alert - t_origin)
                DM_lons.append(DM_cat.anss_lon[i])
                DM_lats.append(DM_cat.anss_lat[i])
                
    
                
                if (DM_cat.E2_mag[i] < E2_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.E2_lat[i],DM_cat.E2_lon[i])
                    E2_dist_err.append(abs(d[0]/1000))
                    E2_mag_err.append((DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_mag_array.append((DM_cat.anss_mag[i], DM_cat.E2_mag[i]))
                    E2_abs_mag_err.append(abs(DM_cat.E2_mag[i] - DM_cat.anss_mag[i]))
                    E2_t_alert = UTCDateTime(DM_cat.E2_alert[i][0:22])
                    E2_t_diff.append(E2_t_alert - t_origin)
                    E2_lons.append(DM_cat.anss_lon[i])
                    E2_lats.append(DM_cat.anss_lat[i])
                    
                if (DM_cat.ON_mag[i] < ON_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.ON_lat[i],DM_cat.ON_lon[i])
                    ON_dist_err.append(abs(d[0]/1000))
                    ON_mag_err.append((DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_abs_mag_err.append(abs(DM_cat.ON_mag[i] - DM_cat.anss_mag[i]))
                    ON_t_alert = UTCDateTime(DM_cat.ON_alert[i][0:22])
                    ON_t_diff.append(ON_t_alert - t_origin)
                    ON_lons.append(DM_cat.anss_lon[i])
                    ON_lats.append(DM_cat.anss_lat[i])
                
                if (DM_cat.VS_mag[i] < VS_min):
                    d = gps2DistAzimuth(DM_cat.anss_lat[i],DM_cat.anss_lon[i],DM_cat.VS_lat[i],DM_cat.VS_lon[i])
                    VS_dist_err.append(abs(d[0]/1000))
                    VS_mag_err.append((DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_abs_mag_err.append(abs(DM_cat.VS_mag[i] - DM_cat.anss_mag[i]))
                    VS_t_alert = UTCDateTime(DM_cat.VS_alert[i][0:22])
                    VS_t_diff.append(VS_t_alert - t_origin)
                    VS_lons.append(DM_cat.anss_lon[i])
                    VS_lats.append(DM_cat.anss_lat[i])
        
 
        
    print('%s %s M>=%s' % (d_dir, date_range, mag))
    
    
    ########## CHOOSE WHETHER TO OUTPUT RESULTS FROM DM, E2, ON, AND/OR VS:
        
    #################### DM ######################   
    #    
    #P.figure()
    #bins = arange(0, 100, 5)
    #n, bins, patches = P.hist(DM_dist_err, bins, histtype='bar')
    #plt.ylabel('Frequency')
    #plt.xlabel('Location Difference (km)')
    ##plt.title('Location Error (Onshore Events Only)')
    #plt.title('DM Location Error - %s, %s M>=%s' % (d_dir, filelist[p][22:-4], mag))
    #P.show()
    #
    #
    #
    #P.figure()
    #mag_bins = arange(-2, 2, 0.25)
    #mag_n, mag_bins, mag_patches = P.hist(DM_mag_err, mag_bins, histtype='bar')
    #plt.ylabel('Frequency')
    #plt.xlabel('DM-ANSS Magnitude Estimate Difference')
    ##plt.title('Magnitude Error (Onshore Events Only)')
    #plt.title('DM Magnitude Error - %s, %s M>=%s' % (d_dir, filelist[p][22:-4], mag))
    #P.show()
    ##
    ##nValues = np.arange(0,max(DM_t_diff))
    #nValues = np.arange(0,22)
    #
    ## setup the normalization and the colormap
    #normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
    #colormap = cm.jet
    #
    ## setup the colorbar
    #scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    #scalarmappaple.set_array(nValues)
    ##
    ##P.figure()
    ##for i in range(0,len(DM_cat.anss_lat)):
    ##    plt.scatter(DM_cat.anss_lon[i],DM_cat.anss_lat[i],c=colormap(normalize((DM_t_diff[i]))))
    ##plt.ylabel('Latitude')
    ##plt.xlabel('Longitude')
    ##plt.title('Time Until First Alert (s)')
    ##plt.colorbar(scalarmappaple)
    ##P.show()
    #
    ## create figure and axes instances
    ##fig = plt.figure(figsize=(8,8))
    ##ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ### create polar stereographic Basemap instance.
    ### BAY AREA:
    ##if d_dir == 'BAY_AREA':
    ##    m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    ##                llcrnrlat=36,urcrnrlat=39.5,\
    ##                llcrnrlon=-124,urcrnrlon=-120,\
    ##                rsphere=6371200.,resolution='l',area_thresh=10000)
    ##            
    ### LA:
    ##elif d_dir == 'LA':
    ##    m = Basemap(projection='stere',lat_0=33,lon_0=-118,\
    ##                llcrnrlat=32.5,urcrnrlat=35,\
    ##                llcrnrlon=-119,urcrnrlon=-116,\
    ##                rsphere=6371200.,resolution='l',area_thresh=10000)
    ##            
    ### CALIFORNIA:
    ##elif d_dir == 'ALL_REGIONS' or d_dir=='CA':
    ##    m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    ##                llcrnrlat=32,urcrnrlat=42,\
    ##                llcrnrlon=-128,urcrnrlon=-113,\
    ##                rsphere=6371200.,resolution='l',area_thresh=10000)
    ###draw coastlines, state and country boundaries, edge of map.
    ###for i in range(0,len(DM_lons)):
    ###    m.scatter(DM_lons[i],DM_lats[i],c=colormap(normalize((DM_t_diff[i]))))
    ##m.drawcoastlines()
    ##m.drawstates()
    ##m.drawcountries()
    ### draw coastlines, country boundaries, fill continents.
    ###m.fillcontinents(lake_color='#99ffff')
    ### draw parallels.
    ##parallels = np.arange(0.,90,1.)
    ##m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    ##m.drawmeridians(np.arange(-130,-112,1),labels=[0,0,0,1])
    ##x, y = m(DM_lons,DM_lats)
    ##for i in range(0,len(x)):
    ##    m.scatter(x[i],y[i],80,c=colormap(normalize((DM_t_diff[i]))),edgecolor='')
    ##plt.colorbar(scalarmappaple)
    ###plt.title('DM Alert Times for Events in %s from %s through %s, M>%s Only' % (d_dir, tstart, tend, mag))
    ##plt.title('DM Alert Times for Events in %s from %s, M>%s Only' % (d_dir, date_range, mag))
    ##P.show()
    #
    #count=0
    #BA_tdiff=[]
    #BA_magerr=[]
    #BA_disterr=[]
    #LA_tdiff=[]
    #LA_magerr=[]
    #LA_disterr=[]
    #notBA_tdiff=[]
    #
    #for i in range(0,len(DM_lons)):
    #    if DM_lats[i]>36.44 and DM_lats[i]<39 and DM_lons[i]>-124 and DM_lons[i]<-120:
    #        #print('Bay Area',DM_cat.anss_origin[i], DM_t_diff[i],mag_err[i])
    #        BA_tdiff.append(DM_t_diff[i])
    #        BA_magerr.append(DM_mag_err[i])
    #        BA_disterr.append(DM_dist_err[i])
    #    elif DM_lats[i]>32.5 and DM_lats[i]<34.2 and DM_lons[i]>-119 and DM_lons[i]<-117:
    #        #print('LA',DM_cat.anss_origin[i], DM_t_diff[i],mag_err[i])
    #        LA_tdiff.append(DM_t_diff[i])
    #        LA_magerr.append(DM_mag_err[i])
    #        LA_disterr.append(DM_dist_err[i])
    #    else:
    #        notBA_tdiff.append(DM_t_diff[i])
    #        
    ######## OUTPUT STATISTICS:
    #print('median(DM_tdiff)',median(DM_t_diff))
    #print('std(DM_tdiff)',std(DM_t_diff))
    #
    #print('Median(DM_dist_err)',median(DM_dist_err))
    #print('Std(DM_dist_err)',std(DM_dist_err))
    #print('Median(DM_abs_mag_err)',median(DM_abs_mag_err))
    #print('Std(DM_abs_mag_err)',std(DM_abs_mag_err))
    #
    #
    ##print('median(DM_BA_tdiff)',median(BA_tdiff))
    ##print('median(DM_LA_tdiff)',median(LA_tdiff))
    ##print('median(DM_notBA_tdiff)',median(notBA_tdiff))
    #    
        
    #################### E2 ######################  
        
    P.figure()
    bins = arange(0, 100, 5)
    n, bins, patches = P.hist(E2_dist_err, bins, histtype='bar')
    plt.ylabel('Number of Events')
    plt.xlabel('Location Difference (km)')
    plt.xlim(0,80)
    #plt.title('Location Error (Onshore Events Only)')
    #plt.title('E2 Location Error - %s, %s M>=%s' % (d_dir, filelist[p][22:-4], mag))
    plt.title('E2 Location Error - %s, %s M>=%s' % (d_dir, date_range, mag))
    P.show()
    
    
    
    P.figure()
    mag_bins = arange(-2, 2, 0.25)
    mag_n, mag_bins, mag_patches = P.hist(E2_mag_err, mag_bins, histtype='bar')
    plt.ylabel('Number of Events')
    plt.xlabel('E2-ANSS Magnitude Estimate Difference')
    #plt.title('Magnitude Error (Onshore Events Only)')
    #plt.title('E2 Magnitude Error - %s, %s M>=%s' % (d_dir, filelist[p][22:-4], mag))
    plt.title('E2 Magnitude Error - %s, %s M>=%s' % (d_dir, date_range, mag))
    P.show()
    #
    #nValues = np.arange(0,max(E2_t_diff))
    nValues = np.arange(0,22)
    
    # setup the normalization and the colormap
    normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
    colormap = cm.jet
    
    # setup the colorbar
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(nValues)
    #
    P.figure()
    for i in range(0,len(E2_lats)):
        plt.scatter(E2_lons[i],E2_lats[i],c=colormap(normalize((E2_t_diff[i]))))
    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.title('Time Until First Alert (s)')
    plt.colorbar(scalarmappaple)
    P.show()
    
        #create figure and axes instances
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # create polar stereographic Basemap instance.
    # BAY AREA:
    if d_dir == 'BAY_AREA':
        m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
                    llcrnrlat=36,urcrnrlat=39.5,\
                    llcrnrlon=-124,urcrnrlon=-120,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
                
    # LA:
    elif d_dir == 'LA':
        m = Basemap(projection='stere',lat_0=33,lon_0=-118,\
                    llcrnrlat=32.5,urcrnrlat=35,\
                    llcrnrlon=-119,urcrnrlon=-116,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
                
    # CALIFORNIA:
    elif (d_dir == 'ALL_REGIONS') or (d_dir=='CA'):
        m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
                    llcrnrlat=32,urcrnrlat=42,\
                    llcrnrlon=-128,urcrnrlon=-113,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    # draw coastlines, country boundaries, fill continents.
    m.fillcontinents(color='lightgrey',lake_color='aqua')
    # draw parallels.
    parallels = np.arange(0.,90,1.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    m.drawmeridians(np.arange(-130,-112,2),labels=[0,0,0,1])
    x, y = m(E2_lons,E2_lats)
    for i in range(0,len(x)):
        m.scatter(x[i],y[i],80,c=colormap(normalize((E2_t_diff[i]))),zorder=100,edgecolor='')
    plt.colorbar(scalarmappaple)
    plt.title('E2 Alert Times for Events in %s from %s, M>%s Only' % (d_dir, date_range, mag))
    P.show()
    
    P.figure()
    for i in range(0,len(E2_mag_array)):
        P.scatter(E2_mag_array[i][1],E2_mag_array[i][0])
    plt.xlabel('ANSS Magnitude')
    plt.ylabel('E2 Magnitude')
    plt.title('E2 and ANSS Magnitudes for Events in %s from %s, M>%s Only' % (d_dir, date_range, mag))
    P.show()
    
    #print('ANSS/ElarmS Mags:', E2_mag_array)
    #print('ANSS/ElarmS Location Difference:', E2_dist_err)
    
    count=0
    BA_tdiff=[]
    BA_magerr=[]
    BA_disterr=[]
    LA_tdiff=[]
    LA_magerr=[]
    LA_disterr=[]
    notBA_tdiff=[]
    
    for i in range(0,len(E2_lons)):
        if E2_lats[i]>36.44 and E2_lats[i]<39 and E2_lons[i]>-124 and E2_lons[i]<-120:
            #print('Bay Area',DM_cat.anss_origin[i], E2_t_diff[i],mag_err[i])
            BA_tdiff.append(E2_t_diff[i])
            BA_magerr.append(E2_mag_err[i])
            BA_disterr.append(E2_dist_err[i])
        elif E2_lats[i]>32.5 and E2_lats[i]<34.2 and E2_lons[i]>-119 and E2_lons[i]<-117:
            #print('LA',DM_cat.anss_origin[i], E2_t_diff[i],mag_err[i])
            LA_tdiff.append(E2_t_diff[i])
            LA_magerr.append(E2_mag_err[i])
            LA_disterr.append(E2_dist_err[i])
        else:
            notBA_tdiff.append(E2_t_diff[i])
            
    ######## OUTPUT STATISTICS:
    
    print('median(E2_tdiff)',median(E2_t_diff))
    print('std(E2_tdiff)',std(E2_t_diff))
    
    print('Median(E2_dist_err)',median(E2_dist_err))
    print('Std(E2_dist_err)',std(E2_dist_err))
    print('Median(E2_abs_mag_err)',median(E2_abs_mag_err))
    print('Std(E2_abs_mag_err)',std(E2_abs_mag_err))
    
    #print('median(E2_BA_tdiff)',median(BA_tdiff))
    #print('median(E2_LA_tdiff)',median(LA_tdiff))
    #print('median(E2_notBA_tdiff)',median(notBA_tdiff))
            
            
            
    ##################### ONSITE ######################    
    #if len(ON_dist_err) >0:    
    #    P.figure()
    #    bins = arange(0, 100, 5)
    #    n, bins, patches = P.hist(ON_dist_err, bins, histtype='bar')
    #    plt.ylabel('Frequency')
    #    plt.xlabel('Location Difference (km)')
    #    plt.title('Onsite Location Error - %s, %s M>=%s' % (d_dir, date_range, mag))
    #    P.show()
    #    
    #    
    #    
    #    P.figure()
    #    mag_bins = arange(-2, 2, 0.25)
    #    mag_n, mag_bins, mag_patches = P.hist(ON_mag_err, mag_bins, histtype='bar')
    #    plt.ylabel('Frequency')
    #    plt.xlabel('ON-ANSS Magnitude Estimate Difference')
    #    plt.title('Onsite Magnitude Error - %s, %s M>=%s' % (d_dir, date_range, mag))
    #    P.show()
    #    #
    #    nValues = np.arange(0,22)
    #    
    #    # setup the normalization and the colormap
    #    normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
    #    colormap = cm.jet
    #    
    #    # setup the colorbar
    #    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    #    scalarmappaple.set_array(nValues)
    #    
    #    P.figure()
    #    for i in range(0,len(ON_lats)):
    #        plt.scatter(ON_lons[i],ON_lats[i],c=colormap(normalize((ON_t_diff[i]))))
    #    plt.ylabel('Latitude')
    #    plt.xlabel('Longitude')
    #    plt.title('Time Until First Alert (s)')
    #    plt.colorbar(scalarmappaple)
    #    P.show()
    #    
    #        #create figure and axes instances
    #    fig = plt.figure(figsize=(8,8))
    #    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    #    # create polar stereographic Basemap instance.
    #    # BAY AREA:
    #    if d_dir == 'BAY_AREA':
    #        m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    #                    llcrnrlat=36,urcrnrlat=39.5,\
    #                    llcrnrlon=-124,urcrnrlon=-120,\
    #                    rsphere=6371200.,resolution='l',area_thresh=10000)
    #                
    #    # LA:
    #    elif d_dir == 'LA':
    #        m = Basemap(projection='stere',lat_0=33,lon_0=-118,\
    #                    llcrnrlat=32.5,urcrnrlat=35,\
    #                    llcrnrlon=-119,urcrnrlon=-116,\
    #                    rsphere=6371200.,resolution='l',area_thresh=10000)
    #                
    #    # CALIFORNIA:
    #    elif (d_dir == 'ALL_REGIONS') or (d_dir=='CA'):
    #        m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    #                    llcrnrlat=32,urcrnrlat=42,\
    #                    llcrnrlon=-128,urcrnrlon=-113,\
    #                    rsphere=6371200.,resolution='l',area_thresh=10000)
    #    m.drawcoastlines()
    #    m.drawstates()
    #    m.drawcountries()
    #    # draw coastlines, country boundaries, fill continents.
    #    # draw parallels.
    #    parallels = np.arange(0.,90,1.)
    #    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #    m.drawmeridians(np.arange(-130,-112,2),labels=[0,0,0,1])
    #    x, y = m(ON_lons,ON_lats)
    #    for i in range(0,len(x)):
    #        m.scatter(x[i],y[i],80,c=colormap(normalize((ON_t_diff[i]))),edgecolor='')
    #    plt.colorbar(scalarmappaple)
    #    plt.title('Onsite Alert Times for Events in %s from %s, M>%s Only' % (d_dir, date_range, mag))
    #    P.show()
    #    
    #    count=0
    #    BA_tdiff=[]
    #    BA_magerr=[]
    #    BA_disterr=[]
    #    LA_tdiff=[]
    #    LA_magerr=[]
    #    LA_disterr=[]
    #    notBA_tdiff=[]
    #    
    #    for i in range(0,len(ON_lons)):
    #        if ON_lats[i]>36.44 and ON_lats[i]<39 and ON_lons[i]>-124 and ON_lons[i]<-120:
    #            BA_tdiff.append(ON_t_diff[i])
    #            BA_magerr.append(ON_mag_err[i])
    #            BA_disterr.append(ON_dist_err[i])
    #        elif ON_lats[i]>32.5 and ON_lats[i]<34.2 and ON_lons[i]>-119 and ON_lons[i]<-117:
    #            LA_tdiff.append(ON_t_diff[i])
    #            LA_magerr.append(ON_mag_err[i])
    #            LA_disterr.append(ON_dist_err[i])
    #        else:
    #            notBA_tdiff.append(ON_t_diff[i])
    #            
    ######## OUTPUT STATISTICS:
    #    print('median(ON_tdiff)',median(ON_t_diff))
    #    print('std(ON_tdiff)',std(ON_t_diff))
    #    
    #    print('Median(ON_dist_err)',median(ON_dist_err))
    #    print('Std(ON_dist_err)',std(ON_dist_err))
    #    print('Median(ON_abs_mag_err)',median(ON_abs_mag_err))
    #    print('Std(ON_abs_mag_err)',std(ON_abs_mag_err))
    #    
    #    #print('median(ON_BA_tdiff)',median(BA_tdiff))
    #    #print('median(ON_LA_tdiff)',median(LA_tdiff))
    #    #print('median(ON_notBA_tdiff)',median(notBA_tdiff))
    #    
    #else:
    #    continue
    #
    #
    
    #################### VS ######################   
    #    
    #P.figure(10)
    #bins = arange(0, 100, 5)
    #n, bins, patches = P.hist(VS_dist_err, bins, histtype='bar')
    #plt.ylabel('Frequency')
    #plt.xlabel('Location Difference (km)')
    ##plt.title('Location Error (Onshore Events Only)')
    #plt.title('VS Location Error')
    #P.show()
    #
    #
    #
    #P.figure(11)
    #mag_bins = arange(-2, 2, 0.25)
    #mag_n, mag_bins, mag_patches = P.hist(VS_mag_err, mag_bins, histtype='bar')
    #plt.ylabel('Frequency')
    #plt.xlabel('VS-ANSS Magnitude Estimate Difference')
    ##plt.title('Magnitude Error (Onshore Events Only)')
    #plt.title('VS Magnitude Error')
    #P.show()
    ##
    ##nValues = np.arange(0,max(VS_t_diff))
    #nValues = np.arange(0,22)
    #
    ## setup the normalization and the colormap
    #normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
    #colormap = cm.jet
    #
    ## setup the colorbar
    #scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    #scalarmappaple.set_array(nValues)
    ##
    ##P.figure(12)
    ##for i in range(0,len(DM_cat.anss_lat)):
    ##    plt.scatter(DM_cat.anss_lon[i],DM_cat.anss_lat[i],c=colormap(normalize((VS_t_diff[i]))))
    ##plt.ylabel('Latitude')
    ##plt.xlabel('Longitude')
    ##plt.title('Time Until First Alert (s)')
    ##plt.colorbar(scalarmappaple)
    ##P.show()
    #
    ## create figure and axes instances
    #fig = plt.figure(figsize=(8,8))
    #ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ## create polar stereographic Basemap instance.
    ## BAY AREA:
    ##m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    ##            llcrnrlat=36,urcrnrlat=39.5,\
    ##            llcrnrlon=-124,urcrnrlon=-120,\
    ##            rsphere=6371200.,resolution='l',area_thresh=10000)
    ## LA:
    #m = Basemap(projection='stere',lat_0=33,lon_0=-118,\
    #            llcrnrlat=32.5,urcrnrlat=35,\
    #            llcrnrlon=-119,urcrnrlon=-116,\
    #            rsphere=6371200.,resolution='l',area_thresh=10000)
    ## CALIFORNIA:
    ##m = Basemap(projection='stere',lat_0=37,lon_0=-120,\
    ##            llcrnrlat=32,urcrnrlat=42,\
    ##            llcrnrlon=-128,urcrnrlon=-113,\
    ##            rsphere=6371200.,resolution='l',area_thresh=10000)
    ## draw coastlines, state and country boundaries, edge of map.
    ##for i in range(0,len(VS_lons)):
    #    #m.scatter(VS_lons[i],VS_lats[i],c=colormap(normalize((VS_t_diff[i]))))
    #m.drawcoastlines()
    #m.drawstates()
    #m.drawcountries()
    ## draw coastlines, country boundaries, fill continents.
    ##m.fillcontinents(lake_color='#99ffff')
    ## draw parallels.
    #parallels = np.arange(0.,90,1.)
    #m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    #m.drawmeridians(np.arange(-130,-112,1),labels=[0,0,0,1])
    #x, y = m(VS_lons,VS_lats)
    #for i in range(0,len(x)):
    #    m.scatter(x[i],y[i],80,c=colormap(normalize((VS_t_diff[i]))),edgecolor='')
    #plt.colorbar(scalarmappaple)
    ##plt.title('VS Alert Times for Events in California from 2014-11-07 through 2015-10-31, M>4 Only')
    ##plt.title('VS Alert Times for Events in California from 2014-11-07 through 2015-10-31, M>3 Only')
    ##plt.title('VS Alert Times for Events in Bay Area from 2014-11-07 through 2015-10-31, M>3 Only')
    #P.show()
    #
    #count=0
    #BA_tdiff=[]
    #BA_magerr=[]
    #BA_disterr=[]
    #LA_tdiff=[]
    #LA_magerr=[]
    #LA_disterr=[]
    #notBA_tdiff=[]
    #
    #for i in range(0,len(VS_lons)):
    #    if VS_lats[i]>36.44 and VS_lats[i]<39 and VS_lons[i]>-124 and VS_lons[i]<-120:
    #        #print('Bay Area',DM_cat.anss_origin[i], VS_t_diff[i],mag_err[i])
    #        BA_tdiff.append(VS_t_diff[i])
    #        BA_magerr.append(VS_mag_err[i])
    #        BA_disterr.append(VS_dist_err[i])
    #    elif VS_lats[i]>32.5 and VS_lats[i]<34.2 and VS_lons[i]>-119 and VS_lons[i]<-117:
    #        #print('LA',DM_cat.anss_origin[i], VS_t_diff[i],mag_err[i])
    #        LA_tdiff.append(VS_t_diff[i])
    #        LA_magerr.append(VS_mag_err[i])
    #        LA_disterr.append(VS_dist_err[i])
    #    else:
    #        notBA_tdiff.append(VS_t_diff[i])
    #        
    ######## OUTPUT STATISTICS:
    ##print('median(VS_tdiff)',median(VS_t_diff))
    ##print('std(VS_tdiff)',std(VS_t_diff))
    ##
    ##print('Median(VS_dist_err)',median(VS_dist_err))
    ##print('Std(VS_dist_err)',std(VS_dist_err))
    ##print('Median(VS_abs_mag_err)',median(VS_abs_mag_err))
    ##print('Std(VS_abs_mag_err)',std(VS_abs_mag_err))
    ##
    ##print('median(VS_BA_tdiff)',median(BA_tdiff))
    ##print('median(VS_LA_tdiff)',median(LA_tdiff))
    ##print('median(VS_notBA_tdiff)',median(notBA_tdiff))
    #
    #########
    
    
    ######## OUTPUT STATISTICS:
    
    print('# Events', event_count)
    print('DM False', DM_false)
    print('DM Missed', DM_missed)
    print('E2 False', E2_false)
    print('E2 Missed', E2_missed)
    #print('ON False', ON_false)
    #print('ON Missed', ON_missed)
    #print('VS False', VS_false)
    #print('VS Missed', VS_missed)
    
    #print('All False', ALL_false)
    #print('All Missed', ALL_missed)
    #print('E2 & ON False', E2_ON_false)
    #print('E2 & ON Missed', E2_ON_missed)
    #print('E2 & VS False', E2_VS_false)
    #print('E2 & VS Missed', E2_VS_missed)
    #print('ON & VS False', ON_VS_false)
    #print('ON & VS Missed', ON_VS_missed)
    p = p+1
    
    