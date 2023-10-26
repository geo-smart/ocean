##################
#
# imports, utility functions
#
##################

import os, sys, time, glob, warnings
from IPython.display import clear_output
warnings.filterwarnings('ignore')
this_dir = os.getcwd()  

from matplotlib import pyplot as plt
from matplotlib import colors as mplcolors
from matplotlib import animation, rc
from matplotlib import dates as mdates

import numpy as np, pandas as pd, xarray as xr
from numpy import datetime64 as dt64, timedelta64 as td64

from ipywidgets import *
from traitlets import dlink
from IPython.display import HTML, Video


##################
#
# parameter configuration
#
##################

# time ranges for midnight and noon profiles, adjusted for UTC
# midn0 - midn1 is a time range for the midnight profile start
# noon0 - noon1 is a time range for the noon profile start
midn0 = td64( 7*60 + 10, 'm')        # 7 hours 10 minutes
midn1 = td64( 7*60 + 34, 'm')        # 7 hours 34 minutes
noon0 = td64(20*60 + 30, 'm')        # 20 hours 30 minutes
noon1 = td64(20*60 + 54, 'm')        # 20 hours 54 minutes 

# global sensor range parameters for charting data: based on osb shallow profiler data

# axis ranges for a variety of sensors
temp_lo,        temp_hi          =     6.5,       11.
sal_lo,         sal_hi           =    31.5,       34.5
do_lo,          do_hi            =    50.0,      300.
chlora_lo,      chlora_hi        =    -0.1,        1.2
bb_lo,          bb_hi            =     0.0007,     0.0020
cdom_lo,        cdom_hi          =     0.6,        1.4

colorT = 'black'
colorS = 'xkcd:blood orange'
colorO = 'xkcd:blue'
colorA = 'xkcd:green'
colorB = 'xkcd:dark cyan'
colorC = 'red'

labelT = 'Temperature'
labelO = 'Oxygen'
labelS = 'Salinity'
labelA = 'Chlor-A'
labelB = 'Backscatter'
labelC = 'CDOM/FDOM'

optionsList = [labelO, labelT, labelS, labelA, labelB, labelC]


########################
#
# Functions and Configuration
#
########################


################
# convenience functions 
################
# abbreviating 'datetime64' and so on
################

def doy(theDatetime): return 1 + int((theDatetime - dt64(str(theDatetime)[0:4] + '-01-01')) / td64(1, 'D'))

def dt64_from_doy(year, doy): return dt64(str(year) + '-01-01') + td64(doy-1, 'D')

def day_of_month_to_string(d): return str(d) if d > 9 else '0' + str(d)

def TimeOverlap(s0, s1, t0, t1):
    '''
    Do two time ranges overlap?
    '''
    if s0 > t0 and s0 < t1: return True
    if s1 > t0 and s1 < t1: return True
    if s0 < t0 and s1 > t1: return True
    return False


#################
# Time series profile metadata load function
#################
# Read in pre-processed profiler metadata for subsequent time-series subsetting.
# Shallow profiler metadata are timestamps for Ascent / Descent / Rest. These are stored 
#   as one-year-duration CSV files in the Profiles subfolder; are read into a Pandas 
#   Dataframe. Columns correspond to ascent start time and so on, as noted in the code.
#################

def ReadProfileMetadata(fnm):
    """
    Profiles are saved by site and year as 12-tuples. Here we read only
    the datetimes (not the indices) so there are only six values. These
    are converted to Timestamps. They correspond to ascend start/end, 
    descend start/end and rest start/end. Timestamps are a bit easier to
    use than datetime64 values, being essentially wrappers around the latter with
    additional utility.
    """
    p = pd.read_csv(fnm, usecols=["1", "3", "5", "7", "9", "11"])
    p.columns=['ascent_start', 'ascent_end', 'descent_start', 'descent_end', 'rest_start', 'rest_end']
    p['ascent_start']  = pd.to_datetime(p['ascent_start'])
    p['ascent_end']    = pd.to_datetime(p['ascent_end'])
    p['descent_start'] = pd.to_datetime(p['descent_start'])
    p['descent_end']   = pd.to_datetime(p['descent_end'])
    p['rest_start']    = pd.to_datetime(p['rest_start'])
    p['rest_end']      = pd.to_datetime(p['rest_end'])
    return p




#######################
# Time series metadata (index range) function
#######################
# Given a time range we want the indices of the profiles within.
#######################
def ProfileListFromTimeWindow(p, t1, t2):
    n = p.shape[0]
    plist = []
    for i in range(n):
        if p['ascent_start'][i] >= t1 and p['ascent_start'][i] <= t2: plist.append(i)
    return plist


def GenerateTimeWindowIndices(p, date0, date1, time0, time1):
    '''
    Given two day boundaries and a time window (UTC) within a day: Return a list
    of indices of profiles that start within both the day and time bounds. This 
    works from the passed dataframe of profile times.
    '''
    nprofiles = len(p)
    pIndices = []
    for i in range(nprofiles):
        a0 = p["ascent_start"][i]
        if a0 >= date0 and a0 <= date1 + td64(1, 'D'):
            delta_t = a0 - dt64(a0.date())
            if delta_t >= time0 and delta_t <= time1: pIndices.append(i)
    return pIndices


def ProfileEvaluation(t0, t1, p):
    '''
    Give this function a time range t0 - t1 and the profile metadata structure; it will
    return how many profiles fall within that time window as well as how many are local
    midnight and local noon 'special' profiles.
    
    Additional: At this time the profile metadata in p is broken up by year of interest and site.
    For example the code above concerns Oregon Slope Base (OSB) and the year 2021. 
    Only profiles through June are available.
    
    This function evaluates profiles within a given time range: How many profiles are there?
    How many 'local noon', how many 'local midnight'? This is a simple way to check profiler 
    operating consistency. This depends in turn on the profiler metadata reliability.
    '''
    global midn0, midn1, noon0, noon1
    
    nTotal = 0
    nMidn = 0
    nNoon = 0
    nNinePerDay = 0

    for i in range(len(p)):
            
        if p["ascent_start"][i] >= t0 and p["ascent_start"][i] <= t1:
            nTotal += 1
            
            if p["descent_end"][i] - p["descent_start"][i] >= td64(60, 'm'):
                
                tProf = p["ascent_start"][i]
                day_time = tProf - dt64(tProf.date())

                if   day_time > midn0 and day_time < midn1: nMidn += 1
                elif day_time > noon0 and day_time < noon1: nNoon += 1
                else: print("found a long descent that did not fit noon or midnight...")
        
    return nTotal, nMidn, nNoon


def GetDiscreteSummaryCastSubset(dsDf, cast, columns):
    '''
    dsDf is a Discrete Summary Dataframe
    cast is a string corresponding to the cast identifier, e.g. 'CTD-001'
    columns is a list of column names to extract from the full Dataframe
    Returns a Dataframe for 'just that cast' and 'just those parameters'
    '''
    return dsDf.loc[(dsDf['cast']==cast)][columns]



def ChartAB(p, xrng, pIdcs, A, Az, Albl, Acolor, B, Bz, Blbl, Bcolor, wid, hgt, \
            z0=-200., z1=0., legA='ascent', legB='ascent'):
    """
    Make a series of charts comparing two types of sensor data, A and B.
    The data are passed in as DataArrays: A and Az are data and z coordinates respectively.
    So A might be P.par (PAR DataArray) and depth Az would be P.z. Both use time as 
    their dimension. Charting is done over a set of passed profile indices pIdcs[].
    The number of profiles charted is constrained: Too many may bog down the kernel.
    """
    global midn0, midn1, noon0, noon1
        
    # if too many charts are requested: Take the first 117 only
    ncharts = len(pIdcs)
    if ncharts > 117: ncharts = 117
    print("Attempting", ncharts, "charts\n")

    # set up the requested number of charts in a vertical column
    fig, axs = plt.subplots(ncharts, 1, figsize=(wid, hgt*ncharts), tight_layout=True)

    # create a list of twin axes, one for each chart
    axstwin0 = [axs[i].twiny() for i in range(ncharts)]

    keyA0, keyA1 = legA + "_start", legA + "_end"
    keyB0, keyB1 = legB + "_start", legB + "_end"
            
    # this index i will range across the dataframe indices for ascent profiles
    for i in range(ncharts):
        
        # Need both a profile index into the profile dataframe p and a chart
        #   index 0, 1, 2, ... These are respectively pIdx and i
        pIdx = pIdcs[i]

        tA0, tA1 = p[keyA0][pIdx], p[keyA1][pIdx]
        tB0, tB1 = p[keyB0][pIdx], p[keyB1][pIdx]
        
        Ax, Ay = A.sel(time=slice(tA0,  tA1)), Az.sel(time=slice(tA0, tA1))
        Bx, By = B.sel(time=slice(tB0,  tB1)), Bz.sel(time=slice(tB0, tB1))
        
        axs[i].plot(Ax, Ay, ms = 4., color=Acolor, mfc=Acolor)
        axstwin0[i].plot(Bx, By, markersize = 4., color=Bcolor, mfc=Bcolor)
        
        # axis ranges
        if i == 0: axs[i].set(title = Albl + ' (' + Acolor + ', lower x-axis) and ' \
                                    + Blbl + ' (' + Bcolor + ', upper x-axis)')

        # Set axis ranges from passed list of pairs xrng[][]
        axs[i].set(     xlim = (xrng[0][0], xrng[0][1]), ylim = (z0, z1))
        axstwin0[i].set(xlim = (xrng[1][0], xrng[1][1]), ylim = (z0, z1))

        # chart timestamp (embellish for noon / midnight)
        #   bug: qualifier will fail if the 'A' data is not ascent-type (non-critical)
        ascent_start_time = 'Start UTC: ' + str(tA0)
        delta_t = tA0-dt64(tA0.date())
        if delta_t > midn0 and delta_t < midn1: ascent_start_time += " MIDNIGHT local"
        if delta_t > noon0 and delta_t < noon1: ascent_start_time += " NOON local"

        xlabel = xrng[0][0] + (xrng[0][1] - xrng[0][0])/2.
        axs[i].text(xlabel, -10., ascent_start_time)
        
    return fig, axs


def CheckHITL(f, t1, t2):
    annot, found_overlap, overlap_rows = pd.read_csv(f), False, []
    for i in range(annot.shape[0]):
        if TimeOverlap(t1, t2, dt64(annot['beginDate'][i]), dt64(annot['endDate'][i])):
            found_overlap = True; overlap_rows.append(i)
    if found_overlap: 
        print('\n', len(overlap_rows), "time overlaps: annotations vs time range of interest. Rows:", overlap_rows)
        for i in range(len(overlap_rows)):
            print('      ', annot['beginDate'][overlap_rows[i]], ' to ', annot['endDate'][overlap_rows[i]])
        print()
    else: print("\nNo overlap between annotation record and time range of interest.\n")



##################
# more parameter configuration
##################
# Load the 2021 Oregon Slope Base profile metadata; and some March 2021 sensor datasets
##################

# Note these are profile times for Axial Base
p = ReadProfileMetadata(os.getcwd()+"/../Profiles/osb2021.csv")

# Some code to test out the above ProfileEvaluation() function
t0, t1 = dt64('2021-03-01'), dt64('2021-04-01')
nDays = (t1 - t0).astype(int)
nTotal, nMidn, nNoon = ProfileEvaluation(t0, t1, p)

print("OOI RCA Oregon Slope Base: Shallow Profiler status report")
print("=========================================================")
print("For March 2021:")
print(nDays, 'days, translates to', nDays*9, 'possible profiles')
print("Actual:")
print(nTotal, 'profiles;', nMidn, 'at local midnight and', nNoon, 'at local noon')


