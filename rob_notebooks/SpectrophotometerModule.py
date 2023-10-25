##################
#
# imports
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

global_lo,    global_hi          =     0.0,        0.3

ba_lo, ba_hi = 0.0, 0.2
oa_lo, oa_hi = 0.15, 0.25

ba28_lo,        ba28_hi          =     global_lo,    global_hi
ba56_lo,        ba56_hi          =     global_lo,    global_hi
oa28_lo,        oa28_hi          =     global_lo,    global_hi
oa56_lo,        oa56_hi          =     global_lo,    global_hi

ba00_lo,        ba00_hi          =     global_lo,    global_hi
ba10_lo,        ba10_hi          =     global_lo,    global_hi
ba20_lo,        ba20_hi          =     global_lo,    global_hi
ba30_lo,        ba30_hi          =     global_lo,    global_hi
ba40_lo,        ba40_hi          =     global_lo,    global_hi
ba50_lo,        ba50_hi          =     global_lo,    global_hi
ba60_lo,        ba60_hi          =     global_lo,    global_hi
ba70_lo,        ba70_hi          =     global_lo,    global_hi
ba80_lo,        ba80_hi          =     global_lo,    global_hi

oa01_lo,        oa01_hi          =     global_lo,    global_hi
oa10_lo,        oa10_hi          =     global_lo,    global_hi
oa20_lo,        oa20_hi          =     global_lo,    global_hi
oa30_lo,        oa30_hi          =     global_lo,    global_hi
oa40_lo,        oa40_hi          =     global_lo,    global_hi
oa50_lo,        oa50_hi          =     global_lo,    global_hi
oa60_lo,        oa60_hi          =     global_lo,    global_hi
oa70_lo,        oa70_hi          =     global_lo,    global_hi
oa80_lo,        oa80_hi          =     global_lo,    global_hi


colorBA28 = 'black'
colorBA56 = 'red'
colorOA28 = 'blue'
colorOA56 = 'gold'


colorBA00 = 'black'
colorBA10 = 'red'
colorBA20 = 'orange'
colorBA30 = 'xkcd:gold'
colorBA40 = 'green'
colorBA50 = 'cyan'
colorBA60 = 'blue'
colorBA70 = 'xkcd:purple blue'
colorBA80 = 'magenta'

colorOA01 = 'black'
colorOA10 = 'red'
colorOA20 = 'orange'
colorOA30 = 'xkcd:gold'
colorOA40 = 'green'
colorOA50 = 'cyan'
colorOA60 = 'blue'
colorOA70 = 'xkcd:purple blue'
colorOA80 = 'magenta'

labelBA28 = 'Beam Attenuation ch 28'
labelBA56 = 'Beam Attenuation ch 56'

labelOA28 = 'Optical Absorption ch 28'
labelOA56 = 'Optical Absorption ch 56'

labelBA00 = 'Beam Attenuation ch 0'
labelBA10 = 'Beam Attenuation ch 10'
labelBA20 = 'Beam Attenuation ch 20'
labelBA30 = 'Beam Attenuation ch 30'
labelBA40 = 'Beam Attenuation ch 40'
labelBA50 = 'Beam Attenuation ch 50'
labelBA60 = 'Beam Attenuation ch 60'
labelBA70 = 'Beam Attenuation ch 70'
labelBA80 = 'Beam Attenuation ch 80'

labelOA01 = 'Optical Absorption ch 1'
labelOA10 = 'Optical Absorption ch 10'
labelOA20 = 'Optical Absorption ch 20'
labelOA30 = 'Optical Absorption ch 30'
labelOA40 = 'Optical Absorption ch 40'
labelOA50 = 'Optical Absorption ch 50'
labelOA60 = 'Optical Absorption ch 60'
labelOA70 = 'Optical Absorption ch 70'
labelOA80 = 'Optical Absorption ch 80'

optionsList = [labelOA28, labelOA56, labelBA28, labelBA56,                                                        \
               labelOA01, labelOA10, labelOA20, labelOA30, labelOA40, labelOA50, labelOA60, labelOA70, labelOA80, \
               labelBA00, labelBA10, labelBA20, labelBA30, labelBA40, labelBA50, labelBA60, labelBA70, labelBA80 ]
               

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



########################
#
# Large Dataset Preparation
#     Commented out when not in use
#
########################

# Received from OOINET: 13 files covering 2021-01-01 through 2021-06-01
# These were renamed 20210101.nc, ..., 20210313.nc, ..., 20210521.nc
# 
###########################
# concatenate into one very large 4.2GB data file
###########################
# ds = []                                                                                 # start with an empty list: to append datasets from the above files
# q=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13']      
# for i in range(len(q)):
#     fnm = '../../data/rca/spectrophotometer/optaa' + q[i] + '.nc'                       # The path '../..' goes to the parent directory, i.e. outside of this repo            
#     ds.append(xr.open_dataset(fnm))                                                     # Append the next file to the Dataset list
#     ds[i] = ds[i].swap_dims({'obs':'time'})                                             # Promote the coordinate 'time' to a dimension of the Dataset
#     ds[i] = ds[i][['optical_absorption', 'beam_attenuation']]                           # Ignoring 30 other data variables: Just retain oa and ba
#     ds[i] = ds[i].rename({'optical_absorption':'oa','beam_attenuation':'ba'})           # Rename those two variables to shortened names
#
# Inspection
# ds[7]                           # examine one of the 13 datasets to make sure it looks correct                                            
# ds[7].time[0:3]
#
# c = xr.concat(ds, 'time')          # concatenates along the time dimension to produce a single Dataset
#
# Inspection
# c
#
# Write to a single file
# c.to_netcdf('../../data/rca/spectrophotometer/concat.nc')
#
# Re-open the result
# c = xr.open_dataset('../../data/rca/spectrophotometer/concat.nc')
#
# Inspect curtain plot 'availability' using depth
# c.depth.plot()
# 
###############################
# Scope of a single source file
###############################
# This work will focus on the March 13 noon profile. From the source dataset:
# s=xr.open_dataset('../../data/rca/optaa/20210313.nc')
# s.time[0] gives 2021-03-13T20:42:16
# s gives Dimensions: 323k observations (obs), 83 wavelength; 6 Coordinates; 32 Data Variables
#     ...many of the Data Variables operate on the obs dimension; so that gets quite large.
# 
###############################
# Create subset file 1
#   2 wavelengths, 3 days: 5 profiles
###############################
# Local repo subset file 1: RepositoryData/rca/optaa/2021-MAR-13_thru_16_chs_28_and_56.nc (4.3MB)
# Metadata: OSB, AC-S, Year 2021, March 13 -- 16, channels 28 and 56 only: oa and ba only
# s=xr.open_dataset('../../data/rca/optaa/20210313.nc')
# s=s.swap_dims({'obs':'time'})
# s=s.sel(time=slice(dt64_from_doy(2021,72), dt64_from_doy(2021,75))).isel(wavelength=[28, 56])
# s=s[['optical_absorption', 'beam_attenuation']]
# s=s.rename({'optical_absorption':'oa', 'beam_attenuation':'ba'})
# s.to_netcdf('../RepositoryData/rca/optaa/2021-MAR-13_thru_16_chs_28_and_56.nc')
#
#################
# Time series metadata load function
#################
# Read in pre-processed profiler metadata for subsequent time-series subsetting.
# Shallow profiler metadata are timestamps for Ascent / Descent / Rest. These are stored 
#   as one-year-duration CSV files in the Profiles subfolder; are read into a Pandas 
#   Dataframe. Columns correspond to ascent start time and so on, as noted in the code.


def ReadProfileMetadata(fnm):
    """
    Profiles are saved by site and year as 12-tuples. Here we read only
    the datetimes (not the indices) so there are only six values. These
    are converted to Timestamps. They correspond to ascend start/end, 
    descend start/end and rest start/end. Timestamps are a bit easier to
    use than datetime64 values, being essentially wrappers around the latter with
    additional utility.
    """
    pDf = pd.read_csv(fnm, usecols=["1", "3", "5", "7", "9", "11"])
    pDf.columns=['ascent_start', 'ascent_end', 'descent_start', 'descent_end', 'rest_start', 'rest_end']
    pDf['ascent_start']  = pd.to_datetime(pDf['ascent_start'])
    pDf['ascent_end']    = pd.to_datetime(pDf['ascent_end'])
    pDf['descent_start'] = pd.to_datetime(pDf['descent_start'])
    pDf['descent_end']   = pd.to_datetime(pDf['descent_end'])
    pDf['rest_start']    = pd.to_datetime(pDf['rest_start'])
    pDf['rest_end']      = pd.to_datetime(pDf['rest_end'])
    return pDf

#######################
# Profile 'index list' generator
#######################
# Given a time range we want the indices of the profiles within
#######################
def GenerateTimeWindowIndices(pDf, date0, date1, time0, time1):
    '''
    Given two day boundaries and a time window (UTC) within a day: Return a list
    of indices of profiles that start within both the day and time bounds. This 
    works from the passed dataframe of profile times.
    '''
    nprofiles = len(pDf)
    pIndices = []
    for i in range(nprofiles):
        a0 = pDf["ascent_start"][i]
        if a0 >= date0 and a0 <= date1 + td64(1, 'D'):
            delta_t = a0 - dt64(a0.date())
            if delta_t >= time0 and delta_t <= time1: pIndices.append(i)
    return pIndices

##################
# Load the 2021 Oregon Slope Base profile metadata; and some March 2021 sensor datasets
##################

# Note these are profile times for Axial Base
p = ReadProfileMetadata(os.getcwd()+"/../Profiles/osb2021.csv")
