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

# global sensor range parameters for charting data: based on osb shallow profiler data

# axis ranges for a variety of sensors
par_lo,         par_hi           =   -10.0,      320.
nitrate_lo,     nitrate_hi       =     0.,        35.
do_lo,          do_hi            =    50.0,      300.
chlora_lo,      chlora_hi        =    -0.1,        1.2
temp_lo,        temp_hi          =     6.5,       11.
sal_lo,         sal_hi           =    31.5,       34.5
bb_lo,          bb_hi            =     0.0007,     0.0020
cdom_lo,        cdom_hi          =     0.6,        1.4
ph_lo,          ph_hi            =     7.6,        8.2
si412_lo,       si412_hi         =     0.0,       80.0
si443_lo,       si443_hi         =     0.0,       80.0
si490_lo,       si490_hi         =     0.0,       80.0
si510_lo,       si510_hi         =     0.0,       80.0
si555_lo,       si555_hi         =     0.0,       80.0
si620_lo,       si620_hi         =     0.0,       15.0
si683_lo,       si683_hi         =     0.0,        6.0
veast_lo,       veast_hi         =    -0.4,        0.4
vnorth_lo,      vnorth_hi        =    -0.4,        0.4
vup_lo,         vup_hi           =    -0.4,        0.4

colorT = 'black'
colorS = 'xkcd:blood orange'
colorO = 'xkcd:blue'
colorA = 'xkcd:green'
colorB = 'xkcd:dark cyan'
colorC = 'red'
colorN = 'xkcd:gold'
colorP = 'magenta'
colorH = 'xkcd:purple blue'

colorTd = 'grey'
colorSd = 'xkcd:yellow orange'
colorOd = 'xkcd:azure'
colorAd = 'xkcd:pale green'
colorBd = 'xkcd:light turquoise'
colorCd = 'xkcd:pinkish'
colorNd = 'xkcd:brownish yellow'
colorPd = 'xkcd:barbie pink'
colorHd = 'xkcd:pastel purple'

labelT = 'Temperature'
labelO = 'Oxygen'
labelS = 'Salinity'
labelA = 'Chlor-A'
labelB = 'Backscatter'
labelC = 'CDOM/FDOM'
labelN = 'Nitrate'
labelP = 'PAR'
labelH = 'pH'

optionsList = [labelO, labelT, labelS, labelA, labelB, labelC, labelN, labelP]


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


###########
# Plot function
###########
# Customized plot to show profiler behavior over one day
###########


    
#################
# Time series metadata load function
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
# Time series metadata (index range) function
#######################
# Given a time range we want the indices of the profiles within.
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




def GetDiscreteSummaryCastSubset(dsDf, cast, columns):
    '''
    dsDf is a Discrete Summary Dataframe
    cast is a string corresponding to the cast identifier, e.g. 'CTD-001'
    columns is a list of column names to extract from the full Dataframe
    Returns a Dataframe for 'just that cast' and 'just those parameters'
    '''
    return dsDf.loc[(dsDf['cast']==cast)][columns]


# Load XArray Datasets from the smaller (intra-repo!) source files

def ReadOSB_March2021_1min():
    data_source = os.getcwd() + '/../RepositoryData/rca/'
    return                                                                         \
        xr.open_dataset(data_source + 'fluor/osb_chlora_march2021_1min.nc'),       \
        xr.open_dataset(data_source + 'fluor/osb_backscatter_march2021_1min.nc'),  \
        xr.open_dataset(data_source + 'fluor/osb_cdom_march2021_1min.nc'),         \
        xr.open_dataset(data_source + 'ctd/osb_temp_march2021_1min.nc'),           \
        xr.open_dataset(data_source + 'ctd/osb_salinity_march2021_1min.nc'),       \
        xr.open_dataset(data_source + 'ctd/osb_doxygen_march2021_1min.nc'),        \
        xr.open_dataset(data_source + 'pH/osb_ph_march2021_1min.nc'),              \
        xr.open_dataset(data_source + 'irrad/osb_spectir_march2021_1min.nc'),      \
        xr.open_dataset(data_source + 'nitrate/osb_nitrate_march2021_1min.nc'),    \
        xr.open_dataset(data_source + 'par/osb_par_march2021_1min.nc'),            \
        xr.open_dataset(data_source + 'current/osb_veast_march2021_1min.nc'),      \
        xr.open_dataset(data_source + 'current/osb_vnorth_march2021_1min.nc'),     \
        xr.open_dataset(data_source + 'current/osb_vup_march2021_1min.nc'),         \
        xr.open_dataset(data_source + 'pCO2/osb_pco2_march2021.nc')


def CompareAscentDescent(p, T, S, O, A, B, C):
    '''Get a sense of variability between ascent and subsequent descent'''
    
    pw = GenerateTimeWindowIndices(p, dt64_from_doy(2021, 65), dt64_from_doy(2021, 65), td64(0, 'h'), td64(24, 'h'))
    ncharts = len(pw)

    fig, axs = plt.subplots(ncharts, 3, figsize=(15, 4*ncharts), tight_layout=True)

    axt0 = [axs[i][0].twiny() for i in range(ncharts)]
    axt1 = [axs[i][1].twiny() for i in range(ncharts)]
    axt2 = [axs[i][2].twiny() for i in range(ncharts)]

    for i in range(pw[0], pw[-1]+1):

        axi = i - pw[0]

        t0, t1, t2 = p["ascent_start"][i], p["ascent_end"][i], p["descent_end"][i]

        Ta = T.sel(time=slice(t0, t1))
        Td = T.sel(time=slice(t1, t2))
        Sa = S.sel(time=slice(t0, t1))
        Sd = S.sel(time=slice(t1, t2))
        Oa = O.sel(time=slice(t0, t1))
        Od = O.sel(time=slice(t1, t2))
        Aa = A.sel(time=slice(t0, t1))
        Ad = A.sel(time=slice(t1, t2))
        Ba = B.sel(time=slice(t0, t1))
        Bd = B.sel(time=slice(t1, t2))
        Ca = C.sel(time=slice(t0, t1))
        Cd = C.sel(time=slice(t1, t2))

        axs[axi][0].plot(Ta.temp, Ta.z, color=colorT, marker='s', ms=4., mfc=colorT)
        axs[axi][0].plot(Td.temp, Td.z, color=colorTd, marker='v', ms=4., mfc=colorTd)
        axt0[axi].plot(Sa.salinity, Sa.z, color=colorS, marker='o', ms=4., mfc=colorS)
        axt0[axi].plot(Sd.salinity, Sd.z, color=colorSd, marker='^', ms=4., mfc=colorSd)

        axs[axi][1].plot(Oa.doxygen, Oa.z, color=colorO, marker='s', ms=4., mfc=colorO)
        axs[axi][1].plot(Od.doxygen, Od.z, color=colorOd, marker='v', ms=4., mfc=colorOd)
        axt1[axi].plot(Aa.chlora, Aa.z, color=colorA, marker='o', ms=4., mfc=colorA)
        axt1[axi].plot(Ad.chlora, Ad.z, color=colorAd, marker='^', ms=4., mfc=colorAd)

        axs[axi][2].plot(Ba.backscatter, Ba.z, color=colorB, marker='s', ms=4., mfc=colorB)
        axs[axi][2].plot(Bd.backscatter, Bd.z, color=colorBd, marker='v', ms=4., mfc=colorBd)
        axt2[axi].plot(Ca.cdom, Ca.z, color=colorC, marker='o', ms=4., mfc=colorC)
        axt2[axi].plot(Cd.cdom, Cd.z, color=colorCd, marker='^', ms=4., mfc=colorCd)

        axs[axi][0].set(ylim=(-200., 0.))
        axs[axi][1].set(ylim=(-200., 0.))
        axs[axi][2].set(ylim=(-200., 0.))

        axs[axi][0].set(xlim=(temp_lo, temp_hi))
        axs[axi][1].set(xlim=(do_lo, do_hi))
        axs[axi][2].set(xlim=(bb_lo, bb_hi))


    axs[0][0].set(title='Temp (black) and Salinity (orange)')
    axs[0][1].set(title='Oxygen (blue) and Chlorophyll (green)')
    axs[0][2].set(title='CDOM (red) and Backscatter (cyan)')

    fig.show()

    # For additional labeling:
    # axs[iC][0].text(7.4, -14, 'S')
    # axs[iC][0].text(10.2, -14, 'T')
    # axs[iC][1].text(170, -30, 'Chl-A')
    # axs[iC][1].text(300, -150, 'DO')
    # axs[iC][2].text(.0007, -20, 'CDOM')
    # axs[iC][2].text(.0013, -75, 'SCATT')  
    
    return


##################
# more parameter configuration
##################
# Load the 2021 Oregon Slope Base profile metadata; and some March 2021 sensor datasets
##################

# profile metadata
p = ReadProfileMetadata(os.getcwd()+"/../Profiles/osb2021.csv")

# sensor data
A, B, C, T, S, O, H, I, N, P, U, V, W, R = ReadOSB_March2021_1min()

# Having loaded the data there are some artifacts to discard in O, T and :
O = O.drop('moles_of_oxygen_per_unit_mass_in_sea_water_profiler_depth_enabled_qc_agg')
T = T.drop('sea_water_temperature_profiler_depth_enabled_qc_agg')
S = S.drop('sea_water_practical_salinity_profiler_depth_enabled_qc_agg')

# ...and then apply the .dropna() to eliminate NaNs

A = A.dropna('time')
B = B.dropna('time')
C = C.dropna('time')
T = T.dropna('time')
S = S.dropna('time')
O = O.dropna('time')
H = H.dropna('time')
I = I.dropna('time')
N = N.dropna('time')
P = P.dropna('time')
U = U.dropna('time')
V = V.dropna('time')
W = W.dropna('time')
R = R.dropna('time')