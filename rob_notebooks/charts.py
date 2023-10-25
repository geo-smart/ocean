import os, sys, time, glob, warnings
from os.path import join as joindir
from IPython.display import clear_output
from matplotlib import pyplot as plt
from matplotlib import colors as mplcolors
import numpy as np, pandas as pd, xarray as xr
from numpy import datetime64 as dt64, timedelta64 as td64

warnings.filterwarnings('ignore')

def doy(theDatetime): return 1 + int((theDatetime - dt64(str(theDatetime)[0:4] + '-01-01')) / td64(1, 'D'))
def dt64_from_doy(year, doy): return dt64(str(year) + '-01-01') + td64(doy-1, 'D')
def day_of_month_to_string(d): return str(d) if d > 9 else '0' + str(d)



def ChartTwoSensors(p, xrng, pidcs, A, Az, Albl, Acolor, Aleg, \
                                    B, Bz, Blbl, Bcolor, Bleg, \
                                    wid, hgt, z0=-200., z1=0.):
    """
    Make a stack of charts with two horizontal axes to compare two sensors A and B.
    The data are in DataArrays: A, Az, B, Bz. pidcs are row indices for the profile
    Dataframe, i.e. a choice of which profiles to plot by index. 
    
    p        pandas Dataframe of indexed profile timestamps
    xrng     list of 2-lists: low-to-high values for the two sensors
    pIdcs    indices within p to use in generating a sequence of paired charts
    A        xarray Dataset: source data of type A (B)
    Az       xarray Dataset: depth data for sensor A (B)
    Albl     string: label for sensor A (B)
    Acolor   string: color for sensor A (B)
    Aleg     'rest', 'ascent' or 'descent'
    wid      width for two charts
    hgt      height for one chart (scaled by number of charts)
    z0, z1   depth range
    """
    
    # empirical values for the day's two longer-duration profiles
    midn0 = td64( 7*60 + 10, 'm')        # 7 hours 10 minutes
    midn1 = td64( 7*60 + 34, 'm')        # 7 hours 34 minutes
    noon0 = td64(20*60 + 30, 'm')        # 20 hours 30 minutes
    noon1 = td64(20*60 + 54, 'm')        # 20 hours 54 minutes 
        
    # limit the number of charts to 100
    ncharts = len(pidcs)
    if ncharts > 100: ncharts = 100
    do_one = True if ncharts == 1 else False
    print("Attempting", ncharts, "charts\n")

    # ncharts x 1: charts in a vertical column
    fig, axs = plt.subplots(ncharts, 1, figsize=(wid, hgt*ncharts), tight_layout=True)
    
    # list of twin axes for sensor B (one for each chart)
    axstwin = axs.twiny() if do_one else [axs[i].twiny() for i in range(ncharts)]

    # profile table p has column headers 'a0z', 'a0t' and so on for r and d
    #   We are interested in the time columns to constrain data selection
    if   Aleg == 'rest':   keyA = ('r0t', 'r1t')
    elif Aleg == 'ascent': keyA = ('a0t', 'a1t')
    else:                  keyA = ('d0t', 'd1t')
    if   Bleg == 'rest':   keyB = ('r0t', 'r1t')
    elif Bleg == 'ascent': keyB = ('a0t', 'a1t')
    else:                  keyB = ('d0t', 'd1t')
  
    # The subsequent code is 'loop over charts: plot each chart, A and B'
    # For this we need both a profile index into the profile dataframe p (from the
    #   passed list pidcs[] *and* a chart index 0, 1, 2, ... These are respectively 
    #   pidx and i.
    for i in range(ncharts):
        
        pidx = pidcs[i]

        tA0, tA1 = p[keyA[0]][pidx], p[keyA[1]][pidx]
        tB0, tB1 = p[keyB[0]][pidx], p[keyB[1]][pidx]
        
        Ax, Ay = A.sel(time=slice(tA0,  tA1)), Az.sel(time=slice(tA0, tA1))
        Bx, By = B.sel(time=slice(tB0,  tB1)), Bz.sel(time=slice(tB0, tB1))
        
        if do_one:
            axs.plot(    Ax, Ay, ms = 4., color=Acolor, mfc=Acolor)
            axstwin.plot(Bx, By, ms = 4., color=Bcolor, mfc=Bcolor)
        else:
            axs[i].plot(    Ax, Ay, ms = 4., color=Acolor, mfc=Acolor)
            axstwin[i].plot(Bx, By, ms = 4., color=Bcolor, mfc=Bcolor)
        
        # axis ranges
        if i == 0: 
            if do_one:
                axs.set(title = Albl + ' (' + Acolor + ', lower x-axis) and ' \
                              + Blbl + ' (' + Bcolor + ', upper x-axis)')
            else:
                axs[i].set(title = Albl + ' (' + Acolor + ', lower x-axis) and ' \
                                 + Blbl + ' (' + Bcolor + ', upper x-axis)')

        # Set axis ranges from passed list of pairs xrng[][]
        if do_one:
            axs.set(    xlim = (xrng[0][0], xrng[0][1]), ylim = (z0, z1))
            axstwin.set(xlim = (xrng[1][0], xrng[1][1]), ylim = (z0, z1))
        else:
            axs[i].set(    xlim = (xrng[0][0], xrng[0][1]), ylim = (z0, z1))
            axstwin[i].set(xlim = (xrng[1][0], xrng[1][1]), ylim = (z0, z1))

        # chart time label
        ascent_start_time = 'Start UTC: ' + str(tA0)
        delta_t = tA0-dt64(tA0.date())
        if delta_t > midn0 and delta_t < midn1: ascent_start_time += " MIDNIGHT local"
        if delta_t > noon0 and delta_t < noon1: ascent_start_time += " NOON local"
        xlabel = xrng[0][0] + 0.2*(xrng[0][1] - xrng[0][0])
        ylabel = -10
        if do_one: axs.text(xlabel, ylabel, ascent_start_time)
        else: axs[i].text(xlabel, ylabel, ascent_start_time)
        
    return fig, axs