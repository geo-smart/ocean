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

print('\nJupyter Notebook running Python {}'.format(sys.version_info[0]))


def ReadProfileMetadata(fnm):
    """
    Profiles are saved in a CSV file as six events per row: Rest start, Rest end, Ascent start,
    Ascent end, Descent start, Descent end. Each event includes a time and depth. There is event
    degeneracy in that for a given row the Rest end is the same event as Ascent start. Likewise
    Ascent end is also Descent start. Descent end is Rest start *in the subsequent row*. The 
    exception is of course in the final row. Event depth is measured as a negative value below
    the zero which is the sea surface.
    """
    df = pd.read_csv(fnm, usecols=["1","2","4","5","7","8","10","11","13","14","16","17"])
    df.columns=['r0t','r0z','r1t','r1z','a0t','a0z','a1t','a1z','d0t','d0z','d1t','d1z']
    df['r0t'] = pd.to_datetime(df['r0t'])
    df['r1t'] = pd.to_datetime(df['r1t'])
    df['a0t'] = pd.to_datetime(df['a0t'])
    df['a1t'] = pd.to_datetime(df['a1t'])
    df['d0t'] = pd.to_datetime(df['d0t'])
    df['d1t'] = pd.to_datetime(df['d1t'])
    return df



sensors = [
['conductivity', 'ctd'], ['density', 'ctd'], ['pressure', 'ctd'], ['salinity', 'ctd'], ['temperature', 'ctd'],
['chlora', 'fluor'], ['bb', 'fluor'], ['fdom', 'fluor'],
['spkir', 'spkir'],
['nitrate', 'nitrate'],
['pco2', 'pco2'],
['do', 'do'],
['par', 'par'],
['ph', 'ph'],
['up', 'vel'], ['east', 'vel'], ['north', 'vel']]

ranges = {
'conductivity':(33.5,36.5),'density':(1022, 1028),'pressure':(0.,200.),'salinity':(31, 35),'temperature':(7, 14),
'chlora':(0.,1.5),'bb':(0.00,0.006),'fdom':(0.5,4.5),
'spkir':(0.0, 15.0),
'nitrate':(0., 35.),
'pco2':(200.0, 1200.0),
'do':(50.0, 400.),
'par':(0.0, 300.),
'ph':(7.6, 8.2),
'up':(-0.4, 0.4),'east':(-0.4, 0.4),'north':(-0.4, 0.4)
}

standard_deviations = {
'conductivity':(0.1, 0.6),'density':(0., .3),'pressure':(0.,10.),'salinity':(.0, .4),'temperature':(.0, .7),
'chlora':(0.0, 0.5),'bb':(0.0,0.003),'fdom':(0.0,0.7),
'spkir':(0.0, .5),
'nitrate':(0., 4.),
'pco2':(0.0, 10.0),
'do':(0.0, 40.),
'par':(0.0, 30.),
'ph':(0., 0.2),
'up':(0., 0.1),'east':(0, 0.1),'north':(0., 0.1)
}

colors = {
'conductivity':'xkcd:maroon','density':'xkcd:brick red','pressure':'xkcd:eggplant','salinity':'cyan','temperature':'red',
'chlora':'green','bb':'xkcd:blood orange','fdom':'xkcd:olive drab',
'spkir':'black',
'nitrate':'black',
'pco2':'black',
'do':'blue',
'par':'red',
'ph':'yellow',
'up':'red','east':'green','north':'blue'
}

sensor_names = {
'conductivity':'Conductivity','density':'Density (kg m-3)','pressure':'Pressure','salinity':'Salinity','temperature':'Temperature (C)',
'chlora':'Chlorophyll-A','bb':'Particulate Backscatter','fdom':'Fluorescent DOM',
'spkir':'Spectral Irradiance',
'nitrate':'Nitrate Concentration',
'pco2':'CO2 Concentration',
'do':'Dissolved Oxygen',
'par':'Photosynthetically Available Radiation',
'ph':'pH',
'up':'Current: Vertical','east':'Current: East','north':'Current: North'
}



def DataFnm(site, instrument, time, sensor): 
    datafnm = './../data/' + site + '_' + instrument + '_' + time + '_' + sensor + '.nc'
    return datafnm
