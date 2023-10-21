import os, sys, time, glob, warnings
from os.path import join as joindir
import numpy as np, pandas as pd, xarray as xr
from numpy import datetime64 as dt64, timedelta64 as td64

warnings.filterwarnings('ignore')


def ReformatDataFile(verbose=False):
    """Read a NetCDF and reformat it, write the result"""
    print('\n\nSpecify input NetCDF data file\n')

    dataLoc         = os.getcwd() + '/../../../data/rca'                  # !!!!!!!!!!!!! adjust for geo-smart level
    osb             = 'OregonSlopeBase'
    oos             = 'OregonOffshore'
    axb             = 'AxialBase'
    sites_list      = [osb, oos, axb]
    plat            = 'platform'
    prof            = 'profiler'
    structures_list = [plat, prof]

    # m = int(input('Site choice: Enter an index 0 1 2 for ' + str(sites_list)))
    # n = int(input('Structure choice: Enter an index 0 1 for ' + str(structures_list)))
    m = 0
    n = 1                                                              # !!!!!!!!!!!!! override option: osb, profiler
    
    resource_folder = joindir(dataLoc, sites_list[m], structures_list[n])
    s = [name for name in os.listdir(resource_folder) if os.path.isdir(joindir(resource_folder, name))]
    print(s)   # This will give PAR, ctd, do, etcetera
    instrument_folder = joindir(resource_folder, s[int(input('Enter index 0, 1, ... to select the instrument: '))])
    l = os.listdir(instrument_folder)
    print(l)
    ds = xr.open_dataset(joindir(instrument_folder, l[int(input('Enter index of NetCDF file to use'))]))

    if verbose: print(ds)
    print()
    print()
    print()

    # user can choose to swap one dimension
    ds_dims      = [i for i in ds.dims]
    ds_coords    = [i for i in ds.coords]
    ds_data_vars = [i for i in ds.data_vars]

    print('\n\nSwap the dimension (e.g. time in place of row or obs)\n')
    print('Dimensions: ' + str(ds_dims))
    print('Coordinates: ' + str(ds_coords))
    print('Data Variables: ' + str(ds_data_vars))

    # old_dim = input('Dimension to swap out (need exact match):')
    old_dim = 'row'                                                    # !!!!!!!!!!!!!!!! hardcoded

    if old_dim in ds_dims:
        # new_dim = input('Coordinate or data variable to swap back in:')
        new_dim = 'time'                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! hardcoded
        
        if new_dim in ds_data_vars or new_dim in ds_coords:
            ds = ds.swap_dims({old_dim:new_dim})

    if verbose: print(ds)

    print('\n\nRename, Drop or Retain Coordinates\n')
    ds_coords = [i for i in ds.coords]
    nC = str(len(ds_coords))
    # print('For ' + nC + ' coordinates: 0 to drop, non-zero string to rename, enter to retain:\n')
    # for c in ds_coords:
    #     print('coord name: ' + c)
    #     a = input('Drop (0), Rename or Keep: ')
    #     if   a == '0': ds = ds.drop(c)
    #     elif len(a):   ds = ds.rename({c:a})
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! skipping this                         

    print('\n\nRename, Drop or Retain Data Variables\n')
    ds_data_vars = [i for i in ds.data_vars]
    nDV = str(len(ds_data_vars))
    print('For ' + nDV + ' data variables: 0 to drop, non-zero string to rename, enter to retain:\n')
    for dv in ds_data_vars:
        print('data variable name: ' + dv)
        a = input('Drop (0), Rename or Keep:')
        if   a == '0': ds = ds.drop(dv)
        elif len(a):   ds = ds.rename({dv:a})


    print('\n\nDrop all Attributes (less any you select)\n')
    ds_attrs_dict = ds.attrs.copy()
    # for k in ds_attrs_dict: print(k + ' '*(40-len(k)) + str(ds_attrs_dict[k]))
    attrs_to_preserve = []
    # print('Enter an attribute to preserve; or just enter by itself when done\n')
    # while True:
    #     s = input('Preserve: ')
    #     if not len(s): break
    #     attrs_to_preserve.append(s)
    for key in ds_attrs_dict: 
        if key not in attrs_to_preserve: ds.attrs.pop(key)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! streamline hardcode
        
    print('\n\nEnsure the new Dimension is sorted (no User action)\n')    
    df   = ds.to_dataframe()
    vals = [xr.DataArray(data=df[c], dims=['time'], coords={'time':df.index}, attrs=ds[c].attrs) for c in df.columns]
    ds  = xr.Dataset(dict(zip(df.columns, vals)), attrs=ds.attrs)

    
    print('\n\nSelect output time window (Format yyyy-mm-dd or enter to use the defaults)\n')
    t0_default, t1_default = '2021-07-01', '2021-08-01'             # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! hardcode
    # t0 = input('start date (' + t0_default + ')')
    # t1 = input('end date   (' + t1_default + ')')
    # if not len(t0): t0 = t0_default
    # if not len(t1): t1 = t1_default
    t0 = dt64(t0_default)
    t1 = dt64(t1_default)                           # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should be (t0), (t1)
    ds = ds.sel(time=slice(t0, t1))
    
    # This code eliminates duplicate-time entries. See data.ipynb for remarks on 
    #   sensor data.
    _, keeper_index = np.unique(ds['time'], return_index=True)
    ds=ds.isel(time=keeper_index)

    print('\n\nHere is the resulting dataset summary view:\n')
    print(ds)
    
    # bug: this plot needs plt.show()-style compulsion; it holds back
    # print("\n\ndepth quick look (assumes a 'z' coord/data variable:\n")
    # ds.z[0:10000].plot()

    outfnm = input('\n\nEnter an output file name. Include the .nc extension (or just enter to skip this): ')
    if len(outfnm): ds.to_netcdf(outfnm)    

    return True



















def ProfCrawler(z, t, verbose = False):
    """
    ProfileCrawler traverses pandas Series s of pressures/depths and matching pandas Series t of times.
    The code is adapted to 1Min per sample and z is increasing up / shallow. (See 'z_direction'.)
    The code returns a set of 6 lists: t0/t1 for ascent, descent and rest intervals
    """
    print(str(z[0]) + ' is initial depth')

    len_z = len(z)
    a0, d0, r0 = [], [], []               # lists for start times: ascents, descents, rests
    
    r0.append((0,t[0],z[0]))

    # Goal is to return 6 lists: r0, r1, a0, a1, d0, d1
    # List entries are triples (i, t, z): Index of, time of, and depth of.
    # The first event is presumed to be a rest r0 to r1. r1 coincides with first ascent a0.
    #   The first possible i is m0; the last is len_z - m1 - 1. len_z is the number of
    #   values in the time series; and we truncate by m0/m1 to permit time window derivatives.
    # Note that z[m0] is m0 minutes ahead of z[0] since samples are assumed to be one per minute.
    #   Likewise z[len_z - m1 - 1] is m1 minutes behind z[len_z - 1].
    # The depth values z[] are taken to be negative down from 0 at the surface.
    # Time-series derivatives are taken in terms of (later) - (earlier) in the normal sense.
    # Conditions for identifying an Ascent start:
    #   - Slope from past to present is less than a threshold (flat)
    #   - Slope from present to future is positive (rising)
    #   - Profiler depth < threshold (eliminates some false positives)
    # And similarly for Descent start and Rest start.
    #   After a detection of ascent start the i search index is bumped forward in time to avoid
    #   subsequent false positives.
    
    m0 = 8
    m1 = 8
    ascent_threshold0 = 0.2
    ascent_threshold1 = .5
    ascent_min_depth = -170.
    ascent_bump_i = 10        # 7 "works" but a bit bigger is maybe no harm
    descent_threshold0 = -0.2
    descent_threshold1 = -0.5
    descent_bump_i = 10
    rest_threshold0 = -0.5                 # -.5, .2, -170 worked pretty well for r0
    rest_threshold1 = 0.2             
    rest_min_depth = -170.
    rest_bump_i = 10

    # This pushes the index forward after an a0 is detected to skip some false positives
    ascent_threshold0 = 0.2
    ascent_threshold1 = .5
    
    i = m0
    while i < len_z - m1:      # i is a candidate index for A/D/R starts
        
        slope0 = (z[i] - z[i-m0])/m0
        slope1 = (z[i+m1] - z[i])/m1
        
        if slope0 <= ascent_threshold0 and slope1 >= ascent_threshold1:
            if z[i] <= ascent_min_depth:
                a0.append((i, t[i], z[i]))
                i += ascent_bump_i
                
        elif slope0 >= descent_threshold0 and slope1 <= descent_threshold1:        
            # (no depth condition) so on to the cases where this is considered a new descent:
            #   No descents found yet OR the most recent descent precedes the most recent ascent
            if (not len(d0)) or (len(a0) and len(d0) and d0[-1][0] < a0[-1][0]):
                d0.append((i, t[i], z[i]))
                i += descent_bump_i
                
        elif slope0 <= rest_threshold0 and abs(slope1) <= rest_threshold1:
            if z[i] <= rest_min_depth:
                if (not len(r0)) or (len(d0) and len(r0) and r0[-1][0] < d0[-1][0]):
                    r0.append((i, t[i], z[i]))
                    i += rest_bump_i

        i += 1

    if verbose: print("there are", len(a0), "ascent starts")
    if verbose: print("there are", len(d0), "descent starts")
    if verbose: print("there are", len(r0), "rest starts")
    
    a1 = d0.copy()             # ascent end = descent start
    d1 = r0[1:].copy()         # first descent ends at start of 2nd rest start
    r1 = a0.copy()             # first rest end = first ascent start

    # The final d1 must still be determined
    i = d0[-1][0] + m0
    while i < len_z - m1:
        slope0 = (z[i] - z[i-m0])/m0
        slope1 = (z[i+m1] - z[i])/m1
        if slope0 <= rest_threshold0 and abs(slope1) <= rest_threshold1:
            if z[i] <= rest_min_depth:
                d1.append((i, t[i], z[i]))
                i = len_z - m1
        i += 1
                
    # redacted: logic check on order of stamp indices
    # Returning lists of tuples: (index, time, depth)
    return a0, a1, d0, d1, r0, r1  



def ProfileGenerator(sourcefnm, s, z, verbose=True):
    """
    Generate a profile CSV file for a site and a time interval implicit in a depth/time NetCDF file 
    Example: result = ProfileWriter('source.nc', 'axb', 'z', True)
      sourcefnm is a NetCDF file containing pressure/depth and time
      s is a site label string
      y0, yN give an inclusive year range
    """
    ds = xr.open_dataset(sourcefnm)

    # if dim not 'time' return False
    
    t0 = ds['time'][0]
    t1 = ds['time'][-1]
    
    a0, a1, d0, d1, r0, r1 = ProfCrawler(ds[z].to_series(), ds['time'].to_series(), True)

    return a0, a1, d0, d1, r0, r1



def ProfileWriter(ofnm, a0, a1, d0, d1, r0, r1):

    print('a0: ' + str(len(a0)) + '    a1: ' + str(len(a1)))
    print('d0: ' + str(len(d0)) + '    d1: ' + str(len(d1)))
    print('r0: ' + str(len(r0)) + '    r1: ' + str(len(r1)))

    if not len(a0) == len(a1):     return False
    if not len(a1) == len(d0):     return False
    if not len(d0) == len(d1):     return False
    if not len(d1) == len(r0):     return False
    if not len(r0) == len(r1):     return False

    profiles = []
    for i in range(len(r0)):
        profiles.append([r0[i][0],r0[i][1],r0[i][2],r1[i][0],r1[i][1],r1[i][2],a0[i][0],a0[i][1],a0[i][2],a1[i][0],a1[i][1],a1[i][2],d0[i][0],d0[i][1],d0[i][2],d1[i][0],d1[i][1],d1[i][2]])

    df = pd.DataFrame(data=np.array([np.array(x) for x in profiles]))
    df.to_csv(ofnm)

    return True