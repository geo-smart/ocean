
# global sensor range parameters for charting data: based on osb shallow profiler data

# axis ranges for a variety of sensors

global_lo,    global_hi          =     0.0,        0.3

ba_lo, ba_hi = 0.0, 0.2
oa_lo, oa_hi = 0.15, 0.25


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

