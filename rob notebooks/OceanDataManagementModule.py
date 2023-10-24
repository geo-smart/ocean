# This code concerns translating large volume datasets into smaller time-range datasets.
#
#   These data files for Oregon Slope Base were pulled from the OOI "data explorer" and stored 
#     in the data directory outside this repository. The data are timestamped at one minute intervals.
#   As each dataset is opened it also has its difficult parameter changed to one more manageable.
#

if False: 
    data_source = '../data/data_explorer_1Min/'

    dsA = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_chlora_1Min.nc')
    dsB = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_backscatter_1Min.nc')
    dsC = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_cdom_1Min.nc')
    dsO = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_doxygen_1Min.nc')
    dsS = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_salinity_1Min.nc')
    dsT = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_temperature_1Min.nc')

    dsA = dsA.rename_vars({"mass_concentration_of_chlorophyll_a_in_sea_water_profiler_depth_enabled":"chlora"})
    dsB = dsB.rename_vars({"flubsct_profiler_depth_enabled":"backscatter"})
    dsC = dsC.rename_vars({"cdomflo_profiler_depth_enabled":"cdom"})
    dsO = dsO.rename_vars({"moles_of_oxygen_per_unit_mass_in_sea_water_profiler_depth_enabled":"doxygen"})
    dsS = dsS.rename_vars({"sea_water_practical_salinity_profiler_depth_enabled":"salinity"})
    dsT = dsT.rename_vars({"sea_water_temperature_profiler_depth_enabled":"temp"})

    # pCO2: no data
    dsG = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_pco2_1Min.nc')
    dsG = dsG.rename_vars({"partial_pressure_of_carbon_dioxide_in_sea_water_profiler_depth_enabled":"pco2"})

    dsH = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_ph_1Min.nc')
    dsI = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_spkir_1Min.nc')
    dsN = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_nitrate_1Min.nc')
    dsP = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_par_1Min.nc')
    dsU = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_veleast_1Min.nc')
    dsV = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_velnorth_1Min.nc')
    dsW = xr.open_dataset(data_source + 'osb/profiler/osb_profiler_velup_1Min.nc')

    dsH = dsH.rename_vars({"sea_water_ph_reported_on_total_scale_profiler_depth_enabled":"ph"})
    dsI = dsI.rename_vars({"spectir_412nm":"si412",  \
                           "spectir_443nm":"si443", \
                           "spectir_490nm":"si490", \
                           "spectir_510nm":"si510", \
                           "spectir_555nm":"si555", \
                           "spectir_620nm":"si620", \
                           "spectir_683nm":"si683"})
    dsN = dsN.rename_vars({"mole_concentration_of_nitrate_in_sea_water_profiler_depth_enabled":"nitrate"})
    dsP = dsP.rename_vars({"downwelling_photosynthetic_photon_flux_in_sea_water_profiler_depth_enabled":"par"})
    dsU = dsU.rename_vars({"eastward_sea_water_velocity_profiler_depth_enabled":"veast"})
    dsV = dsV.rename_vars({"northward_sea_water_velocity_profiler_depth_enabled":"vnorth"})
    dsW = dsW.rename_vars({"upward_sea_water_velocity_profiler_depth_enabled":"vup"})

# Subset 1min data from above to smaller time intervals ranges. Results fit within this repo.
#   month-blocks: All dates use first of the month
#   short variable names to fit everything in one line per sensor type

def WriteMonthLongDataBlocks(s, y, z, m, n, u, v, wd):
    dsT.sel(time=slice(u, v)).to_netcdf(wd + 'ctd/'     + s + '_temp_'        + p + y + '_1min.nc')
    dsS.sel(time=slice(u, v)).to_netcdf(wd + 'ctd/'     + s + '_salinity_'    + p + y + '_1min.nc')
    dsO.sel(time=slice(u, v)).to_netcdf(wd + 'ctd/'     + s + '_doxygen_'     + p + y + '_1min.nc')
    dsA.sel(time=slice(u, v)).to_netcdf(wd + 'fluor/'   + s + '_chlora_'      + p + y + '_1min.nc')
    dsB.sel(time=slice(u, v)).to_netcdf(wd + 'fluor/'   + s + '_backscatter_' + p + y + '_1min.nc')
    dsC.sel(time=slice(u, v)).to_netcdf(wd + 'fluor/'   + s + '_cdom_'        + p + y + '_1min.nc')

    dsH.sel(time=slice(u, v)).to_netcdf(wd + 'pH/'      + s + '_ph_'          + p + y + '_1min.nc')
    dsI.sel(time=slice(u, v)).to_netcdf(wd + 'irrad/'   + s + '_spectir_'     + p + y + '_1min.nc')
    dsN.sel(time=slice(u, v)).to_netcdf(wd + 'nitrate/' + s + '_nitrate_'     + p + y + '_1min.nc')
    dsP.sel(time=slice(u, v)).to_netcdf(wd + 'par/'     + s + '_par_'         + p + y + '_1min.nc')
    dsU.sel(time=slice(u, v)).to_netcdf(wd + 'current/' + s + '_veast_'       + p + y + '_1min.nc')
    dsV.sel(time=slice(u, v)).to_netcdf(wd + 'current/' + s + '_vnorth_'      + p + y + '_1min.nc')
    dsW.sel(time=slice(u, v)).to_netcdf(wd + 'current/' + s + '_vup_'         + p + y + '_1min.nc')

    return

# This March 2021 OSB dataset is used by the subsequent visual 'tour' of sensor types
# s = 'osb'
# y, z = '2021', '2021'
# m, n = '03', '04'
# p = 'march'
# u = dt64(y + '-' + m + '-01')
# v = dt64(z + '-' + n + '-01')
# wd = os.getcwd() + '/RepositoryData/rca/'

# WriteMonthLongDataBlocks(s, y, z, m, n, u, v, wd)

# This 2018 osb segment is for comparison to discrete summary data
# s = 'osb'
# y, z = '2018', '2018'
# m, n = '06', '08'
# p = 'june_july'
# u = dt64(y + '-' + m + '-01')
# v = dt64(z + '-' + n + '-01')
# wd = os.getcwd() + '/RepositoryData/rca/'

# WriteMonthLongDataBlocks(s, y, z, m, n, u, v, wd)




# GLODAP Data Loader
# Requires boto + target directory has write permission
if False:         # disabled once the datasets are loaded into /data/glodap

    glodapTemperatureFnm = data_dir + '/glodap/glodap_temperature.nc'
    glodapSalinityFnm    = data_dir + '/glodap/glodap_salinity.nc'
    glodapOxygenFnm      = data_dir + '/glodap/glodap_oxygen.nc'

    import boto
    from boto.s3.key import Key

    connection = boto.connect_s3(anon=True)
    bucket = connection.get_bucket('fixthisshouldhavesecurebucketnamehere')

    for key in bucket.list(): 
        filename = key.name.encode('utf-8')
        if b'glodap' in filename: 
            if b'salinity.nc' in filename: 
                print ('salinity file is', filename)
                salinityfilename = filename
            if b'temperature.nc' in filename: 
                print ('temperature file is', filename)
                temperaturefilename = filename
            if b'oxygen.nc' in filename: 
                print('oxygen file is', filename)
                oxygenfilename = filename            

    k = Key(bucket)
    k.key = salinityfilename
    k.get_contents_to_filename(glodapSalinityFnm)
    k.key = temperaturefilename
    k.get_contents_to_filename(glodapTemperatureFnm)
    k.key = oxygenfilename
    k.get_contents_to_filename(glodapOxygenFnm)

    print('\ndata load complete for glodap')