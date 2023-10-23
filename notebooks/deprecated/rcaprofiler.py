def ProfileCrawler(s, t, verbose = False):
    """
    ProfileCrawler traverses a passed pandas Series s of pressures and Series t of corresponding times.
    The code was built using data sampled at about 1-minute intervals. Goal: Infer (and return as lists)
        the start and end times for ascent, descent and rest intervals.
    """

    # pandas series, just pressures
    len_s = len(s)
    threshold = 1.
    a0, d0, r0 = [], [], []               # start times for ascents, descents, rests

    for i in range(1, len_s - 5):         # 6 minute window
        
        # catch ascent
        if s[i-1] - s[i]   <= threshold and  \
           s[i]   - s[i+1] >= threshold and  \
           s[i+1] - s[i+2] >= threshold and  \
           s[i+2] - s[i+3] >= threshold and  \
           s[i+3] - s[i+4] >= threshold and  \
           s[i+4] - s[i+5] >= threshold:
            a0.append((i,t[i]))

        # catch descent
        if s[i-1] - s[i]   >= threshold and  \
           s[i]   - s[i+1] <= threshold and  \
           s[i+1] - s[i+2] <= threshold and  \
           s[i+2] - s[i+3] <= threshold and  \
           s[i+3] - s[i+4] <= threshold and  \
           s[i+4] - s[i+5] <= threshold:
            d0.append((i,t[i]))

        # this variant is a little too liberal; false positives ~25%
        # why? Because twice daily there are stops on the way down for pH
        # catch rest
        if i >= 5 and \
           s[i-5] - s[i-4] <= -threshold and  \
           s[i-4] - s[i-3] <= -threshold and  \
           s[i-3] - s[i-2] <= -threshold and  \
           s[i-2] - s[i-1] <= -threshold and  \
           s[i-1] - s[i]   <= -threshold and  \
           s[i]   - s[i+1] >= -threshold:
            r0.append((i,t[i]))

            
    if verbose: print("there are", len(a0), "ascent starts")
    if verbose: print("there are", len(d0), "descent starts")
    if verbose: print("there are", len(r0), "rest starts: now culling extras")

    # keep running the "rest start" list looking for...
    #     ...relative to this particular rest start...
    #     ...a future rest start earlier than the next ascent start
    #         ...if found: delete this particular rest start...
    #             ...and begin again at the top
    # this winnows away false positive rest starts
    while True: 
        profile_counter = 1               
        for i in range(len(r0)-1):
            if profile_counter >= len(a0) - 1: break
            if r0[i+1][1] < a0[profile_counter][1]: 
                r0.remove(r0[i])
                break
            else: profile_counter += 1
        if profile_counter >= len(a0) - 1: break


    if verbose: print("there are", len(a0), "ascent starts")
    if verbose: print("there are", len(d0), "descent starts")
    if verbose: print("there are", len(r0), "rest starts")
                   
    a1 = d0.copy()             # ascent end = descent start
    d1 = r0.copy()             # descent end = rest start
    r1 = a0[1:].copy()         # rest end = next ascent start (cuts off end of year)
       
    # logic check on results
    causal_errors = [0]*3                 # list [0, 0, 0] tracks 3 types of possible error
    smallest = len(a0)
    if len(d0) < smallest: smallest = len(d0)
    if len(r0) < smallest: smallest = len(r0)
    for i in range(smallest):
        if a0[i][0] >= d0[i][0]: causal_errors[0] += 1    # ascent start later than descent start
        if d0[i][0] >= r0[i][0]: causal_errors[1] += 1    # descent start later than rest start
        if a1[i][0] >= d1[i][0]: causal_errors[2] += 1    # ascent end later than descent end (???)
   
    if verbose: print("causal error counts:", causal_errors[0], causal_errors[1], causal_errors[2])
    if verbose: print(len(a0), len(d0), len(r0))

    # Returning lists of tuples: (index, time)
    return a0, a1, d0, d1, r0, r1



def PrintProfileStatistics(a0, a1, d0, d1, r0, r1):
    """ 
    PrintProfileStatistics prints mean and standard deviation for a set of profiles.
    Specifically for a set of Ascents, Descents and Rests. Each passed vector (a0 etc) 
    is a list of tuples. The first element of the tuple is the index of the time in the 
    source data array. The second value is the timestamp for that same element. 
    """
    one_sec = np.timedelta64(1, 's')
    D_asc  = [(dt64(a1[i][1])-dt64(a0[i][1]))/one_sec for i in range(len(a1))]
    D_dsc  = [(dt64(d1[i][1])-dt64(d0[i][1]))/one_sec for i in range(len(d1))]
    D_rst  = [(dt64(r1[i][1])-dt64(r0[i][1]))/one_sec for i in range(len(r1))]

    print('Means, standard deviation for profile phases, in minutes:')
    print('  Ascents:  ', round(np.mean(D_asc)/60., 2), round(np.std(D_asc)/60., 2))
    print('  Descents: ', round(np.mean(D_dsc)/60., 2), round(np.std(D_dsc)/60., 2))
    print('  Rests:    ', round(np.mean(D_rst)/60., 2), round(np.std(D_rst)/60., 2))
    print()
    print('(Recall that two profiles of nine each day have slower, staged descents)')
    print()
    

def PrintProfileEntry(profile):
    """
    A profile is a list of 12 values as six pairs of (index, timestamp) interleaved values
      ascent start:  index, timestamp
      ascent end:    index, timestamp 
      descent start: index, timestamp 
      descent end:   index, timestamp
      rest start:    index, timestamp
      rest end:      index, timestamp
    The indices refer back to the source dataset, likely at 1Min samples. They could be abolished.
    The file that is written from these pre-pends a column of profile indices 0, 1, 2, ...; Hence
    there are 13 columns total.
    """
    print("ascent:  index / start time:", profile[0], profile[1], '       index / end time:', profile[2], profile[3])
    print("descent: index / start time:", profile[4], profile[5], '       index / end time:', profile[6], profile[7])
    print("rest:    index / start time:", profile[8], profile[9], '       index / end time:', profile[10], profile[11])

    
    


    
    
def ProfileWriter(sourcefnm, s, y0, yN, verbose=True):
    """
    Generate Profile CSV files for sites x years. 
    Example: result = ProfileWriter('source.nc', 'axb', 2015, 2021)
      sourcefnm is a NetCDF file containing pressure/depth and time
      s is a list of sites (strings)
      y0, yN give an inclusive year range
    """
    ds = xr.open_dataset(sourcefnm)
    for yr in range(y0, yN+1):
        yrstr   = str(yr)
        yrpostr = str(yr+1)
        dsyr    = ds.sel(time=slice(dt64(yrstr + '-01-01'), dt64(yrpostr + '-01-01')))

        
        # LEFT OFF HERE PARAMETERIZE pressure i.e. dsyr[pname]


        a0, a1, d0, d1, r0, r1 = \
            ProfileCrawler(dsyr.sea_water_pressure_profiler_depth_enabled.to_series(), \
                           dsyr.time.to_series(), True)

        print(len(a0), len(d0), len(r0), 'interval starts')
        print(len(a1), len(d1), len(r1), 'interval ends')


        if len(a0) < 10 or len(a1) < 10 or len(d0) < 10 or len(d1) < 10 or len(r0) < 10 or len(r1) < 10:
            print()
            print('No data: Abandoning this site + year:', site, yrstr)
            print()
        else:

            # we have intervals; do they match? Assume not always. Here is a checking function:
            # CompareShallowProfilerTimestamps(a0, a1, d0, d1, r0, r1)

            ascents, descents, rests = [], [], []
            day_td64                 = pd.to_timedelta(1, unit='D')
            ascent_limit             = pd.to_timedelta(2, unit='H')
            descent_limit            = pd.to_timedelta(2, unit='H')
            rest_limit               = pd.to_timedelta(2, unit='H')
            prior_ascent_start       = a0[0][1] - day_td64
            prior_descent_start      = d0[0][1] - day_td64
            prior_rest_start         = r0[0][1] - day_td64

            end_index = 0     # index into a1
            for i in range(len(a0)):
                all_done = False
                this_start_time = a0[i][1]
                if this_start_time > prior_ascent_start:
                    while a1[end_index][1] <= this_start_time: 
                        end_index += 1
                        if end_index >= len(a1):
                            all_done = True
                            break
                    if all_done: break
                    this_end_time = a1[end_index][1]
                    if this_end_time < this_start_time + ascent_limit:
                        prior_ascent_start = this_start_time
                        ascents.append([a0[i][0], this_start_time, a1[end_index][0], this_end_time])
                if all_done: break

            end_index = 0     # index into d1
            for i in range(len(d0)):
                all_done = False
                this_start_time = d0[i][1]
                if this_start_time > prior_descent_start:
                    while d1[end_index][1] <= this_start_time: 
                        end_index += 1
                        if end_index >= len(d1):
                            all_done = True
                            break
                    if all_done: break
                    this_end_time = d1[end_index][1]
                    if this_end_time < this_start_time + descent_limit:
                        prior_descent_start = this_start_time
                        descents.append([d0[i][0], this_start_time, d1[end_index][0], this_end_time])
                if all_done: break


            end_index = 0     # index into r1
            for i in range(len(r0)):
                all_done = False
                this_start_time = r0[i][1]
                if this_start_time > prior_rest_start:
                    while r1[end_index][1] <= this_start_time: 
                        end_index += 1
                        if end_index >= len(r1):
                            all_done = True
                            break
                    if all_done: break
                    this_end_time = r1[end_index][1]
                    if this_end_time < this_start_time + rest_limit:
                        prior_rest_start = this_start_time
                        rests.append([r0[i][0], this_start_time, r1[end_index][0], this_end_time])
                if all_done: break

            print("found", len(ascents), 'good ascents')
            print("found", len(descents), 'good descents')
            print("found", len(rests), 'good rests')

            # profiles[] will be a list of clean ascend/descend/rest sequences, 12 numbers per sequence
            #   ascend start:  index, timestamp        
            #   ascend end:    index, timestamp      The 'index' refers to the source dataset, typically at "1Min"
            #   descend start: index, timestamp      sampling rate. Note that ascend end = descend start and so on.
            #   descend end:   index, timestamp
            #   rest start:    index, timestamp
            #   rest end:      index, timestamp
            profiles           = []
            descent_index      = 0
            rest_index         = 0

            # This code builds the profiles[] list
            all_done = False
            for i in range(len(ascents)):
                all_done = False
                this_end_ascent_time = ascents[i][3]
                found_matching_descent = False
                while descents[descent_index][1] < this_end_ascent_time:
                    descent_index += 1
                    if descent_index >= len(descents):
                        all_done = True
                        break
                if all_done: break
                if descents[descent_index][1] == ascents[i][3]:
                    this_end_descent_time = descents[descent_index][3]
                    while rests[rest_index][1] < this_end_descent_time:
                        rest_index += 1
                        if rest_index >= len(rests):
                            all_done = True
                            break
                    if all_done: break
                    if rests[rest_index][1] == descents[descent_index][3]:
                        di = descent_index
                        ri = rest_index
                        profiles.append([\
                                         ascents[i][0],   ascents[i][1],   ascents[i][2],   ascents[i][3],   \
                                         descents[di][0], descents[di][1], descents[di][2], descents[di][3], \
                                         rests[ri][0],    rests[ri][1],    rests[ri][2],    rests[ri][3]     \
                                        ])


            # This code removes profiles whose start time is earlier than the prior profile rest end time
            #   This happens when multiple ascend starts are detected for a single actual ascent. It can
            #   result in more than nine profiles per day which is in general unlikely. 
            nTimeSlipsRemoved = 0
            while True: 
                fall_out = True
                for i in range(1, len(profiles)):
                    if profiles[i][1] < profiles[i-1][11]:
                        profiles.remove(profiles[i])
                        nTimeSlipsRemoved += 1
                        fall_out = False
                        break
                if fall_out: break

            # This code looks for and reports on duplicated profile ascent start times
            double_check, fail_index = True, -1
            for i in range(len(profiles)-1):
                if profiles[i][1] == profiles[i+1][1]: 
                    double_check = False
                    fail_index = i
                    break

            if not double_check: PrintProfileEntry(profiles[fail_index])
            else: print('no doubling of profile ascent starts found')

            # This code looks for and reports on non-matching Timestamp sequences:
            #   From ascent to descent and descent to rest.
            double_check, fail_index = True, -1
            for i in range(len(profiles)):
                if profiles[i][3] != profiles[i][5] or profiles[i][7] != profiles[i][9]: 
                    double_check = False
                    fail_index = i
                    break

            # This code compiles a histogram of profiles by doy and it has three faults to be aware of
            #   - Baked in is the assumption that this is at most one year of data
            #   - There is capacity for a leap year with 366 days but it is not explicitly sorted out
            #   - Day of year (doy) usually numbers from 1 but the histogram numbers 
            profile_histogram = [0]*366
            doylist = list(range(366))
            for i in range(len(profiles)):
                profile_histogram[doy(profiles[i][1])-1] += 1

            # This code counts how many days had nine profiles as expected, and how many had more
            #   than nine profiles which is not really possible. So that would indicate false
            #   positives still got through the process here.
            nNines = 0
            nMoreThanNine = 0
            more_than_nine = []
            for i in range(366):
                if profile_histogram[i] == 9: nNines = nNines + 1
                if profile_histogram[i] > 9:  more_than_nine.append(i)


            # print diagnostics from all of the above steps
            print("arrived at", len(profiles), 'good candidate profiles')
            print("after removing", nTimeSlipsRemoved, 'due to time slip error')
            if double_check: print('transitions are self-consistent')
            else: print('double check failed at element', fail_index)
            print('of 365 days,', nNines, 'have nine profiles as desired')
            print('...and', len(more_than_nine), 'had more than nine profiles')

            # If days were found with more than nine profiles: Print some diagnostics
            if len(more_than_nine):
                for i in range(len(more_than_nine)):
                    this_doy = more_than_nine[i] + 1     # convert from Python index to doy 1, 2, 3...
                    print("doy", this_doy, "had more than nine profiles")
                    # print()
                    # print('doy is', this_doy)
                    # print('-------------------')
                    # for j in range(len(profiles)):
                    #     if doy(profiles[j][1]) == this_doy:
                    #         PrintProfileEntry(profiles[j])
                    #         print()

            df = pd.DataFrame(data=np.array([np.array(x) for x in profiles]))
            df.to_csv(os.getcwd() + '/../Profiles/' + site + yrstr + '.csv')

                
                
def ReadProfiles(fnm):
    """
    Profiles are saved by site and year as tuples. Here we read only
    the datetimes; not the corresponding sample indices. There are six 
    values and these are nominally degenerate: ascent_end == descent_start
    and descent_end == rest_start. Read as strings these are converted to 
    Timestamps in the returned dataframe.
    """
    df = pd.read_csv(fnm, usecols=["1", "3", "5", "7", "9", "11"])
    df.columns=['ascent_start', 'ascent_end', 'descent_start', 'descent_end', 'rest_start', 'rest_end']
    df['ascent_start'] = pd.to_datetime(df['ascent_start'])
    df['ascent_end'] = pd.to_datetime(df['ascent_end'])
    df['descent_start'] = pd.to_datetime(df['descent_start'])
    df['descent_end'] = pd.to_datetime(df['descent_end'])
    df['rest_start'] = pd.to_datetime(df['rest_start'])
    df['rest_end'] = pd.to_datetime(df['rest_end'])
    return df