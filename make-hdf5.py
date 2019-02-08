# Copyright 2013-2016 The Salish Sea MEOPAR contributors
# and The University of British Columbia

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Creates forcing HDF5 input files for MOHID based on user input time range
"""

import os
from datetime import datetime, timedelta
from dateutil.parser import parse
import numpy as np
import xarray as xr
import errno
import time
import h5py
from salishsea_tools import utilities
from salishsea_tools import viz_tools
#from scipy.interpolate import griddata

# NEMO input files directory
nemoinput = '/results2/SalishSea/nowcast-green.201806/'

# HRDPS input files directory
hdinput = '/results/forcing/atmospheric/GEM2.5/operational/'

# WW3 input files directory
wwinput = '/opp/wwatch3/nowcast/'

# Output filepath
outpath = '/results2/MIDOSS/forcing/SalishSeaCast/'


## Integer -> String
## consumes time in seconds and outputs a string that gives the time in HH:MM:SS format
def conv_time(time):
    """Give time in HH:MM:SS format.

    :arg time: time elapsed in seconds
    :type integer: :py:class:'int'

    :returns: time elapsed in HH:MM:SS format
    :rtype: :py:class:`str'
    """
    hours = int(time/3600)
    mins = int((time - (hours*3600))/60)
    secs = int((time - (3600 * hours) - (mins *60)))
    return '{}:{}:{}'.format(hours, mins, secs)


def generate_currents_hdf5(timestart, timeend, path, outpath, compression_level = 1):
    """Renerate current forcing HDF5 input file MOHID.

    :arg timestart: date from when to start concatenating
    :type string: :py:class:'str'

    :arg timeend: date at which to stop concatenating
    :type string: :py:class:'str'

    :arg path: path of input files
    :type string: :py:class:'str'

    :arg outpath: path for output files
    :type string: :py:class:'str'

    :arg compression_level: compression level for output file (Integer[1,9])
    :type integer: :py:class:'int'

    :returns: None
    :rtype: :py:class:`NoneType'
    """
    
    # generate list of dates from daterange given
    daterange = [parse(t) for t in [timestart, timeend]]
    U_files = []
    V_files = []
    T_files = []


    # append all filename strings within daterange to lists
    for day in range(np.diff(daterange)[0].days):
        datestamp = daterange[0] + timedelta(days=day)
        datestr1 = datestamp.strftime('%d%b%y').lower()
        datestr2 = datestamp.strftime('%Y%m%d')
        
        # check if file exists. exit if it does not. add path to list if it does.
            # U files
        U_path = f'{path}{datestr1}/SalishSea_1h_{datestr2}_{datestr2}_grid_U.nc'
        if not os.path.exists(U_path):
            print(f'File {U_path} not found. Check Directory and/or Date Range.')
            return
        U_files.append(U_path)
            # V files
        V_path = f'{path}{datestr1}/SalishSea_1h_{datestr2}_{datestr2}_grid_V.nc'
        if not os.path.exists(V_path):
            print(f'File {V_path} not found. Check Directory and/or Date Range.')
            return
        V_files.append(V_path)
            # T files
        T_path = f'{path}{datestr1}/SalishSea_1h_{datestr2}_{datestr2}_grid_T.nc'
        if not os.path.exists(T_path):
            print(f'File {T_path} not found. Check Directory and/or Date Range.')
            return
        T_files.append(T_path)
        
        
    # string: output folder name with date ranges used. end date will be lower by a day than timeend because datasets only go until midnight
    folder = str(datetime(parse(timestart).year, parse(timestart).month, parse(timestart).day).strftime('%d%b%y').lower()) + '-' + str(datetime(parse(timeend).year, parse(timeend).month, parse(timeend).day-1).strftime('%d%b%y').lower())
    # create output directory
    dirname = f'{outpath}hdf5/{folder}/'
    if not os.path.exists(os.path.dirname(dirname)):
        try:
            os.makedirs(os.path.dirname(dirname))
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
                
    # create hdf5 file and create tree structure
    f = h5py.File(f'{dirname}foocurrents.hdf5', 'w')
    times = f.create_group('Time')
    velocity_u = f.create_group('/Results/velocity U')
    velocity_v = f.create_group('/Results/velocity V')
    water_level = f.create_group('/Results/water level')
    
    number_of_files = len(U_files)
    bar = utilities.statusbar('Creating currents forcing file ...')
    for file_index in bar(range(number_of_files)):
        U_raw = xr.open_dataset(U_files[file_index])
        V_raw = xr.open_dataset(V_files[file_index])
        T_raw = xr.open_dataset(T_files[file_index])
        # assume all files have same time_counter markers
        datelist = U_raw.time_counter.values.astype('datetime64[s]').astype(datetime.datetime)
        # unstagger to move U, V to center of grid square
        U  = viz_tools.unstagger_xarray(U_raw.vozocrtx, 'x')
        V  = viz_tools.unstagger_xarray(V_raw.vomecrty, 'y').values[...,:,1:897:,1:397]
        # convert xarrays to numpy arrays and cut off grid edges
        U = U.values[...,:,1:897:,1:397]
        V = V.values[...,:,1:897:,1:397]
        sea_surface = T_raw.sossheig.values[...,:,1:897:,1:397]
        # rotate currents to True North
        current_u, current_v = viz_tools.rotate_vel(U, V)
        # clear memory
        U, V = None, None
        # transpose grid (rotate 90 clockwise)
        current_u = np.transpose(current_u, [0,1,3,2])
        current_v = np.transpose(current_v, [0,1,3,2])
        sea_surface = np.transpose(sea_surface, [0,2,1])
        # flip currents by depth dimension
        current_u = np.flip(current_u, axis = 1)
        current_v = np.flip(current_v, axis = 1)
        # convert nans to 0's and set datatype to float64
        current_u = np.nan_to_num(current_u).astype('float64')
        current_v = np.nan_to_num(current_v).astype('float64')
        sea_surface = np.nan_to_num(sea_surface).astype('float64')
        # make list of time arrays
        datearrays = []
        for date in datelist:
            datearrays.append(np.array([date.year, date.month, date.day, date.hour, date.minute, date.second]).astype('float64'))
        # write u wind values to hdf5
        for i in range(current_u.shape[0]):
            velocity_attr = 'velocity U_' + ((5 - len(str(i + 1))) * '0') + str(i + 1)
            dset = velocity_u.create_dataset(velocity_attr, shape = (40, 396, 896), data = current_u[i],chunks=(40, 396, 896), compression = 'gzip', compression_opts = compression_level)
            metadata = {'FillValue' : np.array([0.]), 'Maximum' : np.array([5.]), 'Minimum' : np.array([-5.]), 'Units' : b'm/s'}
            dset.attrs.update(metadata)
    
        # write v wind values to hdf5
        for i in range(current_v.shape[0]):
            velocity_attr = 'velocity V_' + ((5 - len(str(i + 1))) * '0') + str(i + 1)
            dset = velocity_v.create_dataset(velocity_attr, shape = (40, 396, 896), data = current_v[i],chunks=(40, 396, 896), compression = 'gzip', compression_opts = compression_level)
            metadata = {'FillValue' : np.array([0.]), 'Maximum' : np.array([5.]), 'Minimum' : np.array([-5.]), 'Units' : b'm/s'}
            dset.attrs.update(metadata)
    
        # write  water level values to hdf5

        for i in range(sea_surface.shape[0]):
            level_attr = 'water level_' + ((5 - len(str(i + 1))) * '0') + str(i + 1)
            dset = water_level.create_dataset(level_attr, shape = (396, 896), data = sea_surface[i],chunks=(396, 896), compression = 'gzip', compression_opts = compression_level)
            metadata = {'FillValue' : np.array([0.]), 'Maximum' : np.array([5.]), 'Minimum' : np.array([-5.]), 'Units' : b'm'}
            dset.attrs.update(metadata)
    
        # write time values to hdf5

        for i in range(len(datearrays)):
            time_attr = 'Time_' + ((5 - len(str(i + 1))) * '0') + str(i + 1)
            dset = times.create_dataset(time_attr, shape = (6,), data = datearrays[i],chunks=(6,), compression = 'gzip', compression_opts = compression_level)
            metadata = {'Maximum' : np.array([2016.]), 'Minimum' : np.array([-0.]), 'Units' : b'YYYY/MM/DD HH:MM:SS'} # !!!
            dset.attrs.update(metadata)
            
    f.close()
    return


# input start time
timestart = input('\nEnter the start time in the format |year month day| e.g. 2015 Jan 1:\n--> ')

# input end time. This must be day after required range as upper bound is not inlcuded
timeend = input('\nEnter the end time in the format |year month day| e.g. 2015 Jan 1:\n--> ')

# what parts to run
runs = int(input('\nRun: \n1) All \n2) NEMO ONLY \n3) HRDPS ONLY \n4) WW3 ONLY ?\n--> '))
runsdict = {1: 'All', 2: 'Currents ONLY', 3: 'Winds ONLY', 4: 'WW3 ONLY'}

# user confirmation
ask = input(f'\nProceed with generating input files {runsdict[runs]} from {timestart} to {timeend}?\n--> ')
beganat = time.time()

# run according to user inputs
if ask in ['y', 'yes', 'YES', 'Y']:
    if runs == 2:
        generate_currents_hdf5(timestart, timeend, nemoinput, outpath, compression_level = 1)
    else:
        print('\nAborted')
print('All done\n')
print('Time elapsed: {}'.format(conv_time(time.time()-beganat)))
