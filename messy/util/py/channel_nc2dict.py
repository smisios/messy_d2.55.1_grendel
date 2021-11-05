# -*- coding: utf-8 -*-
"""
Module for offline testing of forpy developments.

get_timesteps
    get number of timesteps in a netcdf file
create_dict
    Read data from a netcdf file and sort them into nested dictionaries and
    lists to mimic online provided fields from forpy.

@author: Bastian Kern, DLR
"""

from netCDF4 import Dataset
import numpy as np
import os.path


def get_timesteps(filepath):
    # open netcdf file
    rootgrp = Dataset(filepath, 'r', format='NETCDF3_64BIT_OFFSET')

    timesteps = rootgrp.dimensions['time'].size

    # close netcdf file
    rootgrp.close()
    return timesteps


def create_dict(filepath, timestep=1):
    # open netcdf file
    rootgrp = Dataset(filepath, 'r', format='NETCDF3_64BIT_OFFSET')
    rootgrp.set_auto_mask(False)
    rootgrp.set_always_mask(False)

    d = {}
    d['source'] = 'netcdf' # to distinguish from on-line transfer
    d['name'] = rootgrp.channel_name
    att_list = rootgrp.ncattrs()

    # channel and global attributes
    catt_dict = {}
    gatt_dict = {}
    for att in att_list:
        if att.startswith('channel'):
            catt_dict[att] = rootgrp.getncattr(att)
        else:
            gatt_dict[att] = rootgrp.getncattr(att)
    # for compatibility with on-line forpy, remove suffix .nc
    catt_dict['channel_file_name'] = os.path.splitext(catt_dict['channel_file_name'])[0]
    d['channel_attributes'] = catt_dict
    d['global_attributes'] = gatt_dict
    # save the channels representation name for later use:
    # there can be variables with different representation in the netcdf file,
    # but only the channels default representation (if set) is written into
    # netcdf.
    try:
        ch_representation = d['global_attributes']['REPRESENTATION']
    except:
        ch_representation = 'UNKNOWN'

    # dimensions
    dim_list = rootgrp.dimensions.keys()
    for dim in dim_list:
        dim_dict = {}
        dim_dict['name'] = dim
        if rootgrp.dimensions[dim].isunlimited():
           dim_dict['length'] = 1
        else:
           dim_dict['length'] = rootgrp.dimensions[dim].size
        dim_dict['dimvars'] = list()
        d[dim] = dim_dict

    # all variables (including dimension variables)
    obj_list = list(rootgrp.variables.keys())
    d['objects'] = list()
    for obj in obj_list:
        obj_dict = {}
        obj_dict['name'] = obj
        # read the object's attributes into list
        oatt_list = rootgrp.variables[obj].ncattrs()
        oatt_dict = {}
        for att in oatt_list:
            oatt_dict[att] = rootgrp.variables[obj].getncattr(att)
        obj_dict['attributes'] = oatt_dict
        # get dimensions of the object
        dimensions = rootgrp.variables[obj].dimensions
        # extract data for one timestep
        if len(dimensions) > 0 and dimensions[0] == 'time':
           ndims = 1 if len(dimensions) == 1 else len(dimensions) - 1
           obj_dict['data'] = np.atleast_1d(rootgrp.variables[obj][timestep-1])
        else:
            obj_dict['data'] = rootgrp.variables[obj][:]
            ndims = len(dimensions)
        # test if object is dimension variable
        # - only one dimension
        # - the dimension variable is available
        # - size of dimension of the object and the dimension variable is same
        if len(dimensions) == 1 and \
            dimensions[0] in d and \
                rootgrp.dimensions[dimensions[0]].size \
                == rootgrp.variables[obj].size:
            d[dimensions[0]]['dimvars'].append(obj_dict)
        # variable
        else:
            # compile a crude representation assumtion
            # not all information is available from the netcdf file
            repr_dict = {}
            repr_dict['rank'] = ndims
            repr_dict['link'] = ''
            axisstr = ''
            # axis string is inverted here, to match with model output
            for i in range(len(dimensions)):
                if 'lon' in dimensions[i]:
                    axisstr = 'X' + axisstr
                elif 'lat' in dimensions[i]:
                    axisstr = 'Y' + axisstr
                elif 'lev' in dimensions[i]:
                    axisstr = 'Z' + axisstr
                elif 'height' in dimensions[i]:
                    axisstr = 'Z' + axisstr
                # if the dimension has a dimension variable
                elif dimensions[i] in rootgrp.variables:
                    # positive attribute available -> vertical dimension
                    if 'positive' in rootgrp.variables[dimensions[i]].ncattrs():
                        axisstr = 'Z' + axisstr
                elif 'time' in dimensions[i]:
                    pass
                else:
                    axisstr = 'N' + axisstr
            repr_dict['axis'] = axisstr
            gdimlen = list()
            for i in range(len(dimensions)):
                if 'time' in dimensions[i]:
                    pass
                else:
                    # if we have to reverse, because python has C-style arrays
                    # insert at 0
                    # if not change this to an append...
                    gdimlen.insert(0, rootgrp.dimensions[dimensions[i]].size)
            repr_dict['gdimlen'] = gdimlen
            repr_dict['name'] = ch_representation
            if ndims == 2:
                if repr_dict['axis'] == 'XY':
                    repr_dict['name'] = 'GP_2D_HORIZONTAL'
            # find the needed dimensions and add diminsion variables
            objdim_list = list()
            for dim in dimensions:
                objdim_list.insert(0, d[dim])
            repr_dict['dimensions'] = objdim_list
            obj_dict['representation'] = repr_dict
            d['objects'].append(obj_dict)
    # close netcdf file
    rootgrp.close()
    return d


def main():
    print('Use this to put data from netcdf in nested dictionary, '
          'list structures.')


if __name__ == "__main__":
    main()
