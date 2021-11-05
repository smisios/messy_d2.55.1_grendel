#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:00:00 2020

@author: Patrick Joeckel, DLR
"""

# for plotting
import matplotlib as mpl
mpl.use('Agg')
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
import numpy as np
# basics
import sys
import os.path
import time
from datetime import datetime, timedelta
import yaml
#
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)
# parallel
from functools import partial
from multiprocessing import Pool, cpu_count

"""
 ################################
 ### GENERAL HELPER FUNCTIONS ###
 ################################
"""

def f_start_timer():
    """
    start a timer
    """
    return time.time()

def f_delta_timer(started):
    """
    stop a timer and print the time interval
    """
    elapsed = time.time() - started
    return elapsed

def f_timer(func):
    """
    decorater for measuring time spent in specific functions
    """
    def helper(*args, **kwargs):
        t0 = f_start_timer()
        func(*args, **kwargs)
        elapsed = f_delta_timer(t0)
        print(f'(elapsed time in {func.__name__}: {elapsed} s)')
    return helper

def f_getn(dict, keys, default=None):
    """
    get value from nested dictionaries according to list of keys
    """
    id_value = dict
    for k in keys:
        id_value = id_value.get(k, None)
        if id_value is None:
            return default
    return id_value

def f_print_nested(d, indent=0):
    for key, value in d.items():
        print('  ' * indent + str(key))
        if isinstance(value, dict):
            f_print_nested(value, indent+1)
        else:
            print('  ' * (indent+1) + str(value))

def f_merge_default(x, y):
    """
    y is the dictionary with defaults
    """
    for key in y:
        if key in x:
            if isinstance(x[key], dict) and isinstance(y[key], dict):
                f_merge_default(x[key], y[key])
        else:
            x[key] = y[key]
    return x

class PrettySafeLoader(yaml.SafeLoader):
    def construct_python_tuple(self, node):
        return tuple(self.construct_sequence(node))
    def construct_python_list(self, node):
        return list(self.construct_sequence(node))

def set_norm(func, args):
    dispatcher = {"mpl.colors.Normalize": mpl.colors.Normalize,
                  "mpl.colors.LogNorm": mpl.colors.LogNorm,
                  "mpl.colors.BoundaryNorm": mpl.colors.BoundaryNorm
              }
    #try:
    return dispatcher[func](**args)
    #except:
    #    return None

"""
 #########################################
 ### CHANNEL SPECIFIC HELPER FUNCTIONS ###
 #########################################
"""

def f_filename_base(dict):
    """
    return file name base
    """
    ch_att = dict["channel_attributes"]
    return ch_att["channel_file_name"]

def f_object_nlu(obj):
    """
    return object name, long name, and units for print title
    """
    n = obj["name"]
    try:
        ln = obj["attributes"]["long_name"]
    except KeyError:
        try:
            ln = obj["attributes"]["standard_name"]
        except KeyError:
            ln = obj["name"]
    if ln == '':
        ln = obj["name"]
    try:
        un = obj["attributes"]["units"]
    except KeyError:
        un = ''
    return n, ln, un

def f_object_repr(obj):
    """
    return representation of object
    """
    # on-line forpy: representation always available
    repr =  obj["representation"]
    name = repr["name"]
    rank = repr["rank"]
    axes = repr["axis"].replace('-','')
    dims = repr["gdimlen"]
    return name, rank, axes, dims

def f_timestamp(d):
    """
    convert time and units to string
    units: "day since YYYY-MM-DD hh:mi:se"
    """
    dimvar = d["dimvars"]
    for v in dimvar:
        if v["name"] == d["name"]:  # time(time)
            a = v["attributes"]
            tunit = a["units"]
            value = v["data"][0]
            pass
    tus = tunit.split(' ')
    tu = tus[0] # days
    #    tus[1] # since
    #              0123456789
    #              12345678910
    t0d = tus[2] # YYYY-MM-DD
    t0t = tus[3] # hr:mi:se
    start = datetime(int(t0d[0:4]),int(t0d[5:7]),int(t0d[8:10]),
                     int(t0t[0:2]),int(t0t[3:5]),int(t0t[6:8]))
    delta = timedelta(value)
    offset = start + delta
    return offset.isoformat(sep=' ', timespec='seconds') + ' UTC'

def f_print_attributes(dict, prefix):
    """
    print attribute dictionary as key, value pairs
    """
    for key in dict:
        print(f'{prefix} A: {key} -> {dict[key]}')

def f_print_dimension_variable(v, prefix):
    """
    print dimension variable
    """
    print(f'{prefix} V: NAME : {v["name"]}')
    att = v["attributes"]
    f_print_attributes(att, prefix + '   ')
    #print(f'{prefix}    DATA : {v["data"]}')

def f_print_dimension(d, prefix):
    """
    print dimension
    """
    print(f'{prefix} D: NAME   : {d["name"]}')
    print(f'{prefix}    LENGTH : {d["length"]}')
    dimvar = d["dimvars"]
    for v in dimvar:
        f_print_dimension_variable(v, prefix + '   ')

def f_print_representation(repr, prefix):
    """
    print representation information
    """
    print(f'{prefix} R: NAME : {repr["name"]}')
    print(f'{prefix} R: RANK : {repr["rank"]}')
    print(f'{prefix} R: LINK : {repr["link"]}')
    print(f'{prefix} R: AXES : {repr["axis"]}')
    print(f'{prefix} R: DIM  : {repr["gdimlen"]}')
    print(f'{prefix} R: DIMENSIONS  :')
    dimensions = repr["dimensions"]
    for d in dimensions:
        f_print_dimension(d, prefix + '   ')

def f_print_object(obj, prefix):
    """
    print object information
    """
    print(f'{prefix} O: {obj["name"]}')
    att = obj["attributes"]
    f_print_attributes(att, prefix + '   ')
    repr = obj["representation"]
    print(f'{prefix}    RERPRESENTATION: ')
    f_print_representation(repr, prefix + '   ')
    print(f'{prefix}    DATA TYPE      : {type(obj["data"])}')
    #print(f'{prefix}    DATA           : {obj["data"]}')

def get_dimvar_zn(obj):
    repr = obj["representation"]
    dimensions = repr["dimensions"]
    for d in dimensions:
        d_name = d["name"]
        dimvars = d["dimvars"]
        for dv in dimvars:
            dv_name = dv["name"]
            if dv_name == d_name:
                break
        try:
            att = dv["attributes"]
        except KeyError:
            return dv
        a = None
        try:
            a = att["long_name"]
        except KeyError:
            pass
        try:
            a = att["standard_name"]
        except KeyError:
            pass
        if a is None:
            return dv
        else:
            if "longitude" in a:
                continue
            if "latitude" in a:
                continue
            return dv

def f_basemodel(dict):
    gl_att = dict["global_attributes"]
    chksum = gl_att["EXEC_CHECKSUM"]
    words = chksum.split()
    bm = os.path.splitext(os.path.basename(words[1]))[0]
    return bm

def set_filename(b, t, c, o):
    """
    constuct filename out of basename (b), time stamp (t), channel name (c),
    and object name (o)
    """
    tf = t.replace(' UTC','').replace(' ','_').replace('-','').replace(':','')
    f = b[:15] + tf + '_' + c + '_' + o
    return f

def get_geo_grid(obj, dict, i):
    """
    extract geographical grid information for object
    Notes: It depends on the basemodel, if the grid axes can be derived
           from the representation (dimension variable) or from a seperate
           object. In the 1st case (representation), the information is
           received from the corresponding dimensions variable attribute
           long_name. In the 2nd case the objects are selected by name
           directly.
    """
    bm = f_basemodel(dict)
    if bm == "cosmo":
        # objects of geographical coordinates
        object_name = ["geolon","geolat"][i]
        # dimension variables' long_name
        long_name = ["rotated longitude","rotated latitude"][i]
    elif bm == "echam5":
        # objects of geographical coordinates
        object_name = [None, None][i]
        # dimension variables' long_name
        long_name = ["longitude","latitude"][i]
    else:
        return None
    ###
    try:
        repr = obj["representation"]
    except KeyError:
        return None
    ###
    # if the representation information is available, search for
    # dimension variable with corresponding attribute long_name
    # the result is the axis in "rotated" coordinates
    dims = repr["dimensions"]
    geo_dvar = None
    for d in dims:
        dimvar = d["dimvars"]
        geo_dvar = next((x for x in dimvar
                         if x["attributes"]["long_name"] == long_name
                     )
                        , None)
        if geo_dvar is not None:
            break
    #
    if geo_dvar is None:
        geo_r = None
    else:
        geo_r = geo_dvar["data"]
    #
    ### next: extract geographical coordinates from objects
    # pre-set rotated with dimension variable (on unrotated grids, they are
    # the same ...
    geo_obj = geo_dvar
    #
    if object_name is not None:
        obj_list = dict["objects"]
        geo_obj = next((x for x in obj_list if x["name"] == object_name), None)
    #
    if geo_obj is None:
        geo = None
    else:
        geo = geo_obj["data"]
    ###
    return geo_r, geo

def get_geo_range(lon, lat):
    return (np.amin(lon),np.amin(lat),np.amax(lon),np.amax(lat))

def init_fig(xyrange, pixels=1024, kml=False):
    """
    returns a matplotlib `fig` and `ax` with approproiate aspect ratio
    xyrange is a tuple with
     - lower left  corner longitude (xmin)
     - lower left  corner latitude  (ymin)
     - upper right corner longitude (xmax)
     - upper right corner latitude  (ymax)
    """
    llcrnrlon = xyrange[0]
    llcrnrlat = xyrange[1]
    urcrnrlon = xyrange[2]
    urcrnrlat = xyrange[3]
    #
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize
    #
    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)
    #
    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    if kml:
        fig = plt.figure(figsize=figsize, frameon=False, dpi=pixels//10)
        ax  = fig.add_axes([0, 0, 1, 1])
    else:
        figsize = (figsize[0] - 1, figsize[1] + 1)
        fig, ax = plt.subplots(figsize=figsize, frameon=True, dpi=pixels//10)
        ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    #
    return fig, ax

def get_object_by_name(dict, name):
    obj_list = dict["objects"]
    for obj in obj_list:
        if obj["name"] == name:
            return obj
    return None

def shift_longitude(data, lon):
    """
    function to convert longitude and data from 0:360 to -180:180
    """
    lon = np.where(lon > 180, lon-360, lon)
    ilon = np.argsort(lon)
    # vector data?
    isvec = isinstance(data, list)
    if isvec:
        d = [q[:,ilon] for q in data]
    else:
        d = data[:,ilon]
    return d, lon[ilon]

def close_gap(data, lon):
    """
    function to close the gap for global data by introducing
    another data point
    """
    # vector data?
    isvec = isinstance(data, list)
    if isvec:
        d_i = []
        for q in data:
            q_i = np.zeros((q.shape[0],q.shape[1]+1))
            q_i[:,0] = q[:,-1]
            q_i[:,1:] = q
            d_i.append(q_i)
    else:
        d_i = np.zeros((data.shape[0],data.shape[1]+1))
        d_i[:,0] = data[:,-1]
        d_i[:,1:] = data
    ###
    l_i = np.zeros(len(lon)+1)
    l_i[0] = lon[-1] - 360
    l_i[1:] = lon
    ###
    return d_i, l_i

def get_plot_object(obj, dict, lswap, lshift, level=None, verbose=False):
    """
    return dictionary with data and annotations for plotting routines
    """
    ch = dict["name"]
    d = {}
    name, lname, units = f_object_nlu(obj)
    # vector data?
    isvec = isinstance(obj["data"], list)
    if level is None:
        d["nztag"] = ''
        d["nzval"] = None
        d["nzuni"] = ''
        d["nztit"] = ''
        d["level"] = None
        d["data"] = obj["data"]
        s = ''
    else:
        d["nztag"] = f'-{level:03}'
        dv = get_dimvar_zn(obj)
        d["nzval"] = f_getn(dv, ["data"], [])[level-1]
        d["nzuni"] = f_getn(dv, ["attributes", "units"], '')
        d["nztit"] = f' (at {d["nzval"]} {d["nzuni"]}) '
        d["level"] = level - 1
        if isvec:
            d["data"] = [q[level-1] for q in obj["data"]]
        else:
            d["data"]  = obj["data"][level-1]
        s = f'@{level}'
    ###
    d["name"] = name
    d["title"] = lname + d["nztit"] + ' [' + units + ']'
    ###
    ### extract horizontal grid information for plotting
    ###
    _, lon = get_geo_grid(obj, dict, 0)
    _, lat = get_geo_grid(obj, dict, 1)
    ###
    ### return on error
    ###
    if lon is None:
        if v:
            print(f' ... WARNING: no longitude information'
                  f' found for OBJECT {name} of CHANNEL {ch}')
        return None
    if lat is None:
        if v:
            print(f' ... WARNING: no latitude  information'
                  f' found for OBJECT {name} of CHANNEL {ch}')
        return None
    d["lat"] = lat
    d["lon"] = lon
    ###
    if lswap:
        if isvec:
            d["data"] = [np.swapaxes(q,0,1) for q in d["data"]]
        else:
            d["data"] = np.swapaxes(d["data"],0,1)
    #
    if lshift:
        d["data"], d["lon"] = shift_longitude(d["data"], d["lon"])
        d["data"], d["lon"] = close_gap(d["data"], d["lon"])
    ###
    ### check for missing values / fill value
    ###
    fill = f_getn(obj, ["attributes", "_FillValue"], None)
    miss = f_getn(obj, ["attributes", "missing_value"], fill)
    if isvec:
        qdat = d["data"]
    else:
        qdat = [d["data"]]
    qd = []
    for q in qdat:
        if miss is not None:
            q = np.ma.masked_invalid(np.ma.masked_values(q, miss))
        else:
            q = np.ma.masked_invalid(q)
        ###
        q.filled(np.nan)
        qd.append(q)
        ###
        if verbose:
            fm = np.ma.count_masked(q, axis=None)
            print(f' ... ... {fm} masked values in OBJECT {name}{s}'
                  f' of CHANNEL {ch}')
            fn = np.isnan(q).sum()
            print(f' ... ... {fn} nan values in OBJECT {name}{s}'
                  f' of CHANNEL {ch}')
        ###
    if isvec:
        d["data"] = qd
    else:
        d["data"] = qd[0]
    ###
    return d

def get_norm(args, data):
    """
    set appropriate arguments, depending on chosen norm
    """
    normargs = {}
    found = False
    normfun = f_getn(args, ["norm"], "mpl.colors.Normalize")
    ###
    ### range
    ###
    vmin = f_getn(args, ["vmin"], data.min())
    vmax = f_getn(args, ["vmax"], data.max())
    #if not vmax > vmin:
    #    vmax = vmin + 42.0
    levels = f_getn(args, ["levels"], 10)
    if isinstance(levels, int):
        bounds = np.linspace(vmin, vmax, levels)
    else:
        ### assuming list here
        bounds = levels
    ###
    if normfun == "mpl.colors.BoundaryNorm":
        cmap = f_getn(args, ["cmap"], plt.rcParams["image.cmap"])
        normargs["ncolors"] = plt.get_cmap(cmap).N
        normargs["boundaries"] = bounds
        normargs["extend"] = f_getn(args, ["extend"], "neither")
        normargs["clip"] = f_getn(args, ["clip"], False)
        found = True
        lognorm = False
    if normfun == "mpl.colors.Normalize":
        normargs["vmin"] = vmin
        normargs["vmax"] = vmax
        normargs["clip"] = f_getn(args, ["clip"], False)
        found = True
        lognorm = False
    if normfun == "mpl.colors.LogNorm":
        normargs["vmin"] = vmin
        normargs["vmax"] = vmax
        normargs["clip"] = f_getn(args, ["clip"], False)
        found = True
        lognorm = True
    if found:
        norm = set_norm(normfun, normargs)
    else:
        norm = None
        lognorm = None
    return norm, lognorm

def get_colorargs(args, data, pf, overlay=False):
    """
    set colorargs depending on pf
    """
    colorargs = {}
    cbarargs = {}
    ### common arguments
    colorargs["cmap"] = f_getn(args, ["cmap"], None)
    colorargs["alpha"] = f_getn(args, ["alpha"], None)
    cbarargs["extend"] = f_getn(args, ["extend"], "neither")
    if overlay:
        colorargs["colors"] = f_getn(args, ["colors"], "gray")
    ### specific arguments
    if pf == "pcolormesh":
        if not overlay:
            colorargs["shading"] = f_getn(args, ["shading"], 'auto')
    if pf == "contourf" or pf == "quiver":
        ### this is for correcting the counter-intuitive behaviour
        ### of countourf (vmin, vmax) have no effect on color bar range
        vmin = f_getn(args, ["vmin"], data.min())
        vmax = f_getn(args, ["vmax"], data.max())
        if not vmax > vmin:
            vmax = vmin + 42.0
        colorargs["vmin"] = vmin
        colorargs["vmax"] = vmax
        levels = f_getn(args, ["levels"], 10)
        if isinstance(levels, int):
            bounds = np.linspace(vmin, vmax, levels)
        else:
            ### assuming a list here
            bounds = levels
        colorargs["levels"] = bounds
        colorargs["extend"] = f_getn(args, ["extend"], "neither")
    ###
    if pf == "quiver":
        colorargs.pop("vmin",None)
        colorargs.pop("vmax",None)
        colorargs.pop("levels",None)
        colorargs.pop("extend",None)
        colorargs.pop("colors",None)
        colorargs["units"]  = f_getn(args, ["units"], "width")
        colorargs["angles"] = f_getn(args, ["angles"], "uv")
        colorargs["scale"]  = f_getn(args, ["scale"], None)
        colorargs["scale_units"] = f_getn(args, ["scale_units"], None)
        if overlay:
            colorargs["color"]  = f_getn(args, ["color"], "black")
    ###
    return colorargs, cbarargs

"""
 ##################################
 ### MAIN ENTRY POINT (LEVEL-0) ###
 ##################################
 - called from fortran/forpy or from wrapper (channel_fpy_wrapper.py)
"""

def fpy_main(dict):
    """
    This is the main entry point from forpy.
    The yaml-configuration file is read for decisions on how to
    proceed further.
    """
    t0 = f_start_timer()
    ### extract channel name
    ch = dict["name"]
    ###
    ### add tuple and list as safe yaml objects
    ###
    PrettySafeLoader.add_constructor(
        u'tag:yaml.org,2002:python/tuple',
        PrettySafeLoader.construct_python_tuple)
    PrettySafeLoader.add_constructor(
        u'tag:yaml.org,2002:python/list',
        PrettySafeLoader.construct_python_list)
    ### read configuration file
    with open('channel.yml') as f:
        #conf = yaml.load(f, Loader=yaml.FullLoader)
        #conf = yaml.load(f, Loader=yaml.SafeLoader)
        conf = yaml.load(f, Loader=PrettySafeLoader)
        f.close()
    ### detect subdirectory setting
    subdir = f_getn(conf, ["global", "subdir"], ".")
    if subdir != ".":
        os.makedirs(subdir, exist_ok=True)
    ### redirect output to log-file
    logfile = f_getn(conf, ["global", "logfile"], None)
    if logfile is not None:
        sys.stdout = open(subdir + '/' + logfile, 'w')
    ### detect verbosity
    verbose = f_getn(conf, ["global", "verbose"], False)
    ### detect function to be called
    func = f_getn(conf, ["channels", ch], None)
    if func is None:
        func = f_getn(conf, ["channels", "default"], fpy_log)
    ###
    ### global configuration
    ###
    gconf = {}
    gconf["subdir"]    = subdir
    gconf["verbose"]   = verbose
    gconf["basemodel"] = f_basemodel(dict)
    gconf["basename"]  = f_filename_base(dict)
    ###
    ### data axes need to be swapped, if input is via forpy, but not
    ### from netcdf file
    try:
        src = dict["source"]
    except KeyError:
        src = "online"
    gconf["lswap"] = src == "online"
    ###
    ### for the ECHAM5 basemodel, the longitudes need to be shifted
    ### from [0,360] to [-180,180] and the "gap" at the dateline needs
    ###  to be closed
    gconf["lshiftlon"] = gconf["basemodel"] == "echam5"
    ###
    ### string with time information
    ###
    gconf["tstring"] = f_timestamp(dict["time"])
    ###
    ###
    ### call function configured for this channel
    if verbose:
        print(f' --------------------------------')
        print(f' --- PYTHON VIA FORPY (START) ---')
        print(f' --------------------------------')
        print(f' - FUNCTION    : fpy_main')
        print(f' - CHANNEL     : {ch}')
        print(f' - OUTPUT DIR  : {subdir}')
        print(f' - BASEMODEL   : {gconf["basemodel"]}')
        print(f'   - SWAP AXES : {gconf["lswap"]}')
        print(f'   - SHIFT LON : {gconf["lshiftlon"]}')
        print(f' - DATE & TIME : {gconf["tstring"]}')
        print(f' - CALLING     : {func}')
    ###
    ### alays call this function
    fpy_def(dict,conf,gconf)
    ###
    globals()[func](dict,conf,gconf)
    ###
    if verbose:
        print(f' --------------------------------')
        print(f' --- PYTHON VIA FORPY (END)   ---')
        print(f' --------------------------------')
    ###
    elapsed = f_delta_timer(t0)
    print(f'(elapsed time in fpy_main: {elapsed} s)')
    ###
    ### close logfile
    ###
    if logfile is not None:
        sys.stdout.close()
    return 0

"""
 ###############################################
 ### LEVEL-1 ROUTINES (for complete channel) ###
 ###############################################
 - called from main entry point (LEVEL-0) as
   configured via yaml-configuration file (channel.yml)
"""

"""
 ###################################
 ### LEVEL-1 ROUTINE FOR TESTING ###
 ### OUTPUT OF INFO TO LOG-FILE  ###
 ###################################
"""

@f_timer
def fpy_log(dict, conf, gconf):
    """
    print contents of channel dictionary
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    bm = gconf["basemodel"]
    fb = gconf["basename"]
    ts = gconf["tstring"]
    ### Note: fpy_log is always "verbose"
    print(f' -----------------------')
    print(f' --- fpy_log (START) ---')
    print(f' -----------------------')
    print(f' - BASEMODEL: {bm}')
    print(f' - CHANNEL  : {ch}')
    #
    # channel specific attributes
    ch_att = dict["channel_attributes"]
    f_print_attributes(ch_att,'             ')
    # global attributes
    gl_att = dict["global_attributes"]
    f_print_attributes(gl_att,'             ')
    # time information
    print(f' - TIME     :')
    f_print_dimension(dict["time"],'             ')
    # this is a list of dictionaries
    print(f' - OBJECT(S): ')
    obj_list = dict["objects"]
    for obj in obj_list:
        f_print_object(obj,'   ')
        ###
        ### extract information for plot
        ###
        ### ... name, long_name, and units
        name, lname, units = f_object_nlu(obj)
        ### ... representation, rank, axes, dimension lengths
        repr, rank, axes, dims  = f_object_repr(obj)
        ###
        title = lname + ' [' + units + ']'
        ###
        ### construct file name
        ###
        fname = set_filename(fb, ts, ch, name)
        ###
        print(f' ... {fname}: {title} -- {ts} ({repr}, {dims})')
        print(f'')
    print(f' -----------------------')
    print(f' --- fpy_log (END)   ---')
    print(f' -----------------------')
    return 0

"""
 ###############################
 ### LEVEL-1 PROLOG ROUTINE  ###
 ### ALWAYS CALLED !!!       ###
 ###############################
"""

@f_timer
def fpy_def(dict, conf, gconf):
    """
    define new objects
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    bm = gconf["basemodel"]
    fb = gconf["basename"]
    ts = gconf["tstring"]
    if v:
        print(f' -----------------------')
        print(f' --- fpy_def (START) ---')
        print(f' -----------------------')
        print(f' - BASEMODEL: {bm}')
        print(f' - CHANNEL  : {ch}')
    ###
    ### function specific configuration
    ###
    fconf = f_getn(conf, ["functions", "fpy_def"], {})
    for key, value in fconf.items():
        ### check if key fits actual channel
        chaname = key.split('::')[0]
        try:
            objname = key.split('::')[1]
        except:
            if v:
                print(f' ... no object name (channelname::objectname),'
                      f' skipping {key}')
            continue
        if not chaname == ch:
            if v:
                print(f' ... CHANNEL is {ch}, skipping {key}')
            continue
        ### check, if object with same name exists already
        o = get_object_by_name(dict, objname)
        if o is not None:
            if v:
                print(f' ... OBJECT {objname} in CHANNEL {ch} exists'
                      f' already, skipping {key}')
            continue
        ### get function, arguments and attributes
        func_dir = f_getn(value, ["function"], None)
        if func_dir is None:
            if v:
                print(f' ... function missing, skipping {key}')
            continue
        func  = next(iter(func_dir))
        fargs = func_dir[func]
        if fargs is None:
            fargs = {"name": objname}
        else:
            fargs["name"] = objname
        atts = f_getn(value, ["attributes"], {})
        obj = globals()[func](dict, conf, gconf, fconf,
                              attributes=atts, **fargs)
        ### append object to list of objects in dictionary
        dict["objects"].append(obj)
    if v:
        print(f' -----------------------')
        print(f' --- fpy_def (END)   ---')
        print(f' -----------------------')
    ###
    return 0

def vector(dict, conf, gconf, fconf, attributes, **kwargs):
    """
    combine two scalar fields to a vector field
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    ### get required arguments
    xobjname = kwargs.pop("x",None)
    yobjname = kwargs.pop("y",None)
    if xobjname is None or yobjname is None:
        if v:
            print(f' ... ... vector: ERROR: x and/or y argument missing')
        return None
    ### retrieve objects
    xobj = get_object_by_name(dict, xobjname)
    if xobj is None:
        if v:
            print(f' ... ... vector: ERROR: x object {xobjname} not available')
        return None
    yobj = get_object_by_name(dict, yobjname)
    if yobj is None:
        if v:
            print(f' ... ... vector: ERROR: y object {yobjname} not available')
        return None
    ### now combine objects
    # - copy x-object (to include all meta-information)
    obj = xobj.copy()
    # - overwrite name and attributes
    obj["name"] = kwargs.pop("name", None)
    for key, value in attributes.items():
        obj["attributes"][key] = value
    ### work on data ...
    x = xobj["data"]
    y = yobj["data"]
    z = (xobj["data"]**2 + yobj["data"]**2)**0.5
    obj["data"] = [x, y, z]
    ###
    return obj

def equivalence(dict, conf, gconf, fconf, attributes, **kwargs):
    """
    define an object as equivalence to another one
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    ### get required arguments
    objname = kwargs.pop("object",None)
    if objname is None:
        if v:
            print(f' ... ... equivalence: ERROR: object argument missing')
        return None
    ### retrieve objects
    xobj = get_object_by_name(dict, objname)
    if xobj is None:
        if v:
            print(f' ... ... equivalence: ERROR: object {objname}'
                  f' not available')
        return None
    ### now combine objects
    # - copy x-object (to include all meta-information)
    obj = xobj.copy()
    # - overwrite name and attributes
    obj["name"] = kwargs.pop("name", None)
    for key, value in attributes.items():
        obj["attributes"][key] = value
    ###
    return obj

"""
 ####################################################
 ### LEVEL-1 ROUTINE FOR PLOTTING HORIZONTAL MAPS ###
 ####################################################
"""

@f_timer
def fpy_map(dict, conf, gconf):
    """
    plot horizontal maps (x,y)
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    bm = gconf["basemodel"]
    fb = gconf["basename"]
    ts = gconf["tstring"]
    sd = gconf["subdir"]
    ###
    ### additional features
    ###
    ### ... function specific configuration
    ###
    fconf = f_getn(conf, ["functions", "fpy_map"], {})
    s     = f_getn(fconf, ["subdirs"], False)
    ipar  = f_getn(fconf, ["parallel"], 0)
    ###
    if v:
        print(f' -----------------------')
        print(f' --- fpy_map (START) ---')
        print(f' -----------------------')
        print(f' - CHANNEL            : {ch}')
        print(f' - OUTPUT DIRECTORY   : {sd}')
        print(f' - BASEMODEL          : {bm}')
        print(f' - SUBDIRECTORIES     : {s}')
        print(f' - DATE AND TIME      : {ts}')
        if ipar <= 0:
            print(f' - MULTIPROCESSING    : OFF')
        else:
            print(f' - MULTIPROCESSING    : {ipar}; {cpu_count()} cores')
    #
    obj_list = dict["objects"]
    if ipar <= 0:
        for obj in obj_list:
            prepare_map_plot(obj, dict, gconf, fconf)
    else:
        n = min(ipar, cpu_count())
        with Pool(n) as pool:
            func = partial(prepare_map_plot, dict=dict,
                           gconf=gconf, fconf=fconf)
            pool.map(func, obj_list)
    if v:
        print(f' -----------------------')
        print(f' --- fpy_map (END)   ---')
        print(f' -----------------------')
    ###
    return 0

def config_complete(conf):
    """
    function to complete configuration internally:
    set std and kml to internal default
    """
    conf_def    = f_getn(conf, ["default"], {})
    conf["std"] = f_merge_default(f_getn(conf, ["std"], {}), conf_def)
    conf["kml"] = f_merge_default(f_getn(conf, ["kml"], {}), conf_def)
    return conf

def prepare_map_plot(obj, dict, gconf, fconf):
    """
    intermediate routine to prepare plot objects
    (allows parallelisation) via multiprocessing
    """
    ### from external dictionary
    ch = dict["name"]
    ### from global configuration
    v  = gconf["verbose"]
    bm = gconf["basemodel"]
    fb = gconf["basename"]
    ts = gconf["tstring"]
    ###
    ### defaults (global and channel specific)
    ###
    deconf    = config_complete(f_getn(fconf, ["default"], {}))
    chconf    = config_complete(f_getn(fconf, [ch], {}))
    fconf_cha = f_merge_default(chconf, deconf)
    ### ... name, long_name, and units
    oname, lname, units = f_object_nlu(obj)
    ### step 1: detect rank and axes and create list of
    ###         potential object names specified in yaml file
    _, rank, axes, dims  = f_object_repr(obj)
    ###
    valid = None
    names = []
    if (axes == 'XY'):
        names.append(ch + '::' + oname)
        one = True
    if (axes == 'XYZ' or axes == 'XYN'):
        for i in range(0, dims[2]):
            names.append(ch + '::' + oname + '@' + f'{i+1}')
        one = False
    if not names:
        if v:
            print(f' ... WARNING: OBJECT {oname} of CHANNEL {ch}'
                  f' has none of the shapes XY, XYZ, or XYN')
        return
    ###
    i = 0
    for name in names:
        i += 1
        # 1st step: check for specific settings cha::ob[@...]
        # 2nd step: merge keys of modes internally
        oconf = config_complete(f_getn(fconf, [name], {}))
        # 3rd step: add missing keys from specific settings cha::obj
        cconf = config_complete(f_getn(fconf,[ch + '::' + oname], {}))
        oconf = f_merge_default(oconf, cconf)
        # 4th step: add missing keys from channel specific setting
        oconf = f_merge_default(oconf, fconf_cha)
        #
        #f_print_nested(oconf)
        #
        skip = f_getn(oconf, ["skip"], None)
        if skip:
            if v:
                print(f' ... WARNING: CHANNEL::OBJECT[@..]'
                      f' {name} skipped'
                      f' according to channel.yml')
            continue
        ###
        if one:
            level = None
        else:
            level = i
        pobj = get_plot_object(obj, dict,
                               gconf["lswap"], gconf["lshiftlon"],
                               level=level, verbose=v)
        ###
        ### check for mode dependent overlays and generate lists
        ### of objects to plot ...
        po = {"std": [pobj], "kml": [pobj]}
        for mode in po:
            overlay = f_getn(oconf, [mode, "overlay"], [])
            ### loop over list
            for ol in overlay:
                ### loop over dictionary
                for od in ol:
                    ### check for level informatin
                    a = od.split('@')
                    olname = a[0]
                    if (len(a) == 1):
                        olevel = None
                    else:
                        olevel = int(a[1])
                ovl_obj = get_object_by_name(dict, olname)
                if ovl_obj is not None:
                    oobj = get_plot_object(ovl_obj, dict,
                                           gconf["lswap"],
                                           gconf["lshiftlon"],
                                           level=olevel, verbose=v)
                    # save full name for plotting routines
                    oobj["name"] = od
                    po[mode].append(oobj)
        ###
        for mode in po:
            skip = f_getn(oconf, [mode, "skip"], False)
            if skip:
                if v:
                    print(f' ... WARNING: CHANNEL::OBJECT[@..]'
                          f' {name} skipped for {mode} mode'
                          f' according to channel.yml')
                continue
            fname = set_filename(fb, ts, ch, oname) + pobj["nztag"] \
                    + '-' + mode
            if v:
                print(f' >>> {fname}: {pobj["title"]}')
                j = 0
                for o in po[mode]:
                    j += 1
                    if j==1:
                        print(f'    ... plotting  {o["title"]}')
                    else:
                        print(f'    ... overlaying {o["title"]}')
            func = "plot_map_" + mode
            globals()[func](po[mode], dict, gconf, fconf, oconf, fname)
            ###


"""
 ##############################################
 ### LEVEL-2 ROUTINES FOR SPECIFIC PURPOSES ###
 ##############################################
 - called from LEVEL-1 routines
"""

#@f_timer
def plot_map_std(obj_list, dict, gconf, fconf, oconf, fname):
    """
    plot horizontal map for one object incl. overlay(s)
    """
    mode = "std"
    ###
    ### create sub-subdirectory on request
    ###
    ch    = dict["name"]
    s     = f_getn(fconf, ["subdirs"], False)
    ###
    nztag = f_getn(obj_list[0], ["nztag"], '')
    ###
    if s:
        sadd = gconf["subdir"] + '/' + ch + '_' + obj_list[0]["name"] \
               + nztag + '-' + mode
        os.makedirs(sadd, exist_ok=True)
        sadd += '/'
    else:
        sadd = ''
    ###
    ### construct name of output file including path
    ###
    outfile = sadd + fname + '.png'
    ###
    ### loop over overlays
    ###
    i = 0
    for obj in obj_list:
        i += 1
        if i == 1:
            overlay = False
            ooconf = oconf[mode]
        else:
            overlay = True
            ooconf = oconf[mode]["overlay"][i-2][obj["name"]]
            if ooconf is None:
                ooconf = {}
        ###
        ### select data from object
        ###
        dx = f_getn(ooconf, ["plotargs", "dx"], 1)
        dy = f_getn(ooconf, ["plotargs", "dy"], 1)
        select_x  = (slice(None, None, dx))
        select_y  = (slice(None, None, dy))
        select_xy = (slice(None, None, dx), slice(None, None, dy))
        ###
        if obj["lon"].ndim == 1:
            lon = obj["lon"][select_x]
        else:
            lon = obj["lon"][select_xy]
        if obj["lat"].ndim == 1:
            lat = obj["lat"][select_y]
        else:
            lat = obj["lat"][select_xy]
        qdata = obj["data"]
        # vector data?
        isvec = isinstance(qdata, list)
        if isvec:
            xdata = qdata[0][select_xy]
            ydata = qdata[1][select_xy]
            data  = qdata[2][select_xy]
        else:
            data = qdata[select_xy]
        ###
        ### select features for plotting from yaml-input
        ###
        # use full range of data as internal default:
        georange  = get_geo_range(obj["lon"], obj["lat"])
        plotrange = f_getn(ooconf, ["plotrange"], georange)
        ###
        plotargs = f_getn(ooconf, ["plotargs"], {}).copy()
        text     = f_getn(ooconf, ["text"], [])
        marker   = f_getn(ooconf, ["marker"], [])
        if not overlay:
            dpi      = f_getn(ooconf, ["dpi"], 100)
        ###
        norm, lognorm = get_norm(plotargs, data)
        if isvec:
            colorargs, cbarargs = get_colorargs(plotargs, data,
                                                'quiver', overlay)
        else:
            colorargs, cbarargs = get_colorargs(plotargs, data,
                                                'contourf', overlay)
        labfmt  = f_getn(plotargs, ["labfmt"], "%f")
        missing = f_getn(plotargs, ["missing"], False)
        ###
        ### ... always add map transformation
        colorargs["transform"] = ccrs.PlateCarree()
        ###
        ### plot
        ###
        if not overlay:
            ###
            fig, ax = init_fig(plotrange, pixels=2048, kml=False)
            #
            # set gridline
            #
            gl = ax.gridlines(draw_labels=True, linewidth=1., color="gray")
            #
            # only at the bottom and at the right hand side
            gl.top_labels = False
            gl.right_labels = False
            #
            # formatting to include °N/°E
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            #
            # size and color
            gl.xlabel_style = {'size': 12, 'color': 'gray'}
            gl.ylabel_style = {'size': 12, 'color': 'gray'}
            #
            #-- add coastlines, country border lines, and grid lines
            ax.coastlines()
            #
            #-- add title
            title = obj["title"] + ' -- ' + gconf["tstring"]
            ax.set_title(title, fontsize=12, fontweight='bold')
            #
            #-- create contour or vector plot
            if isvec:
                cnplot = ax.quiver(lon, lat, xdata, ydata, data, norm=norm,
                                   **colorargs)
            else:
                cnplot = ax.contourf(lon, lat, data, norm=norm, **colorargs)
            #
            ###
            ### add text(s)
            ###
            for t in text:
                #crs = ccrs.PlateCarree()
                ax.text(**t) #, transform = crs._as_mpl_transform(ax))
            ###
            ### add marker(s)
            ###
            for m in marker:
                x = m.pop("x", None)
                y = m.pop("y", None)
                if x is None or y is None:
                    continue
                plt.plot(x, y, **m)
            ### hatch missing values
            ###
            if missing:
                mask = np.ma.getmask(data)
                z1 = np.zeros(data.shape)
                z1[mask] = data.max()
                # use contourf() with proper hatch pattern and alpha value
                cs = ax.contourf(lon, lat, z1, 3,
                                 hatches=['', '....'],  alpha=0.0)
            #cnplot.cmap.set_under('white')
            #cnplot.cmap.set_over('black')
            #cnplot.changed()
            #
            #-- add colorbar
            cbar = plt.colorbar(cnplot, orientation='horizontal',
                                pad=0.08, shrink=0.80, **cbarargs)
            cbarlabel = obj["title"]
            cbar.set_label(cbarlabel, fontsize=12)
            cbar.ax.tick_params(labelsize=12)
            if not lognorm:
                plim = f_getn(plotargs, ["cbarpowerlim"], (-3,3))
                cbar.formatter.set_powerlimits(plim)
                cbar.formatter.set_useMathText(True)
                #
            cbar.update_ticks()
        else:
            ### OVERLAY
            if isvec:
                ovlplot = ax.quiver(lon, lat, xdata, ydata,
                                    **colorargs)
                atxt = 'arrows'
            else:
                ovlplot = ax.contour(lon, lat, data, **colorargs)
                atxt = 'contours'
                cl = plt.clabel(ovlplot, fmt=labfmt)
            cbarlabel = cbarlabel  + '\n' + f'{atxt} ({i-1}): ' + obj["title"]
            cbar.set_label(cbarlabel, fontsize=12)

    ### at end of loop
    #-- save graphic output to PNG file
    fig.savefig(outfile, bbox_inches='tight', dpi=dpi)
    plt.close(fig)
    return 0

#@f_timer
def plot_map_kml(obj_list, dict, gconf, fconf, oconf, fname):
    """
    plot horizontal map for one object incl. overlay(s)
    """
    mode = "kml"
    ###
    ### create sub-subdirectory on request
    ###
    v     = f_getn(gconf, ["verbose"], False)
    ch    = dict["name"]
    s     = f_getn(fconf, ["subdirs"], False)
    nztag = f_getn(obj_list[0], ["nztag"], '')
    ###
    if s:
        sadd = gconf["subdir"] + '/' + ch + '_' + obj_list[0]["name"] \
               + nztag + '-' + mode
        os.makedirs(sadd, exist_ok=True)
        sadd += '/'
    else:
        sadd = ''
    ###
    ### construct name of output file including path
    ###
    outfile = sadd + fname + '.png'
    cbfile  = sadd + 'color_bar.png'
    kmzfile = sadd + fname + '.kmz'
    ###
    ###
    ### loop over overlays
    ###
    ovl_flist = []
    i = 0
    for obj in obj_list:
        i += 1
        if i == 1:
            overlay = False
            ooconf = oconf[mode]
        else:
            overlay = True
            ooconf = oconf[mode]["overlay"][i-2][obj["name"]]
            if ooconf is None:
                ooconf = {}
        #ovlfile = sadd + fname + f'-ovl{i:03}-kml' + '.png'
        ovlfile = sadd + fname + f'-{obj["name"]}-kml' + '.png'
        ovl_flist.append(ovlfile)
        ###
        ### select data from object
        ###
        dx = f_getn(ooconf, ["plotargs", "dx"], 1)
        dy = f_getn(ooconf, ["plotargs", "dy"], 1)
        select_x  = (slice(None, None, dx))
        select_y  = (slice(None, None, dy))
        select_xy = (slice(None, None, dx), slice(None, None, dy))
        ###
        if obj["lon"].ndim == 1:
            lon = obj["lon"][select_x]
        else:
            lon = obj["lon"][select_xy]
        if obj["lat"].ndim == 1:
            lat = obj["lat"][select_y]
        else:
            lat = obj["lat"][select_xy]
        qdata = obj["data"]
        # vector data?
        isvec = isinstance(qdata, list)
        if isvec:
            xdata = qdata[0][select_xy]
            ydata = qdata[1][select_xy]
            data  = qdata[2][select_xy]
        else:
            data = qdata[select_xy]
        ###
        ### select features for plotting from yaml-input
        ###
        # use full range of data as internal default:
        georange  = get_geo_range(obj["lon"], obj["lat"])
        plotrange = f_getn(ooconf, ["plotrange"], georange)
        ###
        plotargs = f_getn(ooconf, ["plotargs"], {}).copy()
        if not overlay:
            dpi      = f_getn(ooconf, ["dpi"], 100)
        ###
        norm, lognorm = get_norm(plotargs, data)
        if isvec:
            colorargs, cbarargs = get_colorargs(plotargs, data,
                                                'quiver', overlay)
        else:
            colorargs, cbarargs = get_colorargs(plotargs, data,
                                                'pcolormesh', overlay)
        labfmt  = f_getn(plotargs, ["labfmt"], "%f")
        missing = f_getn(plotargs, ["missing"], False)
        ###
        ### plot
        ###
        if not overlay:
            ###
            fig, ax = init_fig(plotrange, pixels=1024*10, kml=True)
            #
            if isvec:
                cnplot = ax.quiver(lon, lat, xdata, ydata, data, norm=norm,
                                   **colorargs)
            else:
                #cnplot = ax.pcolormesh(lon, lat, data, norm=norm, **colorargs)
                cnplot = ax.pcolor(lon, lat, data, norm=norm, **colorargs)
            ###
            if missing:
                if v:
                    print(f' ... WARNING: hatching missing values not'
                          f' implemented for kml mode!')
                #mask = np.ma.getmask(data)
                #z1 = np.zeros(data.shape)
                #z1[mask] = data.max()
                ## use pcolor with proper hatch pattern and alpha value
                #plt.pcolor(lon, lat, z1, hatch='...', alpha=0., shading='auto')

            # switch off axes
            ax.set_axis_off()

            #-- save graphic output to PNG file
            fig.savefig(outfile, pad_inches=0.0, dpi=dpi)

            # sperate color-bar
            cbfig = plt.figure(figsize=(2.0, 4.0),
                               facecolor=None, frameon=False)
            cbax = cbfig.add_axes([0.0, 0.05, 0.1, 0.9])
            cb = cbfig.colorbar(cnplot, cax=cbax, **cbarargs)
            #
            cbarlabel = obj["title"]
            cb.set_label(cbarlabel, rotation=-90, color='w', labelpad=20)
            #
            # set colorbar tick color
            cb.ax.yaxis.set_tick_params(color='w')
            cb.ax.xaxis.set_tick_params(color='w')
            # set colorbar edgecolor
            cb.outline.set_edgecolor('w')
            # set colorbar ticklabels
            plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='w')
            plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='w')
            # formatter
            if not lognorm:
                cb.formatter.set_scientific(True)
                plim = f_getn(plotargs, ["cbarpowerlim"], (-2,2))
                cb.formatter.set_powerlimits(plim)
                cb.formatter.set_useMathText(True)
            cb.ax.yaxis.set_offset_position('left')
            cb.ax.yaxis.offsetText.set(color='w')
            ###
            cb.update_ticks()
            ###
            files = [outfile]
            names = [obj["name"]]
            des   = [obj["title"]]
        else:
            ### OVERLAY
            ovl_fig, ovl_ax = init_fig(plotrange, pixels=1024*10, kml=True)
            #
            if isvec:
                ovlplot = ovl_ax.quiver(lon, lat, xdata, ydata, **colorargs)
                atxt = 'arrows'
            else:
                ovlplot = ovl_ax.contour(lon, lat, data, **colorargs)
                atxt = 'contours'
                cl = plt.clabel(ovlplot, fmt=labfmt)

            # switch off axes
            ovl_ax.set_axis_off()
            title = obj["title"]
            cbarlabel = f'{atxt} ({i-1}): ' + title + ' \n' + cbarlabel
            cb.set_label(cbarlabel,
                         rotation=-90,
                         color='w', labelpad=20+(i-1)*8)

            #-- save graphic output to PNG file
            ovl_fig.savefig(ovlfile, pad_inches=0.0, dpi=dpi)
            plt.close(ovl_fig)
            ###
            files.append(ovlfile)
            names.append(obj["name"])
            des.append(title)

    ### at end of loop
    # close figure
    plt.close(fig)

    # finalize color bar
    cbfig.savefig(cbfile, transparent=False, format='png')
    plt.close(cbfig)

    ### timestamp must have format: YYYY-MM-DDThh:mm:ssZ
    timestamp = gconf["tstring"].replace(' UTC','Z').replace(' ','T')
    #print(f'{timestamp}')

    # create kmz file
    make_kml(plotrange, figs=files, colorbar=cbfile,
             kmzfile=kmzfile, name=names, description=des,
             timestamp = timestamp)

    # remote obsolete files contained now in kmz-file
    try:
        os.remove(cbfile)
    except:
        pass
    try:
        os.remove(outfile)
    except:
        pass
    for f in ovl_flist:
        try:
            os.remove(f)
        except:
            pass
    ###
    return 0

def make_kml(pr, figs, colorbar=None, **kw):
    """
    plot range: pr = (llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)

    """
    llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat = pr[0], pr[1], pr[2], pr[3]
    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    timestamp = kw.pop('timestamp',None)
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        #ground.name = kw.pop('name', 'unknown')[draworder]
        ground.name = kw.get('name', ['unknown','unknown'])[draworder-1]
        ground.color = kw.pop('color', '9effffff')
        ground.timestamp.when = timestamp
        ground.atomauthor = kw.pop('author', 'MESSy Consortium')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.get('description','unknown')[draworder-1]
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.name = "colorbar"
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)


def main():
    """
    check, if channel.yml can be read without syntax or format errors
    """
    with open('channel.yml') as f:
        #conf = yaml.load(f, Loader=yaml.FullLoader)
        docs = yaml.load_all(f, Loader=yaml.FullLoader)
        for doc in docs:
            f_print_nested(doc)

if __name__ == "__main__":
    main()
