#!/usr/bin/env python3

"""
program to plot contents of netcdf file using
- channel_nc2dict.py for conversion into nested dictionary structure
- channel.py (ususally used for on-line plotting via forpy)
- channel.yml (imported by channel.py to control everything)

@ author: Patrick Jöckel, DLR
"""

import argparse
import sys
from os.path import dirname
sys.path.append(dirname(__file__))
#print(f'{__file__}')

import channel_nc2dict as c2d
import channel as ch


"""
   define and parse command line paramters
"""
parser = argparse.ArgumentParser(
    description='command line wrapper for testing MESSy forpy plugins')
parser.add_argument('file', metavar='netcdf-file',
                    type=argparse.FileType('r'),
                    help='netcdf file with channel output')
parser.add_argument('-t','--time',metavar='time-step(s)',dest='tsteps',
                    help='time step[s]: start[,stop[,step]]')
args = parser.parse_args()

fname = args.file.name
tsteps = args.tsteps

#print(f'{fname}')
#print(tsteps)

tmax = c2d.get_timesteps(fname)

if tsteps is None:
    start = 1
    stop = tmax
    step = 1
else:
    ts = tsteps.split(",")
    if len(ts) == 1:
        start = int(ts[0])
        stop  = start
        step  = 1
    if len(ts) == 2:
        start = int(ts[0])
        stop  = min(int(ts[1]),tmax)
        step  = 1
    if len(ts) == 3:
        start = int(ts[0])
        stop  = min(int(ts[1]),tmax)
        step  = int(ts[2])

### limit to available time steps
start = min(start, tmax)
stop  = min(stop, tmax)

#print(f'{start} {stop} {step}')
#for i in range(start,stop+1,step):
#    print(f'{i}')

for i in range(start,stop+1,step):
    d = c2d.create_dict(fname,i)
    status = ch.fpy_main(d)
    if status != 0:
        print(f'status = {status}')
        break
