#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2020-01-31 00:05:18 sander>

# rstools
# Rolf Sander, 2018-2019

##############################################################################

import os, sys, shutil
import subprocess
from glob import glob
import re # regexp
from ast import literal_eval
import decimal

HLINE  = '-' * 78
HLINE2 = '*' * 78

##############################################################################

# Florian Strunk, 2019

def isnumber(string):
    try:
        float(string)
        return True
    except ValueError:
        pass
    return False

def isfloat(x):
    return isinstance(literal_eval(x), float)

def isint(x):
    return isinstance(literal_eval(x), int)

# python doesn't support floating point ranges out of the box (?). this is a workaround
def drange(x, y, jump):
    x = decimal.Decimal(x)
    y = decimal.Decimal(y)
    while x < y:
        yield float(x)
        x += decimal.Decimal(jump)

##############################################################################

def corename(path):
    # remove both the directory name and the suffix:
    return os.path.splitext(os.path.basename(path))[0]

##############################################################################
 
def grep_i(searchstring, filewildcard):
    # suffix _i indicates grep option '-i' for case insensitive
    results = []
    allfiles = glob(filewildcard)
    for onefile in allfiles:
        with open(onefile) as f:
            for line in f:
                result = re.search(searchstring, line, re.IGNORECASE)
                if (result):
                    results.append(result.group(0))
    return results
    
##############################################################################

def runcmd(cmd, logfilename, verbosity=1):
    # verbosity = 0 --> no output
    # verbosity = 1 --> print command and logfile
    # verbosity = 2 --> print also HLINEs
    CMDLOGFILE = open(logfilename,'w+', 1)
    if (verbosity>0):
        print()
        if (verbosity>1): print(HLINE)
        print('%s > %s ' % (cmd, os.path.basename(logfilename)), end=' ')
        sys.stdout.flush() # print now, don't wait till cmd has finished
    exitstatus = subprocess.call(
        'time -p '+cmd, stdout=CMDLOGFILE, stderr=CMDLOGFILE, shell=True)
    if (verbosity>0): print('DONE')
    if (verbosity>1): print(HLINE)
    CMDLOGFILE.close()
    if (exitstatus != 0):
        print('\n\n%s' % (HLINE2))
        tail(logfilename, 20)
        print('%s\nsee: %s\nERROR: exitstatus = %d\n%s' % (
            HLINE2, logfilename, exitstatus, HLINE2))
        sys.exit(1)

##############################################################################

def cat(filename):
    with open(filename) as f:
        print(f.read())
    
##############################################################################

# https://gist.github.com/amitsaha/5990310
def tail(filename,lines):
    with open(filename) as f:
        content = f.read().splitlines()
    count = len(content)
    #print 'The file %s has %d lines.' % (filename, count)
    for i in range(max(0,count-int(lines)),count):
        print(content[i])

##############################################################################

# https://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
def rnd(x, sigdig=2):
    from math import log10, floor
    if (x==0.):
        return 0.
    else:
        return round(x, sigdig-int(floor(log10(abs(x))))-1)

##############################################################################

# - not tested yet
# - provide alternative if $TRASH does not exist?

# def rm(filename):
#     import datetime
#     if (os.getenv('TRASH')):
#         # $TRASH exists, move old data to trash directory:
#         trashsubdir = os.getenv(
#             'TRASH') + '/' + datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
#         os.mkdir(trashsubdir)
#         shutil.move(filename, trashsubdir)
#     else:
#         sys.exit('ERROR')

##############################################################################

if __name__ == '__main__':

    # test tail:
    # example usage: ./rstools.py myfile 7
    tail(sys.argv[1], sys.argv[2])

