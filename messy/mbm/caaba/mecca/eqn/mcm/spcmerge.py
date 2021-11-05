#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 12:18:38 sander>

# spcmerge: merge gas.spc with another *.spc file
# Authors: Rolf Sander, Sebastian Tauer, Hartwig Harder (MPI Mainz, 2017)

##############################################################################

import re
import sys

#DEBUG = True
DEBUG = False

if __name__ == '__main__':

    if len(sys.argv)<=1:
        sys.exit('ERROR: Supply the name of an input *.spc file!')
    infilename = sys.argv[1]
    if len(sys.argv)<=2:
        sys.exit('ERROR: Supply the name of an output *.spc file!')
    outfilename = sys.argv[2]
    regexp = re.compile('^ *([A-z0-9_]+) *=')

    # put all species from gas.spc into the gasspcs array:
    GASSPCFILE = open('gas.spc')
    gasspc_data = GASSPCFILE.readlines()
    GASSPCFILE.close()
    gasspcs=[]
    for line in gasspc_data:
        result = regexp.search(line)
        if (DEBUG): print((result==False), line.strip())
        if result: # if not None
            if (DEBUG): print(result.group(1).upper())
            gasspcs.append(result.group(1).upper())
    if (DEBUG): print(gasspcs)
    # check all species in newspc:
    NEWSPCFILE = open(infilename)
    newspc_data = NEWSPCFILE.readlines()
    NEWSPCFILE.close()
    ADDNLFILE = open(outfilename,'w+')

    print('New species in %s:' % (infilename))
    for line in newspc_data:
        line=line.strip()
        result = regexp.search(line)
        if result: # if not None
            newspc = result.group(1)
            if (DEBUG): print('found spc : %s ' % newspc)
            if newspc.upper() not in gasspcs:
                print(newspc, end=' ')
                print(line, file=ADDNLFILE)
    print()
    ADDNLFILE.close()

##############################################################################
