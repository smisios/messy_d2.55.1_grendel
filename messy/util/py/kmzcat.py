#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 06 14:00:00 2021

This script is for concatenating the kmz-outputfiles of one channel/object
(as produced with channel.py/.yml) into one kmz file containing a time series.

@author: Patrick Joeckel, DLR, Apr 2021
"""

import argparse
from pathlib import Path
from zipfile import ZipFile
from pykml import parser as kmlparser
from pykml.factory import KML_ElementMaker as KML
from lxml import etree

"""
   define and parse command line paramters
"""
parser = argparse.ArgumentParser(
    description='python scritp to concatenate kmz files with single time stamp to time series kmz file')

parser.version = '1.0'
parser.add_argument('-v', '--version', action='version')
parser.add_argument('-o','--output',metavar='FILE',dest='outfile',
                    default = 'out.kmz', type=str,
                    help='output file name')
parser.add_argument('-d','--directory',metavar='DIRECTORY',dest='directory',
                    default = '.', type=str,
                    help='directory with kmz files')
parser.add_argument('-n','--nametag',metavar='TAG',dest='nametag',
                    default='', type=str,
                    help='select kmz files with this name-tag in filename')
parser.add_argument('-x','--exclude',metavar='XTAG',dest='extags',
                    action='append', default=[], type=str,
                    help='exclude kmz files with this name-tag in filename (can be used multiple times)')
parser.add_argument('-R', '--remove-original', dest='remorg',
                    action='store_true', default=False,
                        help="remove original kmz files")

args = parser.parse_args()

### always use absolute path
directory = Path(args.directory).resolve()
outfile = Path(args.outfile).resolve()

### always append .kmz
if not outfile.suffix == '.kmz':
    outfile = Path(args.outfile + '.kmz').resolve()

nametag   = args.nametag
extags    = args.extags
remorg    = args.remorg

print(f'DIRECTORY            : {directory}')
print(f'NAME TAG             : {nametag}')
print(f'EXCLUDE TAG(s)       : {extags}')
print(f'REMOVE ORIGINAL FILES: {remorg}')
print(f'OUTPUT FILE          : {outfile}')
#exit()

"""
get list of kmz files
"""

### check, if directory exists
if not directory.is_dir():
    print(f'ERROR: directory {directory} does not exist')
    exit()

### get all files in directory and select .kmz files
flist0 = [x for x in directory.glob('*.kmz') if x.is_file()]
if len(flist0) == 0:
    print(f'ERROR: no *.kmz files found in {directory}')
    exit()

### exclude files with specified name-tags
exclist = [i for i in flist0
           if any(i for j in extags if str(j) in i.name)]
#print(exclist)
#exit()
flist = [i for i in flist0 if i not in exclist]
flist.sort()
if len(flist) == 0:
    print(f'ERROR: no *.kmz files found in {directory} after removal of tag(s) {extags}')
    exit()
#print(flist)
#exit()

"""
initilize new kml-file
"""

### define namespace
nsmap = {"kml":  "http://www.opengis.net/kml/2.2",
         "atom": "http://www.w3.org/2005/Atom",
         "gx":   "http://www.google.com/kml/ext/2.2"}

### open new zipfile
newkmz = ZipFile(outfile, mode='w', allowZip64=True)

newkml = KML.kml()
newdoc = KML.Document()
newkml.append(newdoc)

"""
loop over all kmz files and create new document
"""

i = 0
for kmzfile in flist:
    ### open zip file
    print(f' ... processing {kmzfile}')
    kmz = ZipFile(kmzfile, 'r')
    with kmz.open('doc.kml', 'r') as f:
        root = kmlparser.parse(f).getroot().Document
        i += 1
        ### copy all files, except for doc.kml to new kmz file
        glist = [(s, kmz.read(s)) for s in kmz.namelist()]
        for g in glist:
            fileNameAndPath = g[0]
            actualFile = g[1]
            if fileNameAndPath == 'doc.kml':
                continue
            newkmz.writestr(fileNameAndPath, actualFile)

    ### take these elements from 1st file
    if i == 1:
        ### style information
        style = KML.Style(KML.ListStyle(KML.listItemType('checkHideChildren')),
                          id='hide')
        newdoc.append(style)
        ### first step: camera
        for c in root.findall("kml:Camera",nsmap):
            newdoc.append(c)
        ### second step: ScreenOverlays
        for s in root.findall("kml:ScreenOverlay",nsmap):
            newdoc.append(s)

    ### append GroundOverlays from all files,
    ### but move them into separate folders
    list = root.findall("kml:GroundOverlay",nsmap)
    n = len(list)
    j = 0
    for o in list:
        j += 1
        ###
        ### get name of object
        ###
        name = o.name.text
        #print(f'{name}')
        ###
        ### Note: This is for the time being a hidden feature.
        ###
        if False:
            ###
            ### convert TimeStamp into TimeSpan from t-1 to t
            ###
            ### get time stamp of object
            tstamp = o.TimeStamp.when.text
            if i == 1:
                tstamp_m1 = tstamp
            print(f' ... from {tstamp_m1} to {tstamp}')
            ### delete TimeStamp
            o.remove(o.TimeStamp)
            ### add TimeSpan
            o.TimeSpan = KML.TimeSpan(KML.begin(tstamp_m1),
                                      KML.end(tstamp))
            if j == n:
                tstamp_m1 = tstamp
        ###
        ### check, if new folder is required
        ###
        found = False
        for folder in newdoc.findall("kml:Folder",nsmap):
            fname = folder.name.text
            #print(f'{fname}')
            if fname == name:
                found = True
                print(f' ... ... found folder {fname}')
                break
            #print(f'{found}')
        if not found:
            print(f' ... ... create new folder {name}')
            folder = KML.Folder(id=name)
            folder.name = KML.name(name)
            folder.description = KML.description(o.description.text)
            folder.styleUrl = KML.styleUrl('#hide')
            newdoc.append(folder)
        folder.append(o)

"""
write doc.kml into kmz-file
"""
print(f'WRITING doc.kml to {outfile}')
newkmz.writestr('doc.kml',
                etree.tostring(newkml, pretty_print=True,
                               encoding='utf8', method='xml').decode())
"""
close kmz-file
"""
newkmz.close()

"""
remove original files on request
"""
if remorg:
    for f in flist:
        p = Path(f)
        print(f' ... removing original file {f}')
        p.unlink()

print(f'Done!')
exit()
