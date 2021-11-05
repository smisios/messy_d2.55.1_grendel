#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2020-11-06 18:11:44 sander>

# rstools
# Rolf Sander, 2018-2019

##############################################################################

import os, sys, shutil
assert sys.version_info >= (3, 6)
from glob import glob
import re # regexp
from ast import literal_eval
import decimal

HLINE  = '-' * 78
HLINE2 = '*' * 78

##############################################################################

def evaluate_config_file(configfile, section=None, allow_no_value=False, verbose=False):
    # https://docs.python.org/3.6/library/configparser.html
    import configparser
    config = configparser.ConfigParser(inline_comment_prefixes=('#'),
                                       allow_no_value=allow_no_value)
    # make config case-sensitive:
    # https://stackoverflow.com/questions/1611799/preserve-case-in-configparser
    config.optionxform = str
    # use 'read_file' because 'read' ignores it if ini file doesn't exist:
    with open(configfile) as f:
        config.read_file(f)
    if (section):
        if (verbose):
            print(f'The section [{section}] of the config file {configfile} contains:')
            for key in config[section]:
                print(f'{key:16}= {config[section][key]}')
            print()
        return config[section]
    else:
        return config

##############################################################################
    
def fileselector(basedir, wildcard, text=None):
    pathname = basedir
    while os.path.isdir(pathname):
        dirlist = sorted([os.path.join(pathname,d)+'/' for d in os.listdir(pathname)
                          if os.path.isdir(os.path.join(pathname,d))])
        dirlist += sorted(glob(pathname+'/'+wildcard))
        if (pathname != basedir):
            dirlist = ['..'] + dirlist
        dir_dict = {ind: value for ind, value in enumerate(dirlist)}
        print(f'Directory: {pathname}')
        for key in dir_dict:
            print('%2s) %s' % (key+1, dir_dict[key].replace(basedir+'/','')))
        if (text):
            print(text)
        try:
            selection = dir_dict[int(input())-1]
        except:
            return None
        if (selection=='..'):
            pathname = os.path.realpath(pathname+'/..')
        else:
            pathname = os.path.realpath(selection)
        if (not os.path.isdir(pathname)):
            break
    return pathname

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
 
def grep_i(searchstring, filewildcard, only_matching=False):
    # suffix _i indicates grep option '-i' for case insensitive
    results = []
    allfiles = glob(filewildcard)
    for onefile in allfiles:
        with open(onefile) as f:
            for line in f:
                result = re.search(searchstring, line, re.IGNORECASE)
                if (result):
                    if (only_matching):
                        results.append(result.group(0)) # return only match (grep -o)
                    else:
                        results.append(line.rstrip('\n')) # return complete line
    return results
    
##############################################################################

def runcmd(cmd, logfilename=None, verbosity=1, check=False, pipe=False, env=None):
    import subprocess
    if (logfilename):
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
    else:
        # - verbosity is not implemented for logfilename=None
        # - in python-3.6, universal_newlines=True must be used. Since
        #   python-3.7, text=True is available
        # output.stdout     -> normal output
        # output.stderr     -> error output
        # output.returncode -> exit status
        if (pipe): pipe = subprocess.PIPE
        return subprocess.run(cmd, shell=True, executable='/bin/tcsh', 
            stdout=pipe, stderr=pipe,
            universal_newlines=True, check=check, env=env)

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
#     from datetime import datetime
#     if (os.getenv('TRASH')):
#         # $TRASH exists, move old data to trash directory:
#         trashsubdir = os.getenv(
#             'TRASH') + '/' + datetime.now().strftime('%Y-%m-%d-%H:%M:%S')
#         os.mkdir(trashsubdir)
#         shutil.move(filename, trashsubdir)
#     else:
#         sys.exit('ERROR')

##############################################################################

if __name__ == '__main__':

    # test tail:
    # example usage: ./rstools.py myfile 7
    tail(sys.argv[1], sys.argv[2])

