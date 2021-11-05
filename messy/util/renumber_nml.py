#!/usr/bin/env python
# -*- coding: utf-8 -*- Time-stamp: <2018-11-21 00:37:36 sander>
# Rolf Sander, 2018

__version__ = '1.0'

import os, sys
import re # regexp

HELPTEXT = '''renumber the indices of the array elements in a 
namelist file to avoid overwriting any elements\n
example usage:\n  %s channel.nml ADD_REF''' % (os.path.basename(sys.argv[0]))

def evaluate_command_line_arguments():
    import argparse
    # see also: https://docs.python.org/2/howto/argparse.html
    #           https://docs.python.org/2/library/argparse.html
    parser = argparse.ArgumentParser(description=HELPTEXT,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    parser.add_argument('nmlfile', type=argparse.FileType('r'),
                        help='name of the namelist file, e.g., "channel.nml"')
    parser.add_argument('array', help='name of the array, e.g., "ADD_REF"')
    return parser.parse_args()

if __name__ == '__main__':

    args = evaluate_command_line_arguments()

    tmpfilename = 'tmp.nml'
    nmlfilename = args.nmlfile.name
    TMPFILE = open(tmpfilename,'w+')
    index = 1

    regexp = re.compile(r'(^[ ]*%s)\([0-9 ]*\)' % (args.array))
    for line in args.nmlfile:
        result = regexp.search(line)
        if result is not None:
            # result.group(1) contains the leading spaces and the array name
            replacement = r'%s(%3d)' % (result.group(1), index)
            line0 = line
            index += 1
            line = line.replace(result.group(0), replacement)
            if (line!=line0 and args.verbose):
                print 'BEFORE: %s' % (line0),
                print 'AFTER:  %s' % (line)
        print >> TMPFILE, line,

    args.nmlfile.close()
    TMPFILE.close()

    os.rename(nmlfilename, nmlfilename+'~')
    os.rename(tmpfilename, nmlfilename)
