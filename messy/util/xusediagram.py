#!/usr/bin/env python3

"""
This script creates a diagram of the Fortran USE - dependencies
for a specifid MESSy submodel. It utilizes the graphviz software
(http://www.graphviz.org/) for visualisation.
The script requires the depend.mk (make dependencies), more specifically
echam5/__*/depend.mk and messy/smcl/depend.mk (see input below).
It needs to be called (after configure and gmake) from the MESSy 
$BASEDIR with a submodel name as parameter.
The script will operate in the subdirectory workdir to prevent
results from being zipped.

The original tcsh-script was written by Rolf Sander (MPIC) and has been
adapted and converted to python for common usage by Patrick Joeckel (DLR).

Version 1.1 of 2019-09-12
"""

import argparse
import re
from graphviz import Digraph
#import gv

def f_read_dep(filelist):
    """
    read dependencies from files into dictionary
    """
    # empty dictionary of dependencies
    dep = {}
    for file in filelist:
        # open file
        f = open(file, 'r')
        # read and ignore one header line
        header1 = f.readline()
        # loop over lines and extract variables of interest
        for line in f:
            line = line.strip()
            # separate into columns
            columns = line.split()
            # first column withtout ':' at the end contains name of object
            name = columns[0].strip(':')
            # .o is also removed here, assuming that there is no dot in filename
            name = name.split('.')[0]
            # the rest of the columns contain a list with dependecies
            # - suffixes: .mod, .o, .inc
            # - relative paths
            d = columns[1:]
            # remove relative paths
            da = [i.split('/')[-1] for i in d]
            # remove also suffixes (assuming no dot in basename)
            das = [i.split('.')[0] for i in da]
            # now we build a dictionary of (key, value) pairs with
            #   key: object name
            #   value: list of dependencies
            if name in dep:
                # key already in dictionary: append additional dependencies
                # by converting into set and back into list, we automatically
                # remove doubles, that might have been produced by
                # the dependency maker
                #dep[name].extend(d)
                dep[name] = list(set(dep[name] + das))
            else:
                # key no yet in dictionary: create new entry
                dep[name] = das
                
        # close file
        f.close()
    return(dep)

def f_dep_subset(dep, pattern_list):
    # start with empty dictionary
    sdep = {}
    for key in dep:
        # copy key - value pair, if one of the patterns matches ...
        for pat in pattern_list:
            if pat.search(key):
                # ... either the key
                if key in sdep:
                    sdep[key] = list(set(sdep[key] + dep[key]))
                else:
                    sdep[key] = dep[key]
            else:
                # ... or one of the values
                for l in dep[key]:
                    if pat.search(l):
                        if key in sdep:
                            sdep[key] = list(set(sdep[key] + [l]))
                        else:
                            sdep[key] = [l]
    return sdep

def f_dep2dot(dep, pattern, name, color):
    dot = Digraph(name,
                   graph_attr={'rank': 'same'},
                   node_attr={'shape': 'box',
                              'color': color,
                              'style': 'filled,rounded'} )
    for key in dep:
        if pattern.search(key):
            dot.node(key)
        for l in dep[key]:
            if pattern.search(l):
                dot.node(l)
    return dot

""" 
   main program starts here
"""

"""
   define and parse command line paramters
"""
parser = argparse.ArgumentParser(
    description='create dependency diagram for MESSy submodel')
parser.add_argument('sm', metavar='submodel',
                    help='name of MESSy submodel')
args = parser.parse_args()

sm = args.sm.lower()
SM = args.sm.upper()
#print(sm)
#print(SM)

"""
   read files with dpendncy listings
"""
dep = f_read_dep(['echam5/__LINUX64/depend.mk','messy/smcl/depend.mk'])
#print(dep['master'])

"""
   select subset for specific submodel
"""
#pattern = re.compile("_" + sm + "_*")
print("_" + sm + "_")
pattern = re.compile("_" + sm + "(_.*|$)")
if sm == 'mecca':
    # special case for MECCA / polyMECCA
    padd = re.compile("_" + sm + "[0-9][0-9][0-9]" + "(_.*|$)")
    plist = [pattern, padd]
else:
    plist = [pattern]
dep_sm = f_dep_subset(dep, plist)
#print(dep_sm)
#print(dep_sm['messy_mecca_poly_si'])
#print(dep['messy_mecca'])
#print(dep_sm['messy_mecca'])
#quit()

"""
   convert it into dot graph
   Note: subgraphs need to be defined before edges
   https://graphviz.readthedocs.io/en/stable/manual.html
"""

# declare a new graph
dot = Digraph(SM + 'dependency graph',
              graph_attr = {'clusterrank': 'none',
                            'rankdir': 'RL',
                            'margin': '0.25',
                            'ratio': 'fill',
                            'fontname': 'Helvetica-Bold',
                            'fontsize': '36',
                            'newrank': 'true',
                            'splines': 'ortho'}
)

"""
   add subgraphs (nodes with specific names)   
"""

dot_bmil = f_dep2dot(dep_sm, re.compile('_bi$'), 'cluster_bmil', 'red')
dot.subgraph(dot_bmil)

dot_smil = f_dep2dot(dep_sm, re.compile('_si$'), 'cluster_smil', 'yellow')
dot.subgraph(dot_smil)

dot_e5 = f_dep2dot(dep_sm, re.compile('_e5$'), 'e5', 'purple')
dot.subgraph(dot_e5)

dot_smcl_main = f_dep2dot(dep_sm, re.compile('_main_.*(?<!_bi|_si|_e5)$'),
                          'cluster_main_smcl', 'orange')
dot.subgraph(dot_smcl_main)

dot_smcl =  f_dep2dot(dep_sm, re.compile(r'^((?!_main_).)*$'),
                      'cluster_smcl', 'grey')
dot.subgraph(dot_smcl)

#pat = re.compile(r'^((?!_main_).)*$')
#for key in dep_sm:
#    for l in dep_sm[key]:
#        if pat.search(l):
#            a = "yes"
#        else:
#            a = "no"
#        print(l + ' ' + a)
#
#quit()

"""
   add nodes and edges
"""
for key in dep_sm:
    #dot.node(key)
    for l in dep_sm[key]:
        dot.edge(l,key)

print(dot.source)
dot.render('workdir/' + sm + '.dot', view=True)  

quit()

