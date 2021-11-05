#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 14:08:18 sander>

# rstools
# Rolf Sander, 2018

##############################################################################

# import os, sys, shutil
import graph_tool.all as gt

HLINE  = '-' * 78
HLINE2 = '*' * 78

##############################################################################

def n2v(g, species_name): # n2v = name to vertex
    vertex = gt.find_vertex(g, g.vp.name, species_name)
    if (len(vertex)==1):
        return vertex[0] # return a vertex
    else:
        return -1

##############################################################################
