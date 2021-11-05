#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-04 14:57:43 sander>

# sample points:
# select a file with the sample points (w/o suffix '.nc'):
#samplepointfile = 'skeleton_samplepoints_small'
#samplepointfile = 'skeleton_samplepoints_mbl'
#samplepointfile = 'skeleton_samplepoints_30'
#samplepointfile = 'skeleton_samplepoints_30_loworg'
samplepointfile = 'skeleton_samplepoints_30_lowterp'
#samplepointfile = '' # if empty, create only plots

# drgep:
eps0         = 5E-4 # start value for epsilon_ep
eps_increase = 1.4  # factor for increasing epsilon_ep

# plotting:
plot_delta_skel   = 1 # 0=no plots, 1=all
plot_targets      = 1 # 0=no plots, 1=all
plot_samplepoints = 1 # 0=no plots, 1=all, 2=some
