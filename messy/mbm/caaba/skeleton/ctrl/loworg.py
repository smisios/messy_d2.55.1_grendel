#!/usr/bin/env python
# -*- coding: utf-8 -*- Time-stamp: <2019-05-03 19:44:16 sander>

# loworg: select only sample points where concentrations of organics are low

# sample points:
samplepointfile = 'skeleton_samplepoints_30_loworg'

# drgep:
eps0         = 5E-4 # start value for epsilon_ep
eps_increase = 1.4  # factor for increasing epsilon_ep

# plotting:
plot_delta_skel   = 1 # 0=no plots, 1=all
plot_targets      = 1 # 0=no plots, 1=all
plot_samplepoints = 1 # 0=no plots, 1=all, 2=some
