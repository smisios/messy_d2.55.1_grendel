# -*- Makefile -*-
# ----------------------------------------------
export
SHELL    = sh

# ----------------------------------------------
INSTALLDIR = ../../../bin
OUTPUT = NETCDF

include main.mk

install: all

include depend.mk
