# -*- Makefile -*-

##############################################################################

INSTALLDIR = ../../../bin

##############################################################################

include main.mk

install: all

# list of dependencies (via USE statements)
include depend.mk
