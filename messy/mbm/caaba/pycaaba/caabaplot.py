#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-07 13:49:03 sander>

# caabaplot: plot results from CAABA
# Rolf Sander, 2017

##############################################################################

from math import sin
from viewport import viewport
from netCDF4 import Dataset, num2date
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from cycler import cycler
import os
import sys
from matplotlib.dates import AutoDateFormatter, AutoDateLocator, DateFormatter
from mecca import mecca
# import re # regexp

##############################################################################

import _mecca_spc # created automatically by spc2mpl
spc_names = _mecca_spc.spc_names() # load dictionary
# add family names:
spc_names['NOx'] = r'NO_x'
spc_names['NOy'] = r'NO_y'
spc_names['Clx'] = r'Cl_x'
spc_names['Brx'] = r'Br_x'
spc_names['Ix']  = r'I_x'
spc_names['RGM'] = r'RGM'

def define_family(family, ncvar, tmin, tmax):
    def add_species(pd, units, factor, species):
        if (species in ncvar):
            pd = pd + factor * ncvar[species][tmin:tmax,0,0,0]
            units = ncvar[species].units
        return pd, units
    pd    = 0 # plotdata
    units = False
    if (family=='NOx'): # NOx
        pd, units = add_species(pd, units, 1, 'NO')
        pd, units = add_species(pd, units, 1, 'NO2')
    if (family=='NOy'): # NOy
        pd, units = add_species(pd, units, 1, 'NO')
        pd, units = add_species(pd, units, 1, 'NO2')
        pd, units = add_species(pd, units, 1, 'NO3')
        pd, units = add_species(pd, units, 2, 'N2O5')
        pd, units = add_species(pd, units, 1, 'HONO')
        pd, units = add_species(pd, units, 1, 'HNO3')
        pd, units = add_species(pd, units, 1, 'HNO4')
        pd, units = add_species(pd, units, 1, 'PAN')
        pd, units = add_species(pd, units, 1, 'ClNO2')
        pd, units = add_species(pd, units, 1, 'BrNO2')
        pd, units = add_species(pd, units, 1, 'ClNO3')
        pd, units = add_species(pd, units, 1, 'BrNO3')
    if (family=='Clx'): # total reactive chlorine (excludes HCl)
        pd, units = add_species(pd, units, 1, 'Cl')
        pd, units = add_species(pd, units, 1, 'ClO')
        pd, units = add_species(pd, units, 1, 'HOCl')
        pd, units = add_species(pd, units, 2, 'Cl2O2')
        pd, units = add_species(pd, units, 1, 'ClNO2')
        pd, units = add_species(pd, units, 1, 'ClNO3')
        pd, units = add_species(pd, units, 2, 'Cl2')
        pd, units = add_species(pd, units, 1, 'OClO')
        pd, units = add_species(pd, units, 1, 'BrCl')
        pd, units = add_species(pd, units, 1, 'ICl')
    if (family=='Brx'): # total reactive bromine (excludes HBr)
        pd, units = add_species(pd, units, 1, 'Br')
        pd, units = add_species(pd, units, 1, 'BrO')
        pd, units = add_species(pd, units, 1, 'HOBr')
        pd, units = add_species(pd, units, 1, 'BrNO2')
        pd, units = add_species(pd, units, 1, 'BrNO3')
        pd, units = add_species(pd, units, 2, 'Br2')
        pd, units = add_species(pd, units, 1, 'BrCl')
        pd, units = add_species(pd, units, 1, 'IBr')
    if (family=='Ix'): # total reactive iodine (excludes HI)
        pd, units = add_species(pd, units, 1, 'I')
        pd, units = add_species(pd, units, 1, 'IO')
        pd, units = add_species(pd, units, 1, 'HOI')
        pd, units = add_species(pd, units, 2, 'I2O2')
        pd, units = add_species(pd, units, 1, 'HIO3')
        pd, units = add_species(pd, units, 1, 'INO2')
        pd, units = add_species(pd, units, 1, 'INO3')
        pd, units = add_species(pd, units, 2, 'I2')
        pd, units = add_species(pd, units, 1, 'ICl')
        pd, units = add_species(pd, units, 1, 'IBr')
        pd, units = add_species(pd, units, 1, 'OIO')
    if (family=='RGM'): # reactive gaseous mercury
        pd, units = add_species(pd, units, 1, 'HgO')
        pd, units = add_species(pd, units, 1, 'HgCl')
        pd, units = add_species(pd, units, 1, 'HgCl2')
        pd, units = add_species(pd, units, 1, 'HgBr')
        pd, units = add_species(pd, units, 1, 'HgBr2')
        pd, units = add_species(pd, units, 1, 'ClHgOBr')
        pd, units = add_species(pd, units, 1, 'ClHgBr')
    return pd, units

##############################################################################

class caabaplot(object):
    """ plotting functions for CAABA
    """

    @classmethod
    def plot_0d(cls, modelruns, species, pagetitle, plottitle,
                ncfilename='caaba_mecca.nc', timeformat='', scalefactor=1.,
                tmin=0, tmax=0):
        # print 'modelruns = %s' % (modelruns)
        # print 'species = %s' % (species)
        # print 'pagetitle = %s' % (pagetitle)
        # print 'plottitle = %s' % (plottitle)
        linecolors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']
        ax = viewport.next()
        ax.set_prop_cycle(cycler('color', linecolors))
        baserun = modelruns[0][1]
        if ((viewport.current == 1) and (pagetitle)):
            # on new page, start with legend on a dummy plot:
            for (modelrundir, modelrunname) in modelruns:
                lines = plt.plot([0,0], label=modelrunname)
                if (modelrunname != baserun):
                    plt.setp(lines, linestyle='dotted', linewidth=3)
            plt.axis('off')
            legend = plt.legend(loc='center',
                                mode='expand',
                                fontsize = 'small',
                                title=pagetitle,
                                fancybox=True,
                                shadow=True,
                                borderaxespad=0.)
            plt.setp(legend.get_title(),fontsize='large')
            ax = viewport.next()
            ax.set_prop_cycle(cycler('color', linecolors))
        # loop over all model runs:
        if (tmax>0):
            fulltrange = False # plot only time[tmin:tmax]
        else:
            fulltrange = True # plot full time range
            tmin = 0
        for (modelrundir, modelrunname) in modelruns:
            ncfullfilename = modelrundir+'/'+ncfilename
            ncid = Dataset(ncfullfilename)
            # define time:
            time = ncid.variables['time']
            if (fulltrange):
                tmax = len(time)
            t = num2date(time[tmin:tmax],time.units)
            # plot data:
            units = False
            if (species in ncid.variables):
                # 'species' occurs in current modelrun:
                plotdata = ncid.variables[species][tmin:tmax,0,0,0]
                units = ncid.variables[species].units
            else:
                # check if 'species' refers to a family:
                plotdata, units = define_family(species, ncid.variables, tmin, tmax)
            if (units):
                lines = plt.plot(t, scalefactor*plotdata, label=modelrunname)
                plt.ylabel(units)
                if (modelrunname != baserun):
                    plt.setp(lines, linestyle='dotted', linewidth=2)
            else:
                # create dummy plot in order to cycle to the next color:
                lines = plt.plot(t, 0.*time[tmin:tmax], linewidth=0)
            ncid.close()
        plt.title(plottitle)
        plt.xlabel('time')
        ax.grid(True) # hoizontal and vertical gridlines in plot
        # x-axis:
        ax.xaxis_date() # x-axis is a date
        xtick_locator = AutoDateLocator()
        # define locations of ticks on time axis:
        ax.xaxis.set_major_locator(xtick_locator) # automatic
        #ax.xaxis.set_major_locator(plt.MaxNLocator(5)) # max number of tick intervals
        # define format of ticks on time axis:
        if (timeformat):
            xformatter = DateFormatter(timeformat)
        else:
            xformatter = ticker.FuncFormatter(viewport.timeformat) # function
        #xformatter = AutoDateFormatter(xtick_locator)
        ax.xaxis.set_major_formatter(xformatter)
        # y-axis:
        # adjust yrange:
        # print 'yrange:', plt.ylim()[0], plt.ylim()[1]
        # plt.ylim(plt.ylim()[0] * 0.9, plt.ylim()[1] * 1.1)
        # show all digits instead of using offset for y-axis:
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
        yformatter = ticker.FuncFormatter(
            lambda x, p: viewport.scientificNotation(x))
        ax.yaxis.set_major_formatter(yformatter)
        ax.yaxis.set_major_locator(plt.MaxNLocator(5)) # max number of tick intervals

    @classmethod
    def xxxg(cls, modelruns, plotspecies, pdffile='xxxg',
             pagetitle='Gas-phase', timeformat='', tmin=0, tmax=0):
        print('Plotting these model runs:')
        for modelrun in modelruns:
            print('  %s' % (modelrun))
        viewport.init(4, 4, pdffile+'.pdf', 17, 8) # open pdf
        viewport.newpage()
        print('Plotting these species:')
        for species in plotspecies: # species loop
            print('%s' % (species), end=' ') ; sys.stdout.flush()
            # define plottitle:
            # try to obtain "spc_names[species]", if undefined, use "species":
            plottitle = r'$\sf ' + spc_names.get(species, species) + r'$'
            # mz_rs_20180722+
            # modifying plottitle for aq with regexp is not necessary anymore
            # because it is done in spc2mpl.awk now.
            # regexp = re.compile('_a([0-9][0-9])')
            # result = regexp.search(species)
            # if result is not None: # aq species:
            #     species2 = species.replace(result.group(0),'_a##')
            #     plottitle = r'$\sf ' + spc_names[species2] + r'$'
            #     plottitle = plottitle.replace('\\aq', '(a'+result.group(1)+')')
            # mz_rs_20180722-
            cls.plot_0d(modelruns, species, pagetitle, plottitle,
                        'caaba_mecca.nc', timeformat, tmin=tmin, tmax=tmax)
        print('\nTo show the results, type:\n  qpdfview %s.pdf &\n' % (pdffile))
        viewport.exit() # close pdf

    @classmethod
    def jval(cls, modelruns, plotspecies, pdffile='jval',
             pagetitle='J-values', timeformat='', tmin=0, tmax=0):
        print('Plotting these model runs:')
        for modelrun in modelruns:
            print('  %s' % (modelrun))
        viewport.init(4, 4, pdffile+'.pdf', 17, 8) # open pdf
        viewport.newpage()
        print('Plotting these J-values:')
        for species in plotspecies: # species loop
            print('%s' % (species), end=' ') ; sys.stdout.flush()
            plottitle = species
            cls.plot_0d(modelruns, species, pagetitle, plottitle,
                        'caaba_jval.nc', timeformat, tmin=tmin, tmax=tmax)
        viewport.exit() # close pdf
        print('\nTo show the results, type:\n  qpdfview %s.pdf &\n' % (pdffile))

    @classmethod
    def rxnrates(cls, modelruns, plotrxns, pdffile='rxnrates',
                 pagetitle='Reaction rates', timeformat=''):
        import _rxnrates # created automatically by rxn2mpl
        print('Plotting these model runs:')
        for modelrun in modelruns:
            print('  %s' % (modelrun))
        viewport.init(4, 4, pdffile+'.pdf', 17, 8) # open pdf
        viewport.newpage()
        print('Plotting these rxns:')
        eqn_names = _rxnrates.eqn_names() # load dictionary
        for rxn in plotrxns: # rxn loop
            print('%s' % (rxn), end=' ') ; sys.stdout.flush()
            plottitle = rxn + ': ' + eqn_names[rxn]
            cls.plot_0d(modelruns, 'RR'+rxn, pagetitle, plottitle,
                        'caaba_mecca_rr.nc', timeformat)
        viewport.exit() # close pdf
        print('\nTo show the results, type:\n  qpdfview %s.pdf &\n' % (pdffile))

    @classmethod
    def rxnrates_scaled(cls, modelruns, plotspecies, pdffile='rxnrates_scaled',
                        timeformat=''):
        from _rxnrates_scaled import rxns # created automatically by rxn2mpl
        print('Plotting these model runs:')
        for modelrun in modelruns:
            print('  %s' % (modelrun))
        viewport.init(4, 4, pdffile+'.pdf', 17, 8) # open pdf
        for species in plotspecies:
            print('*** plotting', species)
            viewport.newpage()
            pagetitle='%s prod+loss' % (species)
            for (fct, txt, spc, tex) in rxns:
                if (spc==species):
                    plottitle = '%s: %s (%+d)' % (txt, tex, fct)
                    print(plottitle)
                    cls.plot_0d(modelruns, 'RR'+txt, pagetitle, plottitle,
                                'caaba_mecca_rr.nc', timeformat,
                                scalefactor=fct)
        viewport.newpage()
        viewport.exit() # close pdf
        print('\nTo show the results, type:\n  qpdfview %s.pdf &\n' % (pdffile))

##############################################################################

def makeplots_xxxg(modelruns, option=2):
    mytimeformat=''
    #mytimeformat='%-Hh' # %H:%M:%S https://docs.python.org/2/library/datetime.html
    #option = 1 # plot all species that are included in the nc file from first modelrun
    #option = 2 # plot the species listed in mecca.py up to a certain verbosity
    #option = 3 # plot only selected species
    if (option==1):
        plotspecies = []
        ncid = Dataset(modelruns[0][0]+'/caaba_mecca.nc')
        for var in sorted(ncid.variables):
            if (ncid.variables[var].ndim==4): # exclude lon, lat, lev, time
                plotspecies.append(var)
    elif (option==2):
        plotspecies = mecca.set_species(2)
    elif (option==3):
        plotspecies = ['O3', 'H2O2', 'NO', 'NO2', 'NOx', 'CH4', 'C2H2',
                       'HO12CO3C4', 'APINENE', 'BrO', 'DMS', 'SO2', 'Hg']
    caabaplot.xxxg(modelruns, plotspecies, 'output/xxxg',
                   timeformat=mytimeformat)

##############################################################################

def makeplots_jval(modelruns, option=2):
    mytimeformat=''
    #mytimeformat='%-Hh'
    #option = 1 # plot all species that are included in the nc file from first modelrun
    #option = 2 # plot the species listed in mecca.py up to a certain verbosity
    #option = 3 # plot only selected species
    if (option==1):
        plotspecies = []
        ncid = Dataset(modelruns[0][0]+'/caaba_jval.nc')
        for var in sorted(ncid.variables):
            if (ncid.variables[var].ndim==4): # exclude lon, lat, lev, time
                plotspecies.append(var)
    elif (option==2):
        plotspecies = mecca.set_jvalues(2)
    elif (option==3):
        plotspecies = ['J_O2', 'J_O3P', 'J_O1D', 'J_H2O2', 'J_NO2',
                       'J_NO2O', 'J_NOO2', 'J_N2O5', 'J_HNO3', 'J_HNO4',
                       'J_PAN', 'J_HONO']
    caabaplot.jval(modelruns, plotspecies, 'output/jval',
                   timeformat=mytimeformat)

##############################################################################

def makeplots_rxnrates(option=2):
    #option = 1 # take list of rxns from the first model run to be plotted
    #option = 2 # plot only selected rxns
    if (option==1):
        plotrxns = []
        ncid = Dataset(rxnfilename)
        for rxn in ncid.variables:
            if (rxn[0:2]=='RR'):
                plotrxns.append(rxn[2:])
        ncid.close()
    elif (option==2):
        plotrxns = ['G1000', 'G1001', 'G2100']
    caabaplot.rxnrates(modelruns, plotrxns, 'output/rxnrates')

##############################################################################

def makeplots_rxnrates_scaled():
    plotspecies = ['O3', 'H2O2', 'NO2', 'CH4', 'C2H2']
    #plotspecies = ['O3', 'C2H2', 'H2O2', 'NO2', 'CH4']
    caabaplot.rxnrates_scaled(modelruns, plotspecies, 'output/rxnrates_scaled')

##############################################################################

if __name__ == '__main__':

    #-------------------------------------------------------------------------
    
    # select model runs:
    modelruns = mecca.set_modelruns()
    rxnfilename = modelruns[0][0]+'/caaba_mecca_rr.nc'

    #-------------------------------------------------------------------------
    # plot species:
    #makeplots_xxxg(modelruns,1) # plot all species
    makeplots_xxxg(modelruns)
    #-------------------------------------------------------------------------
    # plot reactions:
    if (os.path.isfile(rxnfilename)):
        makeplots_rxnrates()
        #makeplots_rxnrates_scaled()
    else:
        print('Skipping plots of reaction rates because %s does not exist.' % (rxnfilename))
    #-------------------------------------------------------------------------
    # plot J-values:
    jvalfilename = modelruns[0][0]+'/caaba_jval.nc'
    if (os.path.isfile(jvalfilename)):
        makeplots_jval(modelruns)
    else:
        print('Skipping plots of J values because %s does not exist.' % (jvalfilename))
    #-------------------------------------------------------------------------

##############################################################################
