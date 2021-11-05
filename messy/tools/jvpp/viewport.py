#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2020-11-01 18:37:04 sander>

# viewport: multiple plots on one page
# Rolf Sander, 2016-2019

##############################################################################

import matplotlib.pyplot as plt
import os

class viewport(object):
    """ Define viewports (subfigures)
    """

    # ------------------------------------------------------------------------

    @classmethod
    def init(cls, x, y, pdffile, xsize=None, ysize=None, VERBOSE=False):
        from matplotlib.backends.backend_pdf import PdfPages
        cls.x = x
        cls.y = y
        if (xsize is None):
            cls.xsize = 11.69
        else:
            cls.xsize = xsize
        if (ysize is None):
            cls.ysize = 8.27
        else:
            cls.ysize = ysize
        if (os.path.isfile(pdffile)):
            if VERBOSE: print('deleting old file: %s' % (pdffile))
            os.remove(pdffile)
        cls.pdffile = pdffile
        cls.current = 0 # current subfig number
        cls.pdf_pages = PdfPages(pdffile) # open pdf
        cls.startpage(xsize,ysize)

    @classmethod
    def next(cls): # switch to next subfigure
        if (cls.current == cls.x * cls.y):
            cls.newpage()
        ax = plt.subplot2grid((cls.y, cls.x), (cls.current // cls.x, cls.current % cls.x), 
          rowspan=1, colspan=1)
        cls.current += 1
        #print "%d %d %d" % (cls.current, cls.y, cls.x)
        return ax

    @classmethod
    def startpage(cls, xsize, ysize): # default is A4
        plt.figure(figsize=(xsize,ysize), dpi=100) # start new plot page

    @classmethod
    def finishpage(cls):
        plt.tight_layout()  # fits subplots together
        cls.pdf_pages.savefig() # saves current page and creates new page in pdf
        plt.close()

    @classmethod
    def newpage(cls):
        if (cls.current>0):
            cls.finishpage()
        cls.current = 0
        cls.startpage(cls.xsize,cls.ysize)

    @classmethod
    def exit(cls): # finish creation of pdf
        if (cls.current>0):
            cls.finishpage()
        cls.pdf_pages.close()   # close pdf
        plt.close('all')

    # ------------------------------------------------------------------------

    @staticmethod
    def scientificNotation(value, three=True):
        import numpy as np
        import math
        if (value == 0):
            return '0'
        else:
            e_real = np.log10(np.abs(value))
            e = math.floor(e_real)
            m = np.sign(value) * 10 ** (e_real - e)
            #print '%gE%d --> ' % (m, e),
            if (three):
                shift = e%3
                e -= shift
                m *= 10**(shift)
            #print '%gE%d' % (m, e)
            formatstring = '%gE%d' % (m, e)
            #print formatstring
            return formatstring

    @staticmethod
    # from http://matplotlib.org/api/dates_api.html
    def timeformat(value, pos=None):
        import matplotlib
        value = matplotlib.dates.num2date(value)
        if (pos == 0):
            fmt = '%b %d' # show month only on first tick
        else:
            fmt = '%d'
        formatstring = value.strftime(fmt) # http://strftime.org/
        return formatstring

##############################################################################
