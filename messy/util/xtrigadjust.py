#!/usr/bin/env python3
# coding: utf-8
# alternatively #!/usr/bin/python3

'''
Script to update RG_TRIG() entries in import.nml-files according to
desired start date (and for future releases end date). Only works for
nml files strictly compliant with MESSy input file naming convention:

   <provider>_<version>_<scenario>_<class>_<species>_
   <yyyy[mm][dd]-yyyy[mm][dd]>.nml

Version: 1.0 | Based on script 'xtrigadjust.tcsh' by Patrick Joeckel

Author: Johannes Emmerig, DLR, March 2020

'''

# import modules
from dateutil.relativedelta import relativedelta as rd
import fileinput as fi
import datetime as da
import sys
import os

''' OPEN FILE LOOP '''
while True:
    # input and define file
    file_ = input(
        'Enter setup and file with the following scheme \n> setup/file e.g. STRATOFLY/import_s2000   \n')

    pre_file = '../nml/'
    post_file = '.nml'

    file_ = pre_file + file_.strip() + post_file

    # open file
    try:
        f = open(file_)

    except(FileNotFoundError) as exc:
        print('Please retry and enter an existing file.')

    else:
        break

'''
 Explanation of logic in files (to be found in messy wiki. [__ signifies bold])

     In this example the import starts January 2000 (193,1,204,__1__'),
     is updated once a month (1,'months', 'first',0,) by one step (193,__1__,204,1)
     up to December 2016 (193,1,__204__,1),
     and then reset to January 2016 (**193**,1,204,1).
     The numbers denote the index along the time axis in the netCDF file.
'''

# ask start date
startDate = input('What is the required start date (Format YYYYMM)?\n> ')
startDate = da.date(int(startDate[0:4]), int(startDate[4:6]), 1)
print(startDate)

# ask end date (hardcoded to 21000101 until end date is required in following versions)
#untilDate = input('What is the required end date (Format YYYY-MM)?\n> ')
#untilDate = da.date(int(untilDate.split('-')[0]),int(untilDate.split('-')[1]),1)
untilDate = da.date(2100, 1, 1)

# start, until date check
if startDate >= untilDate:
    print('The end date should be after the start date')
    print('   > Start date: ' + startDate)
    print('   > End date:   ' + untilDate)

# create/overwrite new file
file_new_name = file_[:-4] + '_xtrigadjusted.nml'
f_new = open(file_new_name, 'w')

for line in f:
    '''
    Lines will be changed:
        - that start with RG_TRIG (not !RG_TRIG)
        - that fulfill the date format yyyy[mm][dd]
        -
    '''

    if line.startswith('RG_TRIG('):

        ''' READ LINES '''
        # RG_TRIG(...) = 1,'months','first',0, '...', 1801,1,1812,601, '...',
        update_frequency = line.split(sep=',')[1].strip("'")
        # RG_TRIG(...) = 1,'months','first',0, '...', 1801,1,1812,import_startDate, '...',
        import_startDate = line.split(sep=',')[8]
        # RG_TRIG(...) = 1,'months','first',0, '...', 1801,1,untilDate,601, '...',
        import_untilDate = line.split(sep=',')[7]
        # RG_TRIG(...) = 1,'months','first',0, '...', reset_to_date,1,1812,601, '...',
        reset_to_date = line.split(sep=',')[5]
        # RG_TRIG(...) = 1,'months','first',0, '...', 1801,time_step_size,1812,601, '...',
        time_step_size = line.split(sep=',')[6]

        ''' EXCEPTIONS '''
        # catch missing start and end (date) exceptions
        try:
            nml_file_start = line.split(sep='-')[-2].split('_')[-1]
            nml_file_end = line.split(sep='-')[-1].split('.')[0]
            nml_file_start = str(int(nml_file_start))
            nml_file_end = str(int(nml_file_end))
        except (IndexError, ValueError) as exc:
            f_new.write(line)
            continue

        # to catch exceptions e.g. YYYY without month
        try:
            # use year and month
            nml_file_startDate = da.date(
                int(nml_file_start[:4]), int(nml_file_start[4:]), 1)
            nml_file_endDate = da.date(
                int(nml_file_end[:4]), int(nml_file_end[4:]), 1)
        except BaseException:
            # use year only
            nml_file_startDate = da.date(int(nml_file_start[:4]), 1, 1)
            nml_file_endDate = da.date(int(nml_file_end[:4]), 1, 1)
        if len(nml_file_end) < 4 or len(nml_file_start) < 4:
            f_new.write(line)
            continue

        # lines that do not have to be altered in general
        if update_frequency == 'years':
            f_new.write(line)
            continue

        # lines that do not have to be altered in general as the same year is
        # repeated over and over
        if rd(nml_file_endDate, nml_file_startDate).years == 0:
            f_new.write(line)
            continue

        ''' CALCULATE NEW VALUES '''
        # get difference between dates
        diff_startDate_nml_file_startDate = rd(startDate, nml_file_startDate)
        diff_untilDate_nml_file_startDate = rd(untilDate, nml_file_startDate)
        diff_nml_endDate_startDate = rd(nml_file_endDate, nml_file_startDate)

        # recalculate difference to number of months
        new_import_startDate = diff_startDate_nml_file_startDate.years * \
            12 + diff_startDate_nml_file_startDate.months + 1
        new_import_untilDate = diff_nml_endDate_startDate.years * \
            12 + diff_nml_endDate_startDate.months + 1
        new_reset_to_date = new_import_untilDate + 1 - 12

        # in case there is no data for the input date, the last available year
        # will be set as import_startDate
        if new_import_startDate > new_import_untilDate:
            new_import_startDate = new_reset_to_date

        ''' EXPORT VALUES TO LINE '''
        # add changes in new line (order of cmds is important, first occurence gets replaced)
        newline = line.replace(reset_to_date, str(new_reset_to_date),1)
        newline = newline.replace(import_untilDate, str(new_import_untilDate),1)
        newline = newline.replace(import_startDate, str(new_import_startDate),1)

        # replace old line with new line
        line = newline

    # write lines to new file
    f_new.write(line)

# close files
f_new.close()
f.close()

# Conclusion and potentially helpful commands
print('\nCreated new file ' + file_new_name + ' for you.\n')
print('Further commands:   less ' + file_new_name)
print('                    vimdiff ' + file_new_name + ' ' + file_)
print('                    mv ' + file_new_name + ' ' + file_)
print('\n')

sys.exit()
