#!/usr/bin/tcsh

################################################################################
# ## This script stores the reference structure for a given source file.       #
# ## It uses forcheck to generate the reference structure in both textual and  #
# ## xml format. For the textual format, the reference_structure part of the   #
# ## extensive Forcheck output is extracted.                                   #
# ##                                                                           #
# ## NOTE: Please do not copy this script anywhere. If you need to reference   #
# ##       it from a different location, use a symbolic link. Otherwise it     #
# ##       will not work properly.                                             #
# ##                                                                           #
# ## Usage    : `calltree.tcsh infile outfile'                                 #
# ## Arguments: infile   source file (*.f90) for which to generate reference   #
# ##                     structure                                             #
# ##            outfile  basename of output files in which reference structure #
# ##                     will be stored (outfile.txt and outfile.xml)          #
# ##                                                                           #
# ## Author: Christopher Kaiser, DLR, 2012 (christopher.kaiser@dlr.de)         #
################################################################################

# Check for correct invocation
if ( "$2" == "" ) then
    echo " Usage: `basename $0` <src.f90> <outfile>"
    exit 1
endif

# Check if Forcheck is available
if ( `which forchk` == "" ) then
    echo " Error: Forcheck executable not found."
    exit 2
endif

# Get path to messy/util directory
if ( -l $0 ) then
    set back=`pwd`
    cd `dirname $0`
    set sname=`basename $0`
    set lname=`readlink $sname`
    cd `dirname $lname`
    set p_util=`pwd`
    cd $back
else
    set p_util=`dirname $0`
endif
if ( `echo $p_util | grep "messy/util"` == "" ) then
    echo $p_util
    echo " Error: MESSy root directory not found."
    exit 2
endif

# Get list of required source files
echo " Assembling list of required source files ..."
set srclist=`${p_util}/get_srclist.sh $1`
if ( "$srclist" == "" ) then
    echo " Error: Failed to find source files."
    exit 3
endif
echo " ... done."

# Write header in outfile.txt
echo " Running forchk and extracting reference_structure ..."
echo "## This file contains the reference structure for" > ${2}.txt
echo "## ${1}." >> ${2}.txt
set this=`basename $0`
echo "## It was created by '$this $1 ${2}'." >> ${2}.txt
echo "" >> ${2}.txt

# Run Forcheck and extract reference_structure from output
forchk -l - -shref -refstruc ${2}.xml -ff $srclist | \
    sed -e 's///' -e '/^$/ d' -e '/FORCHECK .*page:/ d' | \
    sed -n '/reference_structure/,/^[^ \n]/ p' | sed '$ d' \
    >> ${2}.txt
set error=$?
# (if Forcheck encountered a fatal error, clean up)
if ( $error == 16 ) then
    rm ${2}.txt
    echo " Error: execution of forchk failed."
    exit $error
endif
echo " ... done."

exit 0
