#!/bin/sh

################################################################################
# ## This script is a supplement to calltree.tcsh. It creates the list of      #
# ## source files to be passed to forchk in calltree.tcsh.                     #
# ##                                                                           #
# ## NOTE: This script will not work properly if it is copied to another       #
# ##       location than messy/util!                                           #
# ##                                                                           #
# ## Author: Christopher Kaiser, DLR, 2012 (christopher.kaiser@dlr.de)         #
################################################################################

################################################################################
# Constants                                                                    #
################################################################################

p_rel=`dirname $0`
back=`pwd`
cd $p_rel
p_util=`pwd`
cd $back
p_messy=`echo $p_util | awk '{sub("/util", "", $1) ; print $1}'`


################################################################################
# Recursive function                                                           #
################################################################################

mk_srclist ()
{
    # Create list of dependencies
    deps=`${p_util}/sfmakedepend.pl --file=- $1 2>&1 1>/dev/null | \
	sed 's/\.\.\..*ule: //g'`
    # Process the dependency list elementwise
    for d in $deps
    do
	# Set filename of dependency (add .f90 to module names, use sfmakedepend
	# output directly otherwise)
	if [ "${d%\.*}" == "$d" ]
	then
	    fname="${d}.f90"
	else
	    fname="$d"
	fi
	# First look for dependency in current working directory;
	# if it is not found there, try the messy root directory
	here=`pwd`
	src_loc=`find $here -name $fname`
	if [ "$src_loc" == "" ]
	then
	    src_tmp=`find ${p_messy} -name $fname`
	else
	    src_tmp=$src_loc
	fi
	# Proceed only if dependency was found
	if [ "$src_tmp" == "" ]
	then
	    echo " Warning: file $fname not found." >&2
	else
	    found=""
	    # Process list of files that match the dependency
	    for s in $src_tmp
	    do
		# Follow any symbolic links and store actual file names
		if [ -L $s ]
		then
		    back=`pwd`
		    cd `dirname $s`
		    srcname=`basename $s`
		    hlp=`readlink $srcname`
		    cd `dirname $hlp`
		    src="`pwd`/`basename $hlp`"
		    cd $back
		else
		    src=$s
		fi
		# If actual file name is not the same as found previously,
		# reset $srclist and exit with error message
		if [ "$found" != "" -a "$found" != $src ]
		then
		    echo " Error: Found multiple files $fname." >&2
		    srclist=""
		    exit 2
		fi
		# Store name of the actual file associated with the dependency
		found=$src
	    done
	    # Add the dependency to $srclist if it is not already present and
	    # check for dependencies of the dependency
	    hit=`echo -e $srclist | grep $src`
	    if [ "$hit" == "" ]
	    then
		srclist="${srclist}${src}\n"
		srcdir=`dirname $src`
		back=`pwd`
		cd $srcdir
		mk_srclist `basename $src`
		cd $back
	    fi
	fi   # found dependency?
    done
}


################################################################################
# Program                                                                      #
################################################################################

# Check for correct invocation
if [ "$1" == "" ]
then
    echo " Usage: `basename $0` <source file>" >&2
    exit 1
fi

# Initialize list of required source files with the first argument
srclist="$1\n"
# Generate the rest of the list
mk_srclist $1

# Output list of required source files (with newlines)
echo -e $srclist

exit 0
