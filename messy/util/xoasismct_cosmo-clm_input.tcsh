#!/bin/tcsh -f

### PREPARING INPUTDIR_OASIS3MCT FOR A SPECIFIC COSMO-CLM-SETUP ###

############################## 
### 0.STEP: SPIN-UP SIMULATION
##############################
### execute a simulation with the new NML_SETUP in the directory WORKDIR with
### the following restrictions in xmessy_mmd: 
### * CLM runs on 1 PE
### * INPUTDIR_OASIS3MCT/${NML_SETUP} does not yet exist
### especially for OASIS_CPL_MODE=AVERAGE (
###            oasis_restart_cosmo.nc and oasis_restart_clm.nc are required)
### * OASIS_LAG_COSMO = OASIS_LAG_CLM = +0
### * EXPORTED in namcouple is replaced by EXPOUT 

set script=`basename $0`
### INIT ###################################################################
set NML_SETUP = ''
set INPUTDIR_OASIS3MCT = ''
set WORKDIR = ''
set OASIS_CPL_MODE = ''

############################################
### 1.STEP: SETTINGS (same as in xmessy_mmd)
############################################

while ($# > 0)
    switch ($1)

      case '-h':
	echo " "
        echo "Usage: $script [-h] -m <mode> -s <setup> -i <directory> -w <directory>"
	echo " "
        echo "  -h  : show this help and exit"
        echo '  -m  : OASIS COUPLIMG MODE (mode = INSTANT or AVERAGE)'
        echo '  -s  : namelist setup'
        echo '  -i  : OASIS3MCT input directory'
        echo '  -w  : working directory'
	echo " "
	exit 1
        breaksw

      case '-m':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <mode> missing\!'
	    echo "Use $script -h for more information\!"
	    exit 1
	else
	    set OASIS_CPL_MODE = "$1"
	    shift
	endif
        breaksw

      case '-s':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <setup> missing\!'
	    echo "Use $script -h for more information\!"
	    exit 1
	else
	    set NML_SETUP = "$1"
	    shift
	endif
        breaksw

      case '-i':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <directory> missing\!'
	    echo "Use $script -h for more information\!"
	    exit 1
	else
	    set INPUTDIR_OASIS3MCT = "$1"
	    shift
	endif
        breaksw

      case '-w':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <directory> missing\!'
	    echo "Use $script -h for more information\!"
	    exit 1
	else
	    set WORKDIR = "$1"
	    shift
	endif
        breaksw
			
      default:
        echo "Unknown option: $1"
	echo "Use $script -h for more information\!"
	exit 1
        breaksw 
    endsw
end

if ("$OASIS_CPL_MODE" == "") then
   echo 'ERROR: -m option must be specified\!'
   echo "Use $script -h for more information\!"
   exit 1
endif
if (("$OASIS_CPL_MODE" != "INSTANT") && ("$OASIS_CPL_MODE" != "AVERAGE")) then
   echo 'ERROR: <mode> must be AVERAGE or INSTANT\!'
   echo "Use $script -h for more information\!"
   exit 1
endif

if ("$NML_SETUP" == "") then
   echo 'ERROR: -s option must be specified\!'
   echo "Use $script -h for more information\!"
   exit 1
endif

if ("$INPUTDIR_OASIS3MCT" == "") then
   echo 'ERROR: -i option must be specified\!'
   echo "Use $script -h for more information\!"
   exit 1
endif

if ("$WORKDIR" == "") then
   echo 'ERROR: -w option must be specified\!'
   echo "Use $script -h for more information\!"
   exit 1
endif

echo '------------------------------------------'
echo '[-m] OASIS_CPL_MODE    : '$OASIS_CPL_MODE
echo '[-s] NML_SETUP         : '$NML_SETUP
echo '[-i] INPUTDIR_OASIS3MCT: '$INPUTDIR_OASIS3MCT
echo '[-w] WORKDIR           : '$WORKDIR
echo '------------------------------------------'

############################################
### 1.STEP: SETTINGS (same as in xmessy_mmd)
############################################
#set NML_SETUP = OASIS/CLM04 
#set INPUTDIR_OASIS3MCT = /work/bb0677/b302010/DATA/
#set WORKDIR = /mnt/lustre01/scratch/b/b302010/tmp19
##set OASIS_CPL_MODE = 'INSTANT'
#set OASIS_CPL_MODE = 'AVERAGE'

#####################################
### 2.STEP: CREATE INPUTDIR_OASIS3MCT
#####################################

mkdir $INPUTDIR_OASIS3MCT/$NML_SETUP
mkdir $INPUTDIR_OASIS3MCT/$NML_SETUP/01
mkdir $INPUTDIR_OASIS3MCT/$NML_SETUP/02

#############################################################################
### 3.STEP: CREATE AND SAVE RESTART-INPUT-FILES IN CASE OF CPL_MODE = AVERAGE
#############################################################################
if ($OASIS_CPL_MODE == 'AVERAGE') then
cd $WORKDIR/01
cdo merge COS*.nc tmp.nc
ncwa -a time tmp.nc oasis_restart_cosmo.nc
rm -f tmp.nc
mv oasis_restart_cosmo.nc  $INPUTDIR_OASIS3MCT/$NML_SETUP/01
cd $WORKDIR/02
cdo merge clm02*.nc tmp.nc
ncwa -a time tmp.nc oasis_restart_clm.nc
rm -f tmp.nc
mv oasis_restart_clm.nc $INPUTDIR_OASIS3MCT/$NML_SETUP/02
endif
#####################################################
### 4.STEP: SAVE NML-SETUP SPECIFIC OASIS-INPUT-FILES
#####################################################
cd $WORKDIR
cp grids.nc areas.nc masks.nc $INPUTDIR_OASIS3MCT/$NML_SETUP
cp 01/rmp_* $INPUTDIR_OASIS3MCT/$NML_SETUP/01
cp 02/rmp_* $INPUTDIR_OASIS3MCT/$NML_SETUP/02

echo "INPUTFILES ARE SUCCESFULLY CREATED IN $INPUTDIR_OASIS3MCT/$NML_SETUP"

exit


