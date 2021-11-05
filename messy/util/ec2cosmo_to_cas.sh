#!/bin/sh -e


# Convert output of the MESSy submodel EC2COSMO (netcdf format) to global 
# casYYYYMMDDhhmmss.nc NetCDF files, which can then be converted with int2clm to# lafYYYYMMDDhh.nc files used to drive the CLM model
# INPUT: MESSy CHANNEL output (netcdf format) 
#              (e.g. files $EXPNAME_YYYYMMDD_HHMM_ec2cosmo.nc)
# OUTPUT: Global NetCDF files casYYYYMMDDhhmmss.nc
#
# usage: ec2cosmo_to_cas.sh <grid-file> <experiment> <start/end year> <month> (<day>)
# usage: ec2cosmo_to_cas.sh <start year> <month> (<day>) <end year> <month> (<day>)
#
# Contact: pantatk@uni-mainz.de
# this script is based on the ECHAM_to_cas.sh script of Andy Dobler and
# Hans-Juergen Panitz rewritten for EMAC by Astrid Kerkweg and optimised by 
# Klaus Pankatz
# Version: 2.0
# Date: 2013-04-18
################################################################################

# FUNCTIONS
################################################################################


function workonfile {

#list of passed variables
SOURCEDIR=$1
EXPNAME=$2
YEAR=$3
MM=$4
DD=$5
HH=$6
MI=$7
WORKDIR=$8

FILE=${SOURCEDIR}/${EXPNAME}${YEAR}${MM}${DD}_${HH}${MI}_ec2cosmo.nc
FILENAME=${EXPNAME}${YEAR}${MM}${DD}_${HH}${MI}_ec2cosmo.nc
#DATESTR=${YEAR}${MM}${DD}${HH}${MI}00
DATESTR=${YEAR}${MM}${DD}${HH}

echo "Check for file $FILE"   
if [ ! -f $FILE ] ; then

    echo "File $FILE not found"
    return 1
fi 

echo "Check for directory $WORKDIR/tmp"
if [ ! -d $WORKDIR/tmp ] ; then

    mkdir $WORKDIR/tmp
fi

echo ... copying file $FILE
cp -f $FILE $WORKDIR/tmp


echo "Creating cas${DATESTR}.nc"

cdo -merge \
    -selname,geosp,st,u,v,aps,sn,wl,FR_LAND -dv2uv ./tmp/$FILENAME \
    -ifthenelse -gec,0.0 -selname,q ./tmp/$FILENAME -selname,q  ./tmp/$FILENAME -gec,0.0 -selname,q ./tmp/$FILENAME \
    -ifthenelse -gec,0.0 -selname,xl ./tmp/$FILENAME -selname,xl ./tmp/$FILENAME -gec,0.0 -selname,xl ./tmp/$FILENAME \
    -ifthenelse -gec,0.0 -selname,xi ./tmp/$FILENAME -selname,xi ./tmp/$FILENAME -gec,0.0 -selname,xi ./tmp/$FILENAME \
    -chname,var1,T_S -ifthenelse -expr,"var1=FR_LAND_MASK;" -selname,FR_LAND_MASK  ./tmp/$FILENAME -expr,"var1=tslm1;" -selname,tslm1  ./tmp/$FILENAME -expr,"var1=seaice*tsi+tsw-seaice*tsw;" -selname,tsi,tsw,seaice ./tmp/$FILENAME \
    -ifthen -expr,"var1=FR_LAND_MASK;" -selname,FR_LAND_MASK  ./tmp/$FILENAME -selname,tsoil  ./tmp/$FILENAME \
    -chname,var2,W_SO_REL -ifthen -expr,"var1=FR_LAND_MASK;" -selname,FR_LAND_MASK ./tmp/$FILENAME -ifthenelse -gtc,1 -setcode,2 -expr,"var2=tsoil-tsoil+ws/wsmx;" -selname,ws,wsmx,tsoil ./tmp/$FILENAME -gtc,1 -setcode,2 -expr,"var2=tsoil-tsoil+ws/wsmx;" -selname,ws,wsmx,tsoil ./tmp/$FILENAME -setcode,2 -expr,"var2=tsoil-tsoil+ws/wsmx;" -selname,ws,wsmx,tsoil ./tmp/$FILENAME \
    -chname,var1,FR_LAND_MASK -settaxis,${YEAR}-${MM}-${DD},${HH}:${MI} -expr,"var1=FR_LAND_MASK;" -selname,FR_LAND_MASK ./tmp/$FILENAME \
    ./tmp/tmp_merge_${DATESTR}.nc

cdo setmissval,-1e20 -remapbil,grid.txt -invertlat -sp2gp ./tmp/tmp_merge_${DATESTR}.nc cas${DATESTR}.nc      


    # rename dimensions
    #ncrename -d ilev,level1 cas${DATESTR}.nc
ncrename \
    -d lev,level \
    -d nhyi,level1 \
    -d belowsf,soil1 \
    -v belowsf,soil1 \
    -v hyai,ak \
    -v hybi,bk \
    -v geosp,FIS \
    -v wl,W_I \
    -v sn,W_SNOW \
    -v u,U \
    -v v,V \
    -v st,T \
    -v aps,PS \
    -v q,QV \
    -v xl,QC \
    -v xi,QI \
    -v tsoil,T_SO cas${DATESTR}.nc

ncatted -O -a standard_name,T,c,c,"air_temperature" \
    -O -a missing_value,T,c,f,-1e20 \
    -O -a standard_name,U,c,c,"grid_eastward_wind" \
    -O -a missing_value,U,c,f,-1e20 \
    -O -a standard_name,V,c,c,"grid_northward_wind" \
    -O -a missing_value,V,c,f,-1e20 \
    -O -a standard_name,PS,c,c,"surface_air_pressure" \
    -O -a missing_value,PS,c,f,-1e20 \
    -O -a standard_name,QV,c,c,"specific_humidity" \
    -O -a missing_value,QV,c,f,-1e20 \
    -O -a standard_name,QC,c,c,"mass_fraction_of_cloud_liquid_water_in_air" \
    -O -a missing_value,QC,c,f,-1e20 \
    -O -a standard_name,QI,c,c,"mass_fraction_of_cloud_ice_in_air" \
    -O -a missing_value,QI,c,f,-1e20 \
    -O -a standard_name,FIS,c,c,"surface_geopotential_height" \
    -O -a missing_value,FIS,c,f,-1e20 \
    -O -a standard_name,T_S,c,c,"surface temperature" \
    -O -a missing_value,T_S,c,f,-1e20 \
    -O -a standard_name,W_SNOW,c,c,"surface_snow_amount" \
    -O -a missing_value,W_SNOW,c,f,-1e20 \
    -O -a standard_name,FR_LAND,c,c,"land_area_fraction" \
    -O -a missing_value,FR_LAND,c,f,-1e20 \
    -O -a standard_name,W_SO_REL,c,c,"relative (scaled to max. field capacity) volume_fraction_of_soil_moisture" \
    -O -a missing_value,W_SO_REL,c,f,-1e20 \
    -O -a standard_name,W_I,c,c,"canopy_water_amount" \
    -O -a missing_value,W_I,c,f,-1e20 cas${DATESTR}.nc \
    -O -a standard_name,T_SO,c,c,"soil_temperature" \
    -O -a missing_value,T_SO,c,f,-1e20 \

# CLEAN UP /tmp

echo "cleaning up"
rm -f $WORKDIR/tmp/tmp_merge_${DATESTR}.nc $WORKDIR/tmp/$FILENAME

echo
echo File cas${DATESTR}.nc finished
return 0
 
}

###########################################################################
#
#  User definition section BEGIN
#
###########################################################################

################################################################################
#
#    Variables in user defintion section
#
# indir:     Directory where to find this script ECHAM_to_cas.sh
# sourcedir: Directory where to find the source data which will be converted to
#            CLM4 compatible format
# targetdir: Directory where the converted data will be stored
#
# head_archive_file: head of name of tar-file; the filename will be extended 
#            by "casyyyymmddhhmiss"
#
################################################################################
################################################################################
#  nco and cdo are needed
################################################################################
################################################################################
# setenv PATH=$PATH:/pool/ia64/cdo/bin/:/pool/ia64/netcdf/nco/bin
################################################################################
###


# Wahl zwischen 14 stelligem Datum und 10 stelligem Datumstring

SOURCEDIR=/miklip/scratch/b302050/DR_g63cn_002
WORKDIR=$SOURCEDIR/int2lm

EXPNAME=DR_g63cn_002___
#EXPNAME=${2}"___"
GRIDFILE=T63_Gaussian.txt

### Define START and END date here in format mm/dd/yy hh:mm:ss

startdate=$(date -u --date="$2/$3/$1 00:00:00" +%s)
echo "Startdate is $startdate or in other words $(date -u -d "1970-01-01 $startdate sec " "+%Y%m%d %H:%M:%S")"

#Sometimes you may want to run the script for just one day
#let enddate=$startdate+86400
enddate=$(date -u --date="$5/$6/$4 00:00:00" +%s) # Zum Beispiel 20060101000000
echo "Enddate is $enddate or in other words $(date -u -d "1970-01-01 $enddate sec " "+%Y%m%d")"

# timestep in seconds
timestep=3600
maxprocs=1


###########################################################################
#
#  User definition section END
#
###########################################################################
#
# CHANGE NOTHING BELOW THIS LINE
###########################################################################

if [ ! -d $WORKDIR ] ; then
    mkdir  ${WORKDIR}
fi


INDIR=`dirname $0`
INDIR=`cd $INDIR; pwd`


if [ ! -f ${INDIR}/grids/${GRIDFILE} ] ; then
    echo "ERROR: grid-file ${GRIDFILE} not available in ${INDIR}/grids"
    exit 1
fi

### Ideally, we should test if the correct grid.txt is given. The cdo command
# won't work otherwise. 

if [ ! -f  ${WORKDIR}/grid.txt ] ; then
cp -f ${INDIR}/grids/${GRIDFILE} ${WORKDIR}/grid.txt
fi

now=$startdate
procs=0
echo "Begin loop at $(date -u -d "1970-01-01 $now sec " "+%Y%m%d %H:%M:%S")"
while [ $now -lt  $enddate ] ; do

echo "Proceeding with next chain element $(date -u -d "1970-01-01 $now sec " "+%Y%m%d %H:%M:%S")"

YEAR=$(date -u -d "1970-01-01 $now sec " "+%Y")
MM=$(date -u -d "1970-01-01 $now sec " "+%m")
DD=$(date -u -d "1970-01-01 $now sec " "+%d")
HH=$(date -u -d "1970-01-01 $now sec " "+%H")
MI=$(date -u -d "1970-01-01 $now sec " "+%M")
SS=$(date -u -d "1970-01-01 $now sec " "+%S")

   ###################################
   # work on file
   ###################################

echo "Work on file ${YEAR}${MM}${DD}${HH}${MI}${SS}"

cd $WORKDIR

workonfile "$SOURCEDIR" "$EXPNAME" "$YEAR" "$MM" "$DD" "$HH" "$MI" "$WORKDIR" &

let procs=procs+1
if [ $procs -gt $maxprocs ] ; then
    wait
    procs=0
fi

let now=$now+$timestep

done  

wait
