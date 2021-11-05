#! /bin/tcsh -f

### This shell-script converts the values along the time axis (name = $2)
### of a netCDF file ($1) into strings of format YYYYMMDD_hhmm. 
### Version 2.1b (without range limited time functions of gawk)
### (C) Patrick Joeckel, DLR, August 2014
 
if (( $# != 2 ) && ($# != 1)) then
    echo ' '
    echo 'Usage:'
    echo " `basename $0` <netCDF-file> [<name of time axis>]"
    echo " Note: name of time axis is autodetected, if not specified."
    echo ' '
    exit 1
endif

set file = $1

if ( $# == 2 ) then
   set tvar = $2
else
  set tvar = `ncdump -h $file | grep -i UNLIMITED | gawk '{print $1}'`
  if ($tvar == "") then
     echo "ERROR: name of time axis could not be detected, possibly due to missing UNLIMITED dimension attribute"
     exit 1
  endif
endif


set tunitstr = `ncdump -h $file | grep $tvar | grep unit | gawk -F '=' '{print $2}' | sed 's|[";]||g'`

# echo $tunitstr

set tunit = `echo $tunitstr | gawk '{print tolower($1)}'`

switch ("$tunit")

 case "second":
 case "seconds":
   set fsec = 1
   breaksw

 case "minute":
 case "minutes":
   set fsec = 60
   breaksw

 case "hour":
 case "hours":
   set fsec = 3600
   breaksw

 case "day":
 case "days":
   set fsec = 86400
   breaksw

 case "month":
 case "months":
  set fsec = 2629800   # 365.25 * 86400 / 12
  breaksw

  case "year":
  case "years":
   set fsec = 31557600 # 365.25 * 86400
   breaksw

 default:
   echo "ERROR: unknown time unit: "$tunitstr
   exit 1
   breaksw
endsw

# echo $fsec

set t0_dy = `echo $tunitstr | gawk '{print $3}' | gawk -F '-' '{printf("%2.2i\n",$3)}'`
set t0_mo = `echo $tunitstr | gawk '{print $3}' | gawk -F '-' '{printf("%2.2i\n",$2)}'`
set t0_yr = `echo $tunitstr | gawk '{print $3}' | gawk -F '-' '{printf("%4.4i\n",$1)}'`
set t0_hr = `echo $tunitstr | gawk '{print $4}' | gawk -F ':' '{printf("%2.2i\n",$1)}'`
set t0_mi = `echo $tunitstr | gawk '{print $4}' | gawk -F ':' '{printf("%2.2i\n",$2)}'`
set t0_se = `echo $tunitstr | gawk '{print $4}' | gawk -F ':' '{printf("%2.2i\n",$3)}'`

# echo $t0_yr $t0_mo $t0_dy $t0_hr $t0_mi $t0_se 

### time axis origin in seconds since 1970-01-01 00:00:00
set t0 = `date "+%s" -u -d "${t0_yr}-${t0_mo}-${t0_dy} ${t0_hr}:${t0_mi}:${t0_se}"`

#echo $t0

set tlist = (`ncdump -v${tvar} $file | gawk 'BEGIN {p=0}; {if ( $1=="data:" ) {p=1;getline}; if (p==1) print}' | sed 's|'$tvar' =||g' | sed 's|[;es}]||g' | tr ',' ' '`)

foreach tval ($tlist)

set tnowsec = `echo $t0 $tval $fsec | gawk '{printf("%i", ($1+$2*$3))}'`

#echo $tnowsec

set tnowstr = `date "+%Y%m%d_%H%M" -u -d "1970-01-01 $tnowsec sec"`

echo -n $tnowstr" "

end 

echo ''

exit 0
