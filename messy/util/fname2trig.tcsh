#!/bin/tcsh -f

# Script to create RG_TRIG entry (e.g., for import.nml)
# based on MESSy input file naming convention
# <provider>_<version>_<scenario>_<class>_<species>_<yyyy[mm][dd]-yyyy[mm][dd]>
#  [.nml,.nc]
#
# Version: 1.0b
# Author: Patrick Joeckel, DLR, March 2015

#set echo verbose

#defaults:
set script=`basename $0`
set filename = ''
set yyyymmdd = ''
set verb  = 0
set n = '@#@'
set m = '@NAME@'
set a = '@ACTION@'

###########################################################################
###########################################################################

while ($# > 0)

   switch ($1)

      case '-h':
        echo ""
        echo "Usage: $script [-h] -f <filename> -s <YYYY[MM[DD]]> [-v] [-n <no>] [-m <name>] [-a <action>]"
        echo ""
        echo "ANALYSE TIME AXIS FROM NETCDF FILE IN MESSy NAMING CONVENTION:"
        echo "<provider>_<version>_<scenario>_<class>_<species>_<yyyy[mm][dd]-yyyy[mm][dd]>[.nml,.nc]"
        echo ""
        echo "Options:"
        echo " -h : show this help and exit"
        echo " -f : file name"
	echo " -s : start date (YYYY[MM[DD]])"
        echo " -v : verbose output (default: no)"
        echo " -n : trigger number"
        echo " -m : trigger name"
        echo " -a : action string"
        echo ""
        shift
        exit 0
        breaksw

      case '-f':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <filename> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            set filename = $1
            shift
        endif
        breaksw

      case '-s':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <YYYY[MM[DD]]> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            set yyyymmdd = $1
            shift
        endif
        breaksw

      case '-n':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <no> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            @ n = $1
            shift
        endif
        breaksw

      case '-m':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <name> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            set m = $1
            shift
        endif
        breaksw

      case '-a':
        shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <action> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            set a = "$1"
            shift
        endif
        breaksw

     case '-v':
        shift
        set verb = 1
        breaksw

      default:
        echo "Unknown option: $1"
        echo "Use $script -h for more information\!"
        exit 1
        breaksw

   endsw

end

if ($filename == '') then
    echo "ERROR: -f option missing\!"
    echo "Use $script -h for more information\!"
    exit 1
endif

if ($yyyymmdd == '') then
    echo "ERROR: -s option missing\!"
    echo "Use $script -h for more information\!"
    exit 1
endif

###########################################################################
###########################################################################

set fname = `basename $filename .nc`
set fname = `basename $fname .nml`

### split at '_' according to convention
# <provider>_<version>_<scenario>_<class>_<species>_<yyyy[mm][dd]-yyyy[mm][dd]>.nc

set prov = `echo $fname | awk -F '_' '{print $1}'`
set vers = `echo $fname | awk -F '_' '{print $2}'`
set scen = `echo $fname | awk -F '_' '{print $3}'`
set clss = `echo $fname | awk -F '_' '{print $4}'`
set spec = `echo $fname | awk -F '_' '{print $5}'`
set tran = `echo $fname | awk -F '_' '{print $6}'`

# for cases where the '-' is an '_' ...
set tran2 = `echo $fname | awk -F '_' '{print $7}'`
if ($tran2 != '') then
   set tran = ${tran}-${tran2}
endif

# special case
if ($tran == const) then
   set tran = 'X-X'
endif

if ($verb == 1) then
   echo '! PROVIDER    : '$prov
   echo '! VERSION     : '$vers
   echo '! SCENARIO    : '$scen
   echo '! CLASS       : '$clss
   echo '! SPECIES     : '$spec
   echo '! TIME RANGE  : '$tran
endif

set t1 = `echo $tran | awk -F '-' '{print $1}'`
set t2 = `echo $tran | awk -F '-' '{print $2}'`

# length of string
set t1l = `echo $t1 | awk '{print length($0)}'`

# echo $t1 $t1l

switch ("$t1l")

 case "1":
   ### X
   set trigunit='years'
   set cycle = '1,1,1,1'
   breaksw

 case "2":
   ### month (MM)
   set trigunit='months'
   set cycle = '1,1,12,$START_MONTH'
   breaksw

 case "4":
   ### year (YYYY)
   set trigunit='years'
   # number of years
   @ ny = ($t2 - $t1) + 1
   # start year
   set sy = `echo $yyyymmdd | awk '{print substr($1,0,4)}'`
   # index of start year
   @ sy = ($sy - $t1) + 1
   if ($sy > $ny) then
      set sy = $ny
   endif
   if ($sy < 1) then
      if ($ny == 1) then
         # cycling the single available year ...
         set sy = $ny
      else
         echo "ERROR: start year out of range ("$sy")\!"
         exit 1
      endif
   endif
   # repeat last year cyclically
   set cycle = "$ny,1,$ny,$sy"
   breaksw

 case "6":
   ### year/month (YYYYMM)
   set trigunit='months'
   # number of months
   set y1 = `echo $t1 | awk '{print substr($1,0,4)}'`
   set m1 = `echo $t1 | awk '{print substr($1,5,2)}'`
   set y2 = `echo $t2 | awk '{print substr($1,0,4)}'`
   set m2 = `echo $t2 | awk '{print substr($1,5,2)}'`
   @ nm = `echo $y1 $m1 $y2 $m2 | awk '{print ($3-$1)*12 + $4 - $2 + 1}'`
   @ nm12 = ($nm - 12) + 1
   # index of start year / month
   @ sy = `echo $yyyymmdd | awk '{print substr($1,0,4)}'`
   @ sm = `echo $yyyymmdd | awk '{print substr($1,5,2)}'`
   if ($sm == 0) then
      echo 'ERROR: start month must be specified (YYYYMM) with -s option\!'
      exit 1
   endif
   if (($sm < 1) || ($sm > 12)) then
      echo 'ERROR: start month out of range [01..12] in -s option ('$sm')\!'
      exit 1
   endif
   @ s = `echo $y1 $m1 $sy $sm | awk '{print ($3-$1)*12 + $4 - $2 + 1}'`
   if ($s > $nm) then
      @ s = ($nm12 + $sm - 1)
   endif
   if ($s < 1) then
      if ($nm == 12) then
         # cycling the single available year
         @ s = $sm
      else
         echo "ERROR: start year/month out of range ("$s")\!"
         exit 1
      endif
   endif
   # repeat last year cyclically
   set cycle = "$nm12,1,$nm,$s"
   breaksw

 case "8":
   ### year/month/day (YYYYMMDD)
   set trigunit='days'
   # years, months, days
   set y1 = `echo $t1 | awk '{print substr($1,0,4)}'`
   set m1 = `echo $t1 | awk '{print substr($1,5,2)}'`
   set d1 = `echo $t1 | awk '{print substr($1,7,2)}'`
   set y2 = `echo $t2 | awk '{print substr($1,0,4)}'`
   @ y2m1 = $y2 - 1
   set m2 = `echo $t2 | awk '{print substr($1,5,2)}'`
   set d2 = `echo $t2 | awk '{print substr($1,7,2)}'`
   #
   @ sy = `echo $yyyymmdd | awk '{print substr($1,0,4)}'`
   @ sm = `echo $yyyymmdd | awk '{print substr($1,5,2)}'`
   @ sd = `echo $yyyymmdd | awk '{print substr($1,7,2)}'`
   if (($sd == 0) || ($sm == 0)) then
      echo 'ERROR: start month and day must be specified (YYYYMMDD) with -s option\!'
      exit 1
   endif
   if (($sm < 1) || ($sm > 12)) then
      echo 'ERROR: start month out of range [01..12] in -s option ('$sm')\!'
      exit 1
   endif
   ### time axis origin in seconds since 1970-01-01 00:00:00
   set tt1   = `date "+%s" -u -d "${y1}-${m1}-${d1} 00:00:00"` 
   set tt2   = `date "+%s" -u -d "${y2}-${m2}-${d2} 00:00:00"` 
   set tt2m1 = `date "+%s" -u -d "${y2m1}-${m2}-${d2} 00:00:00"` 
   set tts   = `date "+%s" -u -d "${sy}-${sm}-${sd} 00:00:00"` 
   # number of days
   @ nd   = `echo $tt1 $tt2   | awk '{print ($2-$1)/86400}'` + 1
   @ ndm1 = `echo $tt1 $tt2m1 | awk '{print ($2-$1)/86400}'` + 2
   # actual day
   @ s = `echo $tt1 $tts | awk '{print ($2-$1)/86400}'` + 1
   # repeat last year cyclically
   if ($s > $nd) then
      @ s = `echo $tt1 $tt2m1 | awk '{print ($2-$1)/86400}'` + 1 + $sd
   endif
   if ($s < 1) then
      echo "ERROR: start year/month/day out of range ("$s")\!"
      exit 1
   endif
   # recycle last year
   set cycle = "$ndm1,1,$nd,$s"
  breaksw

 default:
   echo "ERROR: time axis could not be determined"
   exit 1
   breaksw

endsw

if ($verb == 1) then
   echo '! TRIGGER UNIT: '$trigunit
endif

echo 'RG_TRIG('$n') = '1,\'$trigunit\',\''first'\',0, \'$m\', $cycle, $a,

exit 0
