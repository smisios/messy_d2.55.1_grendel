#!/bin/tcsh -f

# Script to update RG_TRIG entries (e.g., in import.nml) according to
# desired start date. Only works for nml/nc files strictly compliant with
# MESSy input file naming convention:
# <provider>_<version>_<scenario>_<class>_<species>_<yyyy[mm][dd]-yyyy[mm][dd]>
#  [.nml,.nc]
#
# Version: 1.0b
# Author: Patrick Joeckel, DLR, March 2015

#set echo verbose

set nonomatch

if ($# < 2) then
   echo "Usage: `basename $0` <file.nml> <YYYYMMDD>"
   echo "       file.nml: namelist file with RG_TRIG entries"
   echo "       YYYYMMDD: desired start date"
   echo " Note: this script requires fname2trig.tcsh in the same path."
   exit 0
endif

set fname = $1
set yyyymmdd = $2

# 
set mypath = `dirname $0`
set script2 = $mypath/fname2trig.tcsh

# number of lines
set nl = `wc -l $fname | awk '{print $1}'`
#echo $nl

# loop overlines
@ n = 1
while ($n <= $nl)


 set line = "`awk '{if (NR == '$n') print}' $fname`"
# echo "$line"
 set first =  `echo $line | awk '{print toupper($1)}'`
 set ff = `echo $first | sed 's| ||g' | awk -F '(' '{print $1}'`
 if (("$ff" == "\!RG_TRIG") || ("$ff" == "RG_TRIG")) then

    if ("$ff" == "\!RG_TRIG") then
       set pre = '!'
    else
       set pre = ''
    endif

    # last non-blank character must be ','
    set last = `echo "$line" | sed 's| ||g' | awk '{print substr($0,length($0),1)}'`
    if ($last != ',') then
       set line = "${line},"
    endif

    set no = `echo "$line" | awk -F ')' '{print $1}' | awk -F '(' '{print $2}'`
    set name = `echo "$line" | awk -F ',' '{print $5}' | sed 's|'\''||g'`
    set act0 = `echo "$line" | awk -F ',' '{for (i=10;i<NF-1;i++) {printf("%s,",$i)}; NFM1=NF-1;printf("%s\n",$NFM1)}'`
    set act  = `echo "$act0" | sed 's|'\''||g'`
    set nml0  = `echo "$act" | tr ';' '\n' | grep NML`
    set nml  = `echo $nml0 | sed 's|NML=||g'`
#echo "\!OLD: ""$line"
#echo "$act0"
#echo $no $name $act $nml

    echo -n $pre
#   $script2 -v -n $no -m $name -a "$act0" -f $nml -s $yyyymmdd
#    echo $script2 -n $no -m $name -a "$act0" -f $nml -s $yyyymmdd
    $script2 -n $no -m $name -a "$act0" -f $nml -s $yyyymmdd

 else
    echo "$line"
 endif


@ n++
end

exit 0
