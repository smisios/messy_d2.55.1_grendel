#! /bin/tcsh -f

# WITH THIS SCRIPT YOU CAN CHECK IF THE REDIRECTION OF TRACER-OUTPUT
# IN THE channel.nml NAMELIST FILE (which must be available in the local
# directory) IS COMPLETE.
# THE LIST OF TRACERS IS EXTRACTED FROM THE LOG-FILE, YOU NEED TO SPECIFY
# AT THE COMMAND LINE.
# NOTE: THE SCRIPT IS NOT ABLE TO DEAL WITH WILDCARDS IN channel.nml
#
# Author: Patrick Joeckel, DLR, July 2010

if (${#} == 0) then
   echo "USAGE: `basename $0` <log-file>"
   exit 1
endif

set logfile = $1
set chnml   = channel.nml

# check lines to start with
set lstart = (`grep -n "NUMBER OF TRACERS" $logfile | awk -F ':' '{print $1}'`)

# check number of tracers
set notrac = (`grep -n "NUMBER OF TRACERS" $logfile | awk '{print $NF}'`)

# process only half of the entries
set nsets = `echo ${#lstart} | awk '{print $1/2}'`
#set nsets = `echo ${#lstart} | awk '{print $1}'`
#echo $nsets

# loop over tracer sets
@ scnt = 1
while ($scnt <= $nsets)
 if ($notrac[$scnt] != 0) then

    echo "=================================================================="
    echo "TRACER SET $scnt WITH $notrac[$scnt] TRACERS (line $lstart[$scnt])"
    echo "------------------------------------------------------------------"

    @ l1 = $lstart[$scnt] + 20   # add header lines
    @ l2 = $l1 + $notrac[$scnt] - 1

    set trlist = (`awk -v l1=$l1 -v l2=$l2 '{if ((NR>=l1) && (NR<=l2)) {print $3}}' $logfile`)

    # check number of tracers
    set mtrac = ${#trlist}
    if ($mtrac != $notrac[$scnt]) then
       echo "ERROR: NUMBER OF TRACER MISMATCH ($mtrac \!= $notrac[$scnt])"
       exit 1
    endif

    foreach tr ($trlist)

      echo $tr | awk '{printf("%-20s",$1)}'

      set ref = `grep \'$tr\' $chnml`

      if ("$ref" != "") then
        set no = `echo $ref | awk -F '(' '{print $2}' | awk -F ')' '{print $1}'`
        set ch = `echo $ref | awk -F ',' '{print $3}' | sed "s|[']||g"`
        echo $no | awk '{printf("%5i  ",$1)}'
        echo $ch | awk '{printf("%-10s\n",$1)}'
      else
        echo XXXXX | awk '{printf("%5s  \n",$1)}'
      endif

    end

    echo "=================================================================="

 endif
 @ scnt ++
end

exit 0
