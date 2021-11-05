#! /bin/tcsh -f

############################################################################
### SCRIPT TO DETECT CYCLIC DEPENDENCIES IN depend.mk FILES
### Example: 
###   messy/util/cyclic.tcsh echam5/__LINUX/depend.mk mo_geoloc > cyclic.out
###    (!!! this can take a long time !!!)
###   emacs cyclic.out 
###   search for "CYCLE" ...
###   ... and then backwards for decreasing recursion levels ...
###   ... (n-1) (n-2) ... (2) (1)
###
### Author: Patrick Joeckel, MPIC, July 2008
############################################################################

set script = `basename $0`
set spath  = `dirname $0`

if (($# != 2) && ($# != 4)) then
  echo "Usage : $script <dependency file> <module-name>" 
  exit 1
endif

set dpfile = $1
set fname  = $2

if ($3 == "") then
   @ level = 0
else
   @ level = $3
endif

if ($4 == "") then
   set ormod = $2
else
   set ormod = $4
endif

set list = (`awk -F ':' '{if ($1 == "'${fname}.o'") {print $2}}' $dpfile | sed 's|\.o||g' | tr ' ' '\n' | grep -v ".mod" | sort | uniq`)

if ($level == 0) then
   echo ${level}:${fname}: "["$list"]"
endif

if (${#list} > 0) then
  @ level++
  foreach file ($list)

    @ cnt=1
    while ($cnt < $level)
      echo -n "|"
      @ cnt++
    end

    echo "|-> "$file "("$level")"

    if ($file == $ormod) then
       echo " ### CYCLE COMPLETED: $level $file $ormod"
       exit 1
    endif

     ${spath}/${script} $dpfile $file $level $ormod
  end
else
    @ cnt=0
    while ($cnt < $level)
      echo -n "|"
      @ cnt++
    end
    echo "|->#"
endif

exit 0
