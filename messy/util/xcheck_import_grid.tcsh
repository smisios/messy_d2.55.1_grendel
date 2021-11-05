#! /bin/tcsh -f

#set echo verbose

# ---------------------------------------------------------------------
# THIS SCRIPT IS TO CHECK THE INTERNAL CONSISTENCY OF import.nml
# WITH its dependent NAMELIST FILES (import/*.nml)
# AND THE VARIABLE PRESENCE IN THE CORRESPONDING netCDF FILES
# ---------------------------------------------------------------------
# Author: Patrick Joeckel, DLR, May 2011
# ---------------------------------------------------------------------

echo '======================================================================='
echo 'THIS SCRIPT CHECKS THE CONSISTENCY OF YOUR IMPORT_RGT SETUP WITH ALL   '
echo 'DEPENDENT NAMELIST FILES (NCREGRID) AND netCDF FILES'
echo '======================================================================='

if (${#} < 2) then
   echo 'Usage: '`basename $0`' <namelist file> <INPUTDIR_MESSY>'
   echo '        possibe namelist files: import*'
   exit 1
endif
set nmlfile = $1
set INPUTDIR_MESSY = $2

if (! -d import) then
 echo `basename $0`" must be called from within namelist directory."
 exit 1
endif

if (! -r $nmlfile) then
   echo ' *** ERROR: '$nmlfile' is missing or unreadable.'
   exit 1
endif

set bpath = `dirname $0`

# count number of valid lines in namelist file
set nlines = `sed -n '/^\&RGTEVENTS/,/^\//p' $nmlfile | grep -v '^[ ]*!' | grep -i RG_TRIG | wc -l`

# loop over valid lines
@ cnt = 1
while ($cnt <= $nlines)
  echo '----------------------------------------------------------------------'

  # SELECT LINE ACCORDING TO LINE NUMBER
  set line = `sed -n '/^\&RGTEVENTS/,/^\//p' $nmlfile | grep -v '^[ ]*!' | grep -i RG_TRIG | awk '{if (NR=='$cnt') print}'`

  # PARSE LINE
  set no = `echo $line | awk -F '=' '{print $1}' | sed 's|[()]| |g' | awk '{print $2}'`
  set name = `echo $line | awk -F ',' '{print $5}' | sed 's|'\''||g' | sed 's|"||g'`
  set start = `echo $line | awk -F ',' '{print $6}'`
  set step  = `echo $line | awk -F ',' '{print $7}'`
  set stop  = `echo $line | awk -F ',' '{print $8}'`
  set first = `echo $line | awk -F ',' '{print $9}'`
  set act   = `echo $line | awk -F ',' '{for(i=10;i<=NF;i++) {printf("%s,",$i)}; printf("\n")}' | sed 's|"|'\''|g' | awk -F \' '{printf("%s\n", $2)}'`

  # OUTPUT
  echo "NO/CNT/NAME : "$cnt/"RG_TRIG("$no")"/$name
  echo "CYCLE       : "$start $step $stop $first 
  echo "ACTION      : "

  # PARSE ACTION STRING AND SET CORRESPONDING SHELL VARIABLES
  set nblks = `echo $act | tr ';' '\n' | wc -l`
  @ nb = 1
  while ($nb <= $nblks)
    set cmd = `echo $act | tr ';' '\n' | awk '{if (NR=='$nb') print}'` 
    if ("$cmd" != "") then
       echo "              "$cmd
       ### NOTE: THIS IS A TRICK TO SET NML, VAR, etc. IN THIS SHELL ...
       set $cmd
    endif
    @ nb ++
  end

  # (1) CHECK NAMELIST FILE
  echo "NML-FILE    : "$NML
  if (! -r $NML) then
     echo ' *** ERROR: namelist file '$NML' is missing or unreadable!'
     #exit 1
#  else
#     echo "             ... is available."
  endif

  if (${?VAR}) then
     echo "VARIABLE    : "$VAR
  else
     echo 'VARIABLE    : *** WARNING: none explicitly selected!'
  endif

  # (2) CHECK CORRESPONDING netCDF FILE(s)
  set ncfile_l = (`grep -v '^[ ]*!' $NML | grep -i infile | sed 's|,||g' | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's|'\"'||g' | sed 's|$INPUTDIR_MESSY|'$INPUTDIR_MESSY'|g' | sort | uniq`)

  @ n = 0
  foreach ncfile ($ncfile_l)
    @ n ++
    echo "netCDF-FILE : "$ncfile
    if (! -r $ncfile) then
       echo ' *** ERROR: netCDF file '$ncfile' is missing or unreadable!'
       #exit 2
    else
       # (3) OUTPUT TIME INFORMATION
       set timename = `grep -v '^[ ]*!' $NML | grep -ie 'i_timem *=' | awk '{if (NR=='$n') print}' | awk -F '=' '{print $2}' | sed 's|['\'',"]||g'` 
       # echo $timename
       set tlist = (`$bpath/time2yyyymmdd_hhmm.tcsh $ncfile $timename`)
       if ($status == 0) then
         echo "              START: "$tlist[$start]" ("$start")"
         echo "              STEP : "$step
         echo "              STOP : "$tlist[$stop]" ("$stop")"
         if ($first == '$START_MONTH') then
            echo "              FIRST: $first"
         else
            echo "              FIRST: "$tlist[$first]" ("$first")"
         endif
       else
         echo " *** WARNING: time axis could not be determined."
       endif
   endif
  end

  # (4) VARIABLE CONSISTENT ?
  if (${?VAR}) then
     @ found = 0
  endif

  # number of lines with 'var =' ...
  set nvarl = `grep -v '^[ ]*!' $NML | grep -ie 'var *=' | wc -l`
  @ n = 1
  while ($n <= $nvarl) 

    # this is nth line with var = ...
    set varline = `grep -v '^[ ]*!' $NML | grep -ie 'var *=' | awk '{if (NR=='$n') print}' | sed 's|"|'\''|g' | sed 's|"|'\''|g' | awk -F ''\''' '{print $2}' | sed 's| ||g'`
    #echo $varline

    # ... and this is the corresponding netCDF file
    set ncfile = `grep -v '^[ ]*!' $NML | grep -i infile | sed 's|,||g' | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's|'\"'||g' | sed 's|$INPUTDIR_MESSY|'$INPUTDIR_MESSY'|g' | awk '{if (NR=='$n') print}'`

      # parse line var = ...
      set nvars = `echo $varline | tr ';' '\n' | wc -l`
      @ m = 1
      while ($m <= $nvars)
        set rename = `echo $varline | tr ';' '\n' | awk '{if (NR=='$m') print}'`
         if ("$rename" == "") then
            @ m++
            continue
         endif
         echo "   var = ...: "$rename
         set nor = `echo $rename | awk -F '=' '{print NF}'`
         set new = `echo $rename | awk -F '=' '{print $1}' | sed 's|[,:].*$||g'`
         if ($nor == 1) then
            set old = $new
         else
            set old = `echo $rename | awk -F '=' '{print $2}' | sed 's|[,:].*$||g'`
         endif
     
         if (${?VAR}) then
            if ($new == $VAR) then
               @ found ++
            endif
         endif

         if (-r $ncfile) then
            #ncdump -h $ncfile | grep " $old("
            set in = `ncdump -h $ncfile | grep " $old(" | wc -l`
            if ($in == 0) then
               echo ' *** ERROR: Variable '$old' not found in '$ncfile'.'
            else
               echo "              ... "$old" is in "$ncfile
               if (${?VAR}) then
                  if ($VAR == $new) then
                     echo '              OBJECT NAME : '$name'_'$new
                  else
                     echo '              *** WARNING: variable delivered but not used!'
                  endif
               else
                  echo '              OBJECT NAME : '$name'_'$new
               endif
            endif
         else
               echo '             (no object created, netCDF file missing/unreadable)'
         endif
    
         @ m++
       end

    @ n++
  end

   if (${?VAR}) then
     if ($found == 0) then
        echo ' *** WARNING: Variable '$VAR' not explicitly found in '$NML'.'
        # check, if it is in netCDF file(s)
        @ found = 0
        set ncfile_l = (`grep -v '^[ ]*!' $NML | grep -i infile | sed 's|,||g' | awk -F '=' '{print $2}' | sed 's|'\''||g' | sed 's|'\"'||g' | sed 's|$INPUTDIR_MESSY|'$INPUTDIR_MESSY'|g' | sort | uniq`)
        foreach ncfile ($ncfile_l)
           set in = `ncdump -h $ncfile | grep " $VAR(" | wc -l`
           if ($in > 0) then
              echo "             ... "$VAR" is in "$ncfile
              @ found ++
           endif
        end
        if ($found == 0) then
           echo ' *** ERROR: '$VAR' not found in netCDF file(s).'
        else
           echo 'OBJECT NAME: '$name'_'$VAR
        endif
     endif
   
     if ($found == 1) then
        echo "VARIABLE    : "$VAR" has been found in "$NML" (OBJECT "$name"_"$VAR")"
     endif

     if ($found > 1) then
        echo ' *** ERROR: Variable '$VAR' found more than once in '$NML'.'
     endif

     unset VAR
   endif

  echo '----------------------------------------------------------------------'
  @ cnt ++
end

echo '======================================================================='
echo ' FINISHED.'
echo '======================================================================='

exit 0
