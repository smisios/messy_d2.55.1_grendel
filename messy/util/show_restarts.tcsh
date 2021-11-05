#! /bin/tcsh -f

### SHOW INFORMATION ABOUT RESTART DATES
### - TO BE CALLED IN WORKDIR (WITH save/???? - SUBDIRECTORIES)
###
### Author: Patrick Joeckel, DLR, July 2011

#defaults:
set script = `basename $0`
set cha = g3b
set xf  = 0

while ($# > 0)

  switch ($1)

     case '-h':
       echo ""
       echo "Usage: $script [-h] [-c <channel>] [-f]"
       echo ""
       echo "show restart data/time of restart chain elements and cycles"
       echo "(needs to be called from MESSy-workdir)"
       echo ""
       echo "Options:"
       echo "  -h : show this help and exit"
       echo "  -c : channel for inspection (default: $cha)"
       echo "  -f : output in special format for forecast chain"
       echo ""
       shift
       exit 0
       breaksw

     case '-c':
       shift
        if (("$1" == "") || \
            (`echo $1 | awk '{print substr($1,1,1)}'` == "-")) then
            echo 'Error: Argument <channel> missing\!'
            echo "Use $script -h for more information\!"
            exit 1
        else
            set cha = $1
            shift
        endif
        breaksw

     case '-f':
        set xf = 1
        shift
        breaksw

      default:
        echo "Unknown option: $1"
        echo "Use $script -h for more information\!"
        exit 1
        breaksw

   endsw

end

if ($xf == 0) then
   echo
   echo "================================================================="
   echo "CHAIN ELEMENT - CYCLE : RESTART DATE TIME"
   echo "================================================================="
endif

foreach chain (save/????)

   set chain_nr = `basename $chain`

   set cycle_list = (`find $chain -name "restart_????_${cha}.nc" -print | sort`)
    
   @ ncyc = ${#cycle_list}

   @ cnt=1
   while ($cnt <= $ncyc)
      set cyc_nr = `echo ${cycle_list[$cnt]} | awk -F '_' '{print $2}'`
      if ($xf == 0) then
         set dateinfo = `ncdump -h ${cycle_list[$cnt]} | grep restart_date_time | awk -F '=' '{print $2}' | sed 's|;||g'`
      else
         set dateinfo = `ncdump -h ${cycle_list[$cnt]} | grep restart_date_time | awk -F '=' '{print $2}' | sed 's|;||g' | sed 's|\"||g'`
      endif

      if ($cnt == 1) then
         if ($xf == 0) then
           echo -n $chain_nr"          - "$cyc_nr"  : "$dateinfo
         else
           echo -n $chain_nr"            "$cyc_nr"    "$dateinfo
         endif
      else
         if ($xf == 0) then
           echo -n "              - "$cyc_nr"  : "$dateinfo
         else
           echo -n $chain_nr"            "$cyc_nr"    "$dateinfo
         endif
      endif
      if ($cnt == $ncyc) then
         echo " (FILE_DUMP_AND_RESTART)"
         if ($xf == 0) then
            echo "-----------------------------------------------------------------"
         endif
      else
         echo " (FILE_DUMP)"
      endif
   @ cnt++
   end

end

if ($xf == 0) then
   echo "================================================================="
endif

exit 0
