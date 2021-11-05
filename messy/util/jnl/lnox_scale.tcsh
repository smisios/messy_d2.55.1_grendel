#!/bin/tcsh -f

#set echo verbose

set jnlpath = `dirname $0`

### CREATE mc-files
echo '#######################################################################'
echo '### CREATING mc-FILES'
echo '#######################################################################'
foreach inst (`seq -w 1 1 99`)
 if (-d $inst) then
    cd $inst
       if ("$inst" == "01") then
          nc2mc -a -p -x restart -n ECHAM5
       else
          nc2mc -a -p -x restart -n COSMO
       endif
       foreach par (Grewe PaR_T AaP_M AaP_P FinIF)
          nc2mc -a -p -x restart -n lnox_${par}_gp
       end
    cd -
 endif
end
echo '#######################################################################'

### CREATE FLAG FILES
echo '#######################################################################'
echo '### CREATING flag-FILES'
echo '#######################################################################'
foreach inst (`seq -w 2 1 99`)
 if (-d $inst) then
    echo -n GENERATING FLAG-FILE FOR INSTANCE ${inst} ... 
    if (! -e flag_${inst}.nc) then
       echo '(this may take a while)'
       ferret -memsize 1024 -nojnl -script ${jnlpath}/mecon_map2global.jnl ${inst} >&! /dev/null
    else
       echo already available
    endif
 endif
end
echo '#######################################################################'

echo '#######################################################################'
echo "PARAM\tSCALING\t\tINST\tSCALING\tREGIONAL-FF"

foreach inst (`seq -w 2 1 99`)

 if (-d $inst) then

    foreach par (Grewe PaR_T AaP_M AaP_P FinIF)

      set psfile  = lnox_scale_${par}_${inst}.ps
      set logfile = lnox_scale_${par}_${inst}.log

      ferret -memsize 1024 -nojnl -batch lnox_scale_${par}_01.ps \
             -script ${jnlpath}/lnox_scale.jnl \
              01/lnox_${par}_gp.mc \" \" flag_${inst}.nc > lnox_scale_${par}_01.log

      set r_scal_ff = `grep r_scal_ff lnox_scale_${par}_01.log | awk '{print $5}'`
      set ff_ave    = `grep ff_ave lnox_scale_${par}_01.log | awk '{print $4}'`


      ferret -memsize 1024 -nojnl -batch ${psfile} \
             -script ${jnlpath}/lnox_scale.jnl \
              ${inst}/lnox_${par}_gp.mc $ff_ave > ${logfile}

      set r_scal_ff_r = `grep r_scal_ff ${logfile} | awk '{print $5}'`
      set ff_ave_r    = `grep achieve ${logfile} | awk '{print $11}'`

      echo $par"\t"$r_scal_ff"\t"$inst"\t"$r_scal_ff_r"\t"$ff_ave" (G)\t"$ff_ave_r" (R)"

    end

 endif

end

#gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -sOutputFile=combined.pdf *.ps

exit 0
