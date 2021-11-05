#! /bin/tcsh -f

set submodel=$1
set logfile = /dev/null

#  op_pj_20091022+
set echo_style=both    
#  op_pj_20091022-

set batchfile = $2
if ( "$batchfile" != "" ) then
  if ( -f $batchfile ) then
    echo "\nUsing batchfile: $batchfile"
    set batch
    source $batchfile
    if ( "$3" != "" ) then
       set logfile = $3
    endif
  else
    echo "ERROR: $batchfile does not exist"
    exit 1
  endif
endif

if ${?integr} then
   echo "\nFrom ${batchfile}: integr = $integr" | tee -a $logfile 
else
   set defaultintegr = "rosenbrock_mz"
   echo "\nNumerical integrator:"
   cd ./integr
   echo "\nAvailable integrators:"
   set inn      = "0"
   set allfiles = "*.kpp"
   set allfiles = `echo $allfiles | sed 's|\.kpp||g'` # delete suffix .kpp

   foreach i ($allfiles) # list all possibilities
     @ inn=$inn + 1
     echo "$inn) $allfiles[$inn]"
   end
   echo "The Rosenbrock integrators with automatic time-step control"
   echo "(e.g. ros3) are strongly recommended for their ability to cope"
   echo "with stiff sets of differential equations."
   echo "If you choose another integrator, do so at your own risk\!"
   echo "Note: Manual time stepping is switched via the CTRL_KPP namelist in"
   echo "      $submodel.nml."
   echo "Choose an integrator [q=quit, default=$defaultintegr]:"
   set inputstring = "$<"

   if ( "$inputstring" == "q" ) exit 1

   if (($inputstring <= $#allfiles) && ($inputstring >= 1)) then
     set integr = "$allfiles[$inputstring]"
     echo "You selected: $inputstring) $integr"
   else
     set integr = "$defaultintegr"
     echo "Default selection: $integr"
   endif
   cd -
endif

# mz_rs_20150220+

# instead of using the directory def_${submodel}/, kp4.sh will copy the
# integrator file directly into the work directory
# cd ./def_${submodel}
# rm -f integr.kpp
# ln -fs ../integr/$integr.kpp integr.kpp
# cd ..
# mz_rs_20150220-

# OPTIONS:
# scalar    :
# vector    : -v <length>
# deindexing: [-i <0,1,2,3>]

##################################################################
### vector (length) or scalar mode
##################################################################

switch ($integr)

 case "rosenbrock_vec":
    set integr = rosenbrock
    if ${?vlen} then
       echo "\nFrom ${batchfile}: vlen = $vlen"  | tee -a $logfile
    else
       echo "\n Which internal vector length ?"
       set vlen = "$<"
       set number = `echo $vlen | awk '{printf("%i\n",$0)}'`
       if ("$number" != "$vlen") then
         echo "ERROR: Vector length must be integer number: "$vlen
         exit 1
       endif
    endif
    set vopt = "-v $vlen"
    breaksw

 default:
    set vopt = ''
    breaksw

endsw

##################################################################
### indirect indexing
##################################################################

if ${?decomp} then
   echo "\nFrom ${batchfile}: decomp = $decomp" | tee -a $logfile
else
   echo "\nRemove indirect indexing ... ? [0/1/2/3/q, default=0]"
   set decomp = "$<"
endif

switch ($decomp)

 case "q":
    exit 1
    breaksw

 case "1":
 case "2":
 case "3":
    set diopt = "-i $decomp"
    breaksw

 default:
    set diopt = '-i 0'
    breaksw

endsw

echo "starting kp4.sh"
echo "./bin/kp4.sh -m $submodel -s $integr $vopt $diopt"
./bin/kp4.sh -m $submodel -s $integr $vopt $diopt
if ($status != 0) then
   echo "ERROR in kp4.sh / kp4.exe"
   exit 1
endif
echo "kp4.sh finished successfully"

# mz_rs_20150220+
# normally, the submodel path is just the name of the submodel:
set submodelpath = $submodel
# (currently, kp4.tcsh is only used by MECCA. The above definition may
# be used later when, e.g., SCAV is updated to KPP2 and uses KP4 as
# well)
# For MECCA, the core is in mbm/caaba/mecca:
if ( "$submodel" == "mecca" || "$submodel" == "mtchem") then
  set submodelpath = caaba/mecca
endif
# This applies also to polymecca (mecca### with a 3-digit number):
echo $submodel | grep  -q -e '^mecca[0-9][0-9][0-9]$'
if ( $status == 0 ) then
  set submodelpath = caaba/mecca
endif
if ( "$submodel" == "gmxe_aerchem") then
  set submodelpath = gmxe/aerchem
endif
# mz_rs_20150220-

# new code
echo "\nMoving the kp4-generated file messy_${submodel}_kpp.f90 into"
echo "the directory mbm/$submodelpath/smcl"
mv -f messy_${submodel}_kpp.f90 ../../mbm/$submodelpath/smcl/.
# for tracdef.awk
mv -f tmp_workdir/messy_${submodel}_kpp_Parameters.f90 ../../mbm/$submodelpath/tmp_param

# remove temporary files
set tmpfiles = "tmp_workdir/*"

if ${?deltmpkp4} then
   echo "\nFrom ${batchfile}: deltmpkp4 = $deltmpkp4" | tee -a $logfile
else
   echo "\nDo you want to delete the temporary files from kp4"
   ls -1 $tmpfiles
   echo "[y/n/q, default=y]?"
   set deltmpkp4 = "$<"
endif
if ( "$deltmpkp4" == "q" ) exit 1
if ( "$deltmpkp4" != "n" ) then
  if (${?TRASH}) then
    # move to trash directory
    mv -f $tmpfiles $TRASH
  else
    # or delete
    rm $tmpfiles
    rmdir ./tmp_workdir
  endif
endif

exit 0
