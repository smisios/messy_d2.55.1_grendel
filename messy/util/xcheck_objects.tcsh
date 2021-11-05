#! /bin/tcsh -f

#set echo verbose

# ---------------------------------------------------------------------
# THIS SCRIPT IS TO CHECK THE INTERNAL CONSISTENCY OF 
#  - tnudge.nml
#  - offemis.nml
#  - onemis.nml
# WITH import.nml
# ---------------------------------------------------------------------
# Author: Patrick Joeckel, DLR, May 2012
# ---------------------------------------------------------------------

echo '======================================================================='
echo 'THIS SCRIPT CHECKS THE CONSISTENCY OF THE CHANNEL OBJECTS              '
echo 'REQUESTED FROM IMPORT_GRID BY OTHER NAMELIST FILES'
echo '======================================================================='

if (${#} < 3) then
   echo 'Usage: '`basename $0`' <import namelist file> <namelist file> <INPUTDIR_MESSY>'
   echo '        possibe import namelist files  : import*'
   echo '        possibe submodel namelist files: tnudge*, offemis*, onemis*, ch4*, rad*'
   exit 1
endif
set inmlfile = $1
set nmlfile = $2
set INPUTDIR_MESSY = $3

if (! -d import) then
 echo `basename $0`" must be called from within namelist directory."
 exit 1
endif

if (! -r $nmlfile) then
   echo ' *** ERROR: '$nmlfile' is missing or unreadable.'
   exit 1
endif

set bpath = `dirname $0`

## NC and NO: note that the separator in the awk command below is the
##            quotation mark "'" instead of the commy ",". This the 
##            numbers seem odd .. are, however correct.

switch ("`basename $nmlfile`")
  case "tnudge*":
       set F1="TNUDGE_"
#       set NC=3
#       set NO=4
       set NC=6
       set NO=8
       breaksw

  case "offemis*":
       set F1="EMIS_IN"
#       set NC=2
#       set NO=3
       set NO=6
       set NC=4
       breaksw

  case "onemis*":
       set F1="imp_"
#       set NC=1
#       set NO=2
       set NC=2
       set NO=4
       breaksw

  case "ch4*":
       set F1="c_"
#       set NC=1
#       set NO=2
       set NC=2
       set NO=4
       breaksw

  case "rad*":
       set F1="r_inp"
#       set NC=1
#       set NO=2
       set NC=2
       set NO=4
       breaksw

  default:
       echo '*** ERROR: unrecognized namelist file '$nmlfile
       exit 1
       breaksw
endsw

# CREATE LIST OF REQUESTED OBJECTS

#set cha_list = (`grep -v '^[ ]*!' $nmlfile | sed 's|\\!.*$||g' | grep -i '^[ ]*'$F1 | awk -F '=' '{print $2}' | awk -v FS=''',[ ]*''|^''|''$' '{print $'$NC'}' | sed 's|[ ,'\'']||g'`)

#set obj_list = (`grep -v '^[ ]*!' $nmlfile | sed 's|\\!.*$||g' | grep -i '^[ ]*'$F1 | awk -F '=' '{print $2}' | awk -v FS=''',[ ]*''|^''|''$' '{print $'$NO'}' | sed 's|[ ,'\'']||g'`)

set cha_list = (`grep -v '^[ ]*!' $nmlfile | sed 's|\\!.*$||g' | grep -i '^[ ]*'$F1 | awk -F '=' '{print $2}' | awk -F \' '{n=split($'$NO',a,";"); for (i=1;i<=n;i++) print $'$NC'}'`)

set obj_list = (`grep -v '^[ ]*!' $nmlfile | sed 's|\\!.*$||g' | grep -i '^[ ]*'$F1 | awk -F '=' '{print $2}' | awk -F \' '{n=split($'$NO',a,";"); for (i=1;i<=n;i++) print a[i]}'`)

if (${#cha_list} != ${#obj_list}) then
   echo " *** ERROR: Number of channels (${#cha_list}) does not match number of objects (${#obj_list})."
   exit 1
endif

#echo $cha_list
#echo $obj_list
#exit

# CREATE LIST OF OBJECTS DELIVERED FROM IMPORT
echo "Executing"
echo "   $bpath/xcheck_import_grid.tcsh $inmlfile $INPUTDIR_MESSY"
echo "to create list of objects delivered by import."
echo -n "This might take a while ..."

set o_list = (`$bpath/xcheck_import_grid.tcsh $inmlfile $INPUTDIR_MESSY | grep "OBJECT NAME :" | awk -F ':' '{print $2}'`)
#set o_list = (`cat log | grep "OBJECT NAME:" | awk -F ':' '{print $2}'`)

echo " ... done."
echo '-----------------------------------------------------------------------'

@ n = ${#cha_list}
@ cnt = 1
while ($cnt <= $n)

   echo -n "OBJECT $obj_list[$cnt] requested from CHANNEL $cha_list[$cnt] ..."

   switch ("$cha_list[$cnt]")

     case "import_grid":
          @ c2 = 1
          @ found = 0
          while ($c2 <= ${#o_list}) 
            if ("$obj_list[$cnt]" == "$o_list[$c2]") then
               @ found ++
            endif
          @ c2 ++
          end
          if ($found == 0) then
             echo ' *** ERROR: not found in import.'
          endif
          if ($found == 1) then
             echo ' OK.'
          endif
          if ($found > 1) then
             echo ' *** WARNING: found '$found' times.'
          endif

          breaksw

     default:
          echo " *** WARNING: please check manually"
          breaksw

   endsw

@ cnt ++
end

echo '======================================================================='
echo ' FINISHED.'
echo '======================================================================='

exit 0
