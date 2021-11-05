#! /bin/tcsh -f

# ---------------------------------------------------------------------
# Author: Patrick Joeckel, DLR, Nov 2010 
# ---------------------------------------------------------------------

echo '======================================================================='
echo ' THIS SCRIPT CHECKS THE CONSISTENCY OF YOUR IMPORT_TS SETUP:'
echo '  CTRL_TS in import.nml must contain a unique, active entry for all    '
echo '  active import_ts requests of all *.nml'
echo '======================================================================='

if (${#} < 1) then
   echo 'Usage: '`basename $0`' <namelist file>'
   echo '        possibe namelist files: import*'
   exit 1
endif
set nmlfile = $1

foreach file (*.nml)

 if ("$file" == "import.nml") continue
 if ("$file" == "channel.nml") continue
 if ("$file" == "$nmlfile") continue

 set tslist = (` grep 'import_ts' $file | grep -v '^[ ]*!' | awk -F ',' '{print $2}'`)

# if ("$tslist" != "") then
#    echo -n $file": "
#    echo $tslist
# endif

  if ("$tslist" != "") then
     foreach name ($tslist)
       echo '----------------------------------------------------------'
       echo "checking for $name ($file) in CTRL_TS of import.nml:"
       grep -i 'ts(' $nmlfile | grep -v '^[ ]*!' | grep $name
       echo '----------------------------------------------------------'
     end
  endif

end

exit 0
