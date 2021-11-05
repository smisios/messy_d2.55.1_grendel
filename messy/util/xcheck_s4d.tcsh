#! /bin/tcsh -f

# ---------------------------------------------------------------------
# Author: Patrick Joeckel, DLR, Aug 2013 
# ---------------------------------------------------------------------
#set echo verbose

echo '================================================================='
echo ' THIS SCRIPT CHECKS THE INPUT FOR S4D ACCORDING TO YOUR S4D SETUP'
echo '================================================================='

if (${#} < 2) then
   echo 'Usage: '`basename $0`' <namelist file> <INPUTDIR_MESSY>'
   echo '        possibe namelist files: s4d*'
   exit 1
endif

set nmlfile = $1
set INPUTDIR_MESSY = $2

if (! -r $nmlfile) then
   echo ' *** ERROR: '$nmlfile' is missing or unreadable.'
   exit 1
endif

set bpath = `dirname $0`

# count number of valid lines in namelist file
set nlines = `sed -n '/^\&CPL/,/^\//p' $nmlfile | grep -v '^[ ]*!' | grep -i TRACK | wc -l`

# loop over valid lines
@ cnt = 1
while ($cnt <= $nlines)
  echo '----------------------------------------------------------------------'

  # SELECT LINE ACCORDING TO LINE NUMBER
  set line = "`sed -n '/^\&CPL/,/^\//p' $nmlfile | grep -v '^[ ]*!' | grep -i TRACK | awk '{if (NR=='$cnt') print}'`"

  #echo "$line"

#  # PARSE LINE
  set no = `echo "$line" | awk -F '=' '{print $1}' | sed 's|[()]| |g' | awk '{print $2}'`

  set name = `echo "$line" | awk -F '=' '{print $2}' | awk -F ',' '{print $1}' | sed 's|'\''||g'`

  set dpath = `echo "$line" | awk -F '=' '{print $2}' | awk -F ',' '{print $2}' | sed 's|'\''||g' | sed 's|$INPUTDIR_MESSY|'$INPUTDIR_MESSY'|g'`

  set mode = `echo "$line" | awk -F '=' '{print $2}' | awk -F ',' '{print $3}'`

  switch($mode)
   case 0:
     #echo daily
     set modestr = daily
     set flist = (`ls ${dpath}????????.pos`)
     breaksw

   case 1:
     #echo monthly
     set modestr = monthly
     set flist = (`ls ${dpath}??????.pos`)
     breaksw

   default:
     echo ' *** ERROR: Unknown '$mode'; must be 0 (daily) or 1 (monthly).'
     breaksw
  endsw

  set nf = ${#flist}

  echo "TRACK:          " $no
  echo "NAME:           " $name 
  echo "DATAPATH:       " $dpath
  echo "MODE:           " $modestr 
  echo "Number of files:" $nf
  if ($nf <= 0) then
     echo ' *** WARNING: Number of input files is '$nf
  endif

#  echo "List of files:  "   
#   foreach file ($flist)
#     echo "   "$file
#   end

  echo '----------------------------------------------------------------------'
  @ cnt ++
end

echo '======================================================================='
echo ' FINISHED.'
echo '======================================================================='



exit 0
