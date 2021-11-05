#! /bin/tcsh -f

### SCRIPT TO RE-CREATE ALL switch.nml FROM messy/nml/EXAMPLES/switch.nml
### (required for dau-ngrading)
### (C) Patrick Joeckel, DLR, 2011

if (! -d messy/nml) then
   echo "Call this script from BASEDIR of distribution:"
   echo "./messy/util/`basename $0` <option>"
   exit 1
endif

set master = messy/nml/EXAMPLES/switch.nml

if (! -e $master) then
   echo "ERROR: $master does not exist\!"
   exit 1
endif

if ($# < 1) then
   echo "Usage: ./messy/util/`basename $0` <option>"
   echo "       Options:"
   echo "                -t test (use always first to check the changes)"
   echo "                -m modify all switch.nml files by activating"
   echo "                   switches based on $master"
   exit 1
endif

# list of all switch.nml
set swf_list = (`find . -name 'switch.nml' -print`)
#set swf_list = (messy/nml/E5/switch.nml)

# loop over all switch.nml
foreach swf ($swf_list)

  # check, what's ON
  set use_list = (`grep -e '^[ ]*USE_.*=\.\{0,1\}[T,t]' $swf | awk -F '=' '{print $1}'`)

  echo "##############################################################"
  if ("$1" == "-t") then 
     echo "### TESTING $swf"
  endif

  if ("$1" == "-m") then 
     echo "### MODIFYING $swf"
  endif

  cp -f $master ${swf}.new
  foreach use ($use_list)
    # echo $use
    cat ${swf}.new | sed 's|!'$use'=.TRUE.|'$use'=.TRUE.|g' >! ${swf}.tmp
    mv -f ${swf}.tmp ${swf}.new    
  end

  set use_list_new = (`grep -e '^[ ]*USE_.*=\.\{0,1\}[T,t]' ${swf}.new | awk -F '=' '{print $1}'`)

  echo $use_list     | tr ' ' '\n' | sort | uniq >! ${swf}.1.tmp
  echo $use_list_new | tr ' ' '\n' | sort | uniq >! ${swf}.2.tmp

  diff ${swf}.1.tmp ${swf}.2.tmp
  rm -f ${swf}.1.tmp ${swf}.2.tmp

  if ("$1" == "-t") then 
     rm -f ${swf}.new 
  endif

  if ("$1" == "-m") then 
     mv -f ${swf}.new ${swf}
  endif

  echo "##############################################################"

end
