#!/bin/tcsh -f

# list of LEGACY BASE MODELS
set bmslist = (`find . -maxdepth 1 -type l -print | sort`)

# list of TOOLS
set here=`pwd`
cd messy/tools
set toolslist = (`find . -maxdepth 1 -type d -print | sort | sed 's|tag|imtag embudget |g'`)
cd $here

# list of MESSY BASE MODELS
set here=`pwd`
cd messy/mbm
set mbmlist = (`find . -maxdepth 1 -type d -print | sort`)
cd $here


### loop over lists
@ i=1
while ($i <= 3)
  switch ($i)
  case 1: 
         set header = "### LEGACY BASE MODELS ###"
         set list = ($bmslist)
         breaksw
  case 2: 
         set header = "### MESSy TOOLS ###"
         set list = ($toolslist)
         breaksw
  case 3: 
         set header = "### MESSy BASE MODELS ###"
         set list = ($mbmlist)
         breaksw
  endsw

  echo '###------------------------------------------------'
  echo $header
  echo '###------------------------------------------------'

  ### loop over list entries
  @ n=${#list}
  @ j=1
  while ($j <= $n)
    set ex = `basename $list[$j]`.exe
    if ("$ex" == "..exe") then
       @ j++
       continue
    endif
    echo -n $ex | awk '{printf("%26s "),$1}'

    if (-x bin/$ex) then
       echo -n '+' | awk '{printf("%5s\n"),$1}'
    else
       echo -n '-' | awk '{printf("%5s\n"),$1}'
    endif

    @ j++
  end ### loop over list entries

  @ i++
end ### loop over lists

echo '###------------------------------------------------'

exit 0
