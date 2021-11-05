#!/bin/sh -e

DIR=`pwd | sed 's|/.*/||'`
FIL=$DIR'_1r'
if test -r $FIL.zip ; then
  echo "zip file exists! Please (re)move it first:"
  echo "rm $FIL.zip" 
  exit 1
fi
cd data
if test -r $FIL.zip ; then
  echo "zip file exists! Please (re)move it first:"
  echo "rm ../$FIL.zip" 
  exit 1
fi

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -y store symbolic links as the link instead of the referenced file
# -0 store only (do not compress)

zip -ory0 $FIL.zip forcing residui nml ESH_NO rerun* *.rst restart*
mv -f $FIL.zip ..

exit 0
