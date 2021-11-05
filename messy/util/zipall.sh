#!/bin/sh -e

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

DIR=`pwd | sed 's|/.*/||'`
FIL=$DIR
if test -r $FIL.zip ; then
  echo "zip file exists! Please (re)move it first:"
  echo "rm $FIL.zip" 
  exit 1
fi
cd ..
if test -r $FIL.zip ; then
  echo "zip file exists! Please (re)move it first:"
  echo "rm ../$FIL.zip" 
  exit 1
fi

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -y store symbolic links as the link instead of the referenced file

zip -ory $FIL.zip $DIR
mv -f $FIL.zip $DIR/.

exit 0
