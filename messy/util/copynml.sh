#!/bin/sh -e

echo "copying namelist file $1"
if test ! -r $1 ; then
   echo '... namelist file not present'
   exit 1
fi

echo 'cat > $2 << EOF'  >  temporaryfile
echo '! This file was created automatically by copynml.sh, do not edit' \
                        >> temporaryfile
cat $1 | sed 's|!.*||g' \
       | sed 's|( *\([0-9]*\))|(\1)|g' \
       | grep -Ev '^ *$' >> temporaryfile
echo 'EOF' >> temporaryfile

# "." = "source"
. ./temporaryfile
rm -f temporaryfile

#echo '................................................................'
#cat  $2
#echo '................................................................'

exit 0
