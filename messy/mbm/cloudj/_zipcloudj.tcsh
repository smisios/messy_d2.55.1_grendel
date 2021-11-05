#! /bin/tcsh -f
# Time-stamp: <2016-09-27 14:37:40 sander>
# _zipcloudj.tcsh: create a zip file of cloudj code

# if the current directory is only a link, the following cd command
# will jump into the real directory:
cd `/bin/pwd`

set dirname = `basename $PWD`
set zipfile = $PWD/$dirname.zip
if ( -e $zipfile) then
  echo "Renaming old $dirname.zip file to $dirname.zip~"
  mv -f $zipfile $zipfile~
endif

cd ..

# zip options:
# -o make zipfile as old as latest entry
# -r recurse into directories
# -x '...' exclude files
zip -or $zipfile $dirname \
  -x '*~' -x '*.mod' -x '*.exe' -x '*.o' -x '*.a' -x '*.nc' \
  -x '*.log' -x '*.old' -x '*/ferret.jnl' -x '*.zip' -x '*.tar' \
  -x '*.ps' \
  -x '*/tmp/*'

# Symbolic links are now included as the whole files in the zip archive.
# However, some links are internal links (i.e. between directories that
# are both included in the zip file). They must be stored as links. This
# is done by overwriting them in the zip file:

# zip options:
# -y store symbolic links as the link instead of the referenced file
zip -oy $zipfile \
$dirname/messy_cloudj.f90

echo;echo "The zipfile has been created:"
ls -la $zipfile

exit
