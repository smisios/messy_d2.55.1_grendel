#! /bin/tcsh -f
# Time-stamp: <2018-01-09 16:35:49 sander>
# _validate.tcsh: validate Fortran code (TABs and long lines)

set tmpfile = "tmp_file"
unset anyerror

##############################################################################

# check if internal links are okay (and not files):
echo -n "Checking internal links      "
echo -n >  $tmpfile
source _internal_links.tcsh
foreach link ($internal_links)
  if ( -e $link) then
    if ( ! -l $link) then
      #echo "not a link"
      echo $link >> $tmpfile
    endif
  else
    if (("$link" !~ *mecca_kpp*) && ("$link" !~ *caaba_*.nc)) then
      #echo "dead link"
      echo $link >> $tmpfile
    endif
  endif
end
if ( -s $tmpfile != 0 ) then
  echo
  echo "The following should be links but they are not:"
  echo "------------------ start of listing -----------------"
  cat $tmpfile
  echo "------------------  end of listing  -----------------"
  echo "(Note: If _internal_links.tcsh is not up-to-date, run"
  echo "       _update_internal_links.tcsh)"
  set anyerror
else
  echo "OK"
endif

##############################################################################

echo -n "Checking for TABs...         "
grep "	" *.f90 >  $tmpfile
grep "	" *.nml >> $tmpfile
grep "	" *.inc >> $tmpfile
if ( -s $tmpfile != 0 ) then
  echo "The following lines contain TABs:"
  echo "------------------ start of TAB listing -----------------"
  cat $tmpfile
  echo "------------------  end of TAB listing  -----------------"
  set anyerror
else
  echo "OK"
endif

##############################################################################

echo -n "Checking for long lines...   "
# check for any lines that are too long:
# grep '^.\{133,\}' $filelist_f90 > $tmpfile
# check only for non-comments that are too long:
grep '^[^\!]\{133,\}' *.f90 *.nml *.inc > $tmpfile
if ( -s $tmpfile != 0 ) then
  echo "The following lines in f90 files have more than"
  echo "132 characters (comments are not counted here):"
  echo "--------------- start of long line listing --------------"
  cat $tmpfile
  echo "---------------  end of long line listing  --------------"
  set anyerror
else
  echo "OK"
endif

##############################################################################

rm $tmpfile
if (${?anyerror}) exit 1
exit
