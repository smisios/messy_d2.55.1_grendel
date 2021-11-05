# -*- awk -*-
# ## THIS GAWK-SCRIPT SEPARATES THE FORCHECK OUTPUT (.lst) INTO
# ## SEPARATE FILES (.log)
# ## USAGE: gawk -f lst2log.gawk <forcheck-lst-file>
# ## Author: Patrick Joeckel, MPICH, June 2003
#
BEGIN {
  NEWPAGE = 1 ;
  OLDFILE = " ";
  SUMFILE = "fchk_forcheck.log" ;
  LISTFILE = "fchk_filelist.log" ;
  FILE = " " ;
}
#
{
  if ($1 == "")
    {
      NEWPAGE = 1 ;
      getline ;
    } ;
  if (NEWPAGE == 1) {
    if (FILE == " ") {print > SUMFILE};
    getline ;
    if (FILE == " ") {print >> SUMFILE};
    FILE = $NF ;
    sub("\\.f90", "."BASEMODEL".log", FILE) ;
    if (FILE == "" || substr(FILE,1,1) == "-") {FILE = SUMFILE} ;
    NEWPAGE = 0 ;
    if (FILE != OLDFILE)
      {
	if (FILE != SUMFILE) print "! -*- f90 -*-" > FILE ;
	printf("%s\n", FILE) >> LISTFILE ;
	OLDFILE = FILE ;
      }
    getline ;
  }
  print >> FILE ;
}
