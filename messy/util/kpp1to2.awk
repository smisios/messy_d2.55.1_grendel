# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2006-2008
#
# Time-stamp: <2008-07-23 10:19:25 sander>
#
# kpp1to2.awk converts KPP-1 eqn file to KPP-2 syntax
#
# usage:
#   gawk -f kpp1to2.awk infile.eqn
# or:
#   gawk -f kpp1to2.awk infile.rpl
#
# TODO:
# check conversion between "HET" and "PSC"
#
# ----------------------------------------------------------------------------

BEGIN {
  outfile = ARGV[1] "_kpp2"
  print "\nconverting " ARGV[1] " to " outfile "\n"
  printf "" > outfile
}

# ----------------------------------------------------------------------------

{

  currentline = gensub("{#([A-Za-z_0-9]+)}{", "{<\\1> ", "g")
  currentline = gensub("{#([A-Za-z_0-9]+)}", "<\\1>", "g", currentline)
  currentline = gensub("{#([A-Za-z_0-9]+)##}", "<\\1_a##>", "g", currentline)
  currentline = gensub("A##", "Aa##", "g", currentline)
  currentline = gensub("_##", "_a##", "g", currentline)
  currentline = gensub("F95_DECL", "F90_GLOBAL", "g", currentline)
  currentline = gensub("F95_RCONST", "F90_RCONST", "g", currentline)
  currentline = gensub("hline", "myhline", "g", currentline)
  currentline = gensub("([^A-Za-z_0-9])[Kk][Pp][Pp]_", "\\1ind_", "g", currentline)

  # only for *.rpl files (to add a new reaction):
  currentline = gensub("{#}", "<>", "g", currentline)

  printf "%s\n", currentline >> outfile

}

# ----------------------------------------------------------------------------
