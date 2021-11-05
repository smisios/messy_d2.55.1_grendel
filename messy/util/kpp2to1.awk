# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2006-2008
#
# Time-stamp: <2011-02-22 16:50:03 sander>
#
# kpp2to1.awk converts KPP-2 eqn file to KPP-1 syntax
#
# usage:
#   gawk -f kpp2to1.awk infile.eqn
# or:
#   gawk -f kpp2to1.awk infile.rpl
#
# TODO:
# check conversion between "HET" and "PSC"
# further changes necessary for aq rxns:
# (01) -> (as)
# (02) -> (cs)
# _as} -> as}
# _cs} -> cs}
# (01, -> (as,
# (02, -> (cs,

# ----------------------------------------------------------------------------

BEGIN {
  outfile = ARGV[1] "_kpp1"
  print "\nconverting " ARGV[1] " to " outfile "\n"
  printf "" > outfile
}

# ----------------------------------------------------------------------------

{

  # for aqueous phase:
  currentline = gensub("{<([A-Za-z_0-9]+)_a##> ", "{#\\1}{", "g")
  currentline = gensub("<([A-Za-z_0-9]+)_a##>", "{#\\1}", "g", currentline)
  currentline = gensub("a##", "##", "g", currentline)
  # for gas phase:
  currentline = gensub("{<([A-Za-z_0-9]+)> ", "{#\\1}{", "g", currentline)
  currentline = gensub("<([A-Za-z_0-9]+)>", "{#\\1}", "g", currentline)
  currentline = gensub("F90_GLOBAL", "F95_DECL", "g", currentline)
  currentline = gensub("F90_RCONST", "F95_RCONST", "g", currentline)
  currentline = gensub("myhline", "hline", "g", currentline)
  currentline = gensub("([^A-Za-z_0-9])ind_", "\\1KPP_", "g", currentline)

  # only for *.rpl files (to add a new reaction):
  currentline = gensub("<>", "{#}", "g", currentline)

  printf "%s\n", currentline >> outfile

}

# ----------------------------------------------------------------------------
