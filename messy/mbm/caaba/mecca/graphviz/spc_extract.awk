# ----------------------------------------------------------------------------
#
# Author:
#   Rolf Sander, Max-Planck-Institute, Mainz, Germany, 2009
#
# Time-stamp: <2018-12-28 10:32:46 sander>
#
# spc_extract.awk extracts subsets of species from *.spc file
#
# usage: invoke via xgraphvizall
#
# ----------------------------------------------------------------------------

BEGIN {
  #printf "working on %s...\n", ARGV[1]
  logfile    = "spc_extract.log"
  all_spc    = "all.spc"
  C_spc      = "C.spc"
  C1_spc     = "C1.spc"
  C2_spc     = "C2.spc"
  C3_spc     = "C3.spc"
  C4_spc     = "C4.spc"
  C5_spc     = "C5.spc"
  C6_spc     = "C6.spc"
  C7_spc     = "C7.spc"
  C8_spc     = "C8.spc"
  C9_spc     = "C9.spc"
  C10_spc    = "C10.spc"
  N_spc      = "N.spc"
  Cl_spc     = "Cl.spc"
  Br_spc     = "Br.spc"
  I_spc      = "I.spc"
  S_spc      = "S.spc"
  Fe_spc     = "Fe.spc"
  OH_HO2_spc = "OH_HO2.spc"
  dontedit   = "created automatically by spc_extract.awk, DO NOT EDIT!"
  printf "%s\n",    dontedit > logfile
  printf "// %s\n", dontedit > all_spc
  printf "// %s\n", dontedit > C_spc
  printf "// %s\n", dontedit > C1_spc
  printf "// %s\n", dontedit > C2_spc
  printf "// %s\n", dontedit > C3_spc
  printf "// %s\n", dontedit > C4_spc
  printf "// %s\n", dontedit > C5_spc
  printf "// %s\n", dontedit > C6_spc
  printf "// %s\n", dontedit > C7_spc
  printf "// %s\n", dontedit > C8_spc
  printf "// %s\n", dontedit > C9_spc
  printf "// %s\n", dontedit > C10_spc
  printf "// %s\n", dontedit > N_spc
  printf "// %s\n", dontedit > Cl_spc
  printf "// %s\n", dontedit > Br_spc
  printf "// %s\n", dontedit > I_spc
  printf "// %s\n", dontedit > S_spc
  printf "// %s\n", dontedit > Fe_spc
  printf "// %s\n", dontedit > OH_HO2_spc
  # hardcoded files:
  printf "OH\n"  >> OH_HO2_spc
  printf "HO2\n" >> OH_HO2_spc
}

# ----------------------------------------------------------------------------

{
  printf "\nline: |%s|\n", $0 >> logfile
  # delete all comment lines ("//...") from $0:
  gsub("^//.*", "")
  # does current line contain a chemical composition?
  if (match($0, "^ *[A-z0-9_]+ *=[A-z0-9+ ]+;") != 0) {
    # extract the composition "=...;" from current line, and
    # add leading and trailing "+" signs:
    atoms = "+" gensub("^ *[A-z0-9_]+ *=([A-z0-9+ ]+);.*", "\\1", "g", $0) "+"
    printf "atoms = |%s|\n", atoms >> logfile
    # remove all spaces:
    gsub(" +", "", atoms)
    printf "atoms = |%s|\n", atoms >> logfile
    # if suitable, write line to spc files:
    if (match($0, "^RR") == 0) {
      printf "%s\n", $0 >> all_spc
    }
    if (match(atoms, "\\+[0-9]*C\\+") != 0) {
      printf "all C = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C_spc
    }
    if (match(atoms, "\\+[1]?C\\+") != 0) {
      printf "C1    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C1_spc
    }
    if (match(atoms, "\\+2C\\+") != 0) {
      printf "C2    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C2_spc
    }
    if (match(atoms, "\\+3C\\+") != 0) {
      printf "C3    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C3_spc
    }
    if (match(atoms, "\\+4C\\+") != 0) {
      printf "C4    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C4_spc
    }
    if (match(atoms, "\\+5C\\+") != 0) {
      printf "C5    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C5_spc
    }
    if (match(atoms, "\\+6C\\+") != 0) {
      printf "C6    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C6_spc
    }
    if (match(atoms, "\\+7C\\+") != 0) {
      printf "C7    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C7_spc
    }
    if (match(atoms, "\\+8C\\+") != 0) {
      printf "C8    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C8_spc
    }
    if (match(atoms, "\\+9C\\+") != 0) {
      printf "C9    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C9_spc
    }
    if (match(atoms, "\\+10C\\+") != 0) {
      printf "C10   = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> C10_spc
    }
    if (match(atoms, "\\+[0-9]*N\\+") != 0) {
      printf "N     = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> N_spc
    }
    if (match(atoms, "\\+[0-9]*Cl\\+") != 0) {
      printf "Cl    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> Cl_spc
    }
    if (match(atoms, "\\+[0-9]*Br\\+") != 0) {
      printf "Br    = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> Br_spc
    }
    if (match(atoms, "\\+[0-9]*I\\+") != 0) {
      printf "I     = TRUE\n", atoms >> logfile
      printf "%s\n", $0 >> I_spc
    }
    if (match(atoms, "\\+[0-9]*S\\+") != 0) {
       printf "S     = TRUE\n", atoms >> logfile
       printf "%s\n", $0 >> S_spc
    }
    if (match(atoms, "\\+[0-9]*Fe\\+") != 0) {
       printf "Fe    = TRUE\n", atoms >> logfile
       printf "%s\n", $0 >> Fe_spc
    }
  }
}

# ----------------------------------------------------------------------------
