# ----------------------------------------------------------------------------

# Author: Holger Tost, Max-Planck-Institute, Mainz, Germany, 2004
# Time-stamp: <2006-07-18 14:45:38 tost>

# tracmvscav.awk extracts all KPP_XXX from messy_scav_kpp_g_mem.f90 that are not zero
# and inserts these species into several files for calculation of scavenging and
# coupling to gasphase and aerosol scavenging
# the produced files are modified into a module messy_scav_auto_inter.f90 by the
# create_scav_auto script that calls this awk file

# 1) messy_scav_outkpp_e5.inc:    conversion kpp species to tracers
# 2) messy_scav_inkpp_e5.inc:     conversion tracers to kpp species
# 3) messy_scav_l2aercs_e5.inc:   conversion liquid kpp species to coarse mode aerosols
# 4) messy_scav_l2aeras_e5.inc:   conversion liquid kpp species to accumlation mode aerosols
# 5) messy_scav_aer2l_e5.inc:     conversion aerosols of as and cs mode to liquid kpp species
# 6) messy_scav_l2gas_e5.inc:     conversion liquid kpp species to gasphase species
# 7) messy_scav_trac_e5.inc:      check for tracers, their scav behaviour and setting their idt_scav
# 8) messy_scav_idt_mem.inc:      defining idt_scav_"tracer"
# 9) messy_scav_str_kpp.inc:      defining a string array with the names of all used KPP species

# usage:
# gawk -f tracmvscav.awk messy_scav_kpp_g_mem.f90

# Normally, however, tracmvscav.awk is called by xmecca1 and there from create_scav_auto, and not started
# directly from the command line.

# ----------------------------------------------------------------------------

BEGIN {
# define names of output files
conc2mrfile = "messy_scav_outkpp_e5.inc"
mr2concfile = "messy_scav_inkpp_e5.inc"
l2aercsfile = "messy_scav_l2aercs_e5.inc"
l2aerasfile = "messy_scav_l2aeras_e5.inc"
aer2lfile   = "messy_scav_aer2l_e5.inc"
l2gasfile   = "messy_scav_l2gas_e5.inc"
tracfile    = "messy_scav_trac_e5.inc"
idtfile     = "messy_scav_idt_mem.inc"
strfile     = "messy_scav_str_kpp.inc"

# write header line
dontedit = "! -*- f90 -*- this file was created by create_scav_auto, do not edit!"
print dontedit > l2aercsfile
print dontedit > l2aerasfile
print dontedit > conc2mrfile
print dontedit > mr2concfile
print dontedit > tracfile
print dontedit > aer2lfile
print dontedit > l2gasfile
print dontedit > idtfile
print dontedit > strfile
counter = 0

errorstring = ""
}

{
# is current line something like "INTEGER ... KPP_XXX = nnn" with nnn != 0 ?
if (match($0, "INTEGER.*KPP_([A-Za-z0-9_]+) = [^0]", arr) != 0) {
  # split arr (field separator is "_" ) and save basename
  # and subname in array basesub
  hfill  = substr("            ", 1, 12-length(arr[1]))
  split(arr[1], basesub, "_")
  hfill2 = substr("            ", 1, 12-length(basesub[1]))
  counter = counter + 1 
    
  if (basesub[2]=="")  {
  
  # ------------------------------------------------------------------------
    # 2) add tracer to messy_scav_outkpp_e5.inc
    printf "        pmx(idt_scav_%s)%s=  D(KPP_%s)\n",
      arr[1], hfill, arr[1] >> conc2mrfile
    # ------------------------------------------------------------------------
    # 3) add tracer to messy_scav_inkpp_e5.inc
    printf "        D(KPP_%s)%s= pmx(idt_scav_%s)\n",
      arr[1], hfill, arr[1] >> mr2concfile
 
  }
 if (basesub[2]=="l")  {
# l2aercs file: shifting species from l to aerosol and resetting l concentration when successful
   printf "     pmx(idt_scav_%s_cs)%s= pmx(idt_scav_%s_cs)%s+ D(KPP_%s)\n",
     basesub[1], hfill2, basesub[1], hfill2, arr[1] >> l2aercsfile
     printf "     D(KPP_%s)%s          = D(KPP_%s)%s * REAL(1-min(1,idt_scav_%s_cs),dp)  \n",
     arr[1], hfill, arr[1], hfill, basesub[1] >> l2aercsfile
# l2aeras file
      printf "     pmx(idt_scav_%s_as)%s= pmx(idt_scav_%s_as)%s+ D(KPP_%s)\n",
      basesub[1], hfill2, basesub[1], hfill2, arr[1] >> l2aerasfile
      printf "     D(KPP_%s)%s          = D(KPP_%s)%s * REAL(1-min(1,idt_scav_%s_as),dp)  \n",
      arr[1], hfill, arr[1], hfill, basesub[1] >> l2aerasfile
    # ------------------------------------------------------------------------
# aer2l file:  shifting species from aerosol_l to l and resetting aer_l concentration 
      printf "     D(KPP_%s)%s= D(KPP_%s)%s+ aer_field(idt_scav_%s_cs)%s + aer_field(idt_scav_%s_as) \n",
      arr[1], hfill, arr[1], hfill, basesub[1], hfill2, basesub[1] >> aer2lfile
      printf "     aer_field(idt_scav_%s_cs)%s = 0.0_dp \n",
      basesub[1], hfill2 >> aer2lfile
      printf "     aer_field(idt_scav_%s_as)%s = 0.0_dp \n",
      basesub[1], hfill2 >> aer2lfile
  }
# l2gas file:  shifting species from l to gas phase tracers
 if (basesub[2]== "") {
   printf "      pmx(idt_scav_%s)%s = pmx(idt_scav_%s)%s + input_field(KPP_%s_l) \n",
     basesub[1], hfill2, basesub[1], hfill2, basesub[1] >> l2gasfile
    printf "      input_field(KPP_%s_l)%s = input_field(KPP_%s_l)%s * REAL(1-min(1,idt_scav_%s),dp)  \n",
     arr[1], hfill2, arr[1], hfill2, basesub[1] >> l2gasfile  }

 if (basesub[2]=="")  {
     printf "      CALL get_tracer(ierr, GPTRSTR, '%s',idx=idt_scav_%s) \n",
     basesub[1], basesub[1] >> tracfile
     printf "      IF (ierr /=0.or..not. attr(1)%%log_att(idt_scav_%s,lwetdep))%s idt_scav_%s    = 0 \n",
       basesub[1], hfill2, basesub[1] >> tracfile
}
 if (basesub[2]=="l")  {
     printf "      CALL get_tracer(ierr, GPTRSTR, '%s', subname ='cs',idx=idt_scav_%s_cs) \n",
     basesub[1], basesub[1] >> tracfile
     printf "      IF (ierr /=0.or..not. attr(1)%%log_att(idt_scav_%s_cs,lwetdep))%s idt_scav_%s_cs    = 0 \n",
     basesub[1], hfill2, basesub[1] >> tracfile
     printf "      CALL get_tracer(ierr, GPTRSTR, '%s', subname ='as',idx=idt_scav_%s_as) \n",
     basesub[1], basesub[1] >> tracfile
     printf "      IF (ierr /=0.or..not. attr(1)%%log_att(idt_scav_%s_as,lwetdep))%s idt_scav_%s_as    = 0 \n",
     basesub[1], hfill2, basesub[1] >> tracfile
}

 printf "       INTEGER, PUBLIC :: idt_scav_%s%s      = 0   \n",
   arr[1], hfill >> idtfile
   if (basesub[2]=="l")  {
     printf "       INTEGER, PUBLIC :: idt_scav_%s_cs%s   = 0   \n",
       basesub[1], hfill2 >> idtfile
       printf "       INTEGER, PUBLIC :: idt_scav_%s_as%s   = 0   \n",
       basesub[1], hfill2 >> idtfile}

 printf "     str_field_kpp(%s) = '%s' \n",
   counter, arr[1] >> strfile

}
}
# ----------------------------------------------------------------------------

END {
print errorstring
#printf "(warnings about %s are only for ECHAM, not for the box model)\n\n",
#  tracdef
}

# ----------------------------------------------------------------------------
