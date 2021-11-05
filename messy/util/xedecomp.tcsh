#! /bin/tcsh -f

##############################################
### Author: Patrick Joeckel, DLR, Sept 2016
### Use with -h option for explanation.
##############################################

### FOR DEBUGGING ONLY
# set echo verbose

###########################################################################
### TODO: too many options are still displayed (see also qqq)
###
### EXAMPLE USAGE:
### ./messy/util/xedecomp.tcsh -t 42 -l 90 -tpn 36 -n 10
###         #| grep YAXT | grep optimal
###
###########################################################################

############################
### INIT
############################

set nm    = -1
set nlev  = -1
set nodes = -1
set tasks_per_node = -1

set script = `basename $0`

############################
### USER INTEFACE
############################
while ($# > 0)
    switch ($1)

      case '-h':
        echo
        echo "Usage: $script [-h] -t <TRUNCTION> -l <levels> -n <nodes> -tpn <tasks-per-node>"
        echo
        echo "This script provides a listing of all possibilities for the"
        echo "parallel domain decomposition (distributed memory, i.e. MPI)"
        echo "of ECHAM5 and EMAC."
        echo
        echo "Required options:"
	echo " -t  : triangular truncation (without leading T)"
        echo " -l  : number of levels"
        echo " -n  : number of nodes to be used at maximum"
        echo " -tpn: number of MPI tasks to be run in each node"
        echo
        exit 0
        breaksw

      case '-t':
        shift
        set nm = $1
        shift
        breaksw

      case '-l':
        shift
        set nlev = $1
        shift
        breaksw

      case '-n':
        shift
        set nodes = $1
        shift
        breaksw

      case '-tpn':
        shift
        set tasks_per_node = $1
        shift
        breaksw

      default:
        echo "Unknown option: $1"
        echo "Use $script -h for more information\!"
        exit 1
        breaksw
    endsw
end

############################
### CHECK
############################

if ($nm <= 0) then
   echo "Error: truncation missing (-t option)"
   echo "Use $script -h for more information\!"
   exit 1
endif 

if ($nlev <= 0) then
   echo "Error: number of levels missing (-l option)"
   echo "Use $script -h for more information\!"
   exit 1
endif 

if ($nodes <= 0) then
   echo "Error: number of nodes missing (-n option)"
   echo "Use $script -h for more information\!"
   exit 1
endif 

if ($tasks_per_node <= 0) then
   echo "Error: number of tasks_per_node missing (-tpn option)"
   echo "Use $script -h for more information\!"
   exit 1
endif 

############################
### DERIVE SOME PARAMETERS
############################

@ nmp1 = $nm + 1

set nlat = `echo $nm | awk '{m0=int(3*$1/2); if (m0 % 2 == 0) {print m0} else {print m0+1}}'`

set nlath = `echo $nlat | awk '{print $1/2}'`

set nlon = `echo $nlat | awk '{print $1*2}'`

############################
### GENERAL OUTPUT
############################

echo
echo "TRIANGULAR TRUNCATION: "$nm
echo "NUMBER OF LATITUDES:   "$nlat
echo "NUMBER OF LONGITUDES:  "$nlon
echo "NUMBER OF LEVELS:      "$nlev
echo
echo "CONSTRAINTS:"
echo "  1) NPY (nproca) <= nm+1 (="$nmp1")"
echo "  2) NPX (nprocb) <= nlev (="$nlev")"
echo "  3) nlat(="$nlat")/NPY(=nproca) >= 4, if YAXT is NOT used"
echo "  4) nlat(="$nlat")/NPY(=nproca) >= 2"
echo "  5) minimim 2 subdomains for FFSL advection"
echo "     (not checked by this script)"
echo
echo "POSSIBLE DECOMPOSTION FOR NTASK MPI-TASKS:"
echo
echo " NODES NTASKS  NPX  NPY REMARKS"
echo " --------------------------------------------------------------------------"

### loop over nodes
@ ncnt = 1
while ($ncnt <= $nodes)

  set ntasks = `echo $ncnt $tasks_per_node | awk '{print $1*$2}'`

  set npy_list = ()
  set npx_list = ()
  set yxt_list = ()

  @ npy = 1
  # NPY (= NPROCB) MUST BE <= nm+1
  while ($npy <= $nmp1)
   ### NPX == NTASKS / NPY must be integer ...
   set npx = `echo $ntasks $npy | awk '{if ($1%$2 == 0) {print $1/$2} else {print 0};}'`
#qqq+
#   ### ... AND (NPX <= nlev) 
#   if (($npx != 0) && ($npx <= $nlev)) then
   ### ... AND (NPX <= nlev) AND (NPY <= nlat/2)
   if (($npx != 0) && ($npx <= $nlev) && ($npy <= $nlath)) then
#qqq-
      set yxt = `echo $nlat $npy | awk '{r=$1/$2; if (r>=4) {print 0} else {print 1}}'`
      set npy_list = ($npy_list $npy)
      set npx_list = ($npx_list $npx)
      set yxt_list = ($yxt_list $yxt)
   endif
   @ npy ++
  end

  ### check for max. NPY without YAXT
  @ max_no = 0
  @ max_ny = 0
  @ no  = ${#npy_list}
  @ cnt = 1
  while ($cnt <= $no)
       if (($npy_list[$cnt] > $max_ny) && ($yxt_list[$cnt] == 0)) then
          @ max_no = $cnt
          @ max_ny = $npy_list[$cnt]
       endif
  @ cnt ++
  end

  ###########################
  ### LISTING
  ###########################

  @ no  = ${#npy_list}
  @ cnt = 1
  while ($cnt <= $no)
     echo $ncnt $ntasks $npx_list[$cnt] $npy_list[$cnt] | awk '{printf(" %5i %6i %4i %4i",$1,$2,$3,$4)}'
    if ($yxt_list[$cnt] == 1) then
       echo -n ' with YAXT only'
    else
       echo -n '               '
    endif
#    if (($cnt ==  $no) || ($cnt == $max_no)) then
#       echo " (presumably optimal due to max. NPY)"
#    else
       echo
#    endif
  @ cnt ++
  end

@ ncnt++
end
#### end loop over nodes

exit 0
