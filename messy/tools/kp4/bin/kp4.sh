#!/bin/sh -e

###########################################################################
#
#     create_kpp_module
#
#     create code from .f90 sources created by KPP to be used in MESSy
#
#     COPYRIGHT Klaus Ketelsen and MPI-CH   April 2007
#
###########################################################################

#set -eux
set -eu

########################### KPP setup #####################################

KPP_HOME=../../kpp
export KPP_HOME
KPP=$KPP_HOME/bin/kpp
export KPP

########################### working setup #################################

BASE=`pwd`

# default setup
OUTDIR=$BASE
MODE="scalar"
VLEN=1
KEEP="NO"
DE_INDEX=1

# get command line options

while  getopts :m:s:i:v:o:  c       # get options
do case $c in
      m)   SUBMODEL=$OPTARG;;        # which submodel

      s)   KPP_SOLVER=$OPTARG;;      # which kpp solver

      i)   DE_INDEX=$OPTARG;;        # deindexing stage

      v)   MODE="vector"
	   VLEN=$OPTARG;;            # set to vector mode

      o)   OUTDIR=$OPTARG;;          # output directory of generated code

      \?)  echo ${0##*/} "unknown option:" $OPTARG
	 echo "USAGE: ${0##*/} -m submodel -s solver [-i stage -v length -o dir]"
	   exit 1;;
   esac
done

# at least submodel and solver MUST be set
if test "${SUBMODEL:-set}" = "set" ; then
  echo "ERROR: -m option must be specified"
  exit 1
fi
if test "${KPP_SOLVER:-set}" = "set" ; then
  echo "ERROR: -s option must be specified"
  exit 1
fi

# settings
export KPP_SOLVER
PREFIX=messy_${SUBMODEL}_kpp
DEF_PREFIX=${PREFIX}.kpp

case $SUBMODEL in
     mecca*|mtchem)
	 DEFDIR="caaba/mecca"
	 ;;
     scav)
	 DEFDIR="scav/mechanism"
	 ;;
     gmxe_aerchem)
	 DEFDIR="gmxe/aerchem/mechanism"
	 ;;
     *)
	 echo "ERROR: unknown submodel "$SUBMODEL
	 exit 1
	 ;;
esac

echo ${DEFDIR}

# no need for a submodel-specific workdir
WORK=./tmp_workdir

# create or clean working directory
mkdir -p $WORK
rm -rf $WORK/*

# go to working directory
cd $WORK

# depending on kpp version; this might be subject to change

KPP_FILE_LIST="Initialize Integrator LinearAlgebra"
KPP_FILE_LIST="$KPP_FILE_LIST Jacobian Function Rates Util"

KPP_SUBROUTINE_LIST="Initialize"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST INTEGRATE Fun"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Fun_SPLIT"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST CalcStoichNum"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolve KppDecomp WLAMCH WLAMCH_ADD"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Jac_SP Update_RCONST"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST initialize_indexarrays"
KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST initialize_kpp_ctrl"

case $SUBMODEL in
    mecca*|mtchem)
	# from ./messy/tools/kpp/util/UserRateLaws.f90:
	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_arr"
	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_3rd"
# from #INLINE F90_RATES in gas.eqn
#	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_SIV_H2O2"
#	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_limited k_3rd_iupac"
#	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST alpha_AN"
#	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_Op_O2 k_N2_O k_Op_N2"
#	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_RO2_HO2"
	;;
    scav|gmxe_aerchem)
	KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST k_arr"
	;;
esac

if test "$MODE" = "vector" ; then
   # get vector solver
   cp $BASE/templates/${KPP_SOLVER}_vec.f90 ${KPP_SOLVER}.f90
fi

# interface ignore list
KPP_INTERFACE_IGNORE="WAXPY WCOPY"

case $KPP_SOLVER in

    rosenbrock*)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      if test "$MODE" != "vector" ; then
	 KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock"
	 KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FunTemplate JacTemplate"
	 KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Update_SUN"
	 KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST error_output"
      fi;;

#    rosenbrock_mz)
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock FunTemplate"
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST JacTemplate Update_SUN";;

#    rosenbrock)
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST Rosenbrock  FunTemplate"
#      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST JacTemplate Update_SUN";;

    kpp_lsode)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppLsode DLSODE JAC_CHEM"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FUN_CHEM"
      KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE JAC_CHEM KppDecomp"
      KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE KppSolve";;

    kpp_radau5)
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY FUN_CHEM"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST JAC_CHEM SET2ZERO"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST RADAU5 Update_SUN"
      KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST KppSolveCmplx KppDecompCmplx";;

    kpp_sdirk)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SDIRK JAC_CHEM SET2ZERO"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST FUN_CHEM"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE Set2zero SET2ZERO"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE FUN_CHEM";;

    kpp_seulex)
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST WCOPY WSCAL WAXPY"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST ATMSEULEX"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SEULEX_ErrorMsg"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST SEULEX_Integrator FUN_CHEM"
       KPP_SUBROUTINE_LIST="$KPP_SUBROUTINE_LIST JAC_CHEM SEUL"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE SEULEX_Integrator SDIRK"
       KPP_INTERFACE_IGNORE="$KPP_INTERFACE_IGNORE FUN_CHEM SEUL";;

   \?)  echo "THE SELECTED SOLVER IS UNFORTUNATELY NOT YET SUPPORTED BY KP4:"
	echo $KPP_SOLVER
	exit 1;;
esac

KPP_INCLUDE_LIST="Precision Parameters Global JacobianSP Monitor"

# get definition Files

cp ../integr/${KPP_SOLVER}.kpp integr.kpp
#cp ../../../mbm/caaba/mecca/${SUBMODEL}.eqn .
#cp ../../../mbm/caaba/mecca/${SUBMODEL}.spc .
#cp ../../../mbm/caaba/mecca/${PREFIX}.kpp   .
cp ../../../mbm/${DEFDIR}/${SUBMODEL}.eqn .
cp ../../../mbm/${DEFDIR}/${SUBMODEL}.spc .
cp ../../../mbm/${DEFDIR}/${PREFIX}.kpp .

cd $KPP_HOME
KPP_HOME=`pwd`
export KPP_HOME
cd -
KPP_VERSION=`grep KPP_VERSION $KPP_HOME/src/gdata.h | sed 's|.*"\(.*\)".*|\1|g'`
echo "#INLINE F90_GLOBAL"                                          >> integr.kpp
echo "  ! KPP info from xmecca:"                                   >> integr.kpp
echo "  CHARACTER(LEN=*), PUBLIC, PARAMETER :: &"                  >> integr.kpp
echo "    ${SUBMODEL}_spc_file     = '`ls -l ${SUBMODEL}.spc`', &" >> integr.kpp
echo "    ${SUBMODEL}_eqn_file     = '`ls -l ${SUBMODEL}.eqn`', &" >> integr.kpp
echo "    ${SUBMODEL}_spc_file_sum = '`sum ${SUBMODEL}.spc`', &"   >> integr.kpp
echo "    ${SUBMODEL}_eqn_file_sum = '`sum ${SUBMODEL}.eqn`', &"   >> integr.kpp
echo "    kppoption          = '4', &"                             >> integr.kpp
echo "    KPP_HOME           = '$KPP_HOME', &"                     >> integr.kpp
echo "    KPP_version        = '$KPP_VERSION', &"                  >> integr.kpp
echo "    integr             = '$KPP_SOLVER'"                      >> integr.kpp
echo "#ENDINLINE {above lines go to messy_${SUBMODEL}_kpp_global}" >> integr.kpp

# run kpp
echo $KPP $DEF_PREFIX
echo "run kpp"
pwd
$KPP $DEF_PREFIX
# TODO: CHECK KPP OUTPUT STATUS HERE ...

# get templates for C++ program
echo "before header!"
case $SUBMODEL in
    mecca*|mtchem)
	# use specific template for module header
	cp $BASE/templates/module_header .
	;;
    scav)
	# use specific template for module header
	cp $BASE/templates/module_header_scav ./module_header
	;;
    gmxe_aerchem)
	# use specific template for module header
	cp $BASE/templates/module_header_gmxe_aerchem ./module_header
	;;
esac

# CTRL kpp time stepping
cp $BASE/templates/initialize_kpp_ctrl_template.f90 .

# file with subroutine list for c++ program create_kpp_module
for i in $KPP_FILE_LIST
do
  echo ${PREFIX}_${i} >> file_list
done
echo initialize_kpp_ctrl_template >> file_list

# file with subroutine list for c++ program create_kpp_module
for i in $KPP_SUBROUTINE_LIST
do
  echo $i >> subroutine_list
done

# file with include list for c++ program create_kpp_module
for i in $KPP_INCLUDE_LIST
do
  echo ${PREFIX}_${i} >> include_list
done

touch interface_ignore_list
for i in $KPP_INTERFACE_IGNORE
do
  echo $i >> interface_ignore_list
done

# $BASE/src/kp4.exe $PREFIX $MODE $VLEN $DE_INDEX $DE_INDEX_FAST
echo " ... running kp4.exe $PREFIX $MODE $VLEN $DE_INDEX"
$BASE/../../../bin/kp4.exe $PREFIX $MODE $VLEN $DE_INDEX

cp kk_kpp.f90    $OUTDIR/${PREFIX}.f90

exit 0
