#!/bin/tcsh -f

################################################################
### THIS SCRIPT IS TO UPDATE RESTART DIRECTORIES CREATED WITH
### VERSION <= 2.52g FOR USAGE WITH VERSION >= 2.52h:
###   - input: - restart_ECHAM5.nc
###            - restart_g3b.nc
###   - ouput: - restart_e5vdiff.nc
###            - restart_orogw.nc
###            - restart_mlocean.nc
################################################################
### requirements: nco (ncks, ncrename, ncatted)
################################################################
### (C) Patrick Joeckel, DLR, June 2016
################################################################

################################################################
### step-by-step recipe
################################################################
### 1) init_restart -d ... -r ... -c ...
### 2) update WORKDIR in xmessy_mmd
### 3) nml/switch.nml  
###      USE_OROGW=.TRUE. in nml/switch.nml
###      USE_E5VDIFF=.TRUE. in nml/switch.nml
### 4) cp orogw.nml and e5vdiff.nml to nml/. 
### 5) adjust nml/tendency.nml: 
###     gwspect -> gwave
###     vdiff   -> e5vdiff
###     ssodrag -> orogw
### 6) nml/channel.nml (see E5/channel.nml):
###    - set restart-ignore-flag (=T) for e5vdiff, orogw, mlocean:
###      OUT_CHANNEL(.) = 'e5vdiff', $OFT, $OFT, -1, F,T, F,T,F,F,F, F,F, , ,
###      OUT_CHANNEL(.) = 'orogw',   $OFT, $OFT, -1, F,T, F,T,F,F,F, F,F, , ,
###      OUT_CHANNEL(.) = 'mlocan',  $OFT, $OFT, -1, F,T, F,T,F,F,F, F,F, , ,
###    - update for output of e5vdiff, orogw, and mlocean objects
###      (see example list at end of this script)
### 7) cp new executable to $WORKDIR/bin/.
### 8) run this script in $WORKDIR
### 9) submit job
################################################################

# ncdump -h restart_e5vdiff.nc | grep double | awk '{print $2}' | awk -F '(' '{print $1}' | grep -v _x1

### (A) FIND OUT LINK TARGET AND SET FILE NAMES

set target = `ls -la restart_g3b.nc | awk '{print $NF}'`
set tdir   = `dirname $target`
set tfile  = `basename $target`
set tcycle = `echo $tfile | awk -F '_' '{print $2}'`

### SOURCE
set gfile = restart_${tcycle}_g3b.nc
set efile = restart_${tcycle}_ECHAM5.nc

### DESTINATION
set vfile = restart_${tcycle}_e5vdiff.nc
set ofile = restart_${tcycle}_orogw.nc
set mfile = restart_${tcycle}_mlocean.nc
set nfile = restart_${tcycle}_nudg.nc

echo $gfile $efile ' -> '$vfile $ofile $mfile '(update units in '${nfile}')'

set here = `pwd`
cd $tdir

set list = (\
lon \
lat \
lev \
hyam \
hybm \
ahfl \
ahfli \
ahflw \
ahfll \
ahfs \
ahfsi \
ahfsw \
az0hi \
az0hw \
az0hl \
temp2 \
dew2 \
wet_tmp \
ustri \
vstri \
ustrw \
vstrw \
ustrl \
vstrl \
vdis \
ustr \
vstr \
evap \
evapi \
evapw \
evapot_2d  \
evapl_2d \
wind10 \
wind10w \
tke \
tkem \
tkem1 \
az0i \
az0w \
az0l \
t2min \
t2max \
wimax \
)

# temp2_min  <- g3b, t2min
# temp2_max  <- g3b, t2max
# wind10_max <- g3b, wimax
# ahfsl      <- g3b, ahfslac (?)

#set echo verbose

set g3b_list = ()
set e5_list = ()

foreach v ($list)
 set fg3b = `ncdump -h $gfile | grep double | awk '{print $2}' | awk -F '(' '{print $1}' | grep -x $v`

 if ($fg3b == $v) then
    set g3b_list = ($g3b_list $v)
 else
    set fe5 = `ncdump -h $efile | grep double | awk '{print $2}' | awk -F '(' '{print $1}' | grep -x $v`
    if ($fe5 == $v) then
       set e5_list = ($e5_list $v)
    else
       echo $v not found
    endif
 endif
end

# ---------------------
# E5VDIFF

set vstr_g3b = `echo $g3b_list | tr ' ' ','`
set vstr_e5  = `echo $e5_list  | tr ' ' ','`

ncks -O -v ${vstr_g3b} $gfile  $vfile
ncks -A -v ${vstr_e5}  $efile  $vfile

echo rename t2min to temp2_min
ncrename -v t2min,temp2_min $vfile

echo rename t2max to temp2_max
ncrename -v t2max,temp2_max $vfile

echo rename wimax to wind10_max
ncrename -v wimax,wind10_max $vfile

### modify attributes ...
ncatted -a channel_name,global,o,c,"e5vdiff" $vfile

# ---------------------
# OROGW

ncks -O -v ustrgw,vstrgw,vdisgw  $gfile  $ofile
ncks -A -v gworo_du,gworo_dv,gworo_dt,gwlif_du,gwlif_dv,gwlif_dt $efile $ofile
ncatted -a channel_name,global,o,c,"orogw" $ofile

# ---------------------
# MLCOEAN

ncks -O -v amlcorr  $gfile  $mfile
ncatted -a channel_name,global,o,c,"mlocean" $mfile

# ---------------------
# NUDG

if (-r $nfile) then
   ncatted -a units,NAPSFC,o,c,"ln(Pa)" $nfile
   ncatted -a units,NATEMP,o,c,"K" $nfile
   ncatted -a units,NADIV,o,c,"1/s" $nfile
   ncatted -a units,NAVOR,o,c,"1/s" $nfile
endif

# ---------------------

cd $here

ln -s $tdir/$vfile restart_e5vdiff.nc
ln -s $tdir/$ofile restart_orogw.nc
ln -s $tdir/$mfile restart_mlocean.nc
## if present, the nudging restart file ist aleready set by init_restart

exit 0

# -----------------------------------------------------------------------------
# LIST OF RECOMMENDED MODIFICATIONS IN channel.nml
# !!! Numbers in parentheses are arbitrary !!!
# -----------------------------------------------------------------------------
# < OLD
# > NEW
# -----------------------------------------------------------------------------

# > ADD_REF(  )        = 'g3b',     'az0',      'e5vdiff',   '',

# > OUT_CHANNEL(   )   = 'e5vdiff',  $OFT, $OFT, -1, F,T, F,T,F,F,F, F,F, , ,
# > OUT_CHANNEL(   )   = 'orogw',    $OFT, $OFT, -1, F,F, T,F,F,F,F, F,F, , ,
# > OUT_CHANNEL(   )   = 'mlocean',  $OFT, $OFT, -1, F,F, F,T,F,F,F, F,F, , ,


# < OUT_OBJECT(  5)  = 'g3b','wind10w',             F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT(  5)  = 'e5vdiff','wind10w',         F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT(  9)  = 'g3b','az0',                 F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT(  9)  = 'g3b','az0',                 F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 13)  = 'g3b','ustrgw',              F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 14)  = 'g3b','vstrgw',              F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 13)  = 'orogw','ustrgw',            F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 14)  = 'orogw','vstrgw',            F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 22)  = 'g3b','tke',                 F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 23)  = 'g3b','tkem1',               F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 22)  = 'e5vdiff','tke',             F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 23)  = 'e5vdiff','tkem1',           F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 30)  = 'g3b','emter',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 31)  = 'g3b','trsol',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 32)  = 'g3b','emtef',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 33)  = 'g3b','trsof',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 34)  = 'g3b','tkem',                F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 30)  = 'g3b','emter',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 31)  = 'g3b','trsol',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 32)  = 'g3b','emtef',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 33)  = 'g3b','trsof',               F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 34)  = 'e5vdiff','tkem',            F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 42)  = 'g3b','ahfli',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 43)  = 'g3b','ahflw',               F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 44)  = 'g3b','ahfll',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 42)  = 'g3b','ahfli',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 43)  = 'g3b','ahflw',               F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 44)  = 'g3b','ahfll',               F,F, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 46)  = 'g3b','amlcorr',             F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 47)  = 'g3b','amlheatac',           F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 46)  = 'mlocean','amlcorr',         F,F, F,T,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 47)  = 'g3b','amlheatac',           F,T, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 47)  = 'g3b','amlheat',             F,T, F,F,F,F,F, F,F, , ,

# < OUT_OBJECT( 70)  = 'g3b','emtef0',              F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT( 71)  = 'g3b','trsof0',              F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 70)  = 'g3b','emtef0',              F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT( 71)  = 'g3b','trsof0',              F,F, F,F,F,F,F, F,F, , ,

# > OUT_OBJECT( 74)  = 'e5vdiff','ahfs',            F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 75)  = 'e5vdiff','ahfsl',           F,T, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 76)  = 'e5vdiff','ahfsw',           F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 77)  = 'e5vdiff','ahfsi',           F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 78)  = 'e5vdiff','ahfl',            F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 79)  = 'e5vdiff','ahfll',           F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 80)  = 'e5vdiff','ahflw',           F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 81)  = 'e5vdiff','ahfli',           F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 82)  = 'e5vdiff','az0',             F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 83)  = 'e5vdiff','az0l',            F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 84)  = 'e5vdiff','az0w',            F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 85)  = 'e5vdiff','az0i',            F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 86)  = 'e5vdiff','az0h',            F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 87)  = 'e5vdiff','az0hi',           F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 88)  = 'e5vdiff','az0hw',           F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 89)  = 'e5vdiff','az0hl',           F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 90)  = 'e5vdiff','temp2',           F,F, T,T,F,T,T, F,F, , ,
# > OUT_OBJECT( 91)  = 'e5vdiff','dew2',            F,F, T,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 92)  = 'e5vdiff','wet_tmp',         F,F, F,F,F,F,F, F,F, , ,
# > OUT_OBJECT( 93)  = 'e5vdiff','wind10',          F,F, F,T,F,F,T, F,F, , ,
# > OUT_OBJECT( 94)  = 'e5vdiff','evapl_2d',        F,F, F,T,F,F,F, F,F, , ,
# > OUT_OBJECT( 95)  = 'e5vdiff','evapot_2d',       F,F, F,T,F,F,F, F,F, , ,

# < OUT_OBJECT(111)  = 'ECHAM5','az0w',             F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT(112)  = 'ECHAM5','az0i',             F,F, F,F,F,F,F, F,F, , ,
# < OUT_OBJECT(113)  = 'ECHAM5','az0l',             F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT(111)  = 'e5vdiff','az0w',            F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT(112)  = 'e5vdiff','az0i',            F,F, F,F,F,F,F, F,F, , ,
# > !!$OUT_OBJECT(113)  = 's5vdiff','az0l',            F,F, F,F,F,F,F, F,F, , ,

# -----------------------------------------------------------------------------
