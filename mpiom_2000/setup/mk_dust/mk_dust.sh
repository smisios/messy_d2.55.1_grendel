#! /usr/bin/env bash
set +x

GRID="${GRID-TP01}"
POOL="${POOL-/pool/data/MPIOM}"
CDO=${CDO-cdo}

case "$GRID" in
  TP10)
    ie=362
    je=192
    tp='.true.'
    ;;
  TP04)
    ie=802
    je=404
    tp='.true.'
    ;;
  TP01)
    ie=3586
    je=1800
    tp='.true.'
    ;;
  *)
    echo "Grid configuration by this name ($GRID) not supported!" >&2
    exit 1
    ;;
esac

ie2=$(($ie - 2))

# dust deposition is taken from ECHAM5/HAM model run with the Tegen dust scheme
# reference: Stier et al. 2004
# FIXME: this is not currently[1] part of the setup dir at HH
# [1]: current being 2009/04/02, Thomas Jahns
SETUP="${POOL}/setup/SILVIA_ECHAM"


mkdir "$GRID"
cd "$GRID"
#better use remacon, but the n halo needs to be created new (soon available with cdo)
#${CDO} writegrid -selindexbox,1,${ie2},3,${je} -random,/pool/data/MPIOM/${GRID}/${GRID}s.nc ${GRID}s3.nc
#${CDO} remapcon,${GRID}s3.nc -setgrid,t63grid ${SETUP}/dustdep_monthl.nc dustdep_${GRID}.nc

#for the timebeing use bilinear interpolation
"${CDO}" "remapbil,../../anta2nc/${GRID}s.nc" \
  -setgrid,t63grid "${SETUP}/dustdep_monthl.nc" "dustdep_${GRID}.nc"

cat > INPDUST.partab<<EOF
&PARAMETER
  CODE=-1
  NAME=DUST
  LONG_NAME="dust"
  UNITS="kg m-2 yr-1"
/
EOF

y=$(( 365*24*60*60 ))

cdo mulc,$y "dustdep_${GRID}.nc" tt
cdo setpartab,INPDUST.partab tt tt1
cdo sethalo,1,1 tt1 "${GRID}_INPDUST.nc"
cd ..
