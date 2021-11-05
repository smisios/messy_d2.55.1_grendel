#! /bin/sh

dir=`pwd`

echo "checking SMCL:"
cd messy/mbm/caaba/mecca/smcl
kppfiles=`find . -name "*.f90" -print`
cd ../../../../smcl
for kppfile in ${kppfiles}
do
  printf "  %-30s" ${kppfile}
  if test ! -L ${kppfile} ; then
    echo "...re-linked"
    rm -f ${kppfile}
    ln -fs ../mbm/caaba/mecca/smcl/${kppfile} .
  else
    echo "...OK"
  fi
done

cd $dir

echo "checking SMIL:"
cd messy/smil

incfiles="messy_mecca*_c2mr_si.inc messy_mecca*_idt_si.inc messy_mecca*_mr2c_si.inc messy_mecca*_trac_si.inc"

for incfile in ${incfiles}
do
  printf "  %-30s" ${incfile}
  if test ! -L ${incfile} ; then
    echo "...re-linked"
    rm -f ${incfile}
    ln -fs ../mbm/caaba/mecca/${incfile} ${incfile}
  else
    echo "...OK"
  fi
done

cd $dir

exit 0
