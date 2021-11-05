#! /bin/sh -e

sed -n '/^ *F90DEFSMPIOM0/,/^$/ p ' $1 | tr '\n' ' ' | sed 's|\\| |g' | awk '{for(i=3;i<=NF;i++) printf(" %s",$i)}'

exit 0
