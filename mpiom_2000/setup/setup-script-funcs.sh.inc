assert_vars_are_set() {
  for var in "$@" ; do
    if eval test x\$\{${var}+set\} != xset ; then
      echo $var is not set!
      exit 1
    fi
  done
}

