m4_include([m4/acx_prog_fc_v_output.m4])
m4_include([m4/acx_fc_library_ldflags.m4])
m4_if(m4_version_compare(m4_PACKAGE_VERSION,[2.62]),[-1],dnl
  [m4_include([m4/ac_path_progs_feature_check.m4])])
m4_include([m4/acx_lang_c_check_include.m4])
