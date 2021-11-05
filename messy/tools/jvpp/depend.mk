# DO NOT DELETE THIS LINE - used by make depend
jvpp.o: jvpp_mem.o jvpp_step1.o jvpp_step2.o jvpp_step3.o
jvpp_mem.o: messy_cmn_photol_mem.o messy_main_constants_mem.o
jvpp_step1.o: jvpp_mem.o messy_cmn_photol_mem.o messy_main_constants_mem.o
jvpp_step1.o: messy_main_math_spline.o
jvpp_step2.o: jvpp_mem.o messy_cmn_photol_mem.o messy_main_constants_mem.o
jvpp_step3.o: jvpp_mem.o messy_cmn_photol_mem.o messy_main_constants_mem.o
jvpp_step3.o: messy_main_math_lsq.o
messy_main_math_lsq.o: messy_main_constants_mem.o
messy_main_math_spline.o: messy_main_constants_mem.o
