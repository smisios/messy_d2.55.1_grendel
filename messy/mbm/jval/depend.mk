# DO NOT DELETE THIS LINE - used by make depend
jval.o: messy_cmn_photol_mem.o messy_jval.o messy_main_constants_mem.o
jval.o: mo_netcdf.o
messy_jval.o: messy_jval_jvpp.inc
messy_jval.o: messy_cmn_photol_mem.o messy_main_constants_mem.o
messy_jval.o: messy_main_tools.o
messy_main_blather.o: messy_main_constants_mem.o messy_main_tools.o
messy_main_tools.o: messy_main_constants_mem.o
mo_netcdf.o: messy_main_constants_mem.o
jval_column.mod: jval.o
