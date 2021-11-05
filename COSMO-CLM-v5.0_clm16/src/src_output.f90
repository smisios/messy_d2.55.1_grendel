!+ Source Module for writing Grib files
!------------------------------------------------------------------------------

MODULE src_output

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines necessary for writing the result data
!   of the LM. It uses also routines from the module "io_utilities".
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.2        1998/03/30 Ulrich Schaettler
!  Regarding whether digital filtering has been performed.
! 1.7        1998/07/16 Guenther Doms
!  Use of additional global fields to perform output of non-global fields.
! 1.8        1998/08/03 Ulrich Schaettler
!  Use grib parameters from module data_io.f90.
! 1.9        1998/09/16 Guenther Doms
!  Use of parameters 'nincmxt' and 'nincmxu' (replacing 'nincmxn') from
!  data module 'data_runcontrol.f90'.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminate dependency from routine remark.
! 1.14       1998/10/26 Ulrich Schaettler
!  Changed igds to igds_out.
! 1.17       1998/11/17 Ulrich Schaettler
!  Changes for reading and writing ready files and constant variables.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.24       1999/03/01 Guenther Doms
!  Inclusion of the new prognostic 3-D array 'qi'.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO;
!  Subroutine tautsp has been put to module utilities.f90
! 1.30       1999/06/24 Matthias Raschendofer
!  Use 6 additional fields: t, t_g, qv_s, tfm, tfh, tke, gz0.
!  Use 1 additional parameter form module data_runcontrol: ntke
! 1.32       1999/08/24 Guenther Doms
!  Use of postprocessing utilities from new routine 'pp_utilities'
! 1.33       1999/10/14 Guenther Doms
!  Use of postprocessing utility 'caliq' from 'pp_utilities'.
! 1.34       1999/12/10 Ulrich Schaettler
!  Use new Namelist variables and include all module procedures in this file
! 1.38       2000/04/06 Christoph Schraff
!  Correction for mean values when writing analyses (in 'output_grib_data').
! 1.39       2000/05/03 Ulrich Schaettler 
!  Changed names for variables concerned to latitude or longitude.
!  Introduced possibility for database IO.
!  Prepared output for different nests (for interactive nesting).
! 1.40       2000/05/23 Ulrich Schaettler
!  No interpolation to masspoints for u,v for p- and z-interpolation.
! 2.8        2001/07/06 Ulrich Schaettler
!  Eliminated non-necessary variables from the USE-lists;
!  Adapted input to new organization of I/O;
!  Introduced 2D version of p- and z-interpolation
! 2.10       2001/07/24 Ulrich Schaettler
!  Corrected declaration of my_iee in subroutine output_grib_data
! 2.11       2001/09/28 Ulrich Schaettler
!  Eliminated extensive use of LEN_TRIM (which was not performant on e.g. NEC)
!  Corrected a bug when gribbing the pole of the rotation grid (igds_out(21))
! 2.14       2002/02/15 Ulrich Schaettler
!  Modifications for computation of the height of the snow-fall limit
!  Correction in the grib-coding of the upper right corner (igds_out(11))
! 2.17       2002/05/08 Ulrich Schaettler
!  Modifications to perform I/O-communications in irealgrib-format
! 2.19       2002/10/24 Ulrich Schaettler
!  Corrected bugs in writing variables from multi-layer soil model and
!  in calculating the gds-values in case of luvmasspoints=.TRUE.
! 3.5        2003/09/02 Ulrich Schaettler
!  Include output for zenith delay (routine calztd from pp_utilities)
!  Include output for ICW, ICI (integrated cloud water; integrated cloud ice)
!  Corrected output for calrelhum, calomega
! 3.6        2003/12/11 Ulrich Schaettler
!  Adaptations for multi-layer soil model
! 3.7        2004/02/18 Ulrich Schaettler
!  Adaptations for output of synthetic satellite images as GRIB fields;
!  Possibility for specifying additional GRIB tables
!  Renamed phi to rlat
! 3.8        2004/03/23 Ulrich Schaettler
!  Bug-Fix in the treatment of the linked-list for Namelist group /GRIBOUT/
! 3.13       2004/12/03 Ulrich Schaettler
!  Put KIND-parameters for Grib-library to data_parameters;
!  Changed W_ICE to W_SO_ICE (Reinhold Schrodin)
!  Adaptations for output of new variables (graupel scheme, 3D turbulence)
!  Possibility to write 3dimensional variables on p- and z-levels
!                            (Thorsten Reinhardt, Jochen Foerstner)
!  Bug correction for writing gds for u and v in case of p- or z-levels
!                            (Emanuele Zala)
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL; Implemented NL variable nunit_of_time
! 3.16       2005/07/22 Ulrich Schaettler
!  Adapted length of pathname for output-directory
!  Calculate new fields for output of radar images
! 3.18       2006/03/03 Ulrich Schaettler
!  Introduction of writing NetCDF and Restart files
!  Changed treatment of ASCII files for introducing restart possibility
!  Added output variables TQR, TQS, TQG, RELHUM_2M
!  Changes to introduce new type of coding the vertical coordinate parameters
!    (introduced by MeteoSwiss for SLEVE vertical coordinate)
!  Introduction of ldwd_grib_use: if .TRUE., special non-standard Grib Code
!    settings used at DWD are done (setting of ipds(4), tri, reference time
!    for some analysis products)
!  Introduction of possibility to write a subdomain field
!  Determination of grib record length now by return value idims_out(19) from
!    routine grbex1 (instead of idims_out(18)*iwlength)
!  New routines for NetCDF: write_nc_gdefs, write_nc_vdefs
!  New routines smooth_pmsl, smooth_geopot for extreme smoothing pmsl, geopot
!    in mountainous terrain (introduced by Bundeswehr
! 3.19       2006/04/25 Ulrich Schaettler
!  Corrections in the NetCDF output
! 3.21       2006/12/04 Ulrich Schaettler
!  changes in write_nc_gdefs to meet netCDF CF standards (B. Rockel)
!  Correction for specifying soil types for NetCDF output
!  land/sea masks in netCDF output included
!  polgam introduced
!  Save the vertical coordinate parameters in restart files (U. Schaettler)
!  Use nnow for output
!  Adaptations for Ensemble Prediction output (C. Gebhardt)
!  Additional output of several convective indices (D. Leuenberger)
!  Time integrated analysis increment fields introduced (C. Schraff)
! V3_23        2007/03/30 Jochen Foerstner, Michael Baldauf, Ulrich Schaettler
!  Corrected computation of q_sedim (must be done with qi(:,:,:,nnew)
!  Added call to SR calc_ceiling
!  Introduced idbg_level for verbosity of output
! V3_24        2007/04/26 Ulrich Schaettler
!  Eliminated nincmxu, nincmxt and introduced control as for other increments
! V3_25        2007/05/21 Ulrich Schaettler
!  Corrections for writing synthetic satellite images for MSG2
!  Modifications for writing lat/lon values to NetCDF files
! V3_26        2007/05/29 Ulrich Schaettler
!  More debug level output
! V3_27        2007/06/29 Ulrich Schaettler
!  Additional correction for flexible dt in makepds for writing analysis
! V4_1         2007/12/04 Ulrich Schaettler
!  Corrected settings of igds_out for ivctype=3
!  Bug fix for re-initializing rain rates after restart files (Uwe Boehm)
!  Introduced output for SDI (supercell detection indices)
!  Introduced additional clipping of variables, if values are around 0
!  (only for grib-output)
! V4_3         2008/02/25 Ulrich Schaettler
!  Omit grib warnings in case of climate simulations
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
!  Adapted interface of get_timings
! V4_5         2008/09/10 Guenther Zaengl
!  Adaptations for new reference atmosphere
! V4_8         2009/02/16 Ulrich Schaettler, Guenther Zaengl
!  Corrections for NetCDF Input and Grib Output
!  Use noutlevels as maximum number of output levels
!  Adapted interface to SR cal_conv_ind to changes in pp_utilities
!  Adapted GDS encoding and restarts for new reference atmosphere
!  Use p0hl (reference pressure at half levels) for full consistency with
!  new reference atmosphere implementation
!  Bug fix for grib encoding of LM output (affects only sleve coordinate and
!  new reference atmosphere)
!  Add l_ke_in_gds to partly replace ldwd_grib_use
! V4_9         2009/07/16 Ulrich Schaettler, Burkhardt Rockel
!  Corrections in closing YUCHKDAT
!  Vectorization of p_int and z_int routines
!  Adaptations for new reference atmosphere in netCDF output
!  Include ivctype=2 and ivctype=3 in netCDF output
! V4_10        2009/09/11 MCH
!  Computation and output of BRN and HPBL
!  Added compiler directive to use option _on_adb for NEC
! V4_11        2009/11/30 Guenther Zaengl, Lucio Torrisi
!  Adaptations for output using irefatm=3 (const. Brunt-Vaisala frequency)
!  Adaptations for output of additional fluxes
! V4_12        2010/05/11 Michael Baldauf, Ulrich Schaettler
!  Introduced output for potential and (relative) Vorticity
!  Moved subroutine calc_sdi from pp_utilities to here because of formal reasons
!  Renamed t0 to t0_melt because of conflicting names
!  Added more convective indices fields for output; 
!  Introduced unit_of_time=15 for 10 minutes output (Oli Fuhrer)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Ulrich Schaettler, Oliver Fuhrer
!  Eliminated tgrlat, acrlat from interface to SR potential_vorticity_rho
!  Use of nnow instead of itl (for itimelevel) in src_output for some variables
!   VORTIC_U,_V,_W, POT_VORTIC,  T_WATER, and SR caliq, calztd, calomega
! V4_17        2011/02/24 Ulrich Blahak
!  - Adapted interface of exchg_boundaries; 
!  - corrected kzdims(1:20) -> kzdims(1:24);
!  - eliminated my_peri_neigh; 
!  - implemented linear interpolation as alternative to
!    cubic tension splines in SRs z_int() and pint() -- this uses 
!    the new interpolation routine lininterp2D_xinter1D_vec() 
!    from utilities.f90 and may be activated by the new namelist
!    parameter itype_vertint in namelist group gribout (1=cubic spline, 2=linear);
!  - surface values of U, V, W, and T
!    for *linear* z_int()-interpolation now depend on lnosurffluxes_m/h (free-slip b.c or not)
!    FOR NOW: THIS IS NOT DONE FOR P-LEVELS OUTPUT AND Z-LEVELS CUBIC 
!    SPLINES TO PRESERVE "OLD" BEHAVIOUR AND TO AVOID PROBLEMS
!    WITH OVERSHOOTS IN CUBIC SPLINE INTERPOLATION.
!  - surface value of pressure for z_int()-interpolation is now
!    taken to be the surface pressure ps instead of p on the
!    lowest main levels. This changes slightly the results
!    of z-interpolated pressure near the surface.
!  - changed variable name "result"
!    to "results", since "result" is a fortran key word; 
!  - added debug output on processed variables for constant fields output.
! V4_18        2011/05/26 Ulrich Schaettler
!  Introduced conditional compilation for synthetic satellite images
!  Moved NL parameter yform_write to the group(s) /GRIBOUT/ to be able to
!     specify the output format differently for every group. 
!  Adapted NetCDF I/O to deal with 3D external parameter field for sectors of
!   the horizon (for topographical corrections) and its attributes (Anne Roches)
!  Adapted NetCDF I/O to deal with synthetic satellite data (but only MSG)
!  Introduced 4 additional fields for each group of products for syn sat data
!   (Anne Roches et al.)
!  More general if-clauses for SR caliq (Michael Baldauf)
!  Bug fixes in calls to SR calc_Theta_Tppp and potential_vorticity_rho
!   (reported by Jean-Marie Bettems)
!  Exchange of boundaries for AUMFL_S, AVMFL_S, if luvmasspoint=.TRUE.
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for NetCDF and GRIBOUT
!  Check inconsistent RTTOV- and OUTPUT-settings for synthetic satellite images
!     (Robin Faulwetter)
! V4_20        2011/08/31 Ulrich Schaettler
!  Replaced variablename namelist (which is Fortran Keyword) by outblock
!  Bug fix for computing total pressure on highest half level (J. Schmidli)
! V4_21        2011/12/06 Ulrich Blahak
!  Bugfixes p_int for the case itype_vertint=2 (linear vertical interpolation):
!   - monotonically increasing (dummy) p-values are also required below the
!     surface for routine lininterp2D_xinter1D_vec(). 
!   - error in field dimension when calling lininterp2D_xinter1D_vec
!     (-> model crashes) was corrected.
!  Bug in calling sequences of radar_lm_ray: The routine requires
!   hydrometeor densities (kg/m**3), not specific values (kg/kg), so added 
!   necessary multiplications with rho.
!  Initialized variable izerror in SR calc_sdi (Oli Fuhrer)
! V4_23        2012/05/10 Ulrich Schaettler, Oliver Fuhrer, H.-J. Panitz, Ulrich Blahak
!  Use field sqrtg_r_s, dzeta_* from new module grid_metrics_utilities
!  Added computation of total precipitation rate TOT_PR (Oli Fuhrer)
!  The subroutine makegds is used by restart files and must not be embraced by
!    ifdef GRIBDWD
!  CLM:
!  Added support for climatological year with 365 days
!  Introduction of new diagnostic variable for maximum wind speed in 10m height
!  Introduced new field snow_melt
!  Correction of bug (?) related to definition of subregrions
!   for netcdf output (Namelist parameter ydomain = 's')
!  Only for climate mode: reset all necessary precipitation "components" like RAIN_GSP etc.
!   to zero in case that only TOT_PREC is an output variable
!  Comment a few lines in write_nc_gdefs that are not neede any more
!  Changes to allow multi-layer snow model quantities in netCDF output
!  New flag "i" to distinguish between ocean and inland water (lakes) quantities 
!     (for netCDF output only)
!  Write global attributes in netCDF only if they are defined
!  UB:
!  Changed name of l_fi_ps_smooth to l_fi_pmsl_smooth and added l_fi_filter and 
!   l_pmsl_filter in order to be able to independently smooth FI and PMSL
!   with a digital FIR filter, as for all other fields with l_z_filter / l_p_filter.
!   Consequently, eliminated dependency of PMSL-smoothing from l_z_filter.
! V4_24        2012/06/22 Ulrich Schaettler, Hendrik Reich
!  Conditional compilation for GRIBDWD in SR makegds:
!   The vertical coordinate parameters can only be written, if GRIBDWD is set
!   and the Grib library is available, because of packing of REALs to INTEGERs.
!   In case of restart these parameters are not written, but also not needed
!  Adapted length of strings for date variables (HR)
!  Introduced new argument to SR make_fn (lmmss)
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer, Ulrich Blahak
!                         Florian Prill, Hans-Juergen Panitz, Carlos Osuna
!  Replaced qx-variables by using them from the tracer module
!  Implemented second type of gathering 2D fields from all PEs to the one PE that
!   does the (Grib) packing (namelist switch itype_gather) (OF)
!  UB:
!  Introduced output of variables for the 2-moment microphysics (QHAIL, QNXXX,
!   radar reflectivity by subroutine radar_sb_ray()).
!  For z- and p-level interpolated DBZ and hydrometeor quantities (QX, QNXXX),
!   hardwired the linear interpolation instead of cubic splines to
!   prevent undershoots for these highly variable quantities.
!  Added new namelist switch "outblock%loutput_q_densities"
!   for the gribout namelist(s). If set to .true., hydrometeor variables qx, qnx
!   are output in units of kg/m**3 resp. 1/m**3 instead of kg/kg resp. 1/kg.
!   (default: .FALSE. = traditional output).
!  Bugfixes: many changes to make all output fields really consistent
!   to the chosen timelevel of output. Except for restart, this is
!   "nnow" instead of "nnew" since V4.15! Changes are related mainly
!   to rho, qrs, which are now correctly diagnosed before output on
!   both timelevels "itimelevel" and "nnew" and stored in local variables.
!  Adapted call to make_fn according to changes in io_utilities (Uli Schaettler)
!  Adapted computation of reference time, if dt does not fit into the output interval
!  Corrected writing of ready-files at the end of an output step (FP)
!  Adapted interface to write_grib (FP)
!  Introduced nexch_tag for MPI boundary exchange tag to replace ntstep (HJP)
!  In case netcdf async I/O mode is selected, do not write data into netcdf
!   file but send data to I/O PE instead. (CO)
!  Add arguments to output_data subroutine needed to send metadata to asyn I/O PE
!   Move write_nc_gdefs and write_nc_vdefs to netcdf_io.f90
! V4_26        2012/12/06 Ulrich Blahak, Anne Roches, Ulrich Schaettler
!                         Hans-Juergen Panitz, Bojan Skerlak
!  Changed "ytrans_out /= '    ' " to "LEN_TRIM(ytrans_out) > 0" (UB)
!  The new diagnosis of rho and qrs must not be done in case of restart files, 
!    because this can not be re-constructed after a restart
!  Bug fix for itype_gather=2 (AR)
!  Bug fix for calling SR potential_vorticity_rho: the prognostic fields u, v 
!   and w have to be passed with the correct timelevel itl (BS)
!  Write time level ntke for TKE scheme to binary restart file for correct
!   TKE restart. (US)
!  Always construct file name for constant fields with step 0, also for restarts (US)
!  In case of asynchronous I/O print the name of the 'c' file correctly
!   in file YUCHKDAT (HJP)
! V4_27        2013/03/19 Ulrich Schaettler, Astrid Kerkweg
!  Use nmsgchan from data_satellites
!  MESSy interface introduced: get diagnostic output  (AK)
! V4_28        2013/07/12 Ulrich Schaettler
!  Implemented grib_api for writing GRIB(1/2) data
!  Use and set vertical coordinate parameters for output
!  Moved subroutines for setting I/O meta data to new module io_metadata.f90
!  Removed Nest-handling from this subroutine
!  Use subroutines and variables for vertical grid and reference atmospheres
!    from module vgrid_refatm_utils
!  Adapted interface to grib_api routines with special grib_api integer
! V4_29        2013-10-02 Ulrich Schaettler, Astrid Kerkweg (Messy), Ulrich Blahak
!  Corrected placement of ifdef GRIBAPI for gribinit_loop
!  Unification of MESSy interfaces and COSMO Tracer structure
!  For the COSMO-Model only use vcoord and refatm from vgrid_refatm_utils
!  Added check for upper bound of zlev resp. plev in z_int() and p_int() (UB)
! V4_30        2013-11-08 Ulrich Schaettler
!  Renamed ipds to ipds_out to reflect usage for output
! V5_00_clm5   2015-08-24 Ulrich Schaettler, DWD; H.-J. Panitz (HJP), KIT/IMK
!  re-initialization of field ds_real (end of SR output_data)
!   necessary to avoid model crashes in case of usage of large number of processors
! V5_0_clm10   2017-12-20 Ulrich Schaettler, DWD / H.-J. Panitz, IMK/KIT
!  added new Namelist parameter lzint_above_ground to choose z-interpolation to
!  height levels above ground
! V5.0_clm11 2018-10-11 Katherine Oesterried ETHZ / Burkhardt Rockel HZG, netcdf restart files
! V5.0_clm12,  2017-07-14, Ronny Petrik, HZG
!  adapt the method spline2D for handling interpolation to height levels which vary
!   in west-east-direction (i-direction)
!   add the definition of total density as 'RHO_TOT' to output variables
!   treat the wind-sector related statistics 
!   allow rho_tot to be interpolated to z-levels and p-levels
!   allow the TKE to be correctly treated for z-level and p-level outputs
!    (i.e. recalculation of TKE in case of standard turbulence parameterizations)
!   moved the reset of aggregated variables to organize_data (to prevent quickly
!     deletion of variables prior secondary usage)
!   consider the output for 3D variables which use a special third dimension
!     (wind sectors or gust levels)
!   manage to calculate vabs_10m_secav correctly considering wind rose
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Modules used:
#ifdef GRIBAPI
USE grib_api
#endif

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib, & ! KIND-type parameter for real variables in the grib library
    iwlength,  & ! length of an integer word in byte
    intgribf,  & ! KIND-type parameter for fortran files in the grib library
    intgribc,  & ! KIND-type parameter for C files in the grib library
    int_ga       ! integer precision for grib_api: length of message in bytes

USE data_fields, ONLY :  &
    hhl,          & ! geometrical height of model half levels
    hsurf,        & ! height of surface topography
    llandmask,    & ! landpoint mask
    fr_lake,      & ! lake fraction in a grid element [0,1]         (  -  )
    rlat,         & ! geographical latitude
    rain_gsp,     & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp,     & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    grau_gsp,     & ! amount of graupel from grid-scale prec. (sum) (kg/m2)
    hail_gsp,     & ! amount of hail    from grid-scale prec. (sum) (kg/m2)
    rain_con,     & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con,     & ! amount of snow from convective precip. (sum)  (kg/m2)
    snow_melt,    & ! amount of snow melt (sum)                     (kg/m2)
    prr_gsp,      & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp,      & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp,      & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prh_gsp,      & ! precipitation rate of hail,    grid-scale     (kg/m2*s)
    prr_con,      & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con,      & ! precipitation rate of snow, convective        (kg/m2*s)
    pp,           & ! deviation from the reference pressure 
    ps,           & ! surface pressure 
    dp0,          & ! pressure thickness of model layers
    rho0,         & ! base-state density
    rho,          & ! total air density
    p0,           & ! base-state pressure of full model levels
    p0hl,         & ! base-state pressure of half model levels
    t0,           & ! base state temperature
    clc_sgs,      & ! subgrid-scale stratiform cloud cover  
    clc_con,      & ! convective cloud cover
    top_con,      & ! level index of convective cloud top 
    bas_con,      & ! level index of convective cloud base
    pptens          ! pressure tendency

USE data_fields, ONLY :  &
    qrs,          & ! specific precip. water content                (kg/kg)
    crlat,        & ! cosine of transformed latitude
    u,            & ! zonal velocity
    v,            & ! meridional velocity
    w ,           & ! vertical velocity
    t ,           & ! temperature
    t_g,          & ! weighted surface temperature                   (  k  )
    t_2m,         & ! temperature in 2m                             (  K  )
    qv_2m,        & ! specific water vapor content                  (kg/kg)
    w_snow,       & ! water content of snow                         (m H2O)
    p_anai,       & ! deviation from the reference pressure         ( Pa  )
    qv_anai,      & ! specific water vapor content                  (kg/kg)
    qc_anai,      & ! specific cloud water content (via saturation adjustm)
    synme7,       & !
    synmsg,       & !
    fc,           & ! coriolis-parameter                            ( 1/s )
    fccos,        & ! horizontal coriolis-parameter                 ( 1/s )
    acrlat,       & ! 1 / ( crlat * radius of the earth )           ( 1/m )
    tgrlat,       & ! tangens of transformed latitude                 --
    wdirgeo_10m_freq ! frequency distribution of wind directio (  )


USE data_modelconfig,ONLY : &
    czmls,        & ! depth of the main soil layers in m
    czhls,        & ! depth of the half soil layers in m
    msoilgrib,    & ! grib coded depth of main soil levels in centimeters
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    dt,           & ! long time-step
    istartpar,    & ! start index for computations in the parallel program
    jstartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program
    jendpar,      & ! end index for computations in the parallel program
    istart,       & ! start index for the forecast of w, t, qv, qc and pp
    iend,         & ! end index for the forecast of w, t, qv, qc and pp
    jstart,       & ! start index for the forecast of w, t, qv, qc and pp
    jend,         & ! end index for the forecast of w, t, qv, qc and pp
    ie,           & ! number of grid points in zonal direction
    ie_tot,       & ! number of grid points in zonal direction total
    je,           & ! number of grid points in meridional direction
    je_tot,       & ! number of grid points in meridional direction total
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors
    ke,           & ! number of grid points in vertical direction
    ke_tot,       & ! number of grid points in vertical direction total
    ke1,          & ! KE+1
    ke_soil,      & ! number of layers in the multi-layer soil model
    ke_snow,      & ! number of snow layers      !_br 23.01.12
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    klv850,       & !
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    idt_qv,  idt_qc,  idt_qs,   &
    idt_qr,  idt_qi,  idt_qg,   &
    idt_qh,  idt_qnc, idt_qnr,  &
    idt_qni, idt_qns, idt_qnh,  &
    idt_qng


USE data_constants, ONLY : &
    pi,           & ! circle constant
    r_d,          & ! gas constant for dry air
    r_earth,      & ! mean radius of the earth (m)
    g,            & ! gravity acceleration
    rdv,          & ! r_d / r_v
    rvd_m_o,      & ! r_v/r_d - 1
    o_m_rdv,      & ! 1 - r_d/r_v
    cpdr,         & ! 1.0 / cp_d
    t0_melt,      & ! melting temperature
    rho_w,        & ! density of liquid water
    rho_ice,      & ! density of ice          (kg/m^3)
    K_w,          & ! dielectric constant for water
    K_ice,        & ! dielectric constant for ice
    b1,           & !
    b2w,          & !
    b3,           & !
    b4w,          & !
    cp_d,         & ! specific heat capacity of dry air
    lh_v            ! latent heat of condensation

USE data_runcontrol, ONLY : &
    nlastmxu,     & ! last step when vbmax was "nullified"
    nlastmxt,     & ! last step when tmin, tmax were "nullified"
    nnew,         & ! corresponds to ntstep + 1
    nnow,         & ! corresponds to ntstep
    nstart,       & ! first time step of the forecast
    ntstep,       & ! actual time step
    ntke,         & ! actual TKE-time step, corresponds to ntstep
    nvers,        & ! version number of experiment for documentation
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! else with split-explicit scheme (only for l2tls=FALSE!)
    leps,         & ! switch ensemble mode on/off
    iepsmem,      & ! ID of ensemble member (EPS)
    iepstot,      & ! total number ensemble members (EPS)
    iepstyp,      & ! ID of ensemble generation type (EPS)
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_gscp,   & ! type of grid-scale precipitation physics
    lmetr,        & ! lartif_data=.TRUE.:  with metric terms
                    !            =.FALSE.: or without metric terms
    ldfi,         & ! whether digital filtering or not
    lmulti_layer, & ! run multi-layer soil model
    lmulti_snow,  & ! run multi-layer snow model
    nhori,        & ! number of sectors for the horizont array by the topographic
                    ! correction of the radiation
    lradtopo,     & ! if .TRUE., calculate topographic correction of radiation
    lcori_deep,   & ! if =.TRUE.: take cos(phi) coriolis terms into account
    itype_calendar,&! for specifying the calendar used
    psm0,         & ! initial value for mean surface pressure ps
    dsem0,        & ! initial value for mean dry static energy
    msem0,        & ! initial value for mean moist static energy
    kem0,         & ! initial value for mean kinetic energy
    qcm0,         & ! initial value for mean cloudwater content
    yakdat1,      & ! actual date (ydate_ini+ntstep/dt) in the form
                    ! ddmmyyyyhhmmss (day, month, year, hour, min, sec)
    luse_rttov      ! if .true. calculate synthetic satellite images

USE data_runcontrol, ONLY : &
    ltime,        & ! detailed timings of the program are given
    lroutine,     & ! if .TRUE., run an operational forecast
    idbg_level,   & ! to control the verbosity of debug output
    ldebug_io ,   & ! if .TRUE., debug output for I/O
    lprintdeb_all,& ! .TRUE.:  all tasks print debug output
                    ! .FALSE.: only task 0 prints debug output
    luse_rttov,   & ! if rttov-library is used
    lartif_data,  & ! forecast with self-defined artificial data
    lperi_x,      & ! if lartif_data=.TRUE.: periodic boundary conditions in x-dir.
                    ! or with Davies conditions (.FALSE.)
    lperi_y,      & ! if lartif_data=.TRUE.: periodic boundary conditions in y-dir.
                    ! or with Davies conditions (.FALSE.)
    l2dim,        & ! 2 dimensional runs
    nbl_exchg,    & !
    cur_outstep,    & ! current output time step
    cur_outstep_idx,& ! index of current output time step 
    cur_gribout_idx   ! index of current gribout section


USE data_parallel,      ONLY :  &
    lasync_io,      & ! if .TRUE.: the model runs with extra PEs for
                      ! asynchronous IO
    my_cart_id,     & ! rank of this subdomain in the cartesian communicator
    num_compute,    & ! number of compute PEs
    nc_asyn_io,     & ! number of asynchronous I/O PEs (netcdf)
    icomm_cart,     & ! communicator that belongs to the cartesian grid
    imp_reals,      & ! determines the correct REAL type used in the model
                      ! for MPI
    imp_grib,       & ! determines the REAL type for the GRIB library
    imp_integers,   & ! determines the correct INTEGER type used in the model
                      ! for MPI
    nboundlines,    & ! number of boundary lines of the domain for which
                      ! no forecast is computed = overlapping boundary
                      ! lines of the subdomains
    my_cart_neigh,  & ! neighbors of this subdomain in the cartesian grid
    iexch_req,      & ! stores the sends requests for the neighbor-exchange
                      ! that can be used by MPI_WAIT to identify the send
    ldatatypes,     & ! if .TRUE.: use MPI-Datatypes for some communications
    ltime_barrier,  & ! if .TRUE.: use additional barriers for determining the
                      ! load-imbalance
    ncomm_type,     & ! type of communication
    nexch_tag,      & ! tag to be used for MPI boundary exchange
                      !  (in calls to exchg_boundaries)
    sendbuf,        & ! sending buffer for boundary exchange:
                      ! 1-4 are used for sending, 5-8 are used for receiving
    isendbuflen,    & ! length of one column of sendbuf
    nproc, realbuf, intbuf

USE data_io,        ONLY : &
  clen,              & ! length of short name
  nhour_restart,     & ! start-, stop-, inc of writing restart files (tstep)
  nzmxid,            & ! maximum number of NetCDF variabe IDs
  ydate_ini,         & ! start of the forecast 
  ymode_write,       & ! mode for opening the (write) Grib files
  yform_read,        & ! format of the (read) files
  ntrans_out,        & ! Unit Number for writing ready-Files during output
  nuchkdat,          & ! Unit number for checking the I/O data
  yuchkdat,          & ! checking the I/O data
  ytrans_out,        & ! directory for writing ready files
  nsma_stat,         & ! status for soil moisture analysis
  npds,              & ! Dimension for product definition section (pds)
  ngds,              & ! Dimension for grid description section (gds)
  nbms,              & ! Dimension for bit map section (bms)
  nbds,              & ! Dimension for binary data section
  ndsup,             & ! Dimension for dsup
  ndims,             & ! Dimension for idims (contains all dimensions)
  lfd,               & !
  lfa,               & ! Dimension for grib_api message in bytes
  lbm,               & !
  lds,               & !
  inrvert_out,       & ! number of vertical coordinate parameters of output data
  ntrip,             & ! maximum number of timing triples
  noutlevels,        & ! maximum actual existing number of output levels
  nvar,              & ! maximum number of variables in LM variable table
  itype_gather,      & ! switch to determine gather method to use
  max_gribtabs,      & ! maximum number of GRIB tables in LM variable table
  idwdednr,          & ! grib edition number for DWD library
  undefgrib,         & ! value for "undefined" in the grib routines
  undefncdf,         & ! value for "undefined" in the netcdf routines
  undef,             & ! the same as undefgrib but with other KIND-Parameter
  lst_gribtabs,      & ! IDs of GRIB tables use
  nlocaldefnr,       & ! local definition number for GRIB local section
  nactlocdefnr,      & ! to overwrite Namelist parameter with some center default
  nprocess_ini_in,   & ! process gener. identification for initial (analysis)
  nprocess_bd_in,    & ! and for boundary (forecasts) data from input data
  lmmss                ! 10/14 digits date format

USE data_io,        ONLY : &
! Global arrays
  iblock,            & ! array for gribed data
  ymessage,          & ! array for grib2 message (in characters)
  idims_out,         & ! array for all dimensions
  ibmap,             & ! array for
  ipds_out,          & ! product definition section for output
  igds_out,          & ! grid description section
  ibms,              & ! bit map section
  ibds,              & ! binary data section
  dsup,              & ! Parameter for grib routines
  ds_grib,           & ! array for unpacked data
  ds_real,           & ! array for unpacked data
  igrib1_id,         & ! grib1 sample
  igrib2_id,         & ! grib1 sample

! Global types
  pp_nl,             & ! structure for gribout namelist
  var                  ! array for LM variable table

USE data_io,        ONLY : &
  l_ke_in_gds,       & ! explicit GDS entry for number of model levels
  l_ke_in_input,     & ! explicit GDS entry for number of model levels in input data
  lbdclim,           & ! boundary data in climate model     ! PIK  (D.Hauffe)
                       !  (in climate mode also some external parameters have
                       !  to be updated, which are held constant in forecast
                       !  mode; e.g. plant cover, root depth)
  idims_id_out,      & ! array for the IDs of the dimensions of netCDF 
                       ! formatted output
  yncglob_institution,   & ! originating center name
  yncglob_title,         & ! title string for the output
  yncglob_source,        & ! program name and version
  yncglob_project_id,    & ! identification of the project of simulation
  yncglob_experiment_id, & ! identification of the experiment of simulation
  yncglob_contact,       & ! contact e.g. email address
  yncglob_references,    & ! URL, report etc.
  ncglob_realization       ! number of the realization of the experiment


#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
USE data_satellites,          ONLY :  &
    sat_compute, num_sensors, nmsgchan
#endif

!------------------------------------------------------------------------------

USE utilities,                ONLY :  &
    smoother, dfilt4, dfilt8, tautsp, tautsp2D,                     &
    phirot2phi, rlarot2rla, get_utc_date, &
    lininterp2D_xinter1D_vec

USE pp_utilities,             ONLY :  &
    calpmsl,    calprsum,   caltopdc,   calhzero,   calsnowlmt,     &
    calcldepth, calclmod,   calomega,   calrelhum,  caliq,  calztd, &
    radar_lm_ray, cal_conv_ind, calc_ceiling,                       &
    calc_bulk_richardson, calc_pbl_brn,                             &
#ifdef TWOMOM_SB
    radar_sb_ray,                                                   &
#endif
    potential_vorticity_rho

USE grid_metrics_utilities,        ONLY :  &
    sqrtg_r_s,    & ! 1 / square root of G at scalar points       ( 1/m )
    dzeta_dlam,   & ! d zeta / d lambda (for constant phi,    z)
                    ! at the scalar position                      ( 1   )
    dzeta_dphi      ! d zeta / d phi    (for constant lambda, z)
                    ! at the scalar position                      ( 1   )

USE numeric_utilities,             ONLY :  &
    curl, calc_Theta_Tppp, mean_over_box, mean_cov_over_box, vert_avg

USE environment,              ONLY :  &
    model_abort, get_free_unit, release_unit, exchg_boundaries

USE meteo_utilities,  ONLY: calrho

USE parallel_utilities,       ONLY :  &
    gather_values, combine_subarrays, distribute_values, gather_field

USE io_utilities,             ONLY :  &
    open_file, close_file, write_grib, write_gribapi, write_netcdf,        &
    write_restart, make_fn, check_record, compute_grib_intbuffer_length

USE io_metadata,              ONLY :  &
    set_vcoord_refatm_out, make_grib_init, make_grib_grid,                 &
    make_grib_product, makegds

#ifdef NETCDF
USE netcdf_io,            ONLY :  &
    send_asyn_io, nc_orgmdata_length, nc_varmdata_length, &
    write_nc_gdefs,write_nc_vdefs
#endif

USE time_utilities,     ONLY: get_timings, i_computations_O, i_gather_data, &
                                           i_write_data, i_meta_data_w

USE vgrid_refatm_utils, ONLY: vcoord, refatm
USE src_artifdata,      ONLY: lnosurffluxes_m, lnosurffluxes_h

!------------------------------------------------------------------------------
    
USE src_tracer,       ONLY : trcr_get, trcr_errorstr

USE data_tracer,      ONLY : T_ERR_NOTFOUND

!------------------------------------------------------------------------------

#ifdef NETCDF
USE netcdf,           ONLY :   &
  nf90_def_dim,            &
  nf90_def_var,            &
  nf90_enddef,             &
  nf90_redef,              &
  nf90_put_att,            &
  nf90_put_var,            &
  nf90_noerr,              &
  nf90_strerror,           &
  NF90_CHAR,               &
  NF90_DOUBLE,             &
  NF90_FLOAT,              &
  NF90_GLOBAL,             &
  NF90_UNLIMITED
#endif

#ifdef MESSY
!MESSy/BMIL
USE messy_main_channel_bi,    ONLY: L_BM_ORIG_OUTPUT, L_FORCE_calcout
#endif

!==============================================================================

IMPLICIT NONE

!==============================================================================

! string variable to hold grid information
  CHARACTER (LEN=200)     grid_mapping

! for smoothing fi and pmsl, the global hsurf-field is needed
  REAL (KIND=ireals), ALLOCATABLE  :: hsurf_tot(:,:)

INTEGER (KIND=iintegers)             ::          &
  itimelevel, itl              ! for storing the output timelevel

! Some module variables
INTEGER (KIND=iintegers), SAVE  ::     &
  igribid                     ! ID of actual grib message

! pointers for tracers:
REAL (KIND=ireals), POINTER :: &
  qv      (:,:,:)=> NULL() ,&  ! QV at itl
  qc      (:,:,:)=> NULL() ,&  ! QC at itl
  qi      (:,:,:)=> NULL() ,&  ! QI at itl
  qg      (:,:,:)=> NULL() ,&  ! QG at itl
  qr      (:,:,:)=> NULL() ,&  ! QR at itl
  qs      (:,:,:)=> NULL()     ! QS at itl

#ifdef TWOMOM_SB
REAL (KIND=ireals), POINTER :: &
  qh      (:,:,:)=> NULL() ,&  ! QH at itl
  qnc     (:,:,:)=> NULL() ,&  ! QNCLOUD at itl
  qnr     (:,:,:)=> NULL() ,&  ! QNRAIN at itl
  qni     (:,:,:)=> NULL() ,&  ! QNICE at itl
  qns     (:,:,:)=> NULL() ,&  ! QNSNOW at itl
  qng     (:,:,:)=> NULL() ,&  ! QNGRAUPEL at itl
  qnh     (:,:,:)=> NULL()     ! QNHAIL at itl
#endif

REAL (KIND=ireals), ALLOCATABLE, DIMENSION(:,:,:), TARGET :: &
  zrho_itl,       & ! Total density at timelevel itimelevel
  zqrs_itl          ! QRS at timelevel itimelevel

PRIVATE :: qv, qc, qi, qg, qr, qs, zrho_itl, zqrs_itl

#ifdef TWOMOM_SB
PRIVATE :: qh, qnc, qnr, qni, qns, qng, qnh
#endif

#ifdef MESSY
  LOGICAL, SAVE :: l_COSMO_now = .FALSE. ! output required by COSMO
#endif

!==============================================================================
! Module procedures
!==============================================================================

CONTAINS

!==============================================================================
!+ Module procedure in src_output for initializing the output organization
!------------------------------------------------------------------------------

SUBROUTINE init_output (root)

!------------------------------------------------------------------------------
!
! Description:
!  The routine init_output initializes organizational variables of the 
!  output routines dealing with the grib code (such as the dimensions of
!  the different grib code sections). Also the grid description section
!  is initialized (except the location of the lower left grid point, because
!  it depends on the variable (U,V or other)).
!  
! Method:
!
!------------------------------------------------------------------------------

! Pointers with intent(in):
TYPE (pp_nl), POINTER                 :: root

!------------------------------------------------------------------------------
!
! Local scalars:
INTEGER (KIND=iintegers)     :: i,i1,i2,i3,k,n, niostat, izerrstat, itable,  &
                                nuedat, nzrecords, kbot, ktop, izerror,      &
                                iorg_data(3,0:num_compute-1), izdebug,       &
                                izpollat, izpollon, nstatus, modnvers

LOGICAL                      :: lzcheck, lzopen, lzwrite_ended, lzapi1_write,&
                                lzapi2_write

CHARACTER  (LEN=260)         :: yname
CHARACTER  (LEN=25)          :: yzroutine
CHARACTER  (LEN=80)          :: yzerrmsg
CHARACTER  (LEN= 3)          :: yzhead

! Local arrays:
REAL (KIND=ireals)      ::     &
  zvarlev   (ie,je,0:noutlevels) ! variable for vertical interpolation

REAL (KIND=ireals)      ::     &
  pollat_sp, pollon_sp

INTEGER (KIND=iintegers)::     &
  ndate, ntime, nsecond,       & !
  ivar_id(nzmxid)                ! NetCDF-ID of each variable in the list

! Local Pointers:
TYPE (pp_nl), POINTER   :: now

!
!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

yzroutine = 'init_output'
izerrstat = 0_iintegers

! Initialize, whether additional debug output shall be done
IF (ldebug_io) THEN
  IF (lprintdeb_all) THEN
    izdebug = idbg_level
  ELSE
    IF (my_cart_id == 0) THEN
      izdebug = idbg_level
    ELSE
      izdebug = 0
    ENDIF
  ENDIF
ELSE
  izdebug = 0
ENDIF

IF (izdebug > 5) THEN
  WRITE (*,'(A)') '  src_output: starting init_output'
ENDIF

! Set lfd, lds and lbm
lds = ie_tot * je_tot
lbm = 1875
lfd = compute_grib_intbuffer_length(ie_tot, je_tot, iwlength)
lfa = compute_grib_intbuffer_length(ie_tot, je_tot, 1) ! gives length in bytes

! Allocate GRIB arrays
ALLOCATE (iblock(lfd), ibmap(lbm),                 STAT=izerrstat)
ALLOCATE (ds_real(lds), ds_grib(lds), dsup(ndsup), STAT=izerrstat)
ALLOCATE (ymessage(lfa),                           STAT=izerrstat)

! Set vertical coordinate parameters for output
!  (this routine sets inrvert_out and pv_out)
CALL set_vcoord_refatm_out

! Initializations for the grib library
!  moving arraydimensions into idims
!  declaration dimensions
idims_out( 1) = npds
idims_out( 2) = ngds
idims_out( 3) = nbms
idims_out( 4) = nbds
idims_out( 5) = lbm
idims_out( 6) = ndsup
idims_out( 8) = lfd

!  real dimensions
idims_out(11) = 47
idims_out(12) = 25 + inrvert_out
idims_out(13) = 3
idims_out(14) = 5
idims_out(16) = 0

! idims_out(7,15,17) depend on the special namelist group and are set later

! Initialization for ivar_id (should be set also for non netcdf output)
ivar_id(:) = 0_iintegers

!------------------------------------------------------------------------------
! Section 2: Set the grid description section for GRIB1 (DWDLIB, Restart)
!------------------------------------------------------------------------------

CALL makegds

!------------------------------------------------------------------------------
! Section 3: Get grib samples and set constant meta data for grib_api
!------------------------------------------------------------------------------
  
now => root
lzapi1_write = .FALSE.
lzapi2_write = .FALSE.

gribinit_loop: DO

  !----------------------------------------------------------------------------
  ! Section 3.1: Get grib samples
  !----------------------------------------------------------------------------
  
#ifdef GRIBAPI
  IF (izdebug > 5) THEN
    WRITE (*,'(A)') '  src_output: reading grib_api samples for GRIBOUT: ', now%nl_index
  ENDIF

  IF (now%yform_write == 'api1') THEN
    IF (.NOT. lzapi1_write) THEN
      CALL grib_new_from_samples (igrib1_id, 'DWD_rotated_ll_7km_G_grib1', izerrstat)
      IF (izerrstat /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_new_from_sample: for api 1 ', izerrstat
      ENDIF
      lzapi1_write = .TRUE.   ! further api1 GRIBOUT blocks do not need to read sample
    ENDIF

    ! clone this sample to igrib1_sample 
    CALL grib_clone(igrib1_id, now%igribapi_id, izerrstat)
    IF (izerrstat /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone: from sample 1 ', izerrstat
    ENDIF
  ENDIF

  IF (now%yform_write == 'api2') THEN
    IF (.NOT. lzapi2_write) THEN
      CALL grib_new_from_samples (igrib2_id, 'DWD_rotated_ll_7km_G_grib2', izerrstat)
      IF (izerrstat /= GRIB_SUCCESS) THEN
        PRINT *,   ' *** Error in grib_new_from_sample: for api 2 ', izerrstat
      ENDIF
      lzapi2_write = .TRUE.   ! further api2 GRIBOUT blocks do not need to read sample
    ENDIF

    ! clone this sample to igrib1_sample 
    CALL grib_clone(igrib2_id, now%igribapi_id, izerrstat)
    IF (izerrstat /= GRIB_SUCCESS) THEN
      PRINT *,   ' *** Error in grib_clone: from sample 2 ', izerrstat
    ENDIF
  ENDIF
#endif

  !----------------------------------------------------------------------------
  ! Section 3.2: Set constant grib meta data
  !----------------------------------------------------------------------------
  
  IF (izdebug > 5) THEN
    WRITE (*,'(A)') '  src_output: setting constant GRIB meta data'
  ENDIF

  CALL make_grib_init(now)

  IF (ASSOCIATED(now%next)) THEN
    now => now%next
  ELSE
    EXIT gribinit_loop
  ENDIF
ENDDO gribinit_loop

!------------------------------------------------------------------------------
! Section 4: Loop over all GRIBOUT blocks
!------------------------------------------------------------------------------

! All namelist groups are inquired, whether constant fields shall be written.
! If no further output is done in step 0, a ready-file is written (checked
! with lzwrite_ended)

now => root
lzwrite_ended = .FALSE.

gribout_loop: DO

  !----------------------------------------------------------------------------
  ! Section 4.1: Open YUCHKDAT
  !----------------------------------------------------------------------------

  IF ( (now%lcheck .EQV. .TRUE.) .AND.  (my_cart_id == 0) ) THEN
    ! open file YUCHKDAT

    IF (izdebug > 5) THEN
      WRITE (*,'(A)') '  src_output: opening file YUCHKDAT'
    ENDIF

    OPEN(nuchkdat, FILE=yuchkdat, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yzerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
      izerrstat = 2005
      CALL model_abort (my_cart_id, izerrstat, yzerrmsg, yzroutine)
    ENDIF
  ENDIF

  !----------------------------------------------------------------------------
  ! Section 4.2: Write file with constant fields, if required
  !----------------------------------------------------------------------------

  IF (now%lwrite_const .EQV. .TRUE.) THEN

    IF (izdebug > 5) THEN
      WRITE (*,'(A)') '  src_output: writing file with constant fields'
    ENDIF

    ! Create the filename
    ! -------------------

    yzhead = 'lf'//now%ydomain
    ! Construct filename of constant files always for step 0, also in case of restarts
    CALL make_fn (yzhead, ydate_ini, ydate_ini, now%ytunit,'c', 0     , dt, &
                  now%lhour, itype_calendar, now%ydir, yname, lmmss,        &
                  izdebug, izerror)

    ! In case of netcdf, add extension '.nc' to filename
    IF (now%yform_write == 'ncdf') THEN
      yname = yname(1:LEN_TRIM(yname)) // '.nc'
    ENDIF

    ! Add optional suffix to filename
    IF ( LEN_TRIM(now%ysuffix) /= 0 ) THEN
      yname = yname(1:LEN_TRIM(yname)) // TRIM(now%ysuffix)
    ENDIF

#ifdef MESSY
    IF (L_BM_ORIG_OUTPUT) THEN
       l_COSMO_now = .TRUE.
#endif

    ! open the file
    ! -------------

    IF (now%yform_write /= 'ncdf' .OR. nc_asyn_io < 1) THEN

      IF (now%yform_write == 'bina') THEN
        ! get a free unit-number for Fortran OPEN
        CALL get_free_unit (nuedat)
      ENDIF
      CALL open_file (nuedat, yname, ymode_write, now%yform_write, icomm_cart,  &
                      my_cart_id, num_compute, lasync_io, idbg_level,           &
                      yzerrmsg, izerror)
      IF (izerror /= 0) THEN
        CALL model_abort (my_cart_id, 2031, yzerrmsg, yzroutine)
      ENDIF

#ifdef NETCDF
      ! Write global headers for netcdf file
      ! ------------------------------------
      IF (now%yform_write == 'ncdf'  .AND. nc_asyn_io < 1 ) THEN
        CALL write_nc_gdefs (nuedat, now, icomm_cart, num_compute, 'c',  &
                             NF90_FLOAT, -1, yzerrmsg, izerror)
        IF (izerror /= 0) THEN
          CALL model_abort (my_cart_id, 8052, yzerrmsg, yzroutine)
        ENDIF

        CALL write_nc_vdefs (nuedat, now%nyvar_c, now%ilist_c, ivar_id,  &
                             now%luvmasspoint, icomm_cart, num_compute,  &
                             'c', NF90_FLOAT, yzerrmsg, izerror)
        IF (izerror /= 0) THEN
          CALL model_abort (my_cart_id, 8053, yzerrmsg, yzroutine)
        ENDIF
      ENDIF
#endif

    ENDIF

#ifdef MESSY
    ENDIF
#endif

    ! Write the headline in YUCHKDAT for the file with constant variables
    ! -------------------------------------------------------------------

    IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
      WRITE (nuchkdat,'(A,I7)')                                           &
                          'Check the constant data: '
      WRITE (nuchkdat,'(A,A)')                                            &
         '    File:   ',TRIM(yname)
      WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                 &
         '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
      WRITE (nuchkdat,'(A)') '    '
      WRITE (nuchkdat,'(A,A)')                                            &
        '     var       ee    lev         min      ',                     &
        'imin   jmin          max      imax   jmax         mean  '
    ENDIF

    ! loop over all constant variables that should be written
    ! -------------------------------------------------------

    nzrecords = 0
    DO n = 1, now%nyvar_c

      ! location in the variable table
      i1 = now%ilist_c(1,n)
      i2 = now%ilist_c(2,n)
      i3 = now%ilist_c(3,n)

      IF (izdebug >= 5) THEN
        WRITE (*,'(3a,i4,a)') '  src_output: processing ', &
             TRIM(ADJUSTL(var(i1,i2,i3)%name)),' on PE ',my_cart_id,' ...'
      END IF

      SELECT CASE (var(i1,i2,i3)%rank)
      CASE(3)
        kbot = LBOUND(var(i1,i2,i3)%p3,3)
        ktop = UBOUND(var(i1,i2,i3)%p3,3)
        zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p3(:,:,kbot:ktop)

        DO k=kbot,ktop
          nzrecords = nzrecords + 1
          CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,           &
                    zvarlev(1:ie,1:je,k), now, .FALSE., 'c', 0.0_ireals,    &
                    .FALSE., ivar_id(n), n, iorg_data, izdebug)
        ENDDO

      CASE(2)
        IF (now%yvarc(n)(1:LEN_TRIM(now%yvarc(n))) == 'FIS' ) THEN
          zvarlev(1:ie,1:je,1) = hsurf(1:ie,1:je) * g
        ELSE
          zvarlev(1:ie,1:je,1) = var(i1,i2,i3)%p2(1:ie,1:je)
        ENDIF
        nzrecords   = nzrecords+1
        CALL output_data (nuedat, nzrecords, i1,i2,i3, 1, 1,                &
             zvarlev(:,:,1), now, .FALSE., 'c', 0.0_ireals, .FALSE.,        &
             ivar_id(n), n, iorg_data, izdebug)
      END SELECT
    ENDDO

#ifdef MESSY
    IF (L_BM_ORIG_OUTPUT) THEN
#endif

    ! Flush the output buffers and close the file
    ! -------------------------------------------

    CALL output_data (nuedat, -1, -1,-1,-1, -1, -1,                         &
                zvarlev(:,:,1), now, .TRUE., 'c', 0.0_ireals, .FALSE., -1,  &
                -1, iorg_data, izdebug)

    IF( now%yform_write /= 'ncdf' .OR. nc_asyn_io < 1) THEN
      CALL close_file (nuedat, now%yform_write, icomm_cart, my_cart_id,       &
                       num_compute, lasync_io, idbg_level, yzerrmsg, izerror)
      IF (izerror /= 0) THEN
         CALL model_abort (my_cart_id, 2032, yzerrmsg, yzroutine)
      ENDIF
    ENDIF

    IF (now%yform_write == 'bina') THEN
      ! release the unit-number again
      CALL release_unit (nuedat)
    ENDIF

#ifdef MESSY
    ENDIF 
#endif

    ! Indicate that a ready file is needed
    ! ------------------------------------

    lzwrite_ended = .TRUE.

    ! Write a blank line to YUCHKDAT
    IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
      WRITE (nuchkdat,'(A)') '   '
      WRITE (nuchkdat,'(A)') '   '
    ENDIF
    cur_gribout_idx = cur_gribout_idx+1
  ENDIF

  ! close file nuchkdat
  IF ( (now%lcheck) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  IF (ASSOCIATED(now%next)) THEN
    now => now%next
  ELSE
    EXIT gribout_loop
  ENDIF

ENDDO gribout_loop

!------------------------------------------------------------------------------
! Section 5: Check for further output in step 0
!------------------------------------------------------------------------------

IF ( (my_cart_id == 0) .AND. (lzwrite_ended) ) THEN
  now => root

  gribout_loop_2: DO
    IF (now%ngrib(now%nextstep) == 0) THEN
      ! further output will be written: no ready file is necessary here
      lzwrite_ended = .FALSE.
    ENDIF

    ! Add optional suffix to filename
    IF ( LEN_TRIM(now%ysuffix) /= 0 ) THEN
      yname = yname(1:LEN_TRIM(yname)) // TRIM(now%ysuffix)
    ENDIF

    IF (ASSOCIATED(now%next)) THEN
      now => now%next
    ELSE
      EXIT gribout_loop_2
    ENDIF
  ENDDO gribout_loop_2

  ! Write ready-file, if required
  IF (lzwrite_ended) THEN

    ! Write a blank line to YUCHKDAT
    IF (izdebug > 10) &
      WRITE (*,*) "proc ", my_cart_id, ": Enter write_ready."

    IF (LEN_TRIM(ytrans_out) > 0) THEN
      ! Create the filename LMF_forecasttime
      yzhead = 'LMF'
      CALL make_fn (yzhead, yakdat1, ydate_ini,'f',' ', ntstep, dt, .TRUE.,     &
        itype_calendar, ytrans_out, yname, .TRUE., izdebug, izerrstat)

      ! Write the file
      OPEN  (ntrans_out, FILE=yname, FORM='FORMATTED')
      WRITE (ntrans_out, '(A)') 'ready'
      CLOSE (ntrans_out)
    ENDIF

  ENDIF
ENDIF

! Deallocate arrays for IO
DEALLOCATE (iblock, ibmap, ds_real, ds_grib, dsup, ymessage)

IF (izdebug > 5) THEN
  WRITE (*,'(A)') '  src_output: initialization of output ended'
ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE init_output

!==============================================================================
!+ Module procedure in src_output for organizing the output
!------------------------------------------------------------------------------

SUBROUTINE organize_output (outblock, yextension, numlist, ylist, ilist,   &
                            lcout, iout)

#ifdef MESSY
USE messy_main_channel_bi,  ONLY: cosmo_output_list, COSMOOUT       &
                                , LOUTPUT_NOW, js_COSMOm, js_COSMOp &
                                , js_COSMOz, js_COSMOs, js_COSMOc
USE messy_main_tools,       ONLY: int2str
#endif

!------------------------------------------------------------------------------
!
! Description:
!  The routine organize_output is called for every namelist output group and
!  for every of the three output lists (model variables, pressure level 
!  variables and height level variables). In case of pressure level variables 
!  or height level variables the routines p_int and z_int, resp., are called 
!  for vertical interpolation.
!  
!  Parallelization for the output is by layers that should be written
!  to the grib file. Every PE gets a layer and packs it into grib format
!  (in routine output_data).
!
! Method:
!  - Initializations (for the grib library)
!  - Opening the output grib file
!  - Scanning through the list (loop over all variables)
!  - Closing the output grib file
!
! Output files:
!  Output grib files for model variables (without extension), for 
!  pressure level variables (with extension 'p') and for height level
!  variables (with extension 'z').
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
TYPE(pp_nl),              INTENT(IN)     ::    &
  outblock         ! pointer to the namelist group

! Scalar arguments with intent(in):
CHARACTER (LEN=1),        INTENT(IN)     ::    &
  yextension       ! indicates model variables (''), p-('p') or z-levels ('z')

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  numlist          ! number of elements in ylist

CHARACTER (LEN=clen),     INTENT(IN)     ::    &
  ylist(numlist)   ! list of variables for output

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  ilist(3,numlist) ! number of elements in ylist

! The following variables are only meaningful for COSMO/MESSy
! However, for easier code reading they are always defined
LOGICAL                 , INTENT(IN)     ::    &
  lcout              ! output required for COSMO (not necessarily for CHANNEL)
INTEGER (KIND=iintegers), INTENT(IN), OPTIONAL     ::    &
  iout              ! number of output namelist gribout

!------------------------------------------------------------------------------
!
! Local scalars:
  INTEGER (KIND=iintegers) :: i1,i2,i3, i,j,k, n, ktop, kbot, nentry,  &
                              isens, iorg_data(3,0:num_compute-1)
  INTEGER (KIND=iintegers) :: klev, nuedat, niostat, ierrstat, izlen, nzrecords
  CHARACTER (LEN=250)      :: yname
  CHARACTER (LEN= 25)      :: yroutine
  CHARACTER (LEN= 80)      :: yerrmsg
  CHARACTER (LEN= 14)      :: yzdat1, yzdatact
  CHARACTER (LEN= 28)      :: yzdat2


! Local arrays:

REAL (KIND=ireals)   ::     &
#ifndef MESSY
  zvarlev  (ie,je,0:noutlevels),&! variable for vertical interpolation
#endif    
  slev     (0:noutlevels)        ! stores the z- or the p-levels

#ifdef MESSY
  REAL(KIND=ireals), DIMENSION(:,:,:), POINTER :: zvarlev
#endif

REAL (KIND=ireals)   ::     & !
  zenith_t (ie,je),         & ! Arrays for computing the zenith total (dry,
  zenith_w (ie,je),         & ! hydrostatic) delay
  zenith_h (ie,je),         & ! 
  zcape_mu (ie,je),         & ! Arrays for most unstable CAPE
  zcin_mu  (ie,je),         & ! ... and CIN
  zcape_ml (ie,je),         & ! Arrays for mixed layer CAPE
  zcin_ml  (ie,je),         & ! ... and CIN
  zcape_3km(ie,je),         & ! ... and CAPE 3KM
  zlcl_ml  (ie,je),         & ! ... and LCL
  zlfc_ml  (ie,je),         & ! ... and LFC    
  zhelp2d  (ie,je),         & !
  zbrn  (ie,je,ke),         & ! Array for Bulk Richardson Number
  zhelp1(ie,je,ke),         & !
  zhelp2(ie,je,ke),         & !
  zhelp3(ie,je,ke),         & !
  zhelp4(ie,je,ke)            !

REAL (KIND=ireals), PARAMETER   ::     & !
  missing_value = -999.9_ireals

REAL (KIND=ireals)              ::     & !
  zacthour

LOGICAL              ::     &
  lzenith,                  & ! indicates whether to compute zenith delay
  lzrestart,                & ! indicates whether restart-files are written
  lconvind_mu,              & ! indicates whether CAPE_MU/CIN_MU has been computed
  lconvind_ml,              & ! indicates whether CAPE_ML/CIN_ML has been computed
  l_brn,                    & ! indicates whether BULK RICHARDSON NUMBER has been computed
  lwarn,                    & ! to indicate whether a SR issues some warnings
  l_calc_sdi_already_comp     ! prevents from a double computing of Subr. 'calc_sdi'

INTEGER (KIND=iintegers) :: kzdims(24)

INTEGER (KIND=iintegers) :: &
  ivar_id(nzmxid), & ! NetCDF-ID of each variable in the output list
  nzjulianday, izvctype_write, izdebug, izerror

CHARACTER (LEN=3)        :: &
  yzhead           ! characterizes the special kind of the data

#ifdef MESSY
 TYPE(cosmo_output_list),   POINTER :: channeli
 TYPE(cosmo_output_list),   POINTER :: channele
 CHARACTER(LEN=9)                   :: chname
 CHARACTER(LEN=3)                   :: str
 CHARACTER(LEN=1)                   :: tag
#endif

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1a: Initializations
!------------------------------------------------------------------------------

  yroutine(:) = ' '
  yroutine = 'organize_output'
  ierrstat = 0
  izerror  = 0
  lzenith  = .FALSE.
  lconvind_mu  = .FALSE.
  lconvind_ml  = .FALSE.
  l_brn = .FALSE.
  lwarn = .FALSE.
  l_calc_sdi_already_comp = .FALSE.

#ifdef MESSY
  ivar_id(:) = 0.  ! um_ak_20080312
  zbrn(:,:,:) = 0. ! um_gg_20121024
  IF (L_FORCE_calcout) THEN ! output calculations required for channel
     l_COSMO_now = lcout .AND. L_BM_ORIG_OUTPUT
     IF  (.NOT. l_COSMO_now) THEN  ! no COSMO output
        SELECT CASE (yextension)   ! test if output for this channel required
        CASE (' ')
           IF (js_COSMOm(iout) > 0) THEN ! CHANNEL exists
              IF (.NOT. LOUTPUT_NOW(js_COSMOm(iout))) RETURN  ! no output
              ! required for present step => RETURN
           ELSE ! no 'm' variables for this GRIBOUT namelist defined
              RETURN
           ENDIF
        CASE ('p')
           IF (js_COSMOp(iout) > 0) THEN
              IF  (.NOT. LOUTPUT_NOW(js_COSMOp(iout))) RETURN
           ELSE
              RETURN
           ENDIF
        CASE ('z')
           IF (js_COSMOz(iout) > 0) THEN
              IF (.NOT. LOUTPUT_NOW(js_COSMOz(iout))) RETURN
           ELSE
              RETURN
           ENDIF
        CASE ('s')
           IF (js_COSMOs(iout) > 0) THEN
              IF (.NOT. LOUTPUT_NOW(js_COSMOs(iout))) RETURN
           ELSE
              RETURN
           ENDIF
        CASE ('c')
           IF (js_COSMOc(iout) > 0) THEN
              IF (.NOT. LOUTPUT_NOW(js_COSMOc(iout))) RETURN
           ELSE
              RETURN
           ENDIF
        END SELECT
     ELSE ! COSMO output required
     ! NOTHING TODO
     ENDIF
  ELSE ! no output calculations for CHANNEL required
     IF (.NOT. L_BM_ORIG_OUTPUT) RETURN ! no COSMO output required
     l_COSMO_now=.TRUE.
  ENDIF
  IF (numlist == 0) RETURN
  IF (yextension == ' ') THEN
     tag = 'm'
  ELSE
     tag = yextension
  ENDIF
  IF (PRESENT(IOUT)) THEN
     CALL int2str(str,iout)
     chname = 'COSMO'//tag//str
  ELSE
     chname = 'COSMO'//tag
  ENDIF
  channeli => COSMOOUT
  DO
     IF (.NOT. ASSOCIATED(channeli)) THEN
        write (0,*) 'COSMO-OUTPUT not ASSOCIATED '
        CALL model_abort (my_cart_id, 3333, yerrmsg, yroutine)
     ENDIF
     IF (TRIM(chname) == TRIM(channeli%this%label)) EXIT
     channele => channeli
     channeli => channeli%next
  END DO
#endif

  ! Initialize, whether additional debug output shall be done
  IF (ldebug_io) THEN
    IF (lprintdeb_all) THEN
      izdebug = idbg_level
    ELSE
      IF (my_cart_id == 0) THEN
        izdebug = idbg_level
      ELSE
        izdebug = 0
      ENDIF
    ENDIF
  ELSE
    izdebug = 0
  ENDIF

  IF ((yextension == 'o') .OR. (yextension == 'n')) THEN
    lzrestart = .TRUE.
  ELSE
    lzrestart = .FALSE.
  ENDIF

  ! Set lfd, lds and lbm
  lds = ie_tot * je_tot
  lbm = 1875
  lfd = compute_grib_intbuffer_length(ie_tot, je_tot, iwlength)
  lfa = compute_grib_intbuffer_length(ie_tot, je_tot, 1) ! gives length in bytes

  ! Initializations for the grib library
  !  moving arraydimensions into idims
  !  declaration dimensions
  idims_out( 7) = outblock%ie_out_tot * outblock%je_out_tot
  
  !  real dimensions
  idims_out(15) = outblock%ie_out_tot * outblock%je_out_tot
  idims_out(17) = outblock%ie_out_tot * outblock%je_out_tot

  ! Allocate GRIB arrays
  ALLOCATE (iblock(lfd), ibmap(lbm),                    STAT=ierrstat)
  ALLOCATE (ds_real(lds), ds_grib(lds), dsup(ndsup)   , STAT=ierrstat)
  ALLOCATE (ymessage(lfa),                              STAT=ierrstat)

  slev(:)  = 0.0_ireals

  ! gridpoints, simple packing, floating point data
  ibds(2)   = 0
 
  ! nrbit, number of bits
  ibds(5)   = outblock%nrbit

  ! no bitmap
  ibms(3)   = -2

  ! determine the timelevel for output
  IF (lzrestart) THEN
    IF (yextension == 'o') THEN
      IF (.NOT. l2tls) THEN
        itimelevel = nnow   ! for leapfrog
      ELSE
        itimelevel = nnew   ! for Runge-Kutta
      ENDIF
    ELSEIF (yextension == 'n') THEN
      itimelevel = nnew
    ENDIF
  ELSE
    ! use nnow for output (this was nnew before)
    itimelevel = nnow
  ENDIF

  ! Initialization for ivar_id (should be set also for non netcdf output)
  ivar_id(:) = 0_iintegers

!------------------------------------------------------------------------------
!  Section 1b: Gather hsurf field to all PEs, if necessary
!------------------------------------------------------------------------------

  IF (outblock%l_fi_pmsl_smooth) THEN
    ALLOCATE (hsurf_tot(ie_tot,je_tot), STAT=ierrstat)

    IF (num_compute > 1) THEN
      CALL gather_field (hsurf, ie,je, hsurf_tot, ie_tot,je_tot, -1, ierrstat)
    ELSE
      hsurf_tot(:,:) = hsurf(:,:)
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Open the grib file
!------------------------------------------------------------------------------

  ! creating filename
  IF (outblock%lanalysis) THEN
    yzhead = 'la'//outblock%ydomain
  ELSE
    IF (lzrestart) THEN
      yzhead = 'lr'//outblock%ydomain
    ELSE
      yzhead = 'lf'//outblock%ydomain
    ENDIF
  ENDIF

  IF ((outblock%yform_write == 'bina') .AND. (yextension=='o' .OR. yextension=='n')) THEN
    ! The date of the next time step has to be determined to get the proper
    ! file name also for ytunit='d'
    CALL get_utc_date(ntstep+1, ydate_ini, dt, itype_calendar, yzdat1,       &
                      yzdat2, nzjulianday, zacthour)
    CALL make_fn (yzhead, yzdat1, ydate_ini, outblock%ytunit, yextension,    &
                  ntstep+1, dt, outblock%lhour, itype_calendar,              &
                  outblock%ydir, yname, lmmss, izdebug, ierrstat)
  ELSE
    IF (izdebug > 10) THEN
      PRINT *, ' calling make_fn with date/unit: ', yakdat1, ' ', outblock%ytunit
    ENDIF
    CALL make_fn (yzhead, yakdat1, ydate_ini, outblock%ytunit, yextension,   &
                  ntstep, dt, outblock%lhour, itype_calendar,                &
                  outblock%ydir, yname, lmmss, izdebug, ierrstat)
  ENDIF

  ! In case of netcdf, add extension '.nc' to filename
  !BEGIN KOchange
  IF (outblock%yform_write == 'ncdf') THEN
!  IF ((outblock%yform_write == 'ncdf') .AND. (.NOT. (yextension=='o' .OR. yextension=='n'))) THEN
!  IF ((outblock%yform_write == 'ncdf') .AND. (yextension=='o' .OR. yextension=='n')) THEN
  !END KOchange
    yname = yname(1:LEN_TRIM(yname)) // '.nc'
  ENDIF

  ! Add optional suffix to filename
  IF ( LEN_TRIM(outblock%ysuffix) /= 0 ) THEN
    yname = yname(1:LEN_TRIM(yname)) // TRIM(outblock%ysuffix)
  ENDIF

  IF (outblock%yform_write == 'bina') THEN
    ! get a free unit-number for Fortran OPEN
    CALL get_free_unit (nuedat)
  ENDIF

#ifdef MESSY
   IF (l_COSMO_now .AND. L_BM_ORIG_OUTPUT)  THEN
#endif

!  IF (outblock%yform_write /= 'ncdf' .OR. nc_asyn_io < 1 ) THEN
  IF (outblock%yform_write /= 'ncdf' .OR. nc_asyn_io < 1 .OR. (nc_asyn_io > 0 .AND. lzrestart)) THEN
    CALL open_file (nuedat, yname, ymode_write, TRIM(outblock%yform_write),      &
                    icomm_cart, my_cart_id, num_compute, lasync_io, idbg_level,  &
                    yerrmsg, ierrstat)
    IF (ierrstat /= 0) THEN
       CALL model_abort (my_cart_id, 2033, yerrmsg, yroutine)
    ENDIF
  ENDIF
  IF ( (outblock%yform_write == 'bina') .AND. (my_cart_id == 0) .AND. (yextension == 'o')) THEN
    ! write the initial values for the meanvalues and the tke time level
    WRITE (nuedat,IOSTAT=niostat) psm0, dsem0, msem0, kem0, qcm0, ntke
    ! write the vertical coordinate parameters
    IF     (refatm%irefatm == 1) THEN
      izvctype_write = vcoord%ivctype
      IF     (vcoord%ivctype == 1) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%sigm_coord
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%vert_coord
      ENDIF
    ELSEIF (refatm%irefatm == 2) THEN
      izvctype_write = vcoord%ivctype+100
      IF     (vcoord%ivctype == 1) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%sigm_coord
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%vert_coord
      ENDIF
      WRITE (nuedat,IOSTAT=niostat) refatm%delta_t, refatm%h_scal
    ELSEIF (refatm%irefatm == 3) THEN
      izvctype_write = vcoord%ivctype+200
      IF     (vcoord%ivctype == 1) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%sigm_coord
      ELSEIF ( ANY( vcoord%ivctype == (/2,3,4/) ) ) THEN
        WRITE (nuedat,IOSTAT=niostat) izvctype_write, refatm%p0sl, refatm%t0sl, &
              refatm%dt0lp, vcoord%vcflat, vcoord%vert_coord
      ENDIF
      WRITE (nuedat,IOSTAT=niostat) refatm%bvref
    ENDIF
  ENDIF

#ifdef NETCDF
!  IF (outblock%yform_write == 'ncdf' .AND. nc_asyn_io < 1) THEN
  IF (outblock%yform_write == 'ncdf' .AND. (nc_asyn_io < 1 .OR. lzrestart)) THEN
    ! Write global headers for netcdf file
    IF ((yextension == 'o') .OR. (yextension == 'n')) THEN 
      CALL write_nc_gdefs (nuedat, outblock, icomm_cart, num_compute,       &
                           yextension, NF90_DOUBLE, -1, yerrmsg, ierrstat)

    ELSE
      CALL write_nc_gdefs (nuedat, outblock, icomm_cart, num_compute,       &
                           yextension, NF90_FLOAT, -1, yerrmsg, ierrstat)
    ENDIF
    IF (ierrstat /= 0) THEN
      CALL model_abort (my_cart_id, 8052, yerrmsg, yroutine)
    ENDIF
    IF ((yextension == 'o') .OR. (yextension == 'n')) THEN 
      CALL write_nc_vdefs (nuedat, numlist, ilist, ivar_id,                 &
                           outblock%luvmasspoint, icomm_cart, num_compute,  &
                           yextension, NF90_DOUBLE, yerrmsg, ierrstat)
    ELSE 
      CALL write_nc_vdefs (nuedat, numlist, ilist, ivar_id,                 &
                           outblock%luvmasspoint, icomm_cart, num_compute,  &
                           yextension, NF90_FLOAT, yerrmsg, ierrstat)
    ENDIF
    IF (ierrstat /= 0) THEN
      CALL model_abort (my_cart_id, 8053, yerrmsg, yroutine)
    ENDIF
  ENDIF
#endif

#ifdef MESSY
ENDIF
#endif

  ! Write the headline in YUCHKDAT for this file
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    OPEN(nuchkdat, FILE=yuchkdat, FORM=  'FORMATTED', STATUS='UNKNOWN',  &
                   POSITION='APPEND', IOSTAT=niostat)
    IF(niostat /= 0) THEN
      yerrmsg = ' ERROR    *** Error while opening file YUCHKDAT *** '
      ierrstat = 2005
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    ENDIF

    WRITE (nuchkdat,'(A,I7)')                                               &
                        'Check the data in output file for step: ', ntstep
    WRITE (nuchkdat,'(A,A)')                                                &
           '    File:   ',TRIM(yname)
    WRITE (nuchkdat,'(A,I5,A,I5,A,I5)')                                     &
           '    ie_tot =',ie_tot,'   je_tot =',je_tot,'   ke_tot =',ke_tot
    WRITE (nuchkdat,'(A)') '    '
    WRITE (nuchkdat,'(A,A)')                                                &
          '     var       ee    lev         min      ',                     &
          'imin   jmin          max      imax   jmax         mean  '
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Look for output variables in LM variable table
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! 3.1. Retrieve the required microphysics tracers
  !----------------------------------------------------------------------------

  ! Always existing tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = itimelevel, ptr = qv)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = itimelevel, ptr = qc)
  IF (izerror /= 0) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF

  ! Conditionally existing tracers
  CALL trcr_get(izerror, idt_qi, ptr_tlev = itimelevel, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = itimelevel, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg,                       &
                       yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = itimelevel, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = itimelevel, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF

#ifdef TWOMOM_SB
  ! Tracers for the 2-moment scheme
  CALL trcr_get(izerror, idt_qh, ptr_tlev = itimelevel, ptr = qh)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnc, ptr_tlev = itimelevel, ptr = qnc)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnr, ptr_tlev = itimelevel, ptr = qnr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qni, ptr_tlev = itimelevel, ptr = qni)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qns, ptr_tlev = itimelevel, ptr = qns)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qng, ptr_tlev = itimelevel, ptr = qng)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qnh, ptr_tlev = itimelevel, ptr = qnh)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yerrmsg, yroutine)
  ENDIF
#endif

  ! Re-compute qrs and rho for output timelevel itimelevel
  ! (in the long term storage they are only available on timelevel nnew)
  ALLOCATE ( zrho_itl (ie,je,ke), zqrs_itl (ie,je,ke), stat=izerror )
  IF ( izerror /= 0 ) THEN
    CALL model_abort (my_cart_id, 90034, 'Error allocating zrho_itl, zqrs_its', yroutine)
  END IF

  IF (lzrestart) THEN
    zrho_itl(:,:,:) = rho(:,:,:)
    zqrs_itl(:,:,:) = qrs(:,:,:)
  ELSE
    ! Sum up qrs for the output timelevel:
    zqrs_itl = qr
    IF (ASSOCIATED(qi)) THEN
      zqrs_itl = zqrs_itl + qi
    ENDIF
    IF ( ASSOCIATED(qs) ) THEN
      zqrs_itl = zqrs_itl + qs
    END IF
    IF ( ASSOCIATED(qg) ) THEN
      zqrs_itl = zqrs_itl + qg
    END IF

#ifdef TWOMOM_SB
    IF (itype_gscp >= 2000) THEN
      zqrs_itl = zqrs_itl + qh
    ENDIF
#endif

    ! Density on the timelevels itimelevel:
    CALL calrho( t(:,:,:,itimelevel), pp(:,:,:,itimelevel), qv, &
               qc, zqrs_itl, p0, zrho_itl, ie, je, ke, r_d,     &
               rvd_m_o )
  ENDIF

  nzrecords = 0

  ! loop over all variables that should be written and loop over all variables
  ! in the LM variable table until equal elements are found
  write_loop: DO n = 1, numlist

#ifdef MESSY
     zvarlev => channeli%this%vars(n)%ptr(:,:,:,1)
#endif

    ! indices of field in variable table
    i1 = ilist(1,n)
    i2 = ilist(2,n)
    i3 = ilist(3,n)

    izlen = LEN_TRIM(ylist(n))
    IF (ylist(n)(1:izlen) == 'TKE') THEN
      itl = ntke
    ELSE
      itl = itimelevel
    ENDIF

    IF ( izdebug >= 5 ) THEN
      PRINT *, ' src_output: processing ', ylist(n)(1:izlen)
    ENDIF

    SELECT CASE (var(i1,i2,i3)%rank)
      ! pack data depending on the rank

    CASE(4)
      ! vertical interpolation, if necessary
      IF     (yextension == 'p') THEN
        CALL p_int (outblock, i1,i2,i3, n, izdebug, zvarlev(:,:,1:outblock%kepin))
        kbot = 1
        ktop = outblock%kepin
        slev(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      ELSEIF (yextension == 'z') THEN
        CALL z_int (outblock, i1,i2,i3, n, izdebug, zvarlev(:,:,1:outblock%kezin))
        kbot = 1
        ktop = outblock%kezin
        slev(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
      ELSE
        ! calculate additional non-global 4-d fields if required
        IF (ylist(n)(1:izlen) == 'TKE' .AND. (.NOT. lzrestart) ) THEN
          kbot =    1
          ktop = ke+1
          SELECT CASE( itype_turb )
          CASE( 5:8 )
            zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p4(:,:,kbot:ktop,itl)
          CASE default
            zvarlev(:,:,kbot:ktop) = &
              0.5_ireals * (var(i1,i2,i3)%p4(:,:,kbot:ktop,itl))**2
          END SELECT
        ELSEIF ( (ylist(n)(1:izlen) == 'U') .AND.    &
                    (outblock%luvmasspoint) .AND. (.NOT. lzrestart) ) THEN
          ! determine first and last level
          kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          DO k = kbot, ktop
            DO j = 1, je
              DO i = 2, ie
                zvarlev(i,j,k) = 0.5 * (var(i1,i2,i3)%p4(i-1,j,k,itl) + &
                                        var(i1,i2,i3)%p4(i  ,j,k,itl))
              ENDDO
              zvarlev(1,j,k) = zvarlev(2,j,k)
            ENDDO
          ENDDO
        ELSEIF ( (ylist(n)(1:izlen) == 'V') .AND.    &
                    (outblock%luvmasspoint) .AND. (.NOT. lzrestart) ) THEN
          ! determine first and last level
          kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          DO k = kbot, ktop
            DO j = 2, je
              DO i = 1, ie
                zvarlev(i,j,k) = 0.5 * (var(i1,i2,i3)%p4(i,j-1,k,itl) + &
                                        var(i1,i2,i3)%p4(i,j  ,k,itl))
              ENDDO
            ENDDO
            zvarlev(:,1,k) = zvarlev(:,2,k)
          ENDDO
        ELSEIF ((.NOT. lzrestart) .AND.                                 &
                ( ylist(n)(1:izlen) == 'QC' .OR.  ylist(n)(1:izlen) == 'QNCLOUD'   .OR. &
                  ylist(n)(1:izlen) == 'QI' .OR.  ylist(n)(1:izlen) == 'QNICE'     .OR. &
                  ylist(n)(1:izlen) == 'QG' .OR.  ylist(n)(1:izlen) == 'QNGRAUPEL' .OR. &
                  ylist(n)(1:izlen) == 'QR' .OR.  ylist(n)(1:izlen) == 'QNRAIN'    .OR. &
                  ylist(n)(1:izlen) == 'QS' .OR.  ylist(n)(1:izlen) == 'QNSNOW'    .OR. &
                  ylist(n)(1:izlen) == 'QH' .OR.  ylist(n)(1:izlen) == 'QNHAIL') ) THEN
          ! Output of cloud microphysics variables either as densities or mass specific:
          ! determine first and last level
          kbot = LBOUND (var(i1,i2,i3)%p4,3)
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          IF (outblock%loutput_q_densities) THEN
            ! output as densities:
            zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p4(:,:,kbot:ktop,itl) * zrho_itl(:,:,kbot:ktop)
          ELSE
            ! output as mass specific:
            zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p4(:,:,kbot:ktop,itl)
          ENDIF
        ELSE
          ! determine first and last level
          IF ( ((ylist(n)(1:izlen) == 'T_SO') .OR. (ylist(n)(1:izlen) == 'T_SNOW_M')) &
               .AND. (outblock%yform_write == 'ncdf') ) THEN
            IF (lzrestart) THEN
              kbot = 0    !LBOUND (var(i1,i2,i3)%p4,3)+1
            ELSE
              kbot = 1
            ENDIF
          ELSE
#ifdef MESSY
             ! Due to the association of zvarlev to the channel object pointer
             ! only T_SO(1:ktop) can be written out also in grib format
             kbot = 1
#else
             kbot = LBOUND (var(i1,i2,i3)%p4,3)
#endif
          ENDIF
          ktop = UBOUND (var(i1,i2,i3)%p4,3)
          zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p4(:,:,kbot:ktop,itl)
        END IF
      ENDIF

      DO k = kbot, ktop
        nzrecords = nzrecords + 1
        !KOCHANGE T_SO
        IF  ((ylist(n)(1:izlen) == 'T_SO') .AND. lzrestart .AND. &
             (outblock%yform_write == 'ncdf')) THEN
          CALL output_data (nuedat, nzrecords, i1,i2,i3, k+1, ktop,  &
             zvarlev(:,:,k), outblock, .FALSE., yextension, slev(k), &
             lzrestart, ivar_id(n), n, iorg_data, izdebug)
        ELSE
        CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,            &
             zvarlev(1:ie,1:je,k), outblock, .FALSE., yextension, slev(k), &
             lzrestart, ivar_id(n), n, iorg_data, izdebug)
        ENDIF
      ENDDO

    CASE(3)
      IF     (yextension == 's') THEN
#if defined RTTOV7 || defined RTTOV9 || defined RTTOV10
        IF     (ylist(n)(1:izlen) == 'SYNME7') THEN
          ! Look for entry in sat_compute
          nentry = -1
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:8)=='METEOSAT') .AND.    &
                 (sat_compute(isens)%nsat_id        == 7        ) .AND.    &
                 (sat_compute(isens)%ysensor        =='MVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          kbot = 1
          IF (luse_rttov .AND. (nentry > 0)) THEN
             ktop = UBOUND(synme7,3)
             zvarlev(:,:,kbot:ktop) = synme7(:,:,kbot:ktop)
             slev(:) = REAL(nentry)
          ELSE
             ktop = 0
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'SYNMSG') THEN
          ! Look for entry in sat_compute
          nentry = -1
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          kbot = 1
          IF (luse_rttov .AND. (nentry > 0)) THEN
             ktop = UBOUND(synmsg,3)
             zvarlev(:,:,kbot:ktop) = synmsg(:,:,kbot:ktop)
             slev(:) = REAL(nentry)
          ELSE
             ktop = 0
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'MSG_TB') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,1:29:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry,ireals)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_TBC') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,2:30:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry,ireals)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_RAD') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) = synmsg(:,:,3:31:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry,ireals)
        ELSEIF (ylist(n)(1:izlen) == 'MSG_RADC') THEN
          kbot = 1
          ktop = nmsgchan
          zvarlev(:,:,kbot:ktop) =  synmsg(:,:,4:32:4)
          DO isens = 1, num_sensors
            IF ( (sat_compute(isens)%ysatellite(1:3)=='MSG'     ) .AND.    &
                ((sat_compute(isens)%nsat_id        == 1        ) .OR.    &
                 (sat_compute(isens)%nsat_id        == 2        )) .AND.  &
                 (sat_compute(isens)%ysensor        =='SEVIRI'   ) ) THEN
              nentry = isens
            ENDIF
          ENDDO
          slev(:) = REAL(nentry,ireals)
        ENDIF
#endif
      ELSEIF (yextension == 'p') THEN
        ! vertical interpolation on p-levels, if necessary
        CALL p_int (outblock, i1,i2,i3, n, izdebug, zvarlev(:,:,1:outblock%kepin))
        kbot = 1
        ktop = outblock%kepin
        slev(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      ELSEIF (yextension == 'z') THEN
        ! vertical interpolation on z-levels, if necessary
        CALL z_int (outblock, i1,i2,i3, n, izdebug, zvarlev(:,:,1:outblock%kezin))
        kbot = 1
        ktop = outblock%kezin
        slev(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
      ELSE
        !Calculate additional non-global 3-d fields if requird
        IF ( ylist(n)(1:izlen) == 'P' ) THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = p0(:,:,kbot:ktop) + pp(:,:,kbot:ktop,itl)
        ELSEIF (ylist(n)(1:izlen) == 'OMEGA') THEN
          kbot = 1
          ktop = ke
          CALL calomega (zvarlev(:,:,kbot:ktop), pp(:,:,:,nnew),          &
                         pp(:,:,:,itl ), pptens(:,:,:), w(:,:,:,itl),     &
                         rho0(:,:,:), ie, je, ke, dt, g )
        ELSEIF (ylist(n)(1:izlen) == 'CLC') THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = clc_sgs(:,:,kbot:ktop) &
             +  clc_con(:,:,kbot:ktop)*(1.0_ireals -  clc_sgs(:,:,kbot:ktop))
        ELSEIF ( ylist(n)(1:izlen) == 'QRS' ) THEN
          kbot = 1
          ktop = ke
          IF (outblock%loutput_q_densities) THEN
            zvarlev(:,:,kbot:ktop) = zqrs_itl(:,:,kbot:ktop) * zrho_itl(:,:,kbot:ktop)
          ELSE
            zvarlev(:,:,kbot:ktop) = zqrs_itl(:,:,kbot:ktop)
          END IF
        ELSEIF (ylist(n)(1:izlen) == 'Q_SEDIM') THEN
          kbot = 1
          ktop = ke
          IF (ASSOCIATED(qi)) THEN
            zvarlev(:,:,kbot:ktop) = zqrs_itl(:,:,kbot:ktop) - qi(:,:,kbot:ktop)
          ELSE
            zvarlev(:,:,kbot:ktop) = zqrs_itl(:,:,kbot:ktop)
          ENDIF
          IF (outblock%loutput_q_densities) THEN
            zvarlev(:,:,kbot:ktop) = zvarlev(:,:,kbot:ktop) * zrho_itl(:,:,kbot:ktop)
          ENDIF
        ELSEIF (ylist(n)(1:izlen) == 'RELHUM') THEN
          kbot = 1
          ktop = ke
          CALL calrelhum(zvarlev(:,:,kbot:ktop), t(:,:,:,itl), pp(:,:,:,itl),&
                         p0(:,:,:), qv(:,:,:), ie, je, ke,                   &
                         b1, b2w, b3, b4w, rdv, o_m_rdv )
        ELSEIF (ylist(n)(1:izlen) == 'BRN') THEN
          kbot = 1
          ktop = ke
          CALL calc_bulk_richardson(zvarlev(:,:,kbot:ktop),t(:,:,:,itl),   &
                qv(:,:,:), u(:,:,:,itl), v(:,:,:,itl),                     &
                p0(:,:,:)+pp(:,:,:,itl), hsurf(:,:), ps(:,:,itl),          &
                t_2m(:,:), qv_2m(:,:),hhl(:,:,:),                          &
                ie, je, ke, cp_d, r_d,  rvd_m_o, g)
          zbrn(:,:,:) = zvarlev(:,:,kbot:ktop)
          l_brn = .TRUE.
        ELSEIF (ylist(n)(1:izlen) == 'FI_ANAI') THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = p_anai(:,:,kbot:ktop) / zrho_itl(:,:,kbot:ktop)
        ELSEIF (ylist(n)(1:izlen) == 'DBZ') THEN
          kbot = 1
          ktop = ke
          IF (itype_gscp == 3) THEN
            CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
                 klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),           &
                 qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,      &
                 qs(:,:,:)*zrho_itl, z_radar = zvarlev(:,:,kbot:ktop) )
          ELSEIF (itype_gscp == 4) THEN
            CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
                 klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),           &
                 qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,      &
                 qs(:,:,:)*zrho_itl, q_grau  = qg(:,:,:)*zrho_itl,                &
                 z_radar = zvarlev(:,:,kbot:ktop) )
#ifdef TWOMOM_SB
          ELSEIF (itype_gscp >= 100) THEN
            CALL radar_sb_ray (ie, je, ke, pi,                               &
                 klv850, my_cart_id,  t(:,:,:,itl),                          &
                 qc(:,:,:)*zrho_itl,  qr(:,:,:)*zrho_itl,                    &
                 qi(:,:,:)*zrho_itl,  qs(:,:,:)*zrho_itl,                    &
                 qg(:,:,:)*zrho_itl,  qh(:,:,:)*zrho_itl,                    &
                 qnc(:,:,:)*zrho_itl, qnr(:,:,:)*zrho_itl,                   &
                 qni(:,:,:)*zrho_itl, qns(:,:,:)*zrho_itl,                   &
                 qng(:,:,:)*zrho_itl, qnh(:,:,:)*zrho_itl,                   &
                 z_radar = zvarlev(:,:,kbot:ktop) )
#endif
          ENDIF
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_U' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zvarlev(:,:,1:ke), zhelp2(:,:,:), zhelp3(:,:,:))
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_V' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zhelp1(:,:,:), zvarlev(:,:,1:ke), zhelp3(:,:,:))
        ELSEIF ( ylist(n)(1:izlen) == 'VORTIC_W' ) THEN
          kbot = 1
          ktop = ke
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .TRUE., zhelp1(:,:,:), zhelp2(:,:,:), zvarlev(:,:,1:ke))
        ELSEIF ( ylist(n)(1:izlen) == 'POT_VORTIC' ) THEN
          kbot = 1
          ktop = ke
          CALL calc_Theta_Tppp( t(:,:,:,itl ), pp(:,:,:,itl ), p0,           &
                                ie, je, ke, r_d, cp_d, zhelp4)
          CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
                     sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
                     u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
                     .FALSE., zhelp1, zhelp2, zhelp3)

          ! coriolis parameter with cosine only available for deep atmosphere
          IF ( .NOT. lcori_deep ) THEN
            zhelp2d(:,:) = 0.0_ireals
          ELSE
            zhelp2d(:,:) = fccos(:,:)
          ENDIF

          CALL potential_vorticity_rho( ie, je, ke, eddlon, eddlat, r_earth, &
                     fc, zhelp2d, sqrtg_r_s, dzeta_dlam, dzeta_dphi,         &
                     zhelp1, zhelp2, zhelp3, lmetr, zhelp4,                  &
                     u(:,:,:,itl), v(:,:,:,itl), w(:,:,:,itl), zvarlev(:,:,1:ke))
          zvarlev(:,:,1:ke) = zvarlev(:,:,1:ke) / zrho_itl(:,:,1:ke)
        ELSEIF ( ylist(n)(1:izlen) == 'RHO_TOT' ) THEN
          kbot = 1
          ktop = ke
          zvarlev(:,:,kbot:ktop) = zrho_itl(:,:,kbot:ktop)
        ELSE
          kbot = LBOUND(var(i1,i2,i3)%p3,3)
          ktop = UBOUND(var(i1,i2,i3)%p3,3)
          ! decision, if mlf or slf
          IF ( (ktop <= 3) .AND. (var(i1,i2,i3)%levtyp .NE. 333) ) THEN
            ! These are just 2D arrays with a timelevel and no special third dimension
            kbot = 1
            ktop = 1
            zvarlev(:,:,1) = var(i1,i2,i3)%p3(1:ie,1:je,itl)
          ELSE
            zvarlev(:,:,kbot:ktop) = var(i1,i2,i3)%p3(:,:,kbot:ktop)
          ENDIF
        ENDIF
      ENDIF

      ! Distribute the multidimensional fields to the PEs
      DO k=kbot,ktop
        nzrecords = nzrecords + 1
        CALL output_data (nuedat, nzrecords, i1,i2,i3, k, ktop,            &
             zvarlev(1:ie,1:je,k), outblock, .FALSE., yextension, slev(k), &
             lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ENDDO

    CASE(2)
      ! Calculate additional non-global output fields if required
      IF ( ylist(n)(1:izlen) == 'PMSL' ) THEN
        CALL calpmsl( zvarlev(:,:,1), ps(:,:,itl), t(:,:,ke,itl),       &
               rho0(:,:,ke), dp0(:,:,ke), hsurf, ie, je, g, r_d )
      ELSEIF ( ylist(n)(1:izlen) == 'TOT_PREC' ) THEN
        IF (itype_gscp >= 4 .AND. itype_gscp < 2000) THEN
          CALL calprsum( zvarlev(:,:,1), rain_gsp, snow_gsp + grau_gsp, &
                           rain_con, snow_con, ie, je )
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp >= 2000) THEN
          CALL calprsum( zvarlev(:,:,1), rain_gsp, snow_gsp + grau_gsp + hail_gsp, &
                           rain_con, snow_con, ie, je )
#endif
        ELSE
          CALL calprsum( zvarlev(:,:,1), rain_gsp, snow_gsp, rain_con,  &
                         snow_con, ie, je )
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TOT_PR' ) THEN
        IF (itype_gscp >= 4 .AND. itype_gscp < 2000) THEN
          zvarlev(:,:,1) = prr_gsp(:,:) + prs_gsp(:,:) + prg_gsp(:,:) + &
                           prr_con(:,:) + prs_con(:,:)
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp >= 2000) THEN
          zvarlev(:,:,1) = prr_gsp(:,:) + prs_gsp(:,:) + prg_gsp(:,:) + &
                           prh_gsp(:,:) + prr_con(:,:) + prs_con(:,:)
#endif
        ELSE
          zvarlev(:,:,1) = prr_gsp(:,:) + prs_gsp(:,:) +                &
                           prr_con(:,:) + prs_con(:,:)
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'HTOP_DC' ) THEN
        CALL caltopdc( zvarlev(:,:,1), t(:,:,:,itl), p0, pp(:,:,:,itl), &
                       qv(:,:,:), hhl, vcoord%vert_coord, ie, je, ke,   &
                       b1, b2w, b3, b4w, rdv, o_m_rdv, rvd_m_o, cpdr, g )
      ELSEIF ( ylist(n)(1:izlen) == 'FIS' ) THEN
        zvarlev(1:ie,1:je,1) = hsurf(1:ie,1:je) * g
      ELSEIF ( ylist(n)(1:izlen) == 'HTOP_CON'         .OR.      &
               ylist(n)(1:izlen) == 'HTOP_SC' ) THEN
        DO j= 1, je
          DO i = 1, ie 
            zvarlev(i,j,1) = 0.0_ireals
            klev = NINT( top_con(i,j) )
            IF(klev > 0) THEN
              zvarlev(i,j,1) = 0.5*(hhl(i,j,klev)+hhl(i,j,klev+1))
            ENDIF
          ENDDO 
        ENDDO
      ELSEIF ( ylist(n)(1:izlen) == 'HBAS_CON'         .OR.      &
               ylist(n)(1:izlen) == 'HBAS_SC' ) THEN
        DO j= 1, je
          DO i = 1, ie
            zvarlev(i,j,1) = 0.0_ireals
            klev = NINT( bas_con(i,j) )
            IF(klev > 0) THEN
              zvarlev(i,j,1) = hhl(i,j,klev)
            ENDIF
          ENDDO 
        ENDDO
      ELSEIF ( ylist(n)(1:izlen) == 'HZEROCL' ) THEN
        CALL calhzero( zvarlev(:,:,1), t(:,:,:,itl), hhl, vcoord%vert_coord,&
                       ie, je, ke, t0_melt)
      ELSEIF ( ylist(n)(1:izlen) == 'SNOWLMT' ) THEN
        CALL calsnowlmt( zvarlev(:,:,1), t(:,:,:,itl), pp(:,:,:,itl),      &
              p0(:,:,:), qv(:,:,:), hhl, vcoord%vert_coord, ie, je, ke,    &
              t0_melt, 1.3_ireals)
      ELSEIF ( ylist(n)(1:izlen) == 'CLCT_MOD' ) THEN
        CALL calclmod(zvarlev(:,:,1), clc_sgs, clc_con, p0hl, pi, ie, je, ke)
      ELSEIF ( ylist(n)(1:izlen) == 'CLDEPTH' ) THEN
        CALL calcldepth( zvarlev(:,:,1), clc_sgs, clc_con, dp0, ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQV' ) THEN
        CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qv(:,:,: ), ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQC' ) THEN
        CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qc(:,:,: ), ie, je, ke )
      ELSEIF ( ylist(n)(1:izlen) == 'TQR' ) THEN
        IF ( ASSOCIATED(qr) ) THEN
          CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qr(:,:,: ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQS' ) THEN
        IF ( ASSOCIATED(qs) ) THEN
          CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qs(:,:,: ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQG' ) THEN
        IF ( ASSOCIATED(qg) ) THEN
          CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qg(:,:,: ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
#ifdef TWOMOM_SB
      ELSEIF ( ylist(n)(1:izlen) == 'TQH' ) THEN
        IF ( itype_gscp >= 2000 ) THEN
          CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qh(:,:,: ), ie, je, ke )
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
#endif
      ELSEIF ( (ylist(n)(1:izlen) == 'ZTD') .OR.                     &
               (ylist(n)(1:izlen) == 'ZWD') .OR.                     &
               (ylist(n)(1:izlen) == 'ZHD') ) THEN
        IF (.NOT. lzenith) THEN
          CALL calztd( zenith_t(:,:), zenith_w(:,:), zenith_h(:,:),        &
                       zrho_itl, hhl, qv(:,:,: ), ps(:,:,itl ),                 &
                       t(:,:,ke,itl ), hsurf, rlat, pi, ie, je, ke)
          lzenith = .TRUE.
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TWATER' ) THEN
        zvarlev(:,:,2:ke+1) = qv (:,:,1:ke) + qc(:,:,1:ke)
        IF (ASSOCIATED(qi)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qi(:,:,1:ke)
        ENDIF
        IF (ASSOCIATED(qr)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qr(:,:,1:ke)
        ENDIF
        IF (ASSOCIATED(qs)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qs(:,:,1:ke)
        ENDIF
        IF (ASSOCIATED(qg)) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qg(:,:,1:ke)
        ENDIF
#ifdef TWOMOM_SB
         IF (itype_gscp >= 2000) THEN
          zvarlev(:,:,2:ke+1) = zvarlev(:,:,2:ke+1) + qh(:,:,1:ke)
        ENDIF
#endif       
        CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, zvarlev(:,:,2:ke+1), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'TQI') THEN
        IF (ASSOCIATED(qi)) THEN
          CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qi(:,:,:), ie,je,ke)
        ELSE
          zvarlev(:,:,1) = 0.0_ireals
        ENDIF
      ELSEIF ( ylist(n)(1:izlen) == 'TQV_ANAI') THEN
        ! for analysis increments, using rho(nnow) is sufficiently accurate
        CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qv_anai(:,:,:), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'TQC_ANAI') THEN
        CALL caliq( zvarlev(:,:,1), zrho_itl, hhl, qc_anai(:,:,:), ie,je,ke)
      ELSEIF ( ylist(n)(1:izlen) == 'PMSL_ANAI') THEN
        CALL calpmsl( zvarlev(:,:,1), ps(:,:,itl), t(:,:,ke,itl),          &
               rho0(:,:,ke), dp0(:,:,ke), hsurf, ie, je, g, r_d )
        zvarlev(:,:,1) = p_anai(:,:,ke) * zvarlev(:,:,1)                   &
                                        /(p0(:,:,ke) + pp(:,:,ke,itl))
      ELSEIF ( (ylist(n)(1:izlen) == 'AUMFL_S') .AND.     &
               (.NOT. lzrestart) .AND. outblock%luvmasspoint) THEN
        ! Exchange aumfl_s first
        kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
          (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar,                                &
          nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,    &
          20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yerrmsg,         &
          var(i1,i2,i3)%p2(:,:))
        DO j = 1, je
          DO i = 2, ie
            zvarlev(i,j,1) = 0.5 * (var(i1,i2,i3)%p2(i-1,j) +      &
                                    var(i1,i2,i3)%p2(i  ,j))
          ENDDO
          zvarlev(1,j,1) = zvarlev(2,j,1)
        ENDDO
      ELSEIF ( (ylist(n)(1:izlen) == 'AVMFL_S') .AND.     &
               (.NOT. lzrestart) .AND. outblock%luvmasspoint) THEN
        ! Exchange avmfl_s first
        kzdims(1:24)=(/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
        CALL exchg_boundaries                                                &
          (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,&
          ie, je, kzdims, jstartpar, jendpar,                                &
          nbl_exchg, nboundlines, my_cart_neigh, lperi_x, lperi_y, l2dim,    &
          20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yerrmsg,         &
          var(i1,i2,i3)%p2(:,:))
        DO j = 2, je
          DO i = 1, ie
            zvarlev(i,j,1) = 0.5 * (var(i1,i2,i3)%p2(i,j-1) +      &
                                    var(i1,i2,i3)%p2(i,j  ))
          ENDDO
        ENDDO
        zvarlev(:,1,1) = zvarlev(:,2,1)
      ELSEIF (ylist(n)(1:izlen) == 'DBZ_850') THEN
        IF (itype_gscp == 3) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,  &
               klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),            &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,       &
               qs(:,:,:)*zrho_itl, z_radar_850 = zvarlev(:,:,1) )
        ELSEIF (itype_gscp == 4) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,  &
               klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),            &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,       &
               qs(:,:,:)*zrho_itl, q_grau  = qg(:,:,:)*zrho_itl,                 &
               z_radar_850 = zvarlev(:,:,1) )
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp >= 100) THEN
          CALL radar_sb_ray (ie, je, ke, pi,                                 &
               klv850, my_cart_id, t(:,:,:,itl),                             &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl,                       &
               qi(:,:,:)*zrho_itl, qs(:,:,:)*zrho_itl,                       &
               qg(:,:,:)*zrho_itl, qh(:,:,:)*zrho_itl,                       &
               qnc(:,:,:)*zrho_itl, qnr(:,:,:)*zrho_itl,                     &
               qni(:,:,:)*zrho_itl, qns(:,:,:)*zrho_itl,                     &
               qng(:,:,:)*zrho_itl, qnh(:,:,:)*zrho_itl,                     &
               z_radar_850 = zvarlev(:,:,1) )
#endif
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'DBZ_CMAX') THEN
        IF (itype_gscp == 3) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,  &
               klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),            &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,       &
               qs(:,:,:)*zrho_itl, z_radar_cmax = zvarlev(:,:,1) )
        ELSEIF (itype_gscp == 4) THEN
          CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt,  &
               klv850, my_cart_id, itype_gscp, izdebug, t(:,:,:,itl),            &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,       &
               qs(:,:,:)*zrho_itl, q_grau  = qg(:,:,:)*zrho_itl,                 &
               z_radar_cmax = zvarlev(:,:,1) )
#ifdef TWOMOM_SB
        ELSEIF (itype_gscp >= 100) THEN
          CALL radar_sb_ray (ie, je, ke, pi,                                 &
               klv850, my_cart_id, t(:,:,:,itl),                             &
               qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl,                       &
               qi(:,:,:)*zrho_itl, qs(:,:,:)*zrho_itl,                       &
               qg(:,:,:)*zrho_itl, qh(:,:,:)*zrho_itl,                       &
               qnc(:,:,:)*zrho_itl, qnr(:,:,:)*zrho_itl,                     &
               qni(:,:,:)*zrho_itl, qns(:,:,:)*zrho_itl,                     &
               qng(:,:,:)*zrho_itl, qnh(:,:,:)*zrho_itl,                     &
               z_radar_850 = zvarlev(:,:,1) )
#endif
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'SWISS00') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                        &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             swiss00 = zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SWISS12') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                        &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             swiss12 = zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SI') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                        &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             si= zvarlev(:,:,1))
      ELSEIF (ylist(n)(1:izlen) == 'SLI') THEN
        CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                        &
             u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),            &
             p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w, b3, &
             b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,               &
             missing_value, izdebug, lwarn, ierrstat, yerrmsg,             &
             sli= zvarlev(:,:,1))
      ELSEIF ( (ylist(n)(1:izlen) == 'CAPE_MU')  .OR. &
               (ylist(n)(1:izlen) == 'CIN_MU') ) THEN
        IF (.NOT. lconvind_mu) THEN
           CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                     &
                u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),         &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w,  &
                b3, b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,        &
                missing_value, izdebug, lwarn, ierrstat, yerrmsg,          &
                cape_mu=zcape_mu, cin_mu=zcin_mu)
           lconvind_mu = .TRUE.
        ENDIF
      ELSEIF ( (ylist(n)(1:izlen) == 'CAPE_ML')  .OR. &
               (ylist(n)(1:izlen) == 'CIN_ML')   .OR. &
               (ylist(n)(1:izlen) == 'LCL_ML')   .OR. &
               (ylist(n)(1:izlen) == 'LFC_ML')   .OR. &
               (ylist(n)(1:izlen) == 'CAPE_3KM') ) THEN
        IF (.NOT. lconvind_ml) THEN
           CALL cal_conv_ind (t(:,:,:,itl), qv(:,:,:),                     &
                u(:,:,:,itl), v(:,:,:,itl),hsurf(:,:),ps(:,:,itl),         &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:), ie, je, ke, b1, b2w,  &
                b3, b4w, lh_v, cp_d, r_d, rdv, rvd_m_o, o_m_rdv, g,        &
                missing_value, izdebug, lwarn, ierrstat, yerrmsg,          &
                cape_ml=zcape_ml, cin_ml=zcin_ml, lcl_ml=zlcl_ml,          &
                lfc_ml=zlfc_ml, cape_3km=zcape_3km)
           lconvind_ml = .TRUE.
        ENDIF
      ELSEIF (ylist(n)(1:izlen) == 'HPBL') THEN
        IF (.NOT. l_brn) THEN
          kbot = 1
          ktop = ke
#ifndef MESSY
          CALL calc_bulk_richardson(zvarlev(:,:,kbot:ktop),t(:,:,:,itl),  &
               qv(:,:,:), u(:,:,:,itl), v(:,:,:,itl),                     &
               p0(:,:,:)+pp(:,:,:,itl), hsurf(:,:), ps(:,:,itl),          &
               t_2m(:,:), qv_2m(:,:),hhl(:,:,:),                          &
               ie, je, ke, cp_d, r_d,  rvd_m_o, g)
          l_brn = .TRUE.
          zbrn(:,:,:) = zvarlev(:,:,kbot:ktop)
#else
          CALL calc_bulk_richardson(zbrn(:,:,:),t(:,:,:,itl),             &
               qv(:,:,:), u(:,:,:,itl), v(:,:,:,itl),                     &
               p0(:,:,:)+pp(:,:,:,itl), hsurf(:,:), ps(:,:,itl),          &
               t_2m(:,:), qv_2m(:,:),hhl(:,:,:),                          &
               ie, je, ke, cp_d, r_d,  rvd_m_o, g)
          l_brn = .TRUE.
#endif
        ENDIF
        CALL  calc_pbl_brn(t(:,:,:,itl), qv(:,:,:),                        &
                p0(:,:,:)+pp(:,:,:,itl), hhl(:,:,:),hsurf(:,:),            &
                zbrn(:,:,:), ie, je, ke, cp_d, r_d, rvd_m_o, missing_value,&
                zvarlev(:,:,1))
#ifndef MESSY
        zvarlev(:,:,0)=zvarlev(:,:,1)
#endif
      ELSEIF (ylist(n)(1:izlen) == 'CEILING') THEN
        CALL calc_ceiling( zvarlev(:,:,1), clc_sgs, hhl, ie, je, ke )
      ELSEIF ( (ylist(n)(1:izlen) == 'SDI_1') .OR.          &
               (ylist(n)(1:izlen) == 'SDI_2') )  THEN
        IF (.NOT. l_calc_sdi_already_comp ) THEN
          CALL calc_sdi( zvarlev(:,:,1), zvarlev(:,:,2) )
          l_calc_sdi_already_comp = .TRUE.
        END IF
      ELSE
        IF (ASSOCIATED(var(i1,i2,i3)%p2)) THEN
          zvarlev(1:ie,1:je,1) = var(i1,i2,i3)%p2(1:ie,1:je)
        ELSE
          PRINT *, ' *** ERROR: Trying to output unassociated variable: ', &
                     ylist(n)(1:izlen)
        ENDIF
      ENDIF


      nzrecords   = nzrecords+1
      IF     ( ylist(n)(1:izlen) == 'ZTD' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zenith_t(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zenith_t(:,:), &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ELSEIF ( ylist(n)(1:izlen) == 'ZWD' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zenith_W(:,:)
#endif         
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zenith_w(:,:), &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ELSEIF ( ylist(n)(1:izlen) == 'ZHD' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zenith_h(:,:)
#endif         
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zenith_h(:,:), &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_MU' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zcape_mu(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zcape_mu(:,:), &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'CIN_MU' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zcin_mu(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zcin_mu(:,:),  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_ML' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zcape_ml(:,:)
#endif        
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zcape_ml(:,:), &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'CIN_ML' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zcin_ml(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zcin_ml(:,:),  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ELSEIF ( ylist(n)(1:izlen) == 'LCL_ML' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zlcl_ml(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zlcl_ml(:,:),  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'LFC_ML' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zlfc_ml(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zlfc_ml(:,:),  &
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'CAPE_3KM' ) THEN
#ifdef MESSY
         channeli%this%vars(n)%ptr(:,:,1,1) = zcape_3km(:,:)
#endif
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zcape_3km(:,:),&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ELSEIF ( ylist(n)(1:izlen) == 'SDI_1' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zvarlev(:,:,1),&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSEIF ( ylist(n)(1:izlen) == 'SDI_2' ) THEN
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zvarlev(:,:,2),&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)

      ELSE
        CALL output_data (nuedat,nzrecords, i1,i2,i3, 1, 1, zvarlev(:,:,1),&
                          outblock, .FALSE., yextension, slev(1),          &
                          lzrestart, ivar_id(n), n, iorg_data, izdebug)
      ENDIF

    END SELECT

  ENDDO write_loop

!------------------------------------------------------------------------------
! Section 4: Flush output buffers and close grib file
!------------------------------------------------------------------------------

#ifdef MESSY
   IF (l_COSMO_now .AND. L_BM_ORIG_OUTPUT)  THEN
#endif

  CALL output_data (nuedat, -1, -1,-1,-1, -1, -1, zvarlev(:,:,1),              &
                    outblock, .TRUE., yextension, 0.0_ireals, lzrestart, -1,   &
                    -1, iorg_data, izdebug)

!  IF (outblock%yform_write /= 'ncdf' .OR. nc_asyn_io < 1) THEN
  IF (outblock%yform_write == 'ncdf' .AND. (nc_asyn_io < 1 .OR. lzrestart)) THEN
    CALL close_file (nuedat, TRIM(outblock%yform_write), icomm_cart, my_cart_id, &
                     num_compute, lasync_io, idbg_level, yerrmsg, ierrstat)

    IF (ierrstat /= 0) THEN
      CALL model_abort (my_cart_id, 2034, yerrmsg, yroutine)
    ENDIF
  ENDIF

#ifdef MESSY
  ENDIF
#endif

  IF (outblock%yform_write == 'bina') THEN
    ! release the unit-number again
    CALL release_unit (nuedat)
  ENDIF

  ! Write a blank line to YUCHKDAT
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    WRITE (nuchkdat,'(A)') '   '
    WRITE (nuchkdat,'(A)') '   '
  ENDIF

  ! close file nuchkdat
  IF ( (outblock%lcheck) .AND. (my_cart_id == 0) ) THEN
    CLOSE (nuchkdat, STATUS='KEEP')
  ENDIF

  ! Deallocate arrays for IO
  DEALLOCATE (iblock, ibmap, ds_real, ds_grib, dsup, ymessage)

  ! Deallocate arrays for densities
  DEALLOCATE ( zrho_itl, zqrs_itl )

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE organize_output

!==============================================================================
!+ distributes records to PEs and packs them into grib format
!------------------------------------------------------------------------------

SUBROUTINE output_data (nuedat, irec, i1,i2,i3, k, klevels, array2d_real,   &
                        outblock, lflush, yextension, slev, lrestart,       &
                        incdf_var_id, var_index, my_orgdata, idebug)

!------------------------------------------------------------------------------
!
! Description:
!  output_data distributes records to the PEs for packing these into 
!  grib format. First, the records are only gathered from all PEs and
!  every PE stores one record into the variable procarray_xxx. Only if every
!  PE has got a record, the data are packed and written to disk (in the
!  routine "write_xxxx"). If some PEs have got no record because no more
!  records are left, the output buffers (variable procarray_xx) are "flushed".
!
! Method:
!  output_data is called for every record that is processed. The PE that
!  gets a special record, saves the characteristics of this record for the
!  output step later on.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT (IN)    ::    &
  nuedat,  & ! descriptor of the grib file
  irec,    & ! number of record to be processed
  i1,i2,i3,& ! location of the variable in the LM variable table
  k,       & ! number of the actual level
  klevels, & ! number of levels this variable has
  idebug     ! for verbosity debug output

TYPE(pp_nl),              INTENT(IN)     ::    &
  outblock        ! pointer to the namelist group

LOGICAL                 , INTENT (IN)    ::    &
  lflush,  & ! for flushing the output buffers
  lrestart   ! whether restart-files are written or not

CHARACTER (LEN= 1)      , INTENT (IN)    ::    &
  yextension ! to check which output list is processed

REAL (KIND=ireals)      , INTENT (IN)    ::    &
  slev       ! level for vertical interpolated fields
             ! has to be present for p- and z-levels

! Array arguments with intent(inout):
REAL (KIND=ireals)      , INTENT (INOUT) ::    &
  array2d_real (ie,je)      ! values of the variable to be processed

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  incdf_var_id              ! NetCDF-ID of each variable in the output list
                            ! only PE 0 has a reasonable value here

INTEGER (KIND=iintegers), INTENT(IN)     ::    &
  var_index       ! index of the variable within list of variables for
                  ! current gribout

INTEGER (KIND=iintegers), INTENT(INOUT)  ::    &
  my_orgdata(3,0:num_compute-1) ! necessary only for PE 0 to save information

!------------------------------------------------------------------------------
!
! Local variables:

! Local scalars:
INTEGER (KIND=intgribf)         :: iz_ps=1, izerrf
INTEGER (KIND=intgribf)         :: ierrf, tri, my_k_f, my_iee
INTEGER (KIND=iintegers)        :: npe, iz_lfd, izstat, izerror, izlen, i,j, &
                                   irecord_len, implcode, nlastout,          &
                                   ij_out, ij_ful, iz_lfa,                   &
                                   istartlon, jstartlat, iendlon, jendlat,   &
                                   nzjulianday
INTEGER (KIND=iintegers)        :: irec2, npe2, icy, imm, idd, ihh, imi

INTEGER (KIND=int_ga)           :: irecord_lga      ! length of grib record

REAL    (KIND=ireals)           :: zavgfactor, z2, ztgranul, zacthour
REAL    (KIND=ireals)           :: zbias, zgribfactor, array2dreal(ie_max, je_max)
REAL    (KIND=irealgrib)        :: array2d_grib(ie_max,je_max),              &
                                   ds_out(ie_tot*je_tot), undefsub

REAL    (KIND=ireals)           :: ds_api(ie_tot*je_tot)

CHARACTER (LEN=25)              :: yroutine
CHARACTER (LEN=80)              :: yerrmsg
CHARACTER (LEN=clen)            :: yzname 
CHARACTER (LEN= 8)              :: yzdate
CHARACTER (LEN=10)              :: yztime
CHARACTER (LEN=14)              :: yzdatact
CHARACTER (LEN=28)              :: yzdat2


INTEGER (KIND=iintegers), SAVE  :: my_i1, my_i2, my_i3, my_k, my_irec
LOGICAL                 , SAVE  :: loutput
REAL    (KIND=ireals)   , SAVE  :: rmy_slev
INTEGER (KIND=iintegers), SAVE  :: low_irec, high_irec

! Local arrays:

#ifdef NETCDF
INTEGER (KIND=iintegers), SAVE  ::       &
  my_orgdata_asyn( nc_orgmdata_length ), & ! array with metadata of netcdf file
  my_vardata_asyn( nc_varmdata_length  )   ! array with metadata of variable being sent
#endif

REAL (KIND=ireals), ALLOCATABLE, SAVE       :: &
  procarray_real  (:,:,:),  &
  procarray2d_real(:,:,:)

REAL (KIND=irealgrib), ALLOCATABLE, SAVE    :: &
  procarray_grib  (:,:,:),  &
  procarray2d_grib(:,:,:)

INTEGER, DIMENSION(0:num_compute-1) :: sendcnts, sdispls, recvcnts, rdispls

#ifdef GRIBAPI
INTEGER (KIND=kindOfSize) :: ibyte_size_out
#endif

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
!  Section 1: Initializations
!------------------------------------------------------------------------------

  IF (idebug > 15) THEN
    WRITE (*,'(A)') '  src_output: entering output_data'
  ENDIF

  yroutine = 'output_data'
  yzname   = '          '
  izerror  = 0_iintegers
  izerrf   = 0_intgribf
  
  iz_lfd   = INT (lfd, iintegers)
  iz_lfa   = INT (lfa, iintegers)

  ! set undef values
  IF (outblock%yform_write /= 'ncdf') THEN
    undef     = REAL(undefgrib, ireals)
    undefsub  = undefgrib
  ELSE
    undef     = REAL(undefncdf, ireals)
    undefsub  = undefncdf
  ENDIF

  ! security check on itype_gather
  IF (itype_gather /= 1 .AND. itype_gather /= 2) THEN
    WRITE(*,*) 'ERROR: INVALID itype_gather IN output_data '
  ENDIF

  ! allocate buffer memory
  IF (lrestart) THEN
    IF ( .NOT. ALLOCATED(procarray_real) ) THEN
      ALLOCATE(procarray_real(ie_max, je_max, num_compute), STAT=izstat)
      izerror = izerror + izstat
    ENDIF
    IF ( (itype_gather == 2) .AND. (.NOT. ALLOCATED(procarray2d_real)) ) THEN
      ALLOCATE(procarray2d_real(ie_max, je_max, num_compute), STAT=izstat)
      izerror = izerror + izstat
    ENDIF
  ELSE
    IF ( .NOT. ALLOCATED(procarray_grib) ) THEN
      ALLOCATE(procarray_grib(ie_max, je_max, num_compute), STAT=izstat)
      izerror = izerror + izstat
    ENDIF
    IF ( (itype_gather == 2) .AND. (.NOT. ALLOCATED(procarray2d_grib)) ) THEN
      ALLOCATE(procarray2d_grib(ie_max, je_max, num_compute), STAT=izstat)
      izerror = izerror + izstat
    ENDIF
  ENDIF
  IF (izerror /= 0) THEN
  WRITE(*,*) 'ERROR: PROBLEM ALLOCATING MEMORY IN output_data'
  ENDIF

  ! setup house keeping data for all2allv gathering
  IF (irec == 1) THEN
    my_orgdata(:,:) = 0_iintegers
    loutput         = .FALSE.
    low_irec  = 1
    high_irec = 1
  ENDIF

  IF ( low_irec == 0 .AND. .NOT. lflush ) low_irec = irec  ! if just flushed
  IF ( .NOT. lflush ) high_irec = MAX(high_irec,irec)

  IF (.NOT. lrestart) THEN
    IF (outblock%nextstep == 1) THEN
      nlastout = 0
      ! Adaptation, if summation and meanvalues are done between output steps
      ! and the first output step is not 0
      IF (outblock%ngrib( outblock%nextstep) > 0 ) THEN
        nlastout = outblock%ngrib( outblock%nextstep ) -                             &
              (outblock%ngrib( outblock%nextstep+1 )-outblock%ngrib( outblock%nextstep ))
              !US should it be:  nlastout = outblock%ngrib( outblock%nextstep) ???????
      ENDIF
    ELSE
      nlastout = outblock%ngrib( outblock%nextstep-1 )
    ENDIF
  ELSE
    nlastout = outblock%nextstep - nhour_restart(3) * NINT (3600.0_ireals / dt)
  ENDIF

!------------------------------------------------------------------------------
!  Section 2: If this is not a call to only flush the output data,
!             gather the field on the next free PE
!------------------------------------------------------------------------------

  IF (.NOT. lflush) THEN

    IF (idebug > 15) THEN
      WRITE (*,'(A)') '  src_output: gather field on next free task'
    ENDIF

    ! Get the number of the PE to deal with that slice
    npe     = MOD(irec-1,num_compute)

#ifdef NETCDF
    my_orgdata_asyn(1) = cur_outstep_idx
    SELECT CASE ( yextension )
    CASE('c')
      my_orgdata_asyn(2) = 0
    CASE(' ')
      my_orgdata_asyn(2) = 1
    CASE('p')
      my_orgdata_asyn(2) = 2
    CASE('z')
      my_orgdata_asyn(2) = 3
    CASE('s')
      my_orgdata_asyn(2) = 4
    END SELECT

    my_orgdata_asyn(3) = cur_gribout_idx
#endif

    ! Save necessary values for later processing
    IF (my_cart_id == npe ) THEN
      my_irec   = irec
      my_i1     = i1
      my_i2     = i2
      my_i3     = i3
      my_k      = k
      rmy_slev  = slev
      loutput   = .TRUE.
#ifdef NETCDF
      my_vardata_asyn(1) = var_index
      my_vardata_asyn(2) = k
      my_vardata_asyn(3) = klevels
#endif
    ENDIF

#ifdef NETCDF
    ! Processor 0 has to save some organizational data
    IF (my_cart_id == 0) THEN
      my_orgdata(1,npe) = k
      my_orgdata(2,npe) = klevels
      my_orgdata(3,npe) = incdf_var_id
    ENDIF
#endif

    ! scale the fields with time range indicator 3 
    !US-ACHTUNG: this has to be adapted to GRIB2 somehow!!!
    IF ((.NOT. lrestart) .AND. (var(i1,i2,i3)%ntri == 3)) THEN
      IF (ntstep == 0) THEN
        zavgfactor = 1.0_ireals
      ELSE
        IF (lbdclim) THEN
          ! averaging is done between output steps
          zavgfactor = 1.0_ireals / REAL  (ntstep - nlastout, ireals)
        ELSE
          ! averaging is done between beginning of forecast and
          ! actual output step
          IF (l2tls) THEN
            zavgfactor = 1.0_ireals / REAL  (ntstep+1, ireals)
          ELSE
            zavgfactor = 1.0_ireals / REAL  (ntstep, ireals)
          ENDIF
        ENDIF
      ENDIF
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! handle the mean wind speed in single sector classes, only account for 
          !   time period where wind was really active
          IF ( ( TRIM(var(i1,i2,i3)%name) == 'VABS_10M_SECAV' ) .AND. ( wdirgeo_10m_freq(i,j,k) .GT. 0.0_iintegers ) ) THEN
            array2d_real(i,j) = array2d_real(i,j) / wdirgeo_10m_freq(i,j,k)
          ELSE
            array2d_real(i,j) = array2d_real(i,j) * zavgfactor
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! Scale field with factor and bias (only for grib)
    ! transform it to single precision (for grib and Netcdf)
    zbias       = var(i1,i2,i3)%bias
    zgribfactor = var(i1,i2,i3)%factor

    IF (outblock%yform_write == 'grb1' .OR. outblock%yform_write(1:3) == 'api') THEN
      ! Check for undefined values in case of NetCDF-Input and Grib output
      IF (yform_read == 'ncdf') THEN
        ! this is a very pragmatic solution, which might not be satisfactory
        DO j = jstartpar, jendpar
          DO i = istartpar, iendpar
            IF (array2d_real(i,j) == REAL(undefncdf, ireals)) THEN
              array2d_real(i,j) = 0.0_ireals
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! Do an additional clipping, if values are too small
          IF (ABS(array2d_real(i,j)) < 1.0E-15_ireals) THEN
            array2d_grib(i,j) =                                             &
                REAL (((0.0_ireals        + zbias) * zgribfactor), irealgrib)
          ELSE
            array2d_grib(i,j) =                                             &
                REAL (((array2d_real(i,j) + zbias) * zgribfactor), irealgrib)
          ENDIF
        ENDDO
      ENDDO
    ELSEIF (outblock%yform_write == 'ncdf') THEN
      !BEGIN KOCHANGE
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! No scaling for NetCDF
          IF (.NOT.lrestart) THEN
            array2d_grib(i,j) = REAL (array2d_real(i,j), irealgrib)
          ELSE
            ! Do not truncate restart values!
            array2dreal(i,j) = array2d_real(i,j)
          ENDIF
        ENDDO
      ENDDO
      !END KOCHANGE

      ! take care of special fields
      !BEGIN KOCHANGE
      IF (.NOT.lrestart) THEN
        IF (var(i1,i2,i3)%lsm == 'l') THEN
          WHERE (.NOT. llandmask(1:ie,1:je)) array2d_grib(1:ie,1:je) = undefncdf
        ENDIF
        IF (var(i1,i2,i3)%lsm == 's') THEN
          WHERE (llandmask(1:ie,1:je)) array2d_grib(1:ie,1:je) = undefncdf
        ENDIF
        IF (var(i1,i2,i3)%lsm == 'i') THEN
          WHERE (fr_lake(1:ie,1:je) <= 0.5_ireals) array2d_grib(1:ie,1:je) = undefncdf
        ENDIF
        IF (TRIM(var(i1,i2,i3)%name) == 'SNOWLMT' .OR. &
            TRIM(var(i1,i2,i3)%name) == 'HZEROCL' .OR. &
!--br 06.08.19 -->
            TRIM(var(i1,i2,i3)%name) == 'CIN_MU' .OR. &
            TRIM(var(i1,i2,i3)%name) == 'CIN_ML' .OR. &
            TRIM(var(i1,i2,i3)%name) == 'CAPU_MU' .OR. &
            TRIM(var(i1,i2,i3)%name) == 'CAPE_ML') THEN
!--br 06.08.19 <--
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
!--br 06.08.19 -->
!              IF (array2d_grib(i,j) == -999.0_irealgrib) THEN
              IF (array2d_grib(i,j) == -999.0_irealgrib .OR. array2d_grib(i,j) == -999.9_irealgrib) THEN
!--br 06.08.19 <--
                array2d_grib(i,j) = undefncdf
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (TRIM(var(i1,i2,i3)%name) == 'HBAS_CON' .OR. &
            TRIM(var(i1,i2,i3)%name) == 'HTOP_CON') THEN
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              IF (array2d_grib(i,j) == 0.0_irealgrib) THEN
                array2d_grib(i,j) = undefncdf
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (TRIM(var(i1,i2,i3)%name) == 'T_SNOW') THEN
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              IF (w_snow(i,j,nnow) == 0.0_ireals) THEN
                array2d_grib(i,j) = undefncdf
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        IF (TRIM(var(i1,i2,i3)%name) == 'Z0') THEN
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              array2d_grib(i,j) = array2d_grib(i,j) * 0.10197_irealgrib
            ENDDO
          ENDDO
        ENDIF
!<--  rp 2018-05-30      
        !RONNY, take care of wind directions which did not occur
        !IF (TRIM(var(i1,i2,i3)%name) == 'VABS_10M_SECAV' .OR. TRIM(var(i1,i2,i3)%name) == 'VABS_10M_SECMAX') THEN
        ! DO j = jstartpar, jendpar
        !   DO i = istartpar, iendpar
        !      IF (array2d_grib(i,j) == 0.0_irealgrib) THEN
        !        array2d_grib(i,j) = undefncdf
        !      ENDIF
        !    ENDDO
        !  ENDDO
        !ENDIF
        ! END RONNY
        IF ( TRIM(var(i1,i2,i3)%name) == 'WDIRGEO_10M' ) THEN
          DO j = jstartpar, jendpar
            DO i = istartpar, iendpar
              IF ( array2d_grib(i,j) == undefgrib ) THEN
                array2d_grib(i,j) = undefncdf
              ENDIF
            ENDDO
          ENDDO
        ENDIF

      ENDIF !lrestart

    ELSEIF (outblock%yform_write == 'bina') THEN
      DO j = jstartpar, jendpar
        DO i = istartpar, iendpar
          ! no scaling and converting in case of restart files,
          ! because this is not reproducible!
          array2dreal(i,j) = array2d_real(i,j)
        ENDDO
      ENDDO
    ENDIF

    IF ( ldebug_io .AND. (idbg_level >= 10) .AND. (my_cart_id == npe)) THEN
      PRINT *, ' src_output: gathering  ', TRIM(var(my_i1, my_i2, my_i3)%name), my_k, &
               ' to PE ', npe
    ENDIF

    ! Save data if num_compute=1 or for all2allv gathering
    IF (num_compute == 1) THEN
      ! no need to gather data, simply save
      IF (lrestart) THEN
        procarray_real(1:ie,1:je,1) = array2dreal(:,:)
      ELSE
        procarray_grib(:,:,1) = array2d_grib(:,:)
      ENDIF
    ELSE
      ! for all2allv gathering we need to store data
      IF (itype_gather == 2) THEN
        IF (lrestart) THEN
          procarray2d_real(:,:,npe+1) = array2dreal(:,:)
        ELSE
          procarray2d_grib(:,:,npe+1) = array2d_grib(:,:)
        ENDIF
      ENDIF
    ENDIF

  ENDIF

  IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

!------------------------------------------------------------------------------
!  Section 3: Gather full vertical levels onto compute PEs
!------------------------------------------------------------------------------

  IF (num_compute > 1) THEN

    ! Gather the data using MPI_GATHER for each vertical level individually
    IF ( (itype_gather == 1) .AND. (.NOT. lflush) ) THEN
      IF (.NOT. lrestart) THEN
        CALL gather_values (array2d_grib, procarray_grib, ie_max, je_max, &
               num_compute, imp_grib, npe, icomm_cart, yerrmsg, izerror)
      ELSE
        CALL gather_values (array2dreal, procarray_real, ie_max, je_max,     &
               num_compute, imp_reals, npe, icomm_cart, yerrmsg, izerror)
      ENDIF
    ENDIF

    ! Gather the data using MPI_ALL2ALLV for a maximum of num_compute vertical levels
    IF ( (itype_gather == 2) .AND. &
         (lflush .OR. MOD(irec-1,num_compute) == num_compute-1) ) THEN

      IF (low_irec > 0 .AND. high_irec >= low_irec ) THEN

        ! Initialize send counts and receive counts
        sendcnts(:) = 0
        sdispls(:) = 0
        recvcnts(:) = 0
        rdispls(:) = 0

        ! Setup receive counts if this PE should receive a field
        IF (loutput) THEN
          DO npe2 = 0, num_compute-1
            recvcnts(npe2) = ie_max*je_max
            rdispls(npe2) = npe2*ie_max*je_max
          ENDDO
        ENDIF

        ! Setup send counts for all PEs
        DO irec2 = low_irec, high_irec
          npe2 = MOD(irec2-1, num_compute)
          sendcnts(npe2) = ie_max*je_max
          sdispls(npe2) = npe2*ie_max*je_max
        ENDDO

        ! Do all2allv for gather the data on the receiving PEs
        IF (lrestart) THEN
          CALL MPI_ALLTOALLV(procarray2d_real,sendcnts,sdispls,  &
                             imp_reals,procarray_real,recvcnts,  &
                             rdispls, imp_reals,icomm_cart,izerror)
        ELSE
          CALL MPI_ALLTOALLV(procarray2d_grib,sendcnts,sdispls,  &
                             imp_grib,procarray_grib,recvcnts,   &
                             rdispls, imp_grib,icomm_cart,izerror)
        ENDIF

        ! Reset record counters
        low_irec = 0_iintegers
        high_irec = -1_iintegers

      ENDIF

    ENDIF

    IF (ltime) CALL get_timings (i_gather_data, ntstep, dt, izerror)

  ENDIF

!------------------------------------------------------------------------------
!  Section 4: If lflush is .TRUE. or all PEs have gotten data, do the output
!------------------------------------------------------------------------------

  IF ( lflush .OR. MOD(irec-1,num_compute) == num_compute-1) THEN

    IF (idebug > 15) THEN
      WRITE (*,'(A)') '  src_output: flush the buffers'
    ENDIF

   !---------------------------------------------------------------------------
   !  Section 4.1: All PEs that have gotten data must combine the subarrays
   !               and convert data to GRIB
   !---------------------------------------------------------------------------

    IF ( loutput ) THEN

      IF (outblock%yform_write /= 'ncdf') THEN

        IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

        IF (idebug > 15) THEN
          WRITE (*,'(A)') '  src_output: complete grib meta data'
        ENDIF

        IF (outblock%yform_write == 'grb1') THEN
          ! set the namelist-dependent idims_out values
          idims_out( 7) = outblock%ie_out_tot * outblock%je_out_tot
          idims_out(15) = outblock%ie_out_tot * outblock%je_out_tot
          idims_out(17) = outblock%ie_out_tot * outblock%je_out_tot

          ! gridpoints, simple packing, floating point data
          ibds(2)   = 0

          ! nrbit, number of bits
          ibds(5)   = outblock%nrbit

          ! no bitmap
          ibms(3)   = -2
        ENDIF

#ifdef GRIBAPI
        IF (outblock%yform_write(1:3) == 'api') THEN
          CALL grib_clone(outblock%igribapi_id, igribid, izerrf)
          IF (izerrf /= GRIB_SUCCESS) THEN
            PRINT *,   ' *** Error in grib_clone: from outblock sample ', izerrf
          ENDIF
        ENDIF
#endif

        CALL make_grib_grid (outblock, igribid, my_i1,my_i2,my_i3, yextension, idebug)

        ! Determine the reference time
        IF (outblock%lanalysis) THEN
          ! When nudging is active the output fields are treated as analyses
          ! and the reference time is the actual forecast time
          ! In principle we could use yakdat1 to get the entries for grib meta data

          ! But when using a flexible dt, the hour has eventually to be updated
          ! to the nearest output step (which should be a multiple of 0.25 h = 900.0 s

          ztgranul = 900.0_ireals
          z2 = (REAL(ntstep, ireals) * dt) / ztgranul
          IF (ABS(NINT(z2) - z2) > 1E-5_ireals) THEN
            ! determine date again with time step ztgranul and number of steps 
            ! necessary (z2) to reach the same forecast time that wie have now.
            CALL get_utc_date(NINT(z2), ydate_ini, ztgranul, itype_calendar, yzdatact, &
                              yzdat2, nzjulianday, zacthour)
          ELSE
            CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, yzdatact,    &
                              yzdat2, nzjulianday, zacthour)
          ENDIF
        ELSE
          yzdatact(1:14) = ydate_ini(1:14)
        ENDIF

        CALL make_grib_product (outblock, igribid, my_i1,my_i2,my_i3, my_k, yzdatact,      &
                        nlastout, yextension, rmy_slev, lrestart, idebug)

      ENDIF

      IF (ltime) CALL get_timings (i_meta_data_w, ntstep, dt, izerror)

      IF (.NOT. lrestart) THEN
        ! combine the subarrays in the correct order
        CALL combine_subarrays (procarray_grib, ds_grib)

        IF (ltime) CALL get_timings (i_gather_data, ntstep, dt, izerror)

        ! Extra smoothing of fi and pmsl in mountaineous regions:
        IF (outblock%l_fi_pmsl_smooth) THEN
          IF (     (var(my_i1,my_i2,my_i3)%name == 'PMSL     ')             &
              .OR. (var(my_i1,my_i2,my_i3)%name == 'PMSL_ANAI')) THEN
            CALL smooth_pmsl (ds_grib, hsurf_tot, ie_tot, je_tot )
          ENDIF
          IF(var(my_i1,my_i2,my_i3)%name == 'FI      ') THEN
            CALL smooth_geopot (ds_grib, hsurf_tot, ie_tot, je_tot )
          ENDIF
        ENDIF

        ! Possibility for independently smoothing of PMSL, 
        ! decoupled from switch outblock%l_z_filter:
        ! (Note: PMSL is a variable in the group of model-levels!)
        IF (outblock%l_pmsl_filter) THEN
          IF (     (var(my_i1,my_i2,my_i3)%name == 'PMSL     ')             &
              .OR. (var(my_i1,my_i2,my_i3)%name == 'PMSL_ANAI')) THEN
            CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
          ENDIF
        ENDIF

        ! Possibility for independently smoothing of interplated geopotential FI, 
        ! decoupled from general smoothing or not on p- or z-levels:
        IF (outblock%l_fi_filter) THEN
          IF(var(my_i1,my_i2,my_i3)%name == 'FI      ') THEN
            IF ( .NOT.(yextension == 'p' .AND. outblock%l_p_filter) .AND. &
                 .NOT.(yextension == 'z' .AND. outblock%l_z_filter) ) THEN
              CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
              ! .. Otherwise smoothing will be done a few lines below
            END IF
          ENDIF
        ENDIF

        ! Apply a digital smoother for selected fields
        IF(yextension == 'p' .AND. outblock%l_p_filter) THEN
          CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
        ENDIF
        IF(yextension == 'z' .AND. outblock%l_z_filter) THEN
          CALL smoother(ds_grib, ie_tot, je_tot, 4,20 )
        ENDIF

        ! Limit certain variables to the allowed range (again)
        SELECT CASE (var(my_i1,my_i2,my_i3)%name)
        CASE ('RELHUM    ')
          ds_grib(:) = MAX (0.0_irealgrib, MIN(100.0_irealgrib,ds_grib(:)))
        CASE ('QV        ','QC        ','QI        ',    &
              'QR        ','QS        ','QG        ',    &
              'QH        ','QNCLOUD   ','QNICE     ',    &
              'QNRAIN    ','QNSNOW    ','QNGRAUPEL ','QNHAIL    ' )
          ds_grib(:) = MAX (0.0_irealgrib, ds_grib(:))
        END SELECT

        ! Now cut out the proper (sub-)field which was chosen with
        ! slon, slat, elon, elat
        ij_out = 0
        DO j = outblock%j_out_start, outblock%j_out_end
          DO i = outblock%i_out_start, outblock%i_out_end
            ij_out = ij_out+1
            ij_ful = (j-1) * ie_tot + i
            ds_out(ij_out) = ds_grib(ij_ful)
            ds_api(ij_out) = REAL(ds_grib(ij_ful),ireals)
          ENDDO
        ENDDO

        IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

        IF (outblock%yform_write == 'grb1') THEN
#ifdef GRIBDWD
          ! degrib the level
          CALL grbex1(idwdednr, iz_ps, undefgrib, ndims, idims_out, ipds_out, &
                   igds_out, ibms, ibds, ibmap, dsup, ds_out, iblock, ierrf)
          IF (ierrf /= 0) THEN
            yerrmsg = 'error in grbex1'
            CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
          ENDIF

          ! length of GRIB record in bytes
          irecord_len = idims_out(19)
          irecord_lga = INT (idims_out(19), int_ga)
#endif

#ifdef GRIBAPI
      ELSEIF (outblock%yform_write == 'api1') THEN
        ! Set bitmap explicit to "no" again and reduce the size of ibmap
        ! (that seems to be a dummy only)
        CALL grib_set (igribid, 'bitmapPresent',      0)

        ! Set these values again explicit, because they could be overwritten
        ! by input data
        !  ?? ibds (2)     = 0
        CALL grib_set (igribid, 'bitsPerValue',   outblock%nrbit)

        CALL grib_get_message_size(igribid, ibyte_size_out, izerrf)
!ds_out is irealgrib. Or do we need ireals???
        CALL grib_set (igribid, 'values',         ds_api(:), izerrf)
        IF (izerrf /= GRIB_SUCCESS) THEN
          yerrmsg = 'error in grib_set: values'
          CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
        ENDIF

        ! length of GRIB record in bytes
        CALL grib_get_message_size(igribid, ibyte_size_out, izerrf)

        IF (ibyte_size_out <= lfa) THEN
          CALL grib_get(igribid, 'totalLength', irecord_lga)
          CALL grib_copy_message(igribid,  ymessage)
        ELSE
          yerrmsg = 'error with message length: ymessage too small: '
          CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
        ENDIF

      ELSEIF (outblock%yform_write == 'api2') THEN

        ! Set these values explicit
        CALL grib_set (igribid, 'bitsPerValue',   outblock%nrbit)

        CALL grib_get_message_size(igribid, ibyte_size_out, izerrf)
        CALL grib_set (igribid, 'values',         ds_api(:), izerrf)
        IF (izerrf /= GRIB_SUCCESS) THEN
          yerrmsg = 'error in grib_set: values'
          CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
        ENDIF

        ! length of GRIB record in bytes
        CALL grib_get_message_size(igribid, ibyte_size_out, izerrf)

!print *, 'message size:   ', ibyte_size_out, lfa

        IF (ibyte_size_out <= lfa) THEN
          CALL grib_get(igribid, 'totalLength', irecord_lga)
          CALL grib_copy_message(igribid,  ymessage)
        ELSE
          yerrmsg = 'error with message length: ymessage too small: '
          CALL model_abort (my_cart_id, 2022, yerrmsg, yroutine)
        ENDIF
#endif

        ELSE
 
          ! length of netcdf record in words
          irecord_len = outblock%ie_out_tot * outblock%je_out_tot
        ENDIF
        IF (ltime) CALL get_timings (i_meta_data_w, ntstep, dt, izerror)

      ELSE
        ! this is for restart-output
        CALL combine_subarrays (procarray_real, ds_real)
        IF (ltime) CALL get_timings (i_gather_data, ntstep, dt, izerror)

        ! length of restart record in words
        irecord_len = outblock%ie_out_tot * outblock%je_out_tot
      ENDIF

    ELSE

      irecord_len = 0_iintegers
      irecord_lga = 0_int_ga
!
!HJP Begin 2015-08-24
! re-initialization of field ds_real
!  necessary to avoid model crashes in case of usage of large number of
!  processors
      ds_real(:)  = 0.0_ireals
!HJP END   2015-08-24
!

    ENDIF

    !---------------------------------------------------------------------------
    !  Section 4.2: Check data and write output to disk
    !---------------------------------------------------------------------------

    ! check the data, if wanted
    IF (outblock%lcheck .AND. (.NOT. lrestart)) THEN
      my_k_f  = INT (my_k,  intgribf)
      my_iee  = INT (my_i2, intgribf)
      IF (.NOT. loutput) THEN
        my_i1 = 1; my_i2 = 1; my_i3 = 1           ! for safety reasons
      ENDIF
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1,               &
               outblock%i_out_start, outblock%i_out_end,                    &
               outblock%j_out_start, outblock%j_out_end, 1, 1,              &
               undefsub,  var(my_i1,my_i2,my_i3)%name,                      &
               my_iee, my_k_f, loutput, nuchkdat, num_compute,              &
               icomm_cart, my_cart_id, yerrmsg, izerror)
    ELSEIF (lrestart) THEN
    !BEGIN KOCHANGE T_SO
!--br 23.01.19
      IF (.NOT. loutput) THEN
        my_i1 = 1; my_i2 = 1; my_i3 = 1           ! for safety reasons
      ENDIF
      IF ((var(my_i1,my_i2,my_i3)%name == 'T_SO      ') .AND. &
          (outblock%yform_write == 'ncdf')) THEN
        my_k_f = INT(my_k-1, intgribf)
      ELSE
      my_k_f  = INT (my_k,  intgribf)
      ENDIF
      !KOCHANGE END T_SO
!--br 23.01.19 end
      my_iee  = INT (my_i2, intgribf)
!
      ds_grib(:) = REAL (ds_real(:), irealgrib)
      CALL check_record (ds_grib, 1, ie_tot, 1, je_tot, 1, 1, 1, ie_tot, 1, &
               je_tot, 1, 1, undefsub,  var(my_i1,my_i2,my_i3)%name,        &
               my_iee, my_k_f, loutput, nuchkdat, num_compute,              &
               icomm_cart, my_cart_id, yerrmsg, izerror)
    ENDIF

    IF (ltime) CALL get_timings (i_computations_O, ntstep, dt, izerror)

#ifdef MESSY
   IF (l_COSMO_now .AND. L_BM_ORIG_OUTPUT)  THEN
#endif

    SELECT CASE (outblock%yform_write)

    CASE ('grb1')

#ifdef GRIBDWD
      CALL write_grib    (nuedat, iblock, irecord_lga, iz_lfd, icomm_cart, &
                          num_compute, my_cart_id, lflush, lasync_io,      &
                          ltime_barrier, yerrmsg, izerror)
#endif

    CASE ('api1','api2')

#ifdef GRIBAPI
      CALL write_gribapi (nuedat, ymessage, irecord_lga, iz_lfa, num_compute, &
                          icomm_cart, my_cart_id, lflush, lasync_io,          &
                          ltime_barrier, yerrmsg, izerror)
#endif

    CASE ('ncdf')

#ifdef NETCDF
!      IF ( nc_asyn_io > 0 ) THEN
      IF ( nc_asyn_io > 0 .AND. (.NOT. lrestart)) THEN
        CALL send_asyn_io ( ds_out, irecord_len, my_orgdata_asyn,          &
              my_vardata_asyn, outblock,outblock%ie_out_tot,               &
              outblock%je_out_tot, lflush,yerrmsg, izerror )
      ELSE
      !BEGIN KOCHANGE
        IF (lrestart) THEN
        CALL write_netcdf  (nuedat, outblock%ie_out_tot,                     &
                            outblock%je_out_tot, irecord_len, my_orgdata,    &
                            icomm_cart, my_cart_id, num_compute, imp_reals,  &
                            lasync_io, yerrmsg, izerror, data_real=ds_real)
      ELSE
        CALL write_netcdf  (nuedat, outblock%ie_out_tot,                     &
                            outblock%je_out_tot, irecord_len, my_orgdata,    &
                            icomm_cart, my_cart_id, num_compute, imp_grib,   &
                            lasync_io, yerrmsg, izerror, data_grib=ds_out)
        ENDIF
      ENDIF
#endif
    CASE ('bina')

      CALL write_restart (nuedat, ds_real, ie_tot, je_tot, irecord_len,    &
                          ipds_out, npds, igds_out, ngds, icomm_cart,      &
                          my_cart_id, num_compute, imp_reals, lasync_io,   &
                          yerrmsg, izerror)

    END SELECT

    IF (izerror /= 0) THEN
      CALL model_abort (my_cart_id, 2035, yerrmsg, yroutine)
    ENDIF

    IF (ltime) CALL get_timings (i_write_data, ntstep, dt, izerror)

    IF ( ldebug_io .AND. (idbg_level >= 10)) THEN
      PRINT *, ' src_output: saved ', TRIM(var(my_i1, my_i2, my_i3)%name), my_k, ' to disk'
    ENDIF

#ifdef MESSY
  ENDIF
#endif

    ! Reset organizational variables
#ifdef GRIBAPI
    IF (outblock%yform_write(1:3) == 'api') THEN
      CALL grib_release (igribid)
    ENDIF
#endif
    my_orgdata(:,:) = 0_iintegers
    loutput         = .FALSE.

    ! free memory
    IF (lrestart) THEN
      IF (ALLOCATED(procarray_real)) DEALLOCATE(procarray_real)
      IF (ALLOCATED(procarray2d_real)) DEALLOCATE(procarray2d_real)
    ELSE
      IF (ALLOCATED(procarray_grib)) DEALLOCATE(procarray_grib)
      IF (ALLOCATED(procarray2d_grib)) DEALLOCATE(procarray2d_grib)
    ENDIF

ENDIF

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE output_data

!==============================================================================
!==============================================================================
!+ Module procedure in src_output for the p-interpolation
!------------------------------------------------------------------------------

SUBROUTINE p_int (outblock, i1,i2,i3, nlist, idebug, results)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates variables given in the namelist from
!   model leves to pressure levels. The result of the interpolation
!   is given back to the calling procedure in the three dimensional variable
!   "results". 
!
! Method:
!   Column wise interpolation with tension splines.
!   For the interpolation of FI and QV the logarithm of the pressure is
!   used.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
TYPE(pp_nl),              INTENT(IN) ::     &
  outblock       ! pointer to the namelist group

INTEGER (KIND=iintegers), INTENT(IN) ::     &
  i1,i2,i3,      & ! location of the variable to be processed in the LM 
                 ! variable table
  nlist,       & ! location of the variable in the output list
  idebug         ! for debug output

! Array arguments with intent(out):
REAL (KIND=ireals)      , INTENT(OUT)::     &
  results(ie,je,outblock%kepin)

REAL (KIND=ireals) ::       &
  zhelp1(ie,je,ke),         & !
  zhelp2(ie,je,ke)          !

!------------------------------------------------------------------------------
!
! Local parameters
REAL(KIND=ireals),PARAMETER    :: gamma   = 5.5_ireals    ! tension factor
REAL(KIND=ireals),PARAMETER    :: delpchk = 1.0_ireals

! Local scalars:
INTEGER (KIND=iintegers)       :: i, j, k, kint(ie,je)
INTEGER (KIND=iintegers)       :: ierrstat, ierr, izlen
INTEGER (KIND=iintegers)       :: nldim(ie,je)

REAL(KIND=ireals)              :: zt0s, zalnp, zpexp

CHARACTER (LEN=25)             :: yroutine
CHARACTER (LEN=80)             :: yerrmsg
CHARACTER (LEN=clen)           :: yzname

! Local arrays:
REAL(KIND=ireals)              :: fmfl(ie,je,ke)
REAL(KIND=ireals)              :: fexp(ie,ke+4),       &
                                  pexp(ie,ke+4)
REAL(KIND=ireals)              :: zdelp(ie,je)
REAL(KIND=ireals)              :: ztstar(ie,je), zalpha(ie,je), zt0(ie,je)

! Output from tautsp
REAL(KIND=ireals)              :: s_vec    (ie,(ke+4)*6), &
                                  break_vec(ie,(ke+4)*3), &
                                  coef_vec (ie,4,(ke+4)*3)

! Output from spline
!HJP 2017-12-20 Begin
!REAL(KIND=ireals)              :: fpfls(outblock%kepin),  &
!                                 pfls (outblock%kepin)
REAL(KIND=ireals)              :: pfls (ie,outblock%kepin)
!HJP 2017-12-20 End

! Switch to choose between vertical interpolation methods:
! Local value; there is a global namelist parameter itype_vertint
! in each namelist of group GRIBOUT.
INTEGER(KIND=iintegers) :: zitype_vertint

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Set variable that has to be interpolated
!------------------------------------------------------------------------------

yroutine = 'p_int'
yzname   = '          '
yzname(1:LEN_TRIM(outblock%yvarpl(nlist))) =                            &
         outblock%yvarpl(nlist)(1:LEN_TRIM(outblock%yvarpl(nlist)))
izlen = LEN_TRIM(yzname)

! Spline interpolation is the default value, which may be changed for single fields below
! (e.g., if you would like to have splineinterpolation for all fiels but not for QC ...)
zitype_vertint = outblock%itype_vertint

! calculation of u and v on masspoint
! -----------------------------------
#ifdef MESSY
fmfl = 0.0_ireals
#endif

  IF     (yzname(1:izlen) == 'U' ) THEN
    fmfl(2:ie,:,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(2:ie,:,:,itl) &
                                    + var(i1,i2,i3)%p4(1:ie-1,:,:,itl))  
    fmfl(1,:,1:ke)    = fmfl(2,:,1:ke)
  ELSEIF (yzname(1:izlen) == 'V' ) THEN
!CDIR COLLAPSE
    fmfl(:,2:je,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(:,2:je,:,itl) &
                                 + var(i1,i2,i3)%p4(:,1:je-1,:,itl))
    fmfl(:,1,1:ke)    = fmfl(:,2,1:ke)
  ELSEIF (yzname(1:izlen) == 'TKVM' .OR. yzname(1:izlen) == 'TKVH') THEN
!CDIR COLLAPSE
    fmfl(:,:,2:ke)    = var(i1,i2,i3)%p3(:,:,2:ke)
    fmfl(:,:,1)       = fmfl(:,:,2)
  ELSEIF (yzname(1:izlen) == 'FI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    =  0.5_ireals * ( hhl(:,:,1:ke) + hhl(:,:,2:ke+1) ) * g
  ELSEIF (yzname(1:izlen) == 'QRS' ) THEN
    IF (outblock%loutput_q_densities) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    =  zqrs_itl(:,:,:) * zrho_itl(:,:,:)
    ELSE
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    =  zqrs_itl(:,:,:)
    END IF
  ELSEIF (yzname(1:izlen) == 'Q_SEDIM') THEN
    IF (ASSOCIATED(qi)) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke) = zqrs_itl(:,:,:) - qi(:,:,:)
    ELSE
!CDIR COLLAPSE
      fmfl(:,:,1:ke) = zqrs_itl(:,:,:)
    ENDIF
    IF (outblock%loutput_q_densities) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke) = fmfl(:,:,1:ke) * zrho_itl(:,:,:)
    END IF
  ELSEIF (yzname(1:izlen) == 'RELHUM') THEN
    CALL calrelhum(fmfl(:,:,:), t(:,:,:,itl), pp(:,:,:,itl), p0(:,:,:), &
                   qv(:,:,:),ie, je, ke, b1, b2w, b3, b4w, rdv, o_m_rdv)
  ELSEIF (yzname(1:izlen) == 'OMEGA') THEN
    CALL calomega (fmfl(:,:,:), pp(:,:,:,nnew), pp(:,:,:,itl), pptens(:,:,:),&
                   w(:,:,:,itl), rho0(:,:,:), ie, je, ke, dt, g )
  ELSEIF (yzname(1:izlen) == 'VORTIC_W' ) THEN
    CALL curl (ie, je, ke, eddlon, eddlat, r_earth, acrlat, tgrlat,    &
              sqrtg_r_s, dzeta_dlam, dzeta_dphi, lmetr,               &
              u(:,:,:,itl ), v(:,:,:,itl ), w(:,:,:,itl ),            &
             .TRUE., zhelp1(:,:,:), zhelp2(:,:,:), fmfl(:,:,:))
  ELSEIF (yzname(1:izlen) == 'FI_ANAI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p_anai(:,:,1:ke) / zrho_itl(:,:,1:ke)
  ELSEIF (yzname(1:izlen) == 'TKE' ) THEN
    SELECT CASE( itype_turb )
    CASE( 5:8 )
      fmfl(:,:,1:ke) = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    CASE default
      fmfl(:,:,1:ke) = &
          0.5_ireals * (var(i1,i2,i3)%p4(:,:,1:ke,itl))**2
    END SELECT
  ELSEIF (yzname(1:izlen) == 'QC' .OR. yzname(1:izlen) == 'QNCLOUD'   .OR. &
          yzname(1:izlen) == 'QR' .OR. yzname(1:izlen) == 'QNRAIN'    .OR. &
          yzname(1:izlen) == 'QI' .OR. yzname(1:izlen) == 'QNICE'     .OR. &
          yzname(1:izlen) == 'QS' .OR. yzname(1:izlen) == 'QNSNOW'    .OR. &
          yzname(1:izlen) == 'QG' .OR. yzname(1:izlen) == 'QNGRAUPEL' .OR. &
          yzname(1:izlen) == 'QH' .OR. yzname(1:izlen) == 'QNHAIL')         THEN
    IF (outblock%loutput_q_densities) THEN
      ! output as densities:
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl) * zrho_itl(:,:,1:ke)
    ELSE
      ! output as mass specific: 
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    ENDIF
  ELSEIF (yzname(1:izlen) == 'DBZ') THEN
    IF (itype_gscp == 3) THEN
        CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
          klv850, my_cart_id, itype_gscp, idebug, t(:,:,:,itl),               &
          qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,         &
          qs(:,:,:)*zrho_itl, z_radar = fmfl )
    ELSEIF (itype_gscp == 4) THEN
        CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
          klv850, my_cart_id, itype_gscp, idebug, t(:,:,:,itl),               &
          qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,         &
          qs(:,:,:)*zrho_itl, q_grau  = qg(:,:,:)*zrho_itl, z_radar = fmfl )
#ifdef TWOMOM_SB
    ELSEIF (itype_gscp >= 100) THEN
      CALL radar_sb_ray (ie, je, ke, pi,                                 &
           klv850, my_cart_id, t(:,:,:,itl),                             &
           qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl,                       &
           qi(:,:,:)*zrho_itl, qs(:,:,:)*zrho_itl,                       &
           qg(:,:,:)*zrho_itl, qh(:,:,:)*zrho_itl,                       &
           qnc(:,:,:)*zrho_itl, qnr(:,:,:)*zrho_itl,                     &
           qni(:,:,:)*zrho_itl, qns(:,:,:)*zrho_itl,                     &
           qng(:,:,:)*zrho_itl, qnh(:,:,:)*zrho_itl,                     &
           z_radar = fmfl )
#endif
    ENDIF
!CDIR COLLAPSE
    ! Reflectivity interpolation should be done in linear space, not in logarithmic (as alternative, interpolation
    ! in the space of rain rate -- Z^0.6666 --- would also be desireable):
    WHERE (fmfl(:,:,:) >= -90.0)
      fmfl(:,:,:) = 10.0_ireals ** (0.1_ireals * fmfl(:,:,:))
    ELSEWHERE
      fmfl(:,:,:) = 0.0_ireals
    END WHERE
  ELSEIF ( yzname(1:izlen) == 'RHO_TOT' ) THEN
     fmfl(:,:,1:ke) = zrho_itl(:,:,1:ke)
  ELSE
    SELECT CASE( var(i1,i2,i3)%rank )
    CASE(4)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    CASE(3)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p3(:,:,1:ke)
    END SELECT
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Interpolation
!------------------------------------------------------------------------------

! slicewise interpolation
! -----------------------

  ! Hardwired choice of linear interpolation for hydrometeors and DBZ to
  ! prevent spurious overshoots:
  IF ( yzname(1:izlen) == 'QC' .OR. yzname(1:izlen) == 'QR' .OR. &
       yzname(1:izlen) == 'QI' .OR. yzname(1:izlen) == 'QS' .OR. &
       yzname(1:izlen) == 'QG' .OR. yzname(1:izlen) == 'QH' .OR. &
       yzname(1:3)     == 'DBZ'                                  ) THEN
    zitype_vertint = 2
  END IF
#ifdef TWOMOM_SB
  IF ( yzname(1:izlen) == 'QNCLOUD' .OR. yzname(1:izlen) == 'QNRAIN' .OR. &
       yzname(1:izlen) == 'QNICE'   .OR. yzname(1:izlen) == 'QNSNOW' .OR. &
       yzname(1:izlen) == 'QNGRAUPEL'  .OR. yzname(1:izlen) == 'QNHAIL' ) THEN
    zitype_vertint = 2
  END IF
#endif

  DO j=jstartpar,jendpar

    ! Calculation of the pressure at full model levels.
    ! Variables defined on half levels are first averaged to
    ! full levels.
    ! For water vapour (QV) and geopotential height (FI), the 
    ! interpolation on constant pressure levels is logarithmic 
    ! with respect to pressure.

    SELECT CASE(var(i1,i2,i3)%levtyp)
    CASE(110)
      IF (yzname(1:izlen) == 'QV'       .OR.                     &
          yzname(1:izlen) == 'FI'       .OR.                     &
          yzname(1:izlen) == 'QV_ANAI'  .OR.                     &
          yzname(1:izlen) == 'FI_ANAI' ) THEN
        DO k = 1, ke
          DO i = 1, ie
            pexp(i,k) = LOG(p0(i,j,k) + pp(i,j,k,itl))
          ENDDO
        ENDDO
!HJP 2017-12-20 Begin
!       pfls(1:outblock%kepin) = LOG(outblock%plev(1:outblock%kepin))
        DO k = 1, outblock%kepin
          DO i = 1, ie
            pfls(i,k) = LOG(outblock%plev(k))
          ENDDO
        ENDDO
!HJP 2017-12-20 End   
      ELSE
        DO k = 1, ke
          DO i = 1, ie
            pexp(i,k)  = p0(i,j,k) + pp(i,j,k,itl)
          ENDDO
        ENDDO
!HJP 2017-12-20 Begin
!       pfls(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
        DO k = 1, outblock%kepin
          DO i = 1, ie
            pfls(i,k) = outblock%plev(k)
          ENDDO
        ENDDO
!HJP 2017-12-20 End   
      ENDIF
    CASE(109)
      DO k = 2, ke
        DO i = 1, ie
          pexp(i,k) = p0hl(i,j,k)  &
                  + 0.5_ireals*(pp(i,j,k,itl)+pp(i,j,k-1,itl))
        ENDDO
      ENDDO
      DO i = 1, ie
        pexp(i,1) = p0hl(i,j,1) + pp(i,j,1,itl)
             ! best approximation we can get there
             ! before, this was used, but factor 0.5 is wrong
             !                + 0.5_ireals *pp(i,j,1,itl)
      ENDDO
!HJP 2017-12-20 Begin
!     pfls(1:outblock%kepin) = outblock%plev(1:outblock%kepin)
      DO k = 1, outblock%kepin
        DO i = 1, ie
          pfls(i,k) = outblock%plev(k)
        ENDDO
      ENDDO
!HJP 2017-12-20 End   
    CASE DEFAULT
               yerrmsg = 'wrong leveltyp of input field'
               ierrstat = 2004
               CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END SELECT

    ! Check if a plevel for interpolation is above the model top:
    IF ( MINVAL(outblock%plev(1:outblock%kepin)) <= MINVAL(pexp(1:ie,1)) ) THEN
      yerrmsg = 'plev for interpolation above model top!'
      ierrstat = 1004
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END IF

    ! calculate pressure difference of lowest constant pressure level
    ! and surface pressure

    zdelp(:,j) = outblock%plev(outblock%kepin) - ps(:,j,itl)

    ! vertical interpolation colum wise
    ! ---------------------------------

    ! copy slice variable to colum variable and set lowest modellevel
    DO k = 1, ke
      DO i=istartpar,iendpar
        fexp(i,k)  = fmfl(i,j,k)
      ENDDO
    ENDDO
    DO i=istartpar,iendpar
      fexp(i,ke+1)  = fexp(i,ke)
      pexp(i,ke+1)  = ps(i,j,itl)

      IF (zdelp(i,j) > 0.0) THEN
        kint(i,j) = ke + 4
      ELSE
        kint(i,j) = ke + 1
      ENDIF
      nldim(i,j) = (ke+4)*3
    ENDDO

    IF (TRIM(yzname) == 'T') THEN
      ! calculation of the soil temperature (ztstar) from the temperature
      ! of the lowest model level
      DO i=istartpar,iendpar
        ztstar(i,j) = t(i,j,ke,itl) + 0.0065_ireals                       &
                      * 0.5_ireals*(hhl(i,j,ke)-hhl(i,j,ke+1))
        zalpha(i,j) = 0.0065*r_d/g
        zt0   (i,j) = ztstar(i,j)  + 0.0065 * hsurf(i,j)
      ENDDO

      DO i=istartpar,iendpar
        IF (ztstar(i,j) > 315.0 .AND. zt0(i,j) > 315.0) THEN
          ztstar(i,j) = 0.5 * ( 315.0 + ztstar(i,j))
        ELSEIF (ztstar(i,j) < 255.0) THEN
          ztstar(i,j) = 0.5 * ( 255.0 + ztstar(i,j))
        ENDIF
        IF(hsurf(i,j) > 2500_ireals) THEN
          zt0s   = MIN(zt0(i,j),298.0_ireals)
          zalpha(i,j) = r_d * (zt0s - ztstar(i,j))/(hsurf(i,j)*g)
          IF(zt0s < ztstar(i,j)) zalpha(i,j) = 0_ireals
        ELSEIF((hsurf(i,j) <= 2500_ireals)                        &
          .AND.(hsurf(i,j) >= 2000_ireals)) THEN
          zt0s   = 0.002_ireals * (2500_ireals-hsurf(i,j) * zt0(i,j)  &
                   +(hsurf(i,j)-2000_ireals)*MIN(zt0(i,j),298.0_ireals))
          zalpha(i,j) = r_d * (zt0s - ztstar(i,j))/(hsurf(i,j)*g)
          IF(zt0s < ztstar(i,j)) zalpha(i,j) = 0.0_ireals
        ENDIF
      ENDDO
    ENDIF

    IF (TRIM(yzname) == 'FI') THEN
      ! calculation of the soil temperature (ztstar) from the temperature
      ! of the lowest model level
      DO i=istartpar,iendpar
        ztstar(i,j) = t(i,j,ke,itl) + 0.0065_ireals                       &
                      * 0.5_ireals*(hhl(i,j,ke)-hhl(i,j,ke+1))
        zalpha(i,j) = 0.0065*r_d/g
        zt0   (i,j) = ztstar(i,j)  + 0.0065 * hsurf(i,j)
      ENDDO

      DO i=istartpar,iendpar
        IF (ztstar(i,j) <= 290.5 .AND. zt0(i,j) > 290.5) THEN
          zalpha(i,j) = r_d*(290.5 - ztstar(i,j))/(hsurf(i,j)*g)
        ELSEIF (ztstar(i,j) > 290.5 .AND. zt0(i,j) > 290.5) THEN
          ztstar(i,j) = 0.5*(290.5 + ztstar(i,j))
          zalpha(i,j) = 0.0
        ELSEIF (ztstar(i,j) < 255.0) THEN
          ztstar(i,j) = 0.5*(255.0 + ztstar(i,j))
        ENDIF
      ENDDO
    ENDIF


    IF (yzname(1:izlen) == 'T' ) THEN
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          ! extrapolation of the model data if the required pressure level
          ! is below the lowest modellevel
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itl) &
                    + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            zalnp   = zalpha(i,j) * LOG(pexp(i,k)/ps(i,j,itl))
            fexp(i,k) = ztstar(i,j) *                                     &
                      ( 1.0_ireals + zalnp + 0.5*zalnp**2 + 0.1667*zalnp**3)
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (yzname(1:izlen) == 'FI' ) THEN
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          DO k=ke+1,ke+4
            zpexp   = ps(i,j,itl) &
                     + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            zalnp   = zalpha(i,j) * LOG(zpexp/ps(i,j,itl))
            fexp(i,k) = hsurf(i,j) * g                                   &
                          -r_d*ztstar(i,j) * LOG(zpexp/ps(i,j,itl))*     &
                                (1.0_ireals+zalnp                        &
                                +0.5_ireals*zalnp**2                     &
                           +0.1667_ireals*zalnp**3)
            pexp(i,k) = LOG(zpexp)
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO i=istartpar,iendpar
        IF (zdelp(i,j) > 0.0) THEN
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itl) &
                    + (zdelp(i,j) + delpchk)/4.0_ireals* REAL (k-ke, ireals)
            fexp(i,k) = fexp(i,ke)
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    SELECT CASE (zitype_vertint)

    CASE (1)

      ! Spline interpolation over (i,k)-slices
      ! --------------------------------------

      CALL tautsp2D (pexp(:,:), fexp(:,:), kint(:,j) ,ie, istartpar,     &
           iendpar, ke+4, gamma, s_vec, break_vec, coef_vec, nldim, ierr)
      IF (ierr == 2) THEN
        yerrmsg = 'wrong input in tautsp'
        ierrstat = 1004
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in tautsp *** '
        ierrstat = 1005
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      
      CALL spline2D(break_vec,coef_vec,j,nldim,pfls,outblock%kepin,    &
           results, outblock%yvarpl(nlist))
      
    CASE (2)

      ! Linear Interpolation with respect to pressure:
      ! ----------------------------------------------

      ! .. provide monotonically increasing dummy values for pexp(:,ke+1:ke+4)
      !    below the surface, because lininterp2D_xinter1D_vec()
      !    requires that all values of pexp are prescribed:
      DO i=istartpar,iendpar
        IF (zdelp(i,j) <= 0.0) THEN
          DO k=ke+1,ke+4
            pexp(i,k) = ps(i,j,itl) &
                    + delpchk/4.0_ireals* REAL (k-ke, ireals)
            fexp(i,k) = fexp(i,ke)
          ENDDO
        ENDIF
      ENDDO

      ierr = 0
      CALL lininterp2D_xinter1D_vec(pexp(istartpar:iendpar,:), &
           fexp(istartpar:iendpar,:), iendpar-istartpar+1, ke+4, &
!HJP 2017-12-20 Begin
!          pfls(:), results(istartpar:iendpar,j,:), outblock%kepin, ierr)
           pfls(istartpar:iendpar,:), results(istartpar:iendpar,j,:), outblock%kepin, ierr)
!HJP 2017-12-20 End
      
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in lininterp2D_xinter1D_vec *** '
        ierrstat = 1006
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

    CASE default
      
      yerrmsg = ' ERROR *** Wrong value for zitype_vertint *** '
      ierrstat = 1007
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      
    END SELECT

    !.. UB: Correct interpolation artifacts:
    !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$    IF (TRIM(yzname) == 'RELHUM') THEN
!!$      results(istartpar:iendpar,j,:) = MAX (0.0_ireals, MIN(100.0_ireals,results(istartpar:iendpar,j,:)))
!!$    ELSEIF (TRIM(yzname) == 'QC') THEN
!!$      results(istartpar:iendpar,j,:) = MAX (0.0_ireals, results(istartpar:iendpar,j,:) )
!!$    END IF

  ENDDO ! loop over j

  ! In case of radar reflectivity, transform back to log space after interpolation:
  IF (yzname(1:izlen) == 'DBZ') THEN
    WHERE (results(istartpar:iendpar,jstartpar:jendpar,:) >= 1e-9_ireals)
      results(istartpar:iendpar,jstartpar:jendpar,:) = &
           10.0_ireals * LOG10(results(istartpar:iendpar,jstartpar:jendpar,:))
    ELSEWHERE
      results(istartpar:iendpar,jstartpar:jendpar,:) = -99.99
    END WHERE
  END IF

  !.. UB: Set data below the surface to missing values: 
  !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$  SELECT CASE (yzname(1:izlen))
!!$  CASE ('DBZ')
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kepin
!!$          IF (outblock%plev(k) > ps(i,j,itl)) results(i,j,k) = -99.99
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  CASE default
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kepin
!!$          IF (outblock%plev(k) > ps(i,j,itl)) results(i,j,k) = 0.0_ireals
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE p_int

!==============================================================================
!+ Module procedure in src_output for the z-interpolation
!------------------------------------------------------------------------------

SUBROUTINE z_int (outblock, i1,i2,i3, nlist, idebug, results)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine interpolates a number of 3-d variables, given by namelist-
!   input to this routine, from model layers to fixed-height levels (z-levles)
!
! Method:
!   Colum wise interpolation with tension splines.
!   For the interpolation of FI and QV the logarithm of the pressure is
!   used.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
TYPE(pp_nl),              INTENT(IN) ::     &
  outblock       ! pointer to the namelist group

INTEGER (KIND=iintegers), INTENT(IN) ::     &
  i1,i2,i3,      & ! location of the variable to be processed in the LM
                 ! variable table
  nlist,       & ! location of the variable in the output list
  idebug         ! for debug output

! Array arguments with intent(out):
REAL (KIND=ireals)      , INTENT(OUT)::     &
  results(ie,je,outblock%kezin)

!------------------------------------------------------------------------------

! Local parameters
REAL(KIND=ireals),PARAMETER    :: gamma   = 5.5_ireals    ! tension factor

INTEGER (KIND=iintegers)       :: i, j, k, kint(ie)
INTEGER (KIND=iintegers)       :: ierrstat, ierr, izlen
INTEGER (KIND=iintegers)       :: nldim(ie,1)

CHARACTER (LEN=25)             :: yroutine
CHARACTER (LEN=80)             :: yerrmsg
CHARACTER (LEN=clen)           :: yzname

! Local arrays:
REAL(KIND=ireals)              :: fmfl(ie,je,ke+1),  &
                                  zmfl(ie,   ke+1)
REAL(KIND=ireals)              :: fexp(ie,ke+4),  &
                                  zexp(ie,ke+4)

! Output from tautsp
REAL(KIND=ireals)              :: s_vec    (ie,(ke+4)*6), &
                                  break_vec(ie,(ke+4)*3), &
                                  coef_vec (ie,4,(ke+4)*3)

! Output from spline
!HJP 2017-12-20 Begin
!REAL(KIND=ireals)              :: fzfls(outblock%kezin),  &
!                                  zfls (outblock%kezin)
REAL(KIND=ireals)              :: zfls (ie,outblock%kezin)
!HJP 2017-12-20 Begin

! Switch to choose between vertical interpolation methods:
! Local value; there is a global namelist parameter itype_vertint
! in each namelist of group GRIBOUT.
INTEGER(KIND=iintegers) :: zitype_vertint
!
! Switch to choose if values of, U, V, W, and T are extrapolated
! to the surface based on no-slip- and skin-condition
! or if a constant extrapolation is done.
! The former is not so good for cubic spline interpolation
! because of oszillations with unwanted results, so we 
! allow it only for linear interpolation.
LOGICAL :: zlsurfextrapol_noslip

!------------------------------------------------------------------------------
! Section 1: general and grib-io preparations
!------------------------------------------------------------------------------

yroutine = 'z_int'
yzname   = '          '
yzname(1:LEN_TRIM(outblock%yvarzl(nlist))) =                            &
         outblock%yvarzl(nlist)(1:LEN_TRIM(outblock%yvarzl(nlist)))
izlen = LEN_TRIM(yzname)

! The default value is taken from the namelist value, 
! which may be "hard" changed for single fields below
! (e.g., if you would not like to have splineinterpolation for fiels like QC ...)
zitype_vertint = outblock%itype_vertint

! Allow  no-slip- and skin-condition-extrapolation to the
! surface only for linear interpolation:
IF (zitype_vertint == 2) THEN
  zlsurfextrapol_noslip = .TRUE.
ELSE
  zlsurfextrapol_noslip = .FALSE.
END IF

!------------------------------------------------------------------------------
! Section 1: Set variable that has to be interpolated
!------------------------------------------------------------------------------

  ! calculation of u and v on masspoint
  ! -----------------------------------

  IF     (yzname(1:izlen) == 'U' ) THEN
    fmfl(2:ie,:,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(2:ie,:,1:ke,itl)    &
                                    + var(i1,i2,i3)%p4(1:ie-1,:,1:ke,itl))
!CDIR COLLAPSE
    fmfl(1,:,1:ke)    = fmfl(2,:,1:ke)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
    ELSE
      ! For interpolation: Surface velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'V' ) THEN
!CDIR COLLAPSE
    fmfl(:,2:je,1:ke) = 0.5_ireals * (var(i1,i2,i3)%p4(:,2:je,1:ke,itl)    &
                                    + var(i1,i2,i3)%p4(:,1:je-1,1:ke,itl))
    fmfl(:,1,1:ke)    = fmfl(:,2,1:ke)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
    ELSE
      ! For interpolation: Surface velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'VABS_AV' ) THEN
    fmfl(:,:,1:ke)    = var(i1,i2,i3)%p3(:,:,1:ke)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
    ELSE
      ! For interpolation: Surface velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'P' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p0(:,:,1:ke) + pp(:,:,1:ke,itl)
! For interpolation: Surface pressure at the ground!
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = ps(:,:,itl)
  ELSEIF (yzname(1:izlen) == 'W' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    IF (lnosurffluxes_m .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ELSE
      ! For interpolation: Surface vertical velocity = 0
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = 0.0
    END IF
  ELSEIF (yzname(1:izlen) == 'T' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    IF (lnosurffluxes_h .OR. .NOT.zlsurfextrapol_noslip) THEN
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ELSE
      ! For interpolation: Surface temperature = interfacial temperature
!CDIR COLLAPSE
      fmfl(:,:,ke+1)    = t_g(:,:,itl)
    END IF
  ELSEIF (yzname(1:izlen) == 'TKVM' .OR. yzname(1:izlen) == 'TKVH') THEN
!CDIR COLLAPSE
    fmfl(:,:,2:ke)    = var(i1,i2,i3)%p3(:,:,2:ke)
!CDIR COLLAPSE
    fmfl(:,:,1)       = fmfl(:,:,2)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'TKE' ) THEN
    SELECT CASE( itype_turb )
    CASE( 5:8 )
      fmfl(:,:,1:ke1) = var(i1,i2,i3)%p4(:,:,1:ke1,itl)
    CASE default
      fmfl(:,:,1:ke1) = &
          0.5_ireals * (var(i1,i2,i3)%p4(:,:,1:ke1,itl))**2
    END SELECT
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'QRS') THEN
    IF (outblock%loutput_q_densities) THEN
      ! output as densities:
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    = zqrs_itl(:,:,:) * zrho_itl(:,:,:)
    ELSE
      ! output as mass specific:
!CDIR COLLAPSE
      fmfl(:,:,1:ke)    = zqrs_itl(:,:,:)
    END IF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'Q_SEDIM') THEN
    IF (ASSOCIATED(qi)) THEN
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = zqrs_itl(:,:,:) - qi(:,:,:)
    ELSE
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = zqrs_itl(:,:,:)
    ENDIF
    IF (outblock%loutput_q_densities) THEN
!CDIR COLLAPSE
     fmfl(:,:,1:ke)  = fmfl(:,:,1:ke) * zrho_itl(:,:,:)
    END IF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'RELHUM') THEN
    CALL calrelhum(fmfl, t(:,:,:,itl), pp(:,:,:,itl), p0(:,:,:), &
                   qv(:,:,:),ie, je, ke, b1, b2w, b3, b4w, rdv, o_m_rdv)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'OMEGA') THEN
    CALL calomega (fmfl, pp(:,:,:,nnew), pp(:,:,:,itl), pptens(:,:,:), &
                   w(:,:,:,itl), rho0(:,:,:), ie, je, ke, dt, g )
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'FI_ANAI' ) THEN
!CDIR COLLAPSE
    fmfl(:,:,1:ke)    = p_anai(:,:,1:ke) / zrho_itl(:,:,1:ke)
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'QC' .OR. yzname(1:izlen) == 'QNCLOUD'   .OR. &
          yzname(1:izlen) == 'QR' .OR. yzname(1:izlen) == 'QNRAIN'    .OR. &
          yzname(1:izlen) == 'QI' .OR. yzname(1:izlen) == 'QNICE'     .OR. &
          yzname(1:izlen) == 'QS' .OR. yzname(1:izlen) == 'QNSNOW'    .OR. &
          yzname(1:izlen) == 'QG' .OR. yzname(1:izlen) == 'QNGRAUPEL' .OR. &
          yzname(1:izlen) == 'QH' .OR. yzname(1:izlen) == 'QNHAIL')         THEN
    ! Output of cloud microphysics variables either as densities or mass specific:
    IF (outblock%loutput_q_densities) THEN
      ! output as densities:
      fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl) * zrho_itl(:,:,1:ke)
    ELSE
      ! output as mass specific:
      fmfl(:,:,1:ke)    = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    END IF
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
  ELSEIF (yzname(1:izlen) == 'DBZ') THEN
    IF (itype_gscp == 3) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
           klv850, my_cart_id, itype_gscp, idebug, t(:,:,:,itl),            &
           qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,      &
           qs(:,:,:)*zrho_itl, z_radar = fmfl )
    ELSEIF (itype_gscp == 4) THEN
      CALL radar_lm_ray (ie,je,ke, pi, rho_w, rho_ice, K_w, K_ice, t0_melt, &
           klv850, my_cart_id, itype_gscp, idebug, t(:,:,:,itl),            &
           qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl, qi(:,:,:)*zrho_itl,      &
           qs(:,:,:)*zrho_itl, q_grau  = qg(:,:,:)*zrho_itl, z_radar = fmfl )
#ifdef TWOMOM_SB
    ELSEIF (itype_gscp >= 100) THEN
      CALL radar_sb_ray (ie, je, ke, pi,                                 &
           klv850, my_cart_id, t(:,:,:,itl),                             &
           qc(:,:,:)*zrho_itl, qr(:,:,:)*zrho_itl,                       &
           qi(:,:,:)*zrho_itl, qs(:,:,:)*zrho_itl,                       &
           qg(:,:,:)*zrho_itl, qh(:,:,:)*zrho_itl,                       &
           qnc(:,:,:)*zrho_itl, qnr(:,:,:)*zrho_itl,                     &
           qni(:,:,:)*zrho_itl, qns(:,:,:)*zrho_itl,                     &
           qng(:,:,:)*zrho_itl, qnh(:,:,:)*zrho_itl,                     &
           z_radar = fmfl(:,:,1:ke) )
#endif
    ENDIF
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke)
    ! Reflectivity interpolation should be done in linear space, not in logarithmic (as alternative, interpolation
    ! in the space of rain rate -- Z^0.6666 --- would also be desireable):
    WHERE (fmfl(:,:,:) >= -90.0)
      fmfl(:,:,:) = 10.0_ireals ** (0.1_ireals * fmfl(:,:,:))
    ELSEWHERE
      fmfl(:,:,:) = 0.0_ireals
    END WHERE
  ELSEIF ( yzname(1:izlen) == 'RHO_TOT' ) THEN
    fmfl(:,:,1:ke) = zrho_itl(:,:,1:ke)
    fmfl(:,:,ke+1) = fmfl(:,:,ke)
  ELSE
    SELECT CASE( var(i1,i2,i3)%rank )
    CASE(4)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p4(:,:,1:ke,itl)
    CASE(3)
!CDIR COLLAPSE
      fmfl(:,:,1:ke)  = var(i1,i2,i3)%p3(:,:,1:ke)
    END SELECT
!CDIR COLLAPSE
    fmfl(:,:,ke+1)    = fmfl(:,:,ke  )
  ENDIF

!------------------------------------------------------------------------------
! Section 2: Interpolation
!------------------------------------------------------------------------------

  ! Hardwired choice of linear interpolation for hydrometeors and DBZ to
  ! prevent spurious overshoots:
  IF ( yzname(1:izlen) == 'QC' .OR. yzname(1:izlen) == 'QR' .OR. &
       yzname(1:izlen) == 'QI' .OR. yzname(1:izlen) == 'QS' .OR. &
       yzname(1:izlen) == 'QG' .OR. yzname(1:izlen) == 'QH' .OR. &
       yzname(1:3)     == 'DBZ'                                  ) THEN
    zitype_vertint = 2
  END IF
#ifdef TWOMOM_SB
  IF ( yzname(1:izlen) == 'QNCLOUD' .OR. yzname(1:izlen) == 'QNRAIN' .OR. &
       yzname(1:izlen) == 'QNICE'   .OR. yzname(1:izlen) == 'QNSNOW' .OR. &
       yzname(1:izlen) == 'QNGRAUPEL'  .OR. yzname(1:izlen) == 'QNHAIL' ) THEN
    zitype_vertint = 2
  END IF
#endif

  ! slicewise interpolation
  ! -----------------------

  DO j=jstartpar,jendpar

    ! Set the height of the modell levels (zmfl) and the height of
    ! the constant z-levels to interpolate to (zfls).

    SELECT CASE(var(i1,i2,i3)%levtyp)
    CASE(110)
      zmfl(:,1:ke) = 0.5_ireals * ( hhl(:,j,1:ke) + hhl(:,j,2:ke+1) )
      zmfl(:,ke+1) = hsurf(:,j)
!HJP 2017-12-20 Begin
!     zfls(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
!HJP 2017-12-20 End
    CASE(109)
      zmfl(:,1:ke+1) = hhl(:,j,1:ke+1)
!HJP 2017-12-20 Begin
!     zfls(1:outblock%kezin) = outblock%zlev(1:outblock%kezin)
!HJP 2017-12-20 End
    CASE DEFAULT
      yerrmsg = 'wrong leveltyp of input field'
      ierrstat = 2004
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END SELECT

!HJP 2017-12-20 Begin
    ! levels to interpolate to are the same for all levtyps
    IF (outblock%lzint_above_ground) THEN
      ! output to height levels above ground
      DO k = 1, outblock%kezin
        DO i = 1, ie
          zfls(i,k) = hhl(i,j,ke+1) + outblock%zlev(k)
        ENDDO
      ENDDO
    ELSE
      ! output to altitude above mean sea level
      DO k = 1, outblock%kezin
        DO i = 1, ie
          zfls(i,k) = outblock%zlev(k)
        ENDDO
      ENDDO
    ENDIF
!HJP 2017-12-20 End

    ! Check if a zlevel for interpolation is above the model top height:
    IF ( MAXVAL(outblock%zlev(1:outblock%kezin)) >= MAXVAL(zmfl(1:ie,1)) ) THEN
      yerrmsg = 'zlev for interpolation above model top!'
      ierrstat = 1004
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
    END IF

    ! To handle cases where the lowest z-level is below the orography, an
    ! additional data point at z = - 100mf  will be created, assuming a
    ! constant profile.
    ! The interpolated variable value will be set to undefined.

    ! vertical interpolation colum wise
    ! ---------------------------------

    DO i=istartpar,iendpar

      ! copy slice variable to colum variable (in reverse order, because
      ! data points must strictly increase in tautsp) and add an extra
      ! data/function point in any case

      kint (i) = ke + 2
!US   zexp (i,1) = - 100.0_ireals
      zexp (i,1) = MIN (-100.0_ireals, zmfl(i,ke+1) - 20.0_ireals)
      fexp (i,1) = fmfl(i,j,ke+1)
      nldim(i,1) = (ke+4)*3
    ENDDO

    DO k = 1, ke1
      DO i=istartpar,iendpar
        zexp(i,k+1)  = zmfl(i,  ke+2-k)
        fexp(i,k+1)  = fmfl(i,j,ke+2-k)
      ENDDO

    ENDDO ! loop over i


  ! Spline interpolation over total domain or slicewise
  ! ---------------------------------------------------

    SELECT CASE (zitype_vertint)

    CASE (1)

      ! Interpolation by cubic TauT-Splines (Tension Splines):
      ! ------------------------------------------------------

      CALL tautsp2D (zexp(:,:), fexp(:,:), kint(:), ie, istartpar,     &
           iendpar, ke+4, gamma, s_vec, break_vec, coef_vec, nldim, ierr)
      IF (ierr == 2) THEN
        yerrmsg = 'wrong input in tautsp'
        ierrstat = 1004
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF
      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in tautsp *** '
        ierrstat = 1005
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

      CALL spline2D(break_vec,coef_vec,j,nldim,zfls,outblock%kezin,    &
           results,outblock%yvarzl(nlist))

    CASE (2)

      ! Linear Interpolation with respect to height:
      ! --------------------------------------------
      ierr = 0
         CALL lininterp2D_xinter1D_vec(zexp(istartpar:iendpar,1:ke+2), &
             fexp(istartpar:iendpar,1:ke+2), iendpar-istartpar+1, ke+2, &
!HJP 2017-12-20 Begin
!          zfls(:), results(istartpar:iendpar,j,:), outblock%kezin, ierr)
           zfls(istartpar:iendpar,:), results(istartpar:iendpar,j,:), &
           outblock%kezin, ierr)
!HJP 2017-12-20 End  

      IF (ierr /= 0) THEN
        yerrmsg = ' ERROR *** Error in lininterp2D_xinter1D_vec *** '
        ierrstat = 1006
        CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)
      ENDIF

    CASE default

      yerrmsg = ' ERROR *** Wrong value for zitype_vertint *** '
      ierrstat = 1007
      CALL model_abort (my_cart_id, ierrstat, yerrmsg, yroutine)

    END SELECT

     !.. UB: Correct interpolation artifacts:
     !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$     IF (TRIM(yzname) == 'RELHUM') THEN
!!$       results(:,j,:) = MAX (0.0_ireals, MIN(100.0_ireals,results(:,j,:)))
!!$     ELSEIF (TRIM(yzname) == 'QC') THEN
!!$       results(:,j,:) = MAX (0.0_ireals, results(:,j,:) )
!!$     END IF

  ENDDO ! loop over j

  ! In case of radar reflectivity, transform back to log space after interpolation:
  IF (yzname(1:izlen) == 'DBZ') THEN
    WHERE (results(istartpar:iendpar,jstartpar:jendpar,:) >= 1e-9_ireals)
      results(istartpar:iendpar,jstartpar:jendpar,:) = &
           10.0_ireals * LOG10(results(istartpar:iendpar,jstartpar:jendpar,:))
    ELSEWHERE
      results(istartpar:iendpar,jstartpar:jendpar,:) = -99.99
    END WHERE
  END IF

  !.. UB: Set data below the surface to missing values:
  !   (MAY BE COMMENTED IN BY SPECIALIZED USERS, BUT NOT OPERATIONNALY AT DWD!)
!!$  SELECT CASE (yzname(1:izlen))
!!$  CASE ('DBZ')
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kezin
!!$          IF (outblock%zlev(k) < hsurf(i,j)) results(i,j,k) = -99.99
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO
!!$  CASE default
!!$    DO i=istartpar,iendpar
!!$      DO j=jstartpar,jendpar
!!$        DO k=1,outblock%kezin
!!$          IF (outblock%zlev(k) < hsurf(i,j)) results(i,j,k) = 0.0_ireals
!!$        ENDDO
!!$      ENDDO
!!$    ENDDO    
!!$  END SELECT


!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE z_int

!==============================================================================
!+ calculates the function values for given coefficients
!------------------------------------------------------------------------------

SUBROUTINE spline (break, coef, nldim, sptau, spgtau, ng)

!------------------------------------------------------------------------------
!
! Description:
!  spline calculates the values of the interpolation function for given
!  interpolation coefficients and arguments
!
! Method:
!  Spline interpolation.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  nldim, ng    ! Dimensions of the variables
!
! Array arguments with intent(in):
REAL (KIND=ireals),       INTENT(IN)    ::  &
  break(nldim),      & ! arguments for which the function has to be calculated
  coef(4,nldim),     & ! coefficients of the interpolation function
  sptau(ng)            !

!
! Array arguments with intent(out):
REAL (KIND=ireals),       INTENT(OUT)   ::  &
  spgtau(ng)

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)    ::  i, j
REAL (KIND=ireals)          ::  dx, prod

!
!- End of header
!==============================================================================


! Calculate the product (sptau(i)-break(j)) * (sptau(i)-break(j+1))
! and look for the interval where sptau is
  DO i  = 1, ng
    DO j  = 1, nldim-1
      prod   = (sptau(i) - break(j)) * (sptau(i) - break(j+1))

      ! calculate the splines
      IF( prod <= 0.0_ireals ) THEN
        dx = sptau(i) - break(j)
        spgtau(i) = coef(1,j) + dx * (coef(2,j) + dx*0.5*(coef(3,j)     &
                                    + dx/3.0* coef(4,j)))
        EXIT       
      ENDIF
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE SPLINE

!==============================================================================

!option! -pvctl _on_adb
SUBROUTINE spline2D (break_vec, coef_vec, jline, nldim, sptau, ng,    &
                     results, yname)

!------------------------------------------------------------------------------
!
! Description:
!  spline calculates the values of the interpolation function for given
!  interpolation coefficients and arguments
!
! Method:
!  Spline interpolation.
!
!------------------------------------------------------------------------------

! Subroutine / Function arguments
! Scalar arguments with intent(in):
INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  ng, jline            ! Dimensions of variables and line to be processed

INTEGER (KIND=iintegers), INTENT(IN)    ::  &
  nldim(ie,1)

!
! Array arguments with intent(in):
REAL (KIND=ireals),       INTENT(IN)    ::  &
  break_vec(ie,*), & ! arguments for which the function is calculated
  coef_vec (ie,4,*),&! coefficients of the interpolation function
!HJP 2017-12-20 Begin
! sptau(ng)            !
  sptau(ie,ng)       ! every grid point can have different heights
!HJP 2017-12-20 End  


CHARACTER (LEN=clen),     INTENT(IN)    :: yname

! Array arguments with intent(out):
REAL(KIND=ireals),        INTENT(OUT)   :: results(ie,je,ng)

!------------------------------------------------------------------------------

! Local scalars:
INTEGER (KIND=iintegers)    :: i, j, k, maxng, nzl(ie), nzcount, nzelements
REAL (KIND=ireals)          :: dx, prod
LOGICAL                     :: ldone(ie), levelready

!
!- End of header
!==============================================================================

  j = jline
  nzl(:) = 1
  nzelements = iendpar - istartpar + 1

  maxng = MAXVAL (nldim(istartpar:iendpar,1))

  DO k = 1, maxng  !ke+4*3
    nzcount = 0
    ldone(:) = .FALSE.

    DO WHILE (nzcount < nzelements)
      DO i = istartpar, iendpar

        IF (.NOT. ldone(i)) THEN
          ! sptau has to be between break_vec(i,k) and break_vec(i,k+1)

          levelready=(nzl(i) > ng)
          IF (.NOT. levelready) THEN
!HJP 2017-12-20 Begin
!           IF (sptau(nzl(i)) > break_vec(i,k+1)) THEN
            IF (sptau(i,nzl(i)) > break_vec(i,k+1)) THEN
!HJP 2017-12-20 End

              levelready=.TRUE.
            END IF
          END IF

          IF (levelready) THEN
            ! the level nzl(i) is ready
            ldone(i) = .TRUE.
            nzcount    = nzcount + 1
          ELSE
!HJP 2017-12-20 Begin
!           dx = sptau(nzl(i)) - break_vec(i,k)
            dx = sptau(i,nzl(i)) - break_vec(i,k)
!HJP 2017-12-20  End
            results(i,j,nzl(i)) = coef_vec(i,1,k) + dx     *(coef_vec(i,2,k) +   &
                  dx*0.5*(coef_vec(i,3,k) + dx/3.0 * coef_vec(i,4,k)))
            nzl(i) = nzl(i) + 1
          END IF

! if this if is splitted in 2 ifs (as below), it will not vectorize
!          IF (nzl(i) <= ng)  THEN
!            IF (sptau(nzl(i)) <= break_vec(i,k+1)) THEN
!
!              dx = sptau(nzl(i)) - break_vec(i,k)
!              results(i,j,nzl(i)) = coef_vec(i,1,k) + dx     *(coef_vec(i,2,k) +   &
!                    dx*0.5*(coef_vec(i,3,k) + dx/3.0 * coef_vec(i,4,k)))
!              nzl(i) = nzl(i) + 1
!            ELSE
!              ! the level nzl(i) is ready
!              ldone(i) = .TRUE.
!              nzcount    = nzcount + 1
!            ENDIF
!          ELSE
!            ! the level nzl(i) is ready
!            ldone(i) = .TRUE.
!            nzcount    = nzcount + 1
!          ENDIF
        ENDIF

      ENDDO
    ENDDO

  ENDDO

  SELECT CASE (yname)
  CASE ('RELHUM    ')
    DO i = istartpar, iendpar
      results(i,j,1:ng) = MAX (0.0_ireals, MIN(100.0_ireals,results(i,j,1:ng)))
    ENDDO
  CASE ('QV        ','QC        ','QI        ',    &
        'QR        ','QS        ','QG        ',    &
        'QH        ','QNCLOUD   ','QNICE     ',    &
        'QNRAIN    ','QNSNOW    ','QNGRAUPEL ','QNHAIL    ' )
    DO i = istartpar, iendpar
      results(i,j,1:ng) = MAX (0.0_ireals, results(i,j,1:ng) )
    ENDDO
  END SELECT

!------------------------------------------------------------------------------
! End of the subroutine
!------------------------------------------------------------------------------

END SUBROUTINE spline2D

!==============================================================================
!==============================================================================
!+ Module procedure in src_output for the horizontal smoothing
!------------------------------------------------------------------------------

SUBROUTINE smooth_pmsl ( pmsl, hsurf, ie, je )

!------------------------------------------------------------------------------
! special smoothing of pmsl in mountainous terrain
!------------------------------------------------------------------------------

! Parameter list:
INTEGER (KIND=iintegers), INTENT (IN)    :: ie, je
REAL    (KIND=ireals),    INTENT (IN)    :: hsurf(ie,je)

REAL    (KIND=irealgrib), INTENT (INOUT) :: pmsl(ie,je)

! Local Variables
REAL    (KIND=ireals)                    :: zp, hsurf_max, wgt, pmsl_sm(ie,je)
INTEGER (KIND=iintegers)                 :: i, j, ii, jj, n, nmean

!------------------------------------------------------------------------------

! Begin subroutine smooth_pmsl

  PRINT *,' smoothing pmsl over mountainous terrain '

  pmsl_sm(:,:) = REAL (pmsl(:,:), ireals)
  hsurf_max = 750.0_ireals

  DO n = 1,2
    DO j = 1, je
      DO i = 1, ie
        IF (hsurf(i,j) >= hsurf_max) THEN
          IF (hsurf(i,j) >= 1000.0_ireals) THEN
            nmean = 10
          ELSE
            nmean = 5
          ENDIF
          pmsl_sm(i,j) = 0.0_ireals
          zp = 0.0_ireals
          DO jj = j-nmean,j+nmean
            DO ii = i-nmean,i+nmean
              IF ( (ii >= 1) .AND. (ii <= ie) .AND.     &
                   (jj >= 1) .AND. (jj <= je) ) THEN
                IF (hsurf(ii,jj) < hsurf_max) THEN
                  wgt = 2.0_ireals
                ELSE
                  wgt = 1.0_ireals
                ENDIF
                pmsl_sm(i,j) = pmsl_sm(i,j) + wgt * REAL (pmsl(ii,jj), ireals)
                zp = zp + wgt
              ENDIF
            ENDDO
          ENDDO
          pmsl_sm(i,j) = pmsl_sm(i,j) / zp
        ENDIF
      ENDDO
    ENDDO
    pmsl(:,:) = REAL (pmsl_sm(:,:), irealgrib)
  ENDDO

END SUBROUTINE smooth_pmsl

!==============================================================================
!==============================================================================
!+ Module procedure in src_output for the horizontal smoothing
!------------------------------------------------------------------------------

SUBROUTINE smooth_geopot ( geopot, hsurf, ie, je )

!------------------------------------------------------------------------------
! special smoothing of geopot in mountainous terrain
!------------------------------------------------------------------------------

! Parameter list:

INTEGER (KIND=iintegers), INTENT (IN)    :: ie, je
REAL    (KIND=ireals),    INTENT (IN)    :: hsurf(ie,je)
REAL    (KIND=irealgrib), INTENT (INOUT) :: geopot(ie,je)

! Local Variables
REAL    (KIND=ireals)                    :: zp, hsurf_max, wgt, geopot_diff, &
                                            geopot_sm(ie,je)
INTEGER (KIND=iintegers)                 :: i, j, ii, jj, n, nmean

!------------------------------------------------------------------------------

! Begin subroutine smooth_geopot

  PRINT*,' smoothing of geopotential height over mountainous terrain '

  geopot_sm(:,:) = REAL (geopot(:,:), ireals)

  DO n = 1,2
    DO j = 1,je
      DO i = 1,ie
        geopot_diff = g * hsurf(i,j) - REAL(geopot(i,j), ireals)
        IF (geopot_diff > 1000.0_ireals) THEN
          nmean = 10
        ELSE
          nmean = 5
        ENDIF
        geopot_sm(i,j) = 0.0_ireals
        zp = 0.0_ireals
        DO jj = j-nmean,j+nmean
          DO ii = i-nmean,i+nmean
            IF ( (ii >= 1) .AND. (ii <= ie) .AND.     &
                 (jj >= 1) .AND. (jj <= je) ) THEN
              geopot_diff = g * hsurf(ii,jj) - REAL(geopot(ii,jj), ireals)
              wgt = 1.0_ireals
              geopot_sm(i,j) = geopot_sm(i,j) + wgt*REAL(geopot(ii,jj), ireals)
              zp = zp + wgt
            ENDIF
          ENDDO
        ENDDO
        geopot_sm(i,j) = geopot_sm(i,j) / zp
      ENDDO
    ENDDO
    geopot(:,:) = REAL (geopot_sm(:,:), irealgrib)
  ENDDO

END SUBROUTINE smooth_geopot

!==============================================================================
!==============================================================================

SUBROUTINE calc_sdi( sdi_1, sdi_2 )

!------------------------------------------------------------------------------
!
! Description:
!   calculation of the 2 supercell detection indices (SDI)
!
! Method:
!   defined in:
!   Wicker L, J. Kain, S. Weiss and D. Bright, A Brief Description of the
!          Supercell Detection Index, (available from
!   http://www.spc.noaa.gov/exper/Spring_2005/SDI-docs.pdf)
!
!------------------------------------------------------------------------------

REAL    (KIND=ireals),    DIMENSION(1:ie,1:je), INTENT(OUT) :: sdi_1, sdi_2

INTEGER (KIND=iintegers) :: k_center
INTEGER (KIND=iintegers) :: d_idx_x, d_idx_y, d_idx_z
INTEGER (KIND=iintegers) :: i, j, k
REAL    (KIND=ireals)    :: dx, dx_aequ, dy

REAL    (KIND=ireals), ALLOCATABLE, DIMENSION(:,:) :: &
  w_mean,       & ! <w>          mean over box in 'height' k_center
  w2_mean,      & ! <w'w'>       mean over box in 'height' k_center
  zeta_mean,    & ! <zeta>       mean over box in 'height' k_center
  zeta2_mean,   & ! <zeta'zeta'> mean over box in 'height' k_center
  helic_mean      ! <w'zeta'>    mean over box in 'height' k_center

REAL (KIND=ireals), ALLOCATABLE, DIMENSION(:,:,:) :: &
  zeta,            & ! vorticity
  w_s                ! vertical velocity at scalar position

CHARACTER (LEN=80)       :: yzerrmsg
INTEGER (KIND=iintegers) :: izerror

INTEGER (KIND=iintegers) :: kzdims(24)

REAL (KIND=ireals) :: helic_w_corr          ! Correlation coefficient
REAL (KIND=ireals) :: zeta_vert_mean(ie,je)
REAL (KIND=ireals) :: w_crit_SDI2   (ie,je) ! Criterion, if SDI2 = 0 or not

REAL (KIND=ireals) :: EPS = 1.0E-20_ireals

!------------------------------------------------------------------------------

  izerror = 0_iintegers

  ! definition of the integration box:
  !     [ i_center-d_idx_x .. i_center+d_idx_x ]
  !   * [ j_center-d_idx_y .. j_center+d_idx_y ]
  !   * [ k_center-d_idx_z .. k_center+d_idx_z ]
  d_idx_x  = 3
  d_idx_y  = 3
  d_idx_z  = 6

  k_center = 30


  ! Allocations:
  ALLOCATE( zeta       ( 1:ie, 1:je, 1:ke) )
  ALLOCATE( w_s        ( 1:ie, 1:je, 1:ke) )
  ALLOCATE( w_mean     ( 1:ie, 1:je ) )
  ALLOCATE( w2_mean    ( 1:ie, 1:je ) )
  ALLOCATE( zeta_mean  ( 1:ie, 1:je ) )
  ALLOCATE( zeta2_mean ( 1:ie, 1:je ) )
  ALLOCATE( helic_mean ( 1:ie, 1:je ) )

  ! to prevent errors at the boundaries, set some fields to 0:
  sdi_1(:,:)  = 0.0_ireals
  sdi_2(:,:)  = 0.0_ireals
  zeta(:,:,:) = 0.0_ireals
  w_s(:,:,:)  = 0.0_ireals

  ! consistency checks:

  IF ( ( d_idx_x > nboundlines ) .OR. ( d_idx_y > nboundlines ) ) THEN
    yzerrmsg="integration box is too big in horizontal direction!"
    ! if such a big value for d_idx_x or d_idx_y is really needed, then you must increase
    ! nboundlines
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'calc_sdi')
  END IF

  IF ( ( k_center - d_idx_z < 2 ) .OR. ( k_center + d_idx_z > ke ) ) THEN
    yzerrmsg="integration box is too big in vertical direction!"
    CALL model_abort (my_cart_id, 100, yzerrmsg, 'calc_sdi')
  END IF

  ! calculate vorticity and vertical velocity at scalar points:

  dx_aequ = r_earth * (pi/180.0_ireals) * dlon
  dy      = r_earth * (pi/180.0_ireals) * dlat

  DO k = k_center - d_idx_z, k_center + d_idx_z
    DO j = jstart, jend
      dx = dx_aequ * crlat(j,1)
      DO i = istart, iend
        zeta(i,j,k ) = (   ( v(i+1,j,  k,itl) + v(i+1,j-1,k,itl) )         &
          &              - ( v(i-1,j,  k,itl) + v(i-1,j-1,k,itl) ) )       &
          &                        * 0.5_ireals / dx                       &
          &          - (   ( u(i,  j+1,k,itl) + u(i-1,j+1,k,itl) )         &
          &              - ( u(i,  j-1,k,itl) + u(i-1,j-1,k,itl) ) )       &
          &                        * 0.5_ireals / dy

        w_s(i,j,k) = 0.5_ireals * ( w(i,j,k,itl) + w(i,j,k-1,itl) )

      END DO
    END DO
  END DO

  kzdims(1:24)=(/ke,ke,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
  CALL exchg_boundaries                                                    &
    (50+nnew, sendbuf, isendbuflen, imp_reals, icomm_cart, num_compute,    &
     ie, je, kzdims, jstartpar, jendpar, nbl_exchg, nboundlines,           &
     my_cart_neigh, lperi_x, lperi_y, l2dim,                               &
     20000+nexch_tag, .FALSE.   , ncomm_type, izerror, yzerrmsg,           &
     zeta(:,:,:), w_s(:,:,:) )

  ! (exchange of w_s would not be necessary, if it is calculated also at the boundary lines)

  ! --- calculate mean values over the integration box: ------

  CALL mean_over_box( w_s,  w_mean,    k_center, d_idx_x, d_idx_y, d_idx_z,   &
                      crlat, sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_over_box( zeta, zeta_mean, k_center, d_idx_x, d_idx_y, d_idx_z,   &
                      crlat, sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  ! (no exchange needed for w_mean, zeta_mean)

  ! --- calculate covariances over the integration box: -----------

  CALL mean_cov_over_box( w_s,  w_mean,    zeta, zeta_mean, helic_mean,   &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_cov_over_box( w_s,  w_mean,    w_s,  w_mean,    w2_mean,      &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  CALL mean_cov_over_box( zeta, zeta_mean, zeta, zeta_mean, zeta2_mean,   &
                          k_center, d_idx_x, d_idx_y, d_idx_z, crlat,     &
                          sqrtg_r_s, ie, je, ke, istart, iend, jstart, jend)

  ! calculate SDI_1, SDI_2:
  ! call to vectorized version of vert_avg
  CALL vert_avg( zeta_vert_mean,zeta, sqrtg_r_s, ie, je, ke, istart, iend,  &
                 jstart, jend, k_center, d_idx_z)

  ! The meaning of 'w>0' in Wicker et al. is not completely clear, I assume
  ! the following:
  CALL vert_avg( w_crit_SDI2, w_s, sqrtg_r_s, ie, je, ke, istart, iend,     &
                 jstart, jend, k_center, d_idx_z)

  DO j = jstart, jend
    dx = dx_aequ * crlat(j,1)
    DO i = istart, iend

      IF ( ( w2_mean(i,j) > EPS ) .AND. ( zeta2_mean(i,j) > EPS ) ) THEN

        helic_w_corr = helic_mean(i,j) / SQRT(w2_mean(i,j) * zeta2_mean(i,j))

        sdi_1(i,j) = helic_w_corr * zeta_vert_mean(i,j)

        IF ( w_crit_SDI2(i,j) > 0 ) THEN
          sdi_2(i,j) = helic_w_corr * ABS( zeta_vert_mean(i,j) )
        ELSE
          sdi_2(i,j) = 0.0
        END IF

      ELSE
        sdi_1(i,j) = 0.0
        sdi_2(i,j) = 0.0
      END IF

    END DO
  END DO

  DEALLOCATE( zeta, w_s, w_mean, zeta_mean, w2_mean, zeta2_mean, helic_mean )

END SUBROUTINE calc_sdi

!==============================================================================

END MODULE src_output
