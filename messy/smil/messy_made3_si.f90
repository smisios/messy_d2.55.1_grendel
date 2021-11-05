#include "../bmil/messy_main_ppd_bi.inc"

!###############################################################################
!##  NOTE: Comment format `!>' is used for generation of code documentation   ##
!##        from source file via `doxygen'.                                    ##
!##  NOTE2: In order to use doxygen on this file, one of the `update_xte'     ##
!##         subroutine headers has to be temporarily removed!                 ##
!###############################################################################

!> \page interface Couplings and output
!>
!>   For an example setup see \c messy/nml/REF_2000_MADE3 (can be obtained from
!>   the MADE3 maintainer upon request).
!>
!>   \section couplings Couplings
!>
!>     \subsection aeropt AEROPT
!>
!>       The submodel AEROPT calculates aerosol optical depth (AOD), an overall
!>       (i.e., integrated over all modes) aerosol single scattering albedo, and
!>       an overall aerosol asymmetry parameter. These properties can
!>       subsequently be used as an input for the radiation submodel \link rad
!>       RAD \endlink, and they can be output for diagnostic purposes. For the
!>       required calculations, AEROPT uses pre-calculated lookup tables of the
!>       aerosol particle extinction cross section, single scattering albedo,
!>       and asymmetry parameter as a function of the Mie size parameter,
!>       \f$y\f$, and the particle complex refractive index,
!>       \f$\varepsilon\f$. The latter two are computed from the MADE3 aerosol
!>       properties. For each mode, \f$y\f$ is calculated as \f$y=\frac{\pi
!>       D_{\rm p}}{\lambda}\f$, with the median mode diameter \f$D_{\rm p}\f$
!>       and the wavelength \f$\lambda\f$, and \f$\varepsilon\f$ is determined
!>       as the volume average over the refractive indices of the mode's
!>       components. These, again, are tabulated for a number of species as a
!>       function of wavelength.
!>
!>       The above-mentioned tables have to be supplied separately for the
!>       shortwave and for the longwave part of the spectrum. Their paths must
!>       be specified via the CTRL namelist in \c aeropt.nml. For MADE3, only
!>       one set of files is currently available. Hence, unless you create your
!>       own lookup tables, the appropriate line for the AEROPT CTRL namelist
!>       is:
!>       \code
!>       read_lut_sets(n) = m, "$INPUTDIR_MESSY/../EVAL2.3/messy/raw/aeropt/aerorad_sw.nc", "$INPUTDIR_MESSY/../EVAL2.3/messy/raw/aeropt/aerorad_lw.nc",
!>       \endcode
!>       The indices \c n and \c m have to be set in accordance with the rest of
!>       your \c aeropt.nml.
!>
!>       In order to make AEROPT aware of the MADE3 aerosol, a further entry is
!>       required in the CPL namelist in \c aeropt.nml. Use
!>       \code
!>       read_aero_sets(l) = "MADE3" , T, "", "", "", T, F, "made3_gp", "wetradius", "dryradius", "", m, 0.55,,,,,,,,,,,,,,,,
!>       \endcode
!>       here, set the index \c l in accordance with the rest of your \c
!>       aeropt.nml, use the index \c m that you defined in the CTRL namelist
!>       entry, and adapt the "0.55,,,,,,,,,,,,,,,", i.e., up to 16 wavelengths
!>       (in micrometers) at which AOD is calculated for diagnostic purposes,
!>       according to your needs.
!>       \note Leave the rest of the settings as they are, unless you have read
!>         and understood the code where they are used.
!>
!>     \subsection cloud CLOUD
!>
!>       The submodel CLOUD will be able to diagnose CCN from the MADE3 aerosol,
!>       and will eventually also be able to use MADE3 aerosol as an input for
!>       cloud formation. These features are under development and will be
!>       documented here as soon as they have been tested.
!>
!>     \subsection ddep DDEP
!>
!>       The submodel DDEP automatically uses the parameters \c wetradius, \c
!>       densaer, and \c sigma supplied by MADE3 (see \ref made3_gp) to compute
!>       dry deposition of aerosol particles (excl. sedimentation, for which
!>       \ref sedi is used). The only MADE3-specific user setting is in the
!>       CPL namelist in \c ddep.nml: add the full names of the tracers for
!>       which you want to output the dry deposition flux to the variable \c
!>       outflux (separated by ";"). A full name is made up of a base name, an
!>       underscore, and a subname, as explained in \ref tracer_gp.
!>
!>     \subsection mecca MECCA
!>
!>       The submodel MECCA does not require any specific settings for MADE3 in
!>       \c mecca.nml, as the interaction happens via tracers, which are
!>       specified via the \link cpl CPL namelist \endlink in \link namelist
!>       \c made3.nml \endlink.
!>       \note MADE3 has so far only been run (and evaluated) with a simplified
!>         chemistry mechanism. The code (and pdf file) for this \"meccanism\"
!>         can be created via \c xmecca, which can be run from within \c
!>         messy/util/xconfig. The selection of reactions is done via the batch
!>         file \c messy/mbm/caaba/mecca/batch/user_contributed/simple_MADE.bat
!>         (in \c xmecca, first choose no. 3 (\"SUBDIRECTORY:
!>         ./user_contributed\"), then no. 12). A replacement file for certain
!>         reactions is used by this batch file (\c
!>         messy/mbm/caaba/mecca/rpl/mim1-simple_MADE.rpl).
!>
!>     \subsection rad RAD
!>
!>       The submodel RAD will be able to compute the effect of the MADE3
!>       aerosol on solar and terrestrial radiation. This feature is under
!>       development and will be documented here as soon as it has been tested.
!>
!>     \subsection scav SCAV
!>
!>       The SCAV routines for cloud processing of aerosol particles differ for
!>       different aerosol submodels in some cases. Hence, the SCAV code
!>       contains a number of MADE3-specific parts. Nevertheless, a few
!>       parameters can be controlled via the namelist file \c scav.nml.
!>
!>       Specifically, these are:
!>       - \c frac_resnum in the CTRL namelist, the factor applied to the
!>         number of cloud residual aerosol particles upon cloud
!>         evaporation. This factor is applied to all the MADE3 residual modes
!>         (\c as, \c am, \c cs, \c cm; for mode naming see \ref tracer_gp),
!>         rather than only to the hydrophilic mode with the largest particles
!>         (default in SCAV). The value used in the evaluation of EMAC with
!>         MADE3 is
!>         \code
!>         frac_resnum = 0.1,
!>         \endcode
!>       - \c aermod_gp in the CPL namelist has to be set to
!>         \code
!>         aermod_gp = 'made3',
!>         \endcode
!>       - \c made3_params in the CPL namelist should be set to
!>         \code
!>         made3_params = 1,2,3,4,5,6,7,8,9,9,8,'cm',8,4,
!>         \endcode
!>         unless you make changes to the MADE3 and/or SCAV code. The first nine
!>         numbers are the indices of the MADE3 modes in the order \c ks, \c km,
!>         \c ki, \c as, \c am, \c ai, \c cs, \c cm, \c ci (for mode naming see
!>         \ref tracer_gp), followed by the total number of modes. These values
!>         should only be changed in case the MADE3 code is modified. Next, the
!>         index and subname of the mode for dummy tracers are supplied. SCAV
!>         may create such tracers to ensure its chemical mechanism works
!>         properly. These settings are mainly "cosmetic". The last but one
!>         parameter is currently not used, and the last one specifies how many
!>         residual modes shall be considered. This should only be changed if
!>         the MADE3-specific SCAV code is modified.
!>       - \c out_string and \c out_string_aer in the CPL namelist determine for
!>         which species and tracers SCAV creates wet deposition flux channel
!>         objects. For full diagnostic output, and in order to enable the
!>         (post-simulation) calculation of residence times for all species,
!>         they should be set to
!>         \code
!>         out_string = 'H2SO4_l,HSO4m_l;SO4mm_l;HNO3_l;NO3m_l;NH3_l;NH4p_l;HCl_l;Clm_l',
!>         out_string_aer = 'SO4_ks;SO4_km;SO4_ki;SO4_as;SO4_am;SO4_ai;SO4_cs;SO4_cm;SO4_ci;NH4_ks;NH4_km;NH4_ki;NH4_as;NH4_am;NH4_ai;NH4_cs;NH4_cm;NH4_ci;NO3_ks;NO3_km;NO3_ki;NO3_as;NO3_am;NO3_ai;NO3_cs;NO3_cm;NO3_ci;Na_ks;Na_km;Na_ki;Na_as;Na_am;Na_ai;'//&
!>&'Na_cs;Na_cm;Na_ci;Cl_ks;Cl_km;Cl_ki;Cl_as;Cl_am;Cl_ai;Cl_cs;Cl_cm;Cl_ci;POM_ks;POM_km;POM_ki;POM_as;POM_am;POM_ai;POM_cs;POM_cm;POM_ci;BC_km;BC_ki;BC_am;BC_ai;BC_cm;BC_ci;DU_am;DU_ai;DU_cm;DU_ci;N_ks;N_km;N_ki;N_as;N_am;N_ai;N_cs;N_cm;N_ci',
!>         \endcode
!>         (where the latter should actually be on one line)
!>
!>       \note The SCAV&ndash;MADE3 interplay has only been thoroughly tested
!>         with the SCAV settings in \c
!>         messy/nml/REF01_2000_MADE3/scav.nml. Hence, it cannot be assumed that
!>         other settings that work with GMXe, for instance, will also work with
!>         MADE3, even if they should not be submodel dependent.
!>
!>     \subsection sedi SEDI
!>
!>       The submodel SEDI automatically uses the parameters \c wetradius, \c
!>       densaer, and \c sigma supplied by MADE3 (see \ref made3_gp) to compute
!>       sedimentation of aerosol particles. It does not require any
!>       MADE3-specific settings in \c sedi.nml.
!>
!>     \subsection tracer TRACER
!>
!>       Tracer properties are set in messy_made3_si::made3_new_tracer using
!>       subroutines supplied by the submodel TRACER. The relevant settings for
!>       coupling to the base model and to other submodels are those for:
!>       - transport
!>         - \c I_ADVECT = \c ON for all aerosol tracers, except aerosol water,
!>           i.e., all aerosol tracers except aerosol water are advected
!>       - deposition
!>         - \c I_DRYDEP, \c I_SEDI are all set to \c ON, i.e., wet
!>           and dry deposition, as well as sedimentation, are switched on for
!>           all aerosol tracers
!>       - scavenging
!>         - \c I_SCAV = \c ON for all aerosol tracers, i.e., scavenging is
!>           switched on for all aerosol tracers
!>         - \c I_AEROSOL_SOL = \c ON for the soluble and mixed modes, \c OFF
!>           for the insoluble modes, i.e., the latter are considered
!>           hydrophobic, whereas the former are considered hydrophilic
!>         - \c I_AEROSOL_HETICE = \c ON for the mixed and insoluble modes,
!>           \c OFF for the soluble modes
!>
!>   \section output Output
!>
!>     Depending on the settings in your \c channel.nml, you can have the
!>     following output of MADE3 aerosol related quantities. For an example
!>     setup with typical output used for aerosol evaluation see \c
!>     messy/nml/REF01_2000_MADE3/channel.nml.
!>
!>     \note If you choose output of non-instantaneous values in \c channel.nml
!>       the following variable names will have a further extension, like 
!>       \"\c _ave\" for example.
!>
!>     \subsection aeropt_MADE3 *aeropt_MADE3.nc
!>
!>       - \c aot_[ls]w_BNN, \c gamma_sw_BNN, \c omega_sw_BNN contain AOD,
!>         asymmetry parameter, and single scattering albedo per radiation band
!>         (extension \"\c BNN\" for band no. \c NN); some of these can be used
!>         to couple the MADE3 aerosol to the radiation submodel \link rad RAD
!>         \endlink
!>       - \c aot_opt_<mode>_<wavelength>_<species> contain diagnostic AOD,
!>         where \c &lt;mode&gt; (either \"\c MNN\" or \"\c TOT\") corresponds
!>         to the MADE3 mode with index \c NN or the sum over the modes (\"\c
!>         TOT\"), \c &lt;wavelength&gt; is determined by what you request via
!>         the CPL namelist in \link aeropt \c aeropt.nml \endlink and \c
!>         &lt;species&gt; stands for one of the AEROPT particle component
!>         classes (or the sum over all, in which case \c &lt;species&gt; is
!>         \"\c total\"):
!>         - \c bc (black carbon)
!>         - \c du (mineral dust)
!>         - \c h2o (aerosol water)
!>         - \c oc (organic carbon)
!>         - \c ss (sea spray)
!>         - \c waso (other water soluble components)
!>         .
!>         \note These variables contain the AOD \b per \b layer, i.e., for
!>           comparison to satellite measurements, for instance, the values of
!>           all layers have to be summed.
!>
!>     \subsection ddep_gp *ddep_gp.nc
!>
!>       - \c ddepflux_<tracer> contain the dry deposition fluxes of the
!>         \c &lt;tracer&gt;s (see \ref tracer_gp) in units of molecules (or
!>         particles, for number tracers) per square meter per second
!>       - \c ddepfluxsum_<tracer> contain the time integrals of the above,
!>         converted to mass (i.e., [kg m<sup>-2</sup>]) in case of mass tracers
!>
!>     \subsection made3_gp *made3_gp.nc
!>
!>       In this subsection \c &lt;modename&gt; always stands for the
!>       two-character combinations explained in \ref tracer_gp.
!>       - \c wetrad_<modename> contain the modal median particle wet radii
!>       - \c dryrad_<modename> contain the modal median particle dry radii
!>       - \c densaer_<modename> contain the modal particle densities
!>       - \c rhhist_<modename> contain the modal deliquescence histories
!>       - \c sigma contains the MADE3 mode widths
!>       - \c burden_* and \c sink_* should contain data to determine the
!>         conversion lifetime of hydrophobic to hydrophilic black carbon (BC)
!>         and mineral dust (DU), respectively, but these have not been used so
!>         far, so they might not contain the expected data
!>
!>     \subsection offemis *offemis.nc
!>
!>       - \c &lt;source&gt;_<tracer> contain the emissions of \c
!>         &lt;tracer&gt; from \c &lt;source&gt; as defined via \c offemis.nml
!>         in units of molecules per square meter per second or molecules per
!>         cubic meter per second in case of "2D" or "Nx2D" emissions,
!>         respectively; for the different emissions types see \cite
!>         Kerkweg2006; note that no output can be created for "3D" emissions
!>       - \c &lt;source&gt;_<tracer>_vind contain the indices of the vertical
!>         layers to which the above emissions were assigned (only present if
!>         emissions \b are actually assigned to multiple layers)
!>
!>     \subsection scav_gp *scav_gp.nc
!>
!>       This subsection assumes that ice phase "chemistry" is disabled. In case
!>       it is enabled, further output variables will have to be considered,
!>       especially in the calculation of total wet deposition fluxes.
!>
!>       - \c wetflx_<aggregate> contain the sum of the wet deposition fluxes of
!>         the components that make up the \c &lt;aggregate&gt; in the chemical
!>         mechanism of SCAV; these are only available for the \c
!>         &lt;aggregate&gt;s \"\c sulfate\", \"\c nitrate\", and \"\c ammoni\";
!>         for total wet deposition fluxes of the \c &lt;aggregate&gt; the
!>         corresponding \c wetflx_aer_*_&lt;tracer&gt; (see below; e.g., all \c
!>         wetflx_aer_*_NH4_&lt;modename&gt;, where \c &lt;modename&gt; is
!>         defined as in \ref made3_gp) have to be added
!>       - \c wetflx_sum_<aggregate> contain the time integrals of the \c
!>         wetflx_<aggregate>
!>       - \c wetflx_ls_HCl_l, \c wetflx_ls_Clm_l, \c wetflx_cv_HCl_l, and \c
!>         wetflx_cv_Clm_l have to be summed to obtain the chloride wet
!>         deposition flux that corresponds to the \c wetflx_<aggregate> (again,
!>         if you are interested in the total chloride wet deposition flux you
!>         have to add the \c wetflx_aer_*_Cl_&lt;modename&gt;, see below)
!>       - \c wetflx_aer_<cloud>_<tracer>, where \c &lt;cloud&gt; is either
!>         \"\c ls\" or \"\c cv\" for large-scale or convective clouds,
!>         respectively, contain the wet deposition fluxes of the aerosol
!>         tracers requested in the CPL namelist in \link SCAV \c scav.nml
!>         \endlink that did not end up in the corresponding \c
!>         wetflx_<cloud>_<species>_l, which are only created for species
!>         included in the SCAV liquid chemistry mechanism (described above only
!>         for chloride; note that the \c wetflx_<aggregate> described above are
!>         sums of the corresponding \c wetflx_*_&lt;species&gt;_l); hence, if
!>         there are any \c wetflx_*_&lt;species&gt;_l that correspond to a \c
!>         &lt;tracer&gt; base name (see \ref tracer_gp) these have to be added
!>         to the sum of the \c wetflx_aer_*_&lt;basename&gt;_&lt;modename&gt;
!>         over all \c &lt;modename&gt; (defined as in \ref made3_gp) in order
!>         to obtain the total wet deposition flux of the component described by
!>         the \c &lt;basename&gt;; for an example see below
!>       - \c frac_evap_<hydro>_<cloud>_<tracer>_mNN contain the fractions of
!>         each \c &lt;tracer&gt; (see \ref tracer_gp) that are assigned to the
!>         modes with indeces \c NN after evaporation of constituents or
!>         hydrometeors of clouds of type \c &lt;cloud&gt;; cloud constituents
!>         or hydrometeors are identified by \c &lt;hydro&gt;, which can take
!>         the values \"\c snow\" for snow flakes, \"\c rain\" for rain drops,
!>         \"\c iwc\" for in-cloud ice crystals, and \"\c lwc\" for in-cloud
!>         droplets
!>
!>       \par Example:
!>         To obtain the total \c SO4 wet deposition flux, sum up \c
!>         wetflx_sulfate, \c wetflx_aer_ls_SO4_*, and \c wetflx_aer_cv_SO4_*,
!>         which is equivalent to summing \c wetflx_ls_H2SO4_l, \c
!>         wetflx_ls_HSO4m_l, \c wetflx_ls_SO4mm_l, \c wetflx_cv_H2SO4_l, \c
!>         wetflx_cv_HSO4m_l, \c wetflx_cv_SO4mm_l, wetflx_aer_ls_SO4_*, and \c
!>         wetflx_aer_cv_SO4_*.
!>
!>     \subsection sedi_gp *sedi_gp.nc
!>
!>       - \c sediflux_<tracer> contain the sedimentation fluxes of the
!>         \c &lt;tracer&gt;s (see \ref tracer_gp) in units of molecules (or
!>         particles, for number tracers) per square meter per second
!>       - \c loss_<tracer> contain the time integrals of the above, converted
!>         to mass (i.e., [kg m<sup>-2</sup>]) in case of mass tracers
!>
!>     \subsection tracer_gp *tracer_gp.nc
!>
!>       - \c &lt;basename&gt;_<modename> contain the mixing ratios (in [mol
!>         mol<sup>-1</sup>] for mass tracers, in [mol<sup>-1</sup>] for number
!>         tracers) of the aerosol components \c &lt;basename&gt; in the modes
!>         \c &lt;modename&gt;; aerosol tracer \c &lt;basename&gt;s are \"\c
!>         SO4\", \"\c NH4\", \"\c NO3\", \"\c Na\", \"\c Cl\", \"\c POM\", \"\c
!>         BC\", \"\c DU\", and \"\c H2O\" (cf. the \ref index); the \c
!>         &lt;modename&gt;s consist of two letters indicating the size range
!>         and particle class, i.e., first letter
!>         - \"\c k\" for Aitken,
!>         - \"\c a\" for accumulation, and
!>         - \"\c c\" for coarse mode,
!>         .
!>         and second letter
!>         - \"\c s\" for soluble,
!>         - \"\c m\" for mixed, and
!>         - \"\c i\" for insoluble particles,
!>         .
!>         e.g., \"\c SO4_ks\" for sulfate in the soluble Aitken mode.
!>

!> \brief MADE3 interface module.

!> \authors Axel Lauer, DLR Oberpfaffenhofen, 2001-2003 (axel.lauer@dlr.de)
!>   - inclusion of MADE in ECHAM4
!>   - creation of MADE interface to the MESSy framework
!> \authors Valentina Aquila, DLR Oberpfaffenhofen, 2008/09
!>          (valentina.aquila@dlr.de)
!>   - development of MADE-in on the basis of MADE (and under the name MADE)
!> \authors Christopher Kaiser, DLR Oberpfaffenhofen, 2012-2016
!>          (christopher.kaiser@dlr.de)
!>   - adaptation of MADE-in to MESSy2
!>   - development of MADE3 on the basis of MADE-in

!> \version 2.1
!>   - redesign of emissions coupling
!> \version 2.0
!>   - first working version of the MADE3 interface

!> \todo
!>   - A number of simplifications should be possible with the new CBLK (and
!>     other variables') layout, e.g. sfac etc.
!>   - Check FIXMEs.
!>   - Test code to use submodel TENDENCY.

MODULE messy_made3_si

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM, STRLEN_VLONG
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_main_tracer,        ONLY: ON, OFF, STRLEN_FNAME
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY, t_reset_par ! op_mr_20190225
! DEBUG+
!!$  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
! DEBUG-
  USE messy_made3,              ONLY: nspec, nmod, ngas, dim1_cblk, dim2_cblk &
                                    , modstr

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC NULL

  ! Type definitions

  !> Aerosol specs type
  TYPE t_aerspec
     !> Tracer base name
     CHARACTER(LEN=STRLEN_MEDIUM)  :: bname   = ''
     !> Tracer long name
     CHARACTER(LEN=STRLEN_VLONG)   :: lname   = ''
     !> Molar mass [g mol<sup>-1</sup>]
     REAL(dp), POINTER             :: molm    => NULL()
     !> Density [kg m<sup>-3</sup>]
     REAL(dp), POINTER             :: dens    => NULL()
     !> Enable advection
     INTEGER                       :: advec   = ON
  END TYPE t_aerspec

  !> Mode specs type
  TYPE t_mode
     !> Tracer subname
     CHARACTER(LEN=STRLEN_MEDIUM)  :: sname   = ''
     !> Mode long name
     CHARACTER(LEN=STRLEN_VLONG)   :: lname   = ''
     !> Solubility flag
     INTEGER                       :: sol     = ON
     !> Enable higher ice scavenging efficiency
     INTEGER                       :: hice    = ON
  END TYPE t_mode

  !> Gas specs type
  TYPE t_gas
     !> Tracer base name
     CHARACTER(LEN=STRLEN_MEDIUM)  :: bname   = ''
     !> Tracer subname
     CHARACTER(LEN=STRLEN_MEDIUM)  :: sname   = ''
     !> Tracer long name
     CHARACTER(LEN=STRLEN_VLONG)   :: lname   = ''
     !> Molar mass [g mol<sup>-1</sup>]
     REAL(dp), POINTER             :: molm    => NULL()
     !> Effective Henry's law coefficient [mol l<sup>-1</sup> atm<sup>-1</sup>]
     REAL(dp)                      :: henry   = 0.0_dp
     !> Coefficient for dry reaction with surface [-]
     REAL(dp)                      :: dreac   = 0.0_dp
  END TYPE t_gas

  !> Gas coupling type
  TYPE t_gas_cpl
     !> Tracer base name
     CHARACTER(LEN=STRLEN_MEDIUM)  :: bname   = ''
     !> Tracer subname
     CHARACTER(LEN=STRLEN_MEDIUM)  :: sname   = ''
     !> Coupled chemistry submodel (empty for no coupling)
     CHARACTER(LEN=STRLEN_MEDIUM) :: chemmod = ''
  END TYPE t_gas_cpl

  !> Type for emissions assignment via namelist
  TYPE t_emis_cpl
     !> Identifier
     CHARACTER(LEN=STRLEN_MEDIUM)  :: type      = ''
     !> Channel from which to read the flux
     CHARACTER(LEN=STRLEN_CHANNEL) :: inchannel = ''
     !> Object that contains the flux
     CHARACTER(LEN=STRLEN_CHANNEL) :: inobject  = ''
     !> Tracer (full)name to which to assign the flux
     CHARACTER(LEN=STRLEN_VLONG)   :: tracer    = ''
     !> Scaling factor
     REAL(dp)                      :: scale     = 0.0_dp
  END TYPE t_emis_cpl

  !> Emission data type
  TYPE t_emis
     !> Identifier
     CHARACTER(LEN=STRLEN_MEDIUM)              :: type  = ''
     !> Number of tracers associated with this identifier
     INTEGER                                   :: ntrac = 0
     !> Indices of tracers associated with this identifier
     INTEGER, DIMENSION(:), POINTER            :: idt
     !> Emission fluxes assigned to the tracers
     TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER :: flux
     !> Scale factors for the emission fluxes
     REAL(dp), DIMENSION(:), POINTER           :: scale
  END TYPE t_emis

  ! Memory management
  ! ... define 5D pointers (for non-default (channel) memory management)
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: p5_wetrad   => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: p5_dryrad   => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: p5_densaer  => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: p5_rhhist   => NULL()
  ! ... define 4D pointers
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: p4_wetrad   => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: p4_dryrad   => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: p4_densaer  => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: p4_rhhist   => NULL()
! op_ck_20120612+
! op_ck_20120612  No longer required if cloud stuff is done by CLOUD
!!$  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: aernum_4d     => NULL()
!!$  ! ... define cloud droplet number concentration 3D pointer
!!$  REAL(dp), DIMENSION(:,:,:),     POINTER :: cdnc       => NULL()
! op_ck_20120612-
! op_va_20090225+
  ! ... for calculation of lifetime of external mixtures
  REAL(dp), DIMENSION(:,:,:),     POINTER :: p3_BCextIN    => NULL()
  REAL(dp), DIMENSION(:,:,:),     POINTER :: p3_BCextSINK  => NULL()
  REAL(dp), DIMENSION(:,:,:),     POINTER :: p3_DUextIN    => NULL()
  REAL(dp), DIMENSION(:,:,:),     POINTER :: p3_DUextSINK  => NULL()
! op_va_20090225-
! op_cb_20180329+
  ! ... for online DU emissions (ai, ci) output after scaling and threshold
  REAL(dp), DIMENSION(:,:),       POINTER :: p2_DUemfluxai   => NULL()
  REAL(dp), DIMENSION(:,:),       POINTER :: p2_DUemfluxci   => NULL()
! op_cb_20180329-
! op_cb_20180626+
  ! fetched with get_channel_object in init-coupling
  REAL(dp), DIMENSION(:,:), POINTER :: oromea   => NULL()
! op_cb_20180626-
! DEBUG+
!!$  TYPE(PTR_3D_ARRAY), DIMENSION(2,2) :: p3arr_dbg_cblk_1
!!$  TYPE(PTR_3D_ARRAY), DIMENSION(2,2) :: p3arr_dbg_cblk_2
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: p3_dbg_dry_1 => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: p3_dbg_dry_2 => NULL()
! DEBUG-

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

!!$  ! radiation related
!!$  INTEGER :: num_opt_wavelens = 0 ! number of wavelengths for optional
!!$                                  ! diagnostic output

  ! CPL - NAMELIST
  !> Gas phase coupling (NOTE: gas_cpl(5) is ProdH2SO4 here, not SOA!)
  TYPE(t_gas_cpl), DIMENSION(ngas) :: gas_cpl
  ! Emissions coupling
  INTEGER, PARAMETER :: N_MAX_EMIS_CPL = 99
  !> Emissions coupling
  TYPE(t_emis_cpl), DIMENSION(N_MAX_EMIS_CPL) :: emis_cpl
  ! op_cb_20180626+
  !> Enable upper threshold for orography, for DU online emissions
  TYPE(t_reset_par), PUBLIC :: orogr_thr_4DUemis = t_reset_par(.FALSE.,1.e4_dp)
  ! op_cb_20180626-
! op_va_20090225+ !TEST on ext.BC lifetime
  !> Enable calculation of externally mixed BC and DU lifetimes
  LOGICAL :: l_bctime = .FALSE.
! op_va_20090225- !TEST on ext.BC lifetime
  !> Names of unnecessary tracers
  CHARACTER(LEN=STRLEN_FNAME) :: notrac(nspec*nmod)

! op_cb_20180626+
!!$  NAMELIST /CPL/ gas_cpl, emis_cpl, l_bctime, notrac
  NAMELIST /CPL/ gas_cpl, emis_cpl, orogr_thr_4DUemis, l_bctime, notrac
! op_cb_20180626-

  ! Available emission types
  INTEGER, PARAMETER :: N_EMIS_AVAIL = 4
  CHARACTER(STRLEN_MEDIUM), DIMENSION(N_EMIS_AVAIL), PARAMETER :: emis_avail = &
       (/'SS  ', 'DU  ', 'OPOM', 'SOA '/)
  ! Emissions coupled
  LOGICAL, DIMENSION(N_EMIS_AVAIL) :: l_emis
  ! Number of coupled emission types
  INTEGER :: n_emis

  ! Internal index for ProdH2SO4 tracer
  INTEGER, PARAMETER :: i_ph2so4 = 5

  ! Aerosol species
  TYPE(t_aerspec)    :: aerspec(nspec)
  ! Mode specs
  TYPE(t_mode)       :: mode(nmod)
  ! Gas phase species
  TYPE(t_gas)        :: gspec(ngas)
  ! Tracer indices (NOTE: itrac(5,gas) is ProdH2SO4 here, not SOA!)
  INTEGER            :: itrac(dim1_cblk,dim2_cblk)

  ! Pointer to geometric standard deviations (i.e., mode widths)
  REAL(dp), DIMENSION(:), POINTER :: p1_sigma
! op_ck_20120612+
! op_ck_20120612  No longer required if cloud stuff is done by CLOUD
!!$  REAL(dp), POINTER               :: actscheme_str
! op_ck_20120612-

  ! Emissions data
  TYPE(t_emis), DIMENSION(:), ALLOCATABLE :: emis

  ! SUBROUTINES
  ! Initializes MADE3 parameters
  PUBLIC  :: made3_initialize
  ! Creates MADE3 tracers
  PUBLIC  :: made3_new_tracer
  ! Creates the MADE3 channel (incl. objects)
  PUBLIC  :: made3_init_memory
  ! Checks and initializes requested couplings
  PUBLIC  :: made3_init_coupling
  ! Distributes online SS and DU emissions among species and modes
  PUBLIC  :: made3_vdiff
  ! Entry point to MADE3 core routines
  PUBLIC  :: made3_physc
  ! Deallocates memory allocated in #made3_init_memory
  PUBLIC  :: made3_free_memory   ! deallocate radius field

  ! Note: the following lines are required only for creation of the
  ! documentation via doxygen
  PRIVATE :: made3_read_nml_cpl
  PRIVATE :: echam_to_made3
  PRIVATE :: made3_to_echam
  PRIVATE :: update_xte

CONTAINS

!-------------------------------------------------------------------------------

  !> \brief Initializes MADE3 parameters
  !> \details Reads CTRL and CPL namelists, broadcasts results, calls
  !>   messy_made3::made3_initialize_core, and sets further MADE3 constants that
  !>   are used in the interface.

  SUBROUTINE made3_initialize

    ! External parameters and subroutines
    USE messy_main_constants_mem, ONLY: STRLEN_SHORT
    USE messy_main_tools,         ONLY: find_next_free_unit
    USE messy_main_mpi_bi,        ONLY: p_bcast, p_parallel_io, p_io
    USE messy_main_blather_bi,    ONLY: error_bi
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,   ONLY: mtend_get_handle, mtend_register &
         , mtend_id_tracer
#endif
    USE messy_main_tracer,        ONLY: OFF
    ! op_cb_20171107: add epsilon_max
    USE messy_made3,              ONLY: sigma, alpha, diff, epsilon_max       &
         , rset_nucsize, l_eqsam, l_coag, l_cond, l_nuc, l_rename             &
         , nspec, ngas, mw, rho, i_mwaero, i_mwgas                            &
         , akn, akns, sooti, acc, accs, sootj, cor, cors, sootc               &
         , i_so4, i_nh4, i_no3, i_ss, i_cl, i_pom, i_bc, i_bctag, i_du, i_h2o &
         , i_h2so4, i_nh3, i_hno3, i_hcl, i_soa                               &
         , made3_initialize_core, made3_read_nml

    IMPLICIT NONE
    INTRINSIC TRIM
  
    ! Local
    CHARACTER(LEN=*), PARAMETER :: substr = 'made3_initialize'
    INTEGER                     :: iou        ! I/O unit
    INTEGER                     :: status     ! error status flag
    INTEGER                     :: js, jg, je ! loop indices

    ! Read MADE3 CTRL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100, 200)
       CALL made3_read_nml(status, iou)
       IF (status /= 0)  CALL error_bi('error in made3_read_nml', substr)
    END IF

    ! Broadcast results
    CALL p_bcast(sigma(:),       p_io)
    CALL p_bcast(alpha(:,:),     p_io)
    CALL p_bcast(diff(:),        p_io)
    CALL p_bcast(epsilon_max,    p_io)  ! op_cb_20171107
    CALL p_bcast(rset_nucsize%l, p_io)
    CALL p_bcast(rset_nucsize%v, p_io)
    CALL p_bcast(l_eqsam,        p_io)
    CALL p_bcast(l_coag,         p_io)
    CALL p_bcast(l_cond,         p_io)
    CALL p_bcast(l_nuc,          p_io)
    CALL p_bcast(l_rename,       p_io)

    ! Set default gas tracer names
    gas_cpl(i_h2so4)%bname  = 'H2SO4'
    gas_cpl(i_nh3)%bname    = 'NH3'
    gas_cpl(i_hno3)%bname   = 'HNO3'
    gas_cpl(i_hcl)%bname    = 'HCl'
    gas_cpl(i_ph2so4)%bname = 'ProdH2SO4'

    ! Set gas names
    gspec(i_h2so4)%bname = 'H2SO4'
    gspec(i_h2so4)%lname = 'Sulfuric acid (gas)'
    gspec(i_nh3)%bname   = 'NH3'
    gspec(i_nh3)%lname   = 'Ammonia (gas)'
    gspec(i_hno3)%bname  = 'HNO3'
    gspec(i_hno3)%lname  = 'Nitric acid (gas)'
    gspec(i_hcl)%bname   = 'HCl'
    gspec(i_hcl)%lname   = 'Hydrochloric acid (gas)'
    gspec(i_soa)%bname   = 'SOApre'
    gspec(i_soa)%lname   = 'SOA precursors (gas)'

    ! Set gas parameters
    DO jg = 1, ngas
       gspec(jg)%molm => mw(jg,i_mwgas)
    END DO
    gspec(i_h2so4)%henry = 8.7e11_dp
    gspec(i_nh3)%henry   = 58.0_dp   ! from MECCA1
    gspec(i_hno3)%henry  = 1.e4_dp   ! from MECCA1
    gspec(i_hno3)%dreac  = 1.0_dp    ! from MECCA1
    gspec(i_hcl)%henry   = 1.e14_dp  ! from MECCA1
    gspec(i_hcl)%dreac   = 1.0_dp    ! from MECCA1

    ! Do not exclude any tracers by default
    notrac = ''

    ! Read MADE3 CPL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL made3_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('error in made3_read_nml_cpl', substr)
    END IF

    ! Broadcast results
    ! ... chemistry
    CALL p_bcast(gas_cpl(i_h2so4)%bname,    p_io)
    CALL p_bcast(gas_cpl(i_h2so4)%sname,    p_io)
    CALL p_bcast(gas_cpl(i_h2so4)%chemmod,  p_io)
    CALL p_bcast(gas_cpl(i_nh3)%bname,      p_io)
    CALL p_bcast(gas_cpl(i_nh3)%sname,      p_io)
    CALL p_bcast(gas_cpl(i_nh3)%chemmod,    p_io)
    CALL p_bcast(gas_cpl(i_hno3)%bname,     p_io)
    CALL p_bcast(gas_cpl(i_hno3)%sname,     p_io)
    CALL p_bcast(gas_cpl(i_hno3)%chemmod,   p_io)
    CALL p_bcast(gas_cpl(i_hcl)%bname,      p_io)
    CALL p_bcast(gas_cpl(i_hcl)%sname,      p_io)
    CALL p_bcast(gas_cpl(i_hcl)%chemmod,    p_io)
    CALL p_bcast(gas_cpl(i_ph2so4)%bname,   p_io)
    CALL p_bcast(gas_cpl(i_ph2so4)%sname,   p_io)
    CALL p_bcast(gas_cpl(i_ph2so4)%chemmod, p_io)
    ! ... emissions
    DO je = 1, N_MAX_EMIS_CPL
       CALL p_bcast(emis_cpl%type,      p_io)
       CALL p_bcast(emis_cpl%inchannel, p_io)
       CALL p_bcast(emis_cpl%inobject,  p_io)
       CALL p_bcast(emis_cpl%tracer,    p_io)
       CALL p_bcast(emis_cpl%scale,     p_io)
    END DO
    ! ... others
    CALL p_bcast(orogr_thr_4DUemis%l, p_io) ! op_cb_20180626
    CALL p_bcast(orogr_thr_4DUemis%v, p_io) ! op_cb_20180626
! op_va_20090225+
    CALL p_bcast(l_bctime,       p_io)
! op_va_20090225-
    CALL p_bcast(notrac(:),      p_io)

    ! Set MADE3 parameters (MUST be called AFTER `made3_read_nml'!!!)
    CALL made3_initialize_core

    ! Set aerosol species names
    aerspec(i_so4)%bname   = 'SO4'
    aerspec(i_so4)%lname   = 'Sulfate'
    aerspec(i_nh4)%bname   = 'NH4'
    aerspec(i_nh4)%lname   = 'Ammonium'
    aerspec(i_no3)%bname   = 'NO3'
    aerspec(i_no3)%lname   = 'Nitrate'
    aerspec(i_ss)%bname    = 'Na'
    aerspec(i_ss)%lname    = 'Sea spray components w/o Cl'
    aerspec(i_cl)%bname    = 'Cl'
    aerspec(i_cl)%lname    = 'Chloride'
    aerspec(i_pom)%bname   = 'POM'
    aerspec(i_pom)%lname   = 'particulate org. matter'
    aerspec(i_bc)%bname    = 'BC'
    aerspec(i_bc)%lname    = 'Black carbon'
    aerspec(i_bctag)%bname = 'BCtag'
    aerspec(i_bctag)%lname = 'Tagged black carbon'
    aerspec(i_du)%bname    = 'DU'
    aerspec(i_du)%lname    = 'Mineral dust'
    aerspec(i_h2o)%bname   = 'H2O'
    aerspec(i_h2o)%lname   = 'Aerosol water'

    ! Set aerosol species parameters
    DO js = 1, nspec
       aerspec(js)%molm => mw(js,i_mwaero)
       aerspec(js)%dens => rho(js)
    END DO
! op_ck_20131128 FIXME: Why is advection switched off, but convection, vertical
!                       diffusion, and removal are switched on? H2O is diagnosed
!                       from the other tracers anyway...
! op_ck_20140612        Nobody really knows a reason, so switching on advection
!                       as well should be tested in the future.
    aerspec(i_h2o)%advec = OFF

    ! Set mode names
    mode(akn)%sname   = 'ks'
    mode(akn)%lname   = 'soluble Aitken mode'
    mode(akns)%sname  = 'km'
    mode(akns)%lname  = 'mixed Aitken mode'
    mode(sooti)%sname = 'ki'
    mode(sooti)%lname = 'insoluble Aitken mode'
    mode(acc)%sname   = 'as'
    mode(acc)%lname   = 'soluble accumulation mode'
    mode(accs)%sname  = 'am'
    mode(accs)%lname  = 'mixed accumulation mode'
    mode(sootj)%sname = 'ai'
    mode(sootj)%lname = 'insoluble accumulation mode'
    mode(cor)%sname   = 'cs'
    mode(cor)%lname   = 'soluble coarse mode'
    mode(cors)%sname  = 'cm'
    mode(cors)%lname  = 'mixed coarse mode'
    mode(sootc)%sname = 'ci'
    mode(sootc)%lname = 'insoluble coarse mode'

    ! Set mode parameters
    mode(akn)%hice    = OFF
    mode(acc)%hice    = OFF
    mode(cor)%hice    = OFF
    mode(sooti)%sol   = OFF
    mode(sootj)%sol   = OFF
    mode(sootc)%sol   = OFF

    ! Set default tracer indices to -1 to recognize unset tracers
    itrac         = -1

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    ! op_ck_20160810+
    ! register all tracers (workaround)
    CALL mtend_register(my_handle, mtend_id_tracer)
    ! op_ck_20160810-
#endif

  END SUBROUTINE made3_initialize

!-------------------------------------------------------------------------------

  !> \brief Creates MADE3 tracers
  !> \details Creates tracers for MADE3 aerosol components and numbers (and
  !>   dummy gas phase tracers in case coupling to the gas phase is switched
  !>   off), and sets their attributes.

  SUBROUTINE made3_new_tracer

    ! External parameters and subroutines
    USE messy_main_constants_mem,   ONLY: STRLEN_MEDIUM, STRLEN_VLONG
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer,          ONLY: new_tracer, set_tracer, get_tracer  &
         , AIR, AEROSOL, OFF, ON, AMOUNTFRACTION, NUMBERDENSITY             &
         , I_ADVECT, I_DRYDEP, I_SEDI, I_SCAV                               &
         , I_AEROSOL_MODE, I_AEROSOL_SOL, I_AEROSOL_HETICE                  &
         , R_MOLARMASS, R_PSS  , R_DRYREAC_SF, R_AEROSOL_DENSITY            &
         , S_AEROSOL_MODEL
    USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi
    USE messy_made3,                ONLY: ngas, i_num, gas

    IMPLICIT NONE
    INTRINSIC TRIM, SIZE

    ! Default settings for new tracers:
    !   see smcl/messy_main_tracer.f90

    ! Local
    LOGICAL                      :: dont = .FALSE.
    INTEGER                      :: status, idt, jg, jm, js, jt
    CHARACTER(LEN=*), PARAMETER  :: substr   = 'made3_new_tracer'
    CHARACTER(LEN=STRLEN_MEDIUM) :: bname    = ''
    CHARACTER(LEN=STRLEN_MEDIUM) :: sname    = ''
    CHARACTER(LEN=STRLEN_VLONG)  :: lname    = ''
    CHARACTER(LEN=STRLEN_VLONG)  :: halt_msg = ''


    CALL start_message_bi(modstr, 'check/define tracers', substr)

    ! Scavenging flags:
    !   I_SCAV           = OFF ---> don't calculate scavenging at all
    !                    = ON  ---> calculate impact scavenging
    !   I_AEROSOL_SOL    = ON  ---> calculate nucleation scavenging
    !                               (only if I_SCAV == ON)
    !   I_AEROSOL_HETICE = ON  ---> higher ice scav. coefficient because of
    !                               heterogeneous ice scav.

    ! Define aerosol tracers
    loop_species: DO js = 1, nspec

       ! Set tracer base and long names
       bname = TRIM(aerspec(js)%bname)
       lname = TRIM(aerspec(js)%lname)//' mixing ratio'

       loop_modes: DO jm = 1, nmod

          ! Set tracer subname and halt message
          sname    = TRIM(mode(jm)%sname)
          halt_msg = TRIM(substr)//' ('//TRIM(bname)//'_'//TRIM(sname)//')'

          ! Do not create unnecessary tracers
          dont = .FALSE.
          DO jt = 1, SIZE(notrac)
             IF (TRIM(notrac(jt)) == TRIM(bname)//'_'//TRIM(sname)) THEN
                dont = .TRUE.
                EXIT
             END IF
          END DO
          IF (dont) CYCLE

          ! Create tracer and set attributes
          CALL new_tracer(status, GPTRSTR, bname, modstr                       &
               , idx=idt, subname=sname                                        &
               , longname=TRIM(lname)//' - '//TRIM(mode(jm)%lname)             &
               , unit='mol mol-1', medium=AEROSOL, quantity=AMOUNTFRACTION)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_ADVECT, i=aerspec(js)%advec)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=jm)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=mode(jm)%sol)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE  &
               , i=mode(jm)%hice)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=aerspec(js)%molm)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY &
               , r=aerspec(js)%dens)
          CALL tracer_halt(halt_msg, status)
          CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
          CALL tracer_halt(halt_msg, status)

          ! Store tracer index
          itrac(js,jm) = idt

       END DO loop_modes

    END DO loop_species

    ! Define number tracers

    loop2_modes: DO jm = 1, nmod

       ! Set tracer subname and halt message
       sname = TRIM(mode(jm)%sname)
       halt_msg = TRIM(substr)//' (N_'//TRIM(sname)//')'

       CALL new_tracer(status, GPTRSTR, 'N', modstr                           &
            , idx=idt, subname=sname                                          &
            , longname='Aerosol number mixing ratio - '//TRIM(mode(jm)%lname) &
         , unit='mol-1', medium=AEROSOL, quantity=NUMBERDENSITY)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_SEDI, i=ON)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_MODE, i=jm)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_SOL, i=mode(jm)%sol)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, I_AEROSOL_HETICE, i=mode(jm)%hice)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=1.0_dp)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, R_AEROSOL_DENSITY, r=1.0_dp)
       CALL tracer_halt(halt_msg, status)
       CALL set_tracer(status, GPTRSTR, idt, S_AEROSOL_MODEL, s=modstr)
       CALL tracer_halt(halt_msg, status)

       ! Store tracer index
       itrac(i_num,jm) = idt

    END DO loop2_modes

    ! Define tracers for gas species that are not coupled to chemistry
    loop_gases: DO jg = 1, ngas

       IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE

       IF (TRIM(gas_cpl(jg)%chemmod) == '') THEN

          ! Check if tracer has already been defined
          CALL get_tracer(status, GPTRSTR, TRIM(gas_cpl(jg)%bname), idx=idt)
          IF (status .NE. 0) THEN
             bname = TRIM(gspec(jg)%bname)
          ELSE
             itrac(jg,gas) = idt
             CYCLE
          END IF

          lname = TRIM(gspec(jg)%lname)//' mixing ratio'

          CALL new_tracer(status, GPTRSTR, TRIM(bname), TRIM(modstr)   &
               , idx=idt, longname='MADE3 dummy '//TRIM(lname) &
               , unit='mol mol-1', medium=AIR)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)
          CALL set_tracer(status, GPTRSTR, idt, I_DRYDEP, i=ON)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)
          CALL set_tracer(status, GPTRSTR, idt, I_SCAV, i=ON)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)
          CALL set_tracer(status, GPTRSTR, idt, R_MOLARMASS, r=gspec(jg)%molm)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)
          CALL set_tracer(status, GPTRSTR, idt, R_PSS  , r=gspec(jg)%henry)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)
          CALL set_tracer(status, GPTRSTR, idt, R_DRYREAC_SF, r=gspec(jg)%dreac)
          CALL tracer_halt(TRIM(substr)//' ('//TRIM(bname)//')', status)

          itrac(jg,gas) = idt

       END IF

    END DO loop_gases

    CALL end_message_bi(modstr,'check/define tracers', substr)

  END SUBROUTINE made3_new_tracer

!-------------------------------------------------------------------------------

  !> \brief Creates the MADE3 channel (incl. objects)
  !> \details Creates new dimension for modes and representations that use this
  !>   dimension; creates the MADE3 channel (\c made3_gp) and the objects for
  !>   - wet radius (\c wetradius, \c wetrad_xx),
  !>   - dry radius (\c dryradius, \c dryrad_xx),
  !>   - particle density (\c densaer, \c densaer_xx),
  !>   - deliquescence history (\c rhhist, \c rhhist_xx),
  !>   - BC and DU aging related variables, and
  !>   - modal standard deviations (\c sigma).
  !>   .
  !>   The `xx' stands for the mode name abbreviation as described under \ref
  !>   interface. The subroutine also sets the channel object `sigma' to the
  !>   values from ::messy_made3.

  SUBROUTINE made3_init_memory

    ! External subroutines and parameters
    USE messy_main_constants_mem,      ONLY: STRLEN_MEDIUM
    USE messy_main_tracer,             ONLY: STRLEN_FNAME
    USE messy_main_blather_bi,         ONLY: start_message_bi, end_message_bi
    USE messy_main_grid_def_mem_bi,    ONLY: nproma, nlev, ngpblks
    USE messy_main_timer,              ONLY: lstart
    USE messy_main_channel,            ONLY: new_channel, new_channel_object &
         , new_attribute, STRLEN_CHANNEL
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO        &
         , set_representation_decomp, IRANK, PIOTYPE_COL, repr_def_axes
    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: GP_3D_MID, DC_GP                &
         , DIMID_LON, DIMID_LAT, DIMID_LEV, DC_BC                            &
! op_cb_20180405+
         , GP_2D_HORIZONTAL                                                  &
! op_cb_20180405-
         , gp_nseg, gp_start, gp_cnt, gp_meml, gp_memu                       &
         , SCALAR
    USE messy_made3,                   ONLY: nmod, sigma
! DEBUG+
!!$         , dbg_cblk_1, dbg_cblk_2, dbg_dry_1, dbg_dry_2
! DEBUG-

    IMPLICIT NONE
    INTRINSIC TRIM

    ! Local
    CHARACTER(len=*), PARAMETER           :: substr = ''
    CHARACTER(len=STRLEN_CHANNEL)         :: cname  = ''
    CHARACTER(len=STRLEN_MEDIUM)          :: bname  = ''
    CHARACTER(len=STRLEN_FNAME)           :: oname  = ''
    INTEGER                               :: status
    INTEGER                               :: jm
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: mem => NULL()
    INTEGER                               :: DIMID_NMODE
    INTEGER                               :: REPR_MADE3_4D_NMOD
    INTEGER                               :: REPR_MADE3_1D
!!$
!!$    ! Parallel I/O decomposition
!!$    INTEGER                          :: nseg = 0
!!$    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
!!$    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
!!$    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
!!$    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()


    CALL start_message_bi(modstr, 'Initialize memory', substr)
    
    ALLOCATE(p5_wetrad(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    p5_wetrad(:,:,:,:,:) = 0._dp
    ALLOCATE(p5_dryrad(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    p5_dryrad(:,:,:,:,:) = 0._dp
    ALLOCATE(p5_densaer(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    p5_densaer(:,:,:,:,:) = 0._dp
    ALLOCATE(p5_rhhist(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    p5_rhhist(:,:,:,:,:) = 0._dp

    ! --- 1) Construct the MADE3 channel ---------------------------------------

    ! Create new dimension: modes
    CALL new_dimension(status, DIMID_NMODE, 'MADE3_NMODE', nmod)
    CALL channel_halt(substr, status)

    ! Create new representations

    ! ... 4D with lon, lat, lev, and mode
    CALL new_representation(status, REPR_MADE3_4D_NMOD, &
         'REPR_MADE3_4D_NMOD'    &
         , rank = 4, link = 'xxxx', dctype = DC_GP               &
         , dimension_ids = (/ &
            _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NMODE) /) &
         , ldimlen       = (/ &
            _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_   &        ! E: 3,1,4,2
                            , _IY_XYZN_, _IZ_XYZN_ /)       & ! C: 3,1,2,4
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
         )
    CALL channel_halt(substr, status)

    ! ... 1D for modes
    CALL new_representation(status, REPR_MADE3_1D  &
         , 'REPR_MADE3_1D'                         &
         , rank = 1, link = 'x---', dctype = DC_BC &
         , dimension_ids = (/ DIMID_NMODE /)       &
         , ldimlen       = (/ AUTO /)              &
         , axis = 'N---'                           &
         )
    CALL channel_halt(substr, status)

!!$    ! Set decompositions for parallel I/O
!!$! op_ck_20131203 NOTE: This has not been checked! (It was probably copied from
!!$!                      messy/echam5/smil/messy_gmxe_e5.f90)
!!$
!!$    ! ... REPR_MADE3_4D_NMOD
!!$    nseg = gp_nseg
!!$    ALLOCATE(start(nseg,IRANK))
!!$    ALLOCATE(cnt(nseg,IRANK))
!!$    ALLOCATE(meml(nseg,IRANK))
!!$    ALLOCATE(memu(nseg,IRANK))
!!$    
!!$    start(:,:) = gp_start(:,:)
!!$    cnt(:,:)   = gp_cnt(:,:)
!!$    meml(:,:)  = gp_meml(:,:)
!!$    memu(:,:)  = gp_memu(:,:)
!!$    
!!$    cnt(:,3)   = nmod
!!$    memu(:,3)  = nmod
!!$    
!!$    CALL set_representation_decomp(status, REPR_MADE3_4D_NMOD &
!!$         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
!!$    CALL channel_halt(substr, status)
!!$    
!!$    DEALLOCATE(start) ; NULLIFY(start)
!!$    DEALLOCATE(cnt)   ; NULLIFY(cnt)
!!$    DEALLOCATE(meml)  ; NULLIFY(meml)
!!$    DEALLOCATE(memu)  ; NULLIFY(memu)
!!$
!!$    ! ... REPR_MADE3_1D
!!$    nseg = gp_nseg
!!$    ALLOCATE(start(nseg,IRANK))
!!$    ALLOCATE(cnt(nseg,IRANK))
!!$    ALLOCATE(meml(nseg,IRANK))
!!$    ALLOCATE(memu(nseg,IRANK))
!!$    
!!$    start(:,:) = 1
!!$    cnt(:,:)   = 1
!!$    meml(:,:)  = 1
!!$    memu(:,:)  = 1
!!$    
!!$    cnt(:,1)   = nmod
!!$    memu(:,1)  = nmod
!!$    
!!$    CALL set_representation_decomp(status, REPR_MADE3_1D &
!!$         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
!!$    CALL channel_halt(substr, status)
!!$    
!!$    DEALLOCATE(start) ; NULLIFY(start)
!!$    DEALLOCATE(cnt)   ; NULLIFY(cnt)
!!$    DEALLOCATE(meml)  ; NULLIFY(meml)
!!$    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! Create MADE3 channel
    cname = TRIM(modstr)//'_gp'
    CALL new_channel(status, TRIM(cname), reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    ! --- 2) Create channel objects --------------------------------------------

    ! Wet radii
    
    mem => p5_wetrad(:,:,:,:,1)
    CALL new_channel_object(status, TRIM(cname), 'wetradius' &
         , p4=p4_wetrad, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE3_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'wetradius', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'wetradius' &
         , 'long_name', c='Wet particle radius')
    CALL channel_halt(substr, status)

    bname = 'wetrad'
! op_ck_20131203 FIXME: Do we need these objects? And is GP_3D_MID the right
!                       representation for them?
    loop_modes: DO jm = 1, nmod

       mem => p5_wetrad(_RI_XYZN_(:,:,:,jm),:)
       oname = TRIM(bname)//'_'//TRIM(mode(jm)%sname)
       CALL new_channel_object(status, TRIM(cname), TRIM(oname), mem=mem)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Wet particle radius - '//TRIM(mode(jm)%lname))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units', c='m')
       CALL channel_halt(substr, status)

    END DO loop_modes

    ! Dry radii
    
    mem => p5_dryrad(:,:,:,:,1)
    CALL new_channel_object(status, TRIM(cname), 'dryradius' &
         , p4=p4_dryrad, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE3_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'dryradius', 'units', c='m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'dryradius' &
         , 'long_name', c='Dry particle radius')
    CALL channel_halt(substr, status)

    bname = 'dryrad'
! op_ck_20131203 FIXME: Do we need these objects? And is GP_3D_MID the right
!                       representation for them?
    loop2_modes: DO jm = 1, nmod

       mem => p5_dryrad(_RI_XYZN_(:,:,:,jm),:)
       oname = TRIM(bname)//'_'//TRIM(mode(jm)%sname)
       CALL new_channel_object(status, TRIM(cname), TRIM(oname), mem=mem)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Dry particle radius - '//TRIM(mode(jm)%lname))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units', c='m')
       CALL channel_halt(substr, status)

    END DO loop2_modes

    ! Particle density

    bname = 'densaer'
    mem => p5_densaer(:,:,:,:,1)
    CALL new_channel_object(status, TRIM(cname), TRIM(bname) &
         , p4=p4_densaer, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE3_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(bname), 'units', c='kg m-3')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(bname) &
         , 'long_name', c='Particle density')
    CALL channel_halt(substr, status)

    loop3_modes: DO jm = 1, nmod

       mem => p5_densaer(_RI_XYZN_(:,:,:,jm),:)
       oname = TRIM(bname)//'_'//TRIM(mode(jm)%sname)
       CALL new_channel_object(status, TRIM(cname), TRIM(oname), mem=mem)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Particle density - '//TRIM(mode(jm)%lname))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units', c='kg m-3')
       CALL channel_halt(substr, status)

    END DO loop3_modes

    ! Deliquescence history

    bname = 'rhhist'
    mem => p5_rhhist(:,:,:,:,1)
    CALL new_channel_object(status, TRIM(cname), TRIM(bname) &
         , p4=p4_rhhist, mem=mem, lrestreq=.TRUE., reprid=REPR_MADE3_4D_NMOD)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(bname), 'units', c='')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(bname) &
         , 'long_name', c='Deliquescence history')
    CALL channel_halt(substr, status)

! op_ck_20131203 FIXME: Do we need these objects? And is GP_3D_MID the right
!                       representation for them?
    loop4_modes: DO jm = 1, nmod

       mem => p5_rhhist(_RI_XYZN_(:,:,:,jm),:)
       oname = TRIM(bname)//'_'//TRIM(mode(jm)%sname)
       CALL new_channel_object(status, TRIM(cname), TRIM(oname), mem=mem)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Deliquescence history - '//TRIM(mode(jm)%lname))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units', c='')
       CALL channel_halt(substr, status)

    END DO loop4_modes

    ! (Fixed) standard deviations of log-normal distributions

    CALL new_channel_object(status, TRIM(cname), 'sigma' &
         , p1=p1_sigma, lrestreq=.TRUE., reprid=REPR_MADE3_1D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'sigma' &
         , 'long_name', c='Standard deviations of MADE3 modes')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), 'sigma', 'units', c='')
    CALL channel_halt(substr, status)

! op_va_20090225+
    IF (l_bctime) THEN ! accumulated element (lav=true)
                       ! op_ck_20120604       ^^^^^^^^ <- deprecated,
                       ! op_ck_20120604                   now set in channel.nml

       oname = 'burden_BCext'
       CALL new_channel_object(status, TRIM(cname), TRIM(oname)     &
            , p3=p3_BCextIN, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units' &
            , c='mol mol-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Incoming externally mixed BC')
       CALL channel_halt(substr, status)

       oname = 'sink_BCext'
       CALL new_channel_object(status, TRIM(cname), TRIM(oname) &
            , p3=p3_BCextSINK, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units' &
            , c='mol mol-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Sink of externally mixed BC (BCin-BCout)')
       CALL channel_halt(substr, status)

       oname = 'burden_DUext'
       CALL new_channel_object(status, TRIM(cname), TRIM(oname) &
            , p3=p3_DUextIN, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units' &
            , c='mol mol-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Incoming externally mixed DU')
       CALL channel_halt(substr, status)

       oname = 'sink_DUext'
       CALL new_channel_object(status, TRIM(cname), TRIM(oname) &
            , p3=p3_DUextSINK, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname), 'units' &
            , c='mol mol-1')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, TRIM(cname), TRIM(oname) &
            , 'long_name', c='Sink of externally mixed DU (DUin-DUout)')
       CALL channel_halt(substr, status)
   
    ENDIF

! op_cb_20180329+
    ! Channel objects for DU online emissions (after scaling and threshold
    ! calculation)
    oname = 'DU_emflux_ai'
    CALL new_channel_object(status, TRIM(cname), TRIM(oname) &
         , p2=p2_DUemfluxai, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(oname), &
           'long_name', c='Dust emission flux ai')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(oname), &
         'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    oname = 'DU_emflux_ci'
    CALL new_channel_object(status, TRIM(cname), TRIM(oname) &
         , p2=p2_DUemfluxci, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(oname), &
           'long_name', c='Dust emission flux ci')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, TRIM(cname), TRIM(oname), &
         'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)
! op_cb_20180329-

! op_va_20090225+
! DEBUG+
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_1_cl_as', p3=p3arr_dbg_cblk_1(1,1)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_1_cl_cs', p3=p3arr_dbg_cblk_1(1,2)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_1_num_as', p3=p3arr_dbg_cblk_1(2,1)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_1_num_cs', p3=p3arr_dbg_cblk_1(2,2)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_2_cl_as', p3=p3arr_dbg_cblk_2(1,1)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_2_cl_cs', p3=p3arr_dbg_cblk_2(1,2)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_2_num_as', p3=p3arr_dbg_cblk_2(2,1)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_cblk_2_num_cs', p3=p3arr_dbg_cblk_2(2,2)%ptr)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_dry_1', p3=p3_dbg_dry_1)
!!$    CALL new_channel_object(status, TRIM(cname), 'dbg_dry_2', p3=p3_dbg_dry_2)
!!$    ALLOCATE(dbg_cblk_1(2,2,nproma))
!!$    ALLOCATE(dbg_cblk_2(2,2,nproma))
!!$    ALLOCATE(dbg_dry_1(nproma))
!!$    ALLOCATE(dbg_dry_2(nproma))
! DEBUG-

    ! --- 3) Initialization of channel object values
    p1_sigma = sigma

!!$    IF (l_cloudphysics) THEN
!!$
!!$! op_ck_20120612+
!!$! op_ck_20120612  Activation scheme will be set via CLOUD namelist. In this case
!!$!                 SCALAR can be removed from the USE statements.
!!$       ! number of selected aerosol activation scheme (for Lohmann cloud
!!$       ! module, ECHAM5 version)
!!$
!!$       CALL new_channel_object(status, modstr, 'actscheme' &
!!$            , p0=actscheme_str, reprid=SCALAR)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'actscheme' &
!!$            , 'long_name', c='selected aerosol activation scheme')
!!$       CALL channel_halt(substr, status)
!!$
!!$       actscheme_str = REAL(act_scheme)
!!$
! op_ck_20120612  CDNC will be computed by CLOUD
!!$       ! number concentration of activated aerosols (CDNC)
!!$       CALL new_channel_object(status, modstr, 'aercdnc' &
!!$            , p3=cdnc, lrestreq=.TRUE.)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'aercdnc', 'units', c='1/m3')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'aercdnc' &
!!$            , 'long_name', c='activated aerosols (CDNC)')
!!$       CALL channel_halt(substr, status)
!!$
! op_ck_20120612  CLOUD takes aerosol number concentrations from tracers
!!$       CALL new_channel_object(status, modstr, 'aernum' &
!!$            , p4=aernum_4d, reprid=REPR_MADE3_4D_NMOD, lrestreq=.TRUE.)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'aernum', 'units' &
!!$            , c='1/m3')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr, 'aernum' &
!!$            , 'long_name', c='aerosol number concentration')
!!$       CALL channel_halt(substr, status)
!!$
! op_ck_20120612  No longer required if cloud stuff is done by CLOUD. In this
!                 case, lstart can also be removed from the USE statements.
!!$! new channel objects will be reset to 0.0 automatically
!!$!       IF (lstart) THEN
!!$!          cdnc(:,:,:)       = 0.0_dp
!!$!       END IF
!!$! op_ck_20120612-
!!$
!!$    END IF ! if l_cloudphysics

    CALL end_message_bi(modstr,'Initialize memory', substr)

  END SUBROUTINE made3_init_memory

!-------------------------------------------------------------------------------

  !> \brief Checks and initializes requested couplings
  !> \details Checks if tracers and channel objects requested for couplings
  !>   (i.e., chemistry and emissions) are available and stores the
  !>   corresponding tracer indices or associates the corresponding pointers,
  !>   respectively.  Writes relevant coupling information to standard output.

  SUBROUTINE made3_init_coupling

    ! External parameters and subroutines
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_channel,          ONLY: STRLEN_CHANNEL, get_channel_object &
                                         , get_attribute
! op_cb_20180627+
    USE messy_main_channel_error_bi, ONLY: channel_halt
! op_cb_20180627-
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , error_bi, info_bi
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, ti_gp
    USE messy_main_tracer,           ONLY: STRLEN_FNAME, get_tracer &
                                         , NUMBERDENSITY
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tools,            ONLY: int2str
    ! op_ck_20160810+
    ! op_ck_20160810  Registration is too late here!
!!$#ifdef MESSYTENDENCY
!!$    USE messy_main_tendency_bi,   ONLY: mtend_register, mtend_id_tracer
!!$    USE messy_made3,              ONLY: ngas, gas, i_h2so4, dim1_cblk, dim2_cblk
!!$#else
    ! op_ck_20160810-
    USE messy_made3,                 ONLY: ngas, gas, i_h2so4
    ! op_ck_20160810+
!!$#endif
    ! op_ck_20160810-

    IMPLICIT NONE
    INTRINSIC TRIM, SIZE, ADJUSTL

    CHARACTER(LEN=*), PARAMETER      :: substr = 'made3_init_coupling'
    CHARACTER(LEN=STRLEN_MEDIUM)     :: chemmod = ''
    CHARACTER(LEN=STRLEN_FNAME)      :: tname = ''
    CHARACTER(LEN=STRLEN_MEDIUM)     :: units
    CHARACTER(LEN=2)                 :: nn
    INTEGER                          :: status               ! error status flag
    INTEGER                          :: idt = 0              ! tracer index
    INTEGER                          :: jg, je, je2, jt, jt2 ! loop indices
    INTEGER, DIMENSION(N_EMIS_AVAIL) :: ntrac_emis
    LOGICAL                          :: avail
    ! op_ck_20160810+
    ! op_ck_20160810  Registration is too late here!
!!$#ifdef MESSYTENDENCY
!!$    INTEGER                          :: jm, js
!!$#endif
    ! op_ck_20160810-

    CALL start_message_bi(modstr,'Check couplings', substr)

! op_cb_20180626+
    ! -------------
    ! Get additional fields for sub-grid scale orography (SSO)
    ! op_pj_20160617: these objects can currently not be moved (from g3b)
    !                 to here, because they are initialised in ioinitial.f90
    CALL get_channel_object(status,'g3b','oromea',p2=oromea)
    CALL channel_halt(substr, status)
! op_cb_20180626-

    ! Check couplings with chemistry submodel

    DO jg = 1, ngas

       IF (TRIM(gas_cpl(jg)%chemmod) /= '') THEN

          ! Build tracer name for output (log)
          IF (TRIM(gas_cpl(jg)%sname) /= '') THEN
             tname = TRIM(gas_cpl(jg)%bname)//'_'//TRIM(gas_cpl(jg)%sname)
          ELSE
             tname = TRIM(gas_cpl(jg)%bname)
          END IF

          ! Check if tracer exists and is provided by requested submodel
          CALL info_bi('Checking for gas phase tracer '//TRIM(tname)//' ...' &
               , substr)
          CALL get_tracer(status, GPTRSTR, TRIM(gas_cpl(jg)%bname) &
               , subname=TRIM(gas_cpl(jg)%sname), idx=idt, submodel=chemmod)
          CALL tracer_halt(substr, status)
          IF (TRIM(gas_cpl(jg)%chemmod) /= TRIM(chemmod))      &
               CALL error_bi('Gas phase tracer '//TRIM(tname)//&
               ' not provided by submodel '//TRIM(gas_cpl(jg)%chemmod) &
               , substr)
          CALL info_bi('... OK', substr)

          ! Set tracer index
          itrac(jg,gas) = idt

       END IF

    END DO

    ! Gas phase production of H2SO4(g) (SO2 + OH ---> H2SO4)
    IF (TRIM(gas_cpl(i_ph2so4)%bname) == '') THEN
       CALL info_bi('No gas phase production of H2SO4 specified, '//&
            'using H2SO4 tendency instead!', substr)
       itrac(i_ph2so4,gas) = itrac(i_h2so4,gas)
    ENDIF

    ! op_ck_20160810+
    ! op_ck_20160810  Registration is too late here!
!!$#ifdef MESSYTENDENCY
!!$    ! Indices of all tracers that MADE3 modifies should be set by now, so
!!$    ! register them in TENDENCY now.
!!$    DO jm = 1, dim2_cblk
!!$       DO js = 1, dim1_cblk
!!$          IF (itrac(js,jm) .NE. -1) &
!!$               CALL mtend_register(my_handle, mtend_id_tracer, idt=itrac(js,jm))
!!$       END DO
!!$    END DO
!!$#endif
    ! op_ck_20160810-

    ! Check couplings with online emissions submodel(s)

    l_emis(:) = .FALSE.
    ntrac_emis(:) = 0
    n_emis = 0

    ! Check which emission types are used and store number of tracers associated
    ! with each
    DO je = 1, N_MAX_EMIS_CPL
       avail = .FALSE.
       IF (emis_cpl(je)%type.NE.'') THEN
          DO je2 = 1, N_EMIS_AVAIL
             IF (TRIM(emis_avail(je2)).EQ.TRIM(emis_cpl(je)%type)) THEN
                avail = .TRUE.
                l_emis(je2) = .TRUE.
                IF (TRIM(emis_cpl(je)%tracer).NE.'') THEN
                   ntrac_emis(je2) = ntrac_emis(je2) + 1
                END IF
             END IF
          END DO
          IF (.NOT.avail) THEN
             CALL error_bi('Emission type '//TRIM(emis_cpl(je)%type)//&
                  ' not available!', substr)
          END IF
       END IF
    END DO

    ! Determine number of coupled emission types
    DO je = 1, N_EMIS_AVAIL
       IF (l_emis(je)) THEN
          n_emis = n_emis + 1
       END IF
    END DO

! op_ck_20151002  FIXME: Move allocation (and the above) to init_memory?
    ! Allocate the array emis and its elements
    ALLOCATE(emis(n_emis))
    je = 1
    DO je2 = 1, N_EMIS_AVAIL
       IF (l_emis(je2)) THEN
          emis(je)%type  = emis_avail(je2)
          emis(je)%ntrac = ntrac_emis(je2)
          ALLOCATE(emis(je)%idt(ntrac_emis(je2)))
          emis(je)%idt(:) = -1
          IF (ntrac_emis(je2).EQ.0) THEN
             ALLOCATE(emis(je)%flux(1))
             NULLIFY(emis(je)%flux(1)%ptr)
             ALLOCATE(emis(je)%scale(1))
          ELSE
             ALLOCATE(emis(je)%flux(ntrac_emis(je2)))
             DO jt = 1, ntrac_emis(je2)
                NULLIFY(emis(je)%flux(jt)%ptr)
             END DO
             ALLOCATE(emis(je)%scale(ntrac_emis(je2)))
          END IF
          emis(je)%scale(:) = 1._dp
          je = je + 1
       END IF
    END DO

    ! Fill the elements of array emis
    emis_loop: DO je2 = 1, n_emis

       jt = 1

       emis_cpl_loop: DO je = 1, N_MAX_EMIS_CPL

          type_exists: IF (TRIM(emis_cpl(je)%type).EQ.TRIM(emis(je2)%type)) THEN

             ! ... fluxes
             CALL info_bi('Checking for '//TRIM(emis_cpl(je)%type)//&
                  ' emissions ...', substr)
             CALL get_channel_object(status, TRIM(emis_cpl(je)%inchannel) &
                  , TRIM(emis_cpl(je)%inobject), p2=emis(je2)%flux(jt)%ptr)
             IF (status /= 0) THEN
                CALL error_bi('Could not obtain emissions flux from '//&
                     TRIM(emis_cpl(je)%inobject)//' in '//&
                     TRIM(emis_cpl(je)%inchannel)//'!', substr)
             ELSE
                CALL info_bi('... OK', substr)
             END IF

             tracer_specified: IF (emis(je2)%ntrac.GT.0) THEN

                ! ... tracer indices
                DO jt2 = 1, SIZE(ti_gp)
                   IF (TRIM(ti_gp(jt2)%tp%ident%fullname).EQ.&
                        TRIM(emis_cpl(je)%tracer)) THEN
                      emis(je2)%idt(jt) = jt2
                      EXIT
                   END IF
                END DO
                IF (emis(je2)%idt(jt).EQ.-1) THEN
                   CALL int2str(nn, je)
                   CALL error_bi('Tracer '//TRIM(emis_cpl(je)%tracer)//&
                        ' specified in emis_cpl('//nn//') not found!' &
                        , substr)
                END IF

                ! ... units
! op_mr_20160420+
! Skip units check for online dust: the code spots wrong units for number
! emissions since these are derived from mass emissions in made3.nml
!!$                IF (TRIM(emis_cpl(je)%inchannel).NE.'import_grid') THEN
                IF (TRIM(emis_cpl(je)%inchannel).NE.'import_grid' .AND. &
                   .NOT.(TRIM(emis_cpl(je)%type).EQ.'DU' .AND. &
                         TRIM(emis_cpl(je)%inchannel).EQ.'onemis')) THEN
! op_mr_20160420-

                   CALL get_attribute(status, TRIM(emis_cpl(je)%inchannel) &
                        , TRIM(emis_cpl(je)%inobject), 'units', c=units)
                   IF (status == 0) THEN
                      IF (ti_gp(emis(je2)%idt(jt))%tp%ident%quantity == &
                           NUMBERDENSITY) THEN   ! number emissions
                         IF (TRIM(units).NE.'m-2 s-1') THEN
                            CALL error_bi('Units of '//&
                                 TRIM(emis_cpl(je)%inobject)//' in '//&
                                 TRIM(emis_cpl(je)%inchannel)//&
                                 ' must be [m-2 s-1]!', substr)
                         END IF
                      ELSE   ! mass emissions
                         IF ((TRIM(units).NE.'kg m-2 s-1').AND.&
                              (TRIM(units).NE.'kg (POC)  m-2 s-1')) THEN
                            CALL error_bi('Units of '//&
                                 TRIM(emis_cpl(je)%inobject)//' in '//&
                                 TRIM(emis_cpl(je)%inchannel)//&
                                 ' must be [kg m-2 s-1]!', substr)
                         END IF
                      END IF
                   ELSE
                      CALL info_bi('Could not determine units of '//&
                           TRIM(emis_cpl(je)%inobject)//' in '//&
                           TRIM(emis_cpl(je)%inchannel)//&
                           '! Assuming [kg m-2 s-1] or [m-2 s-1].', substr)
                   END IF
                END IF
              
             END IF tracer_specified

             ! ... scaling factors
             emis(je2)%scale(jt) = emis_cpl(je)%scale
             jt = jt + 1

          END IF type_exists

       END DO emis_cpl_loop

    END DO emis_loop

    IF (p_parallel_io) THEN
       ! Write info to standard output
       WRITE(*,*) ''
       WRITE(*,*) '--------------------------------------------------------------'
       WRITE(*,*) '--------------------------------------------------------------'
       WRITE(*,*) '---    Settings for aerosol submodel MADE3 (couplings)     ---'
       WRITE(*,*) '---                                                        ---'
       WRITE(*,*) '--- Gas phase coupling:                                    ---'
       DO jg = 1, ngas
          CALL get_tracer(status, GPTRSTR, itrac(jg,gas), fullname=tname &
               , submodel=chemmod)
          CALL tracer_halt(substr, status)
          IF (jg .EQ. i_ph2so4) THEN
             WRITE(*,'(1X,A4,A28,A12,A6,A8,A4)') '--- ', &
                  'H2SO4(g) production: ', TRIM(tname), ' from ', &
                  TRIM(chemmod), ' ---'
          ELSE
             WRITE(*,'(1X,A4,A26,A2,A12,A6,A8,A4)') '--- ', &
                  TRIM(gspec(jg)%lname), ': ', TRIM(tname), ' from ', &
                  TRIM(chemmod), ' ---'
          END IF
       END DO
       WRITE(*,*) '---                                                        ---'
       IF (.NOT.ANY(l_emis)) THEN
          WRITE(*,*) '--- Online (distributed) emissions: none                   ---'
       ELSE
          WRITE(*,*) '--- Online (distributed) emissions:                        ---'
          DO je = 1, n_emis
             WRITE(*,'(1X,A8,A50,A4)') "---     ", TRIM(emis(je)%type) &
                  , ' ---'
          END DO
       END IF
       WRITE(*,*) '--------------------------------------------------------------'
       WRITE(*,*) '--------------------------------------------------------------'
       WRITE(*,*) ''
    END IF

    CALL end_message_bi(modstr,'Check couplings', substr)
    
  END SUBROUTINE made3_init_coupling
  
!-------------------------------------------------------------------------------

  !> \brief Distributes online SS and DU emissions among species and modes
  !> \details Assigns online sea spray emissions to the soluble accumulation and
  !>   coarse modes and distributes them among the Na and Cl tracers; assigns
  !>   online dust emissions to the insoluble accumulation and coarse modes and
  !>   calculates the number fluxes from the mass fluxes.

  SUBROUTINE made3_vdiff

    ! External parameters and subroutines
    USE messy_main_tracer_mem_bi, ONLY: ti_gp
    USE messy_main_grid_def_mem_bi,ONLY: kproma, jrow
#ifdef ECHAM5
    USE messy_main_data_bi,       ONLY: pxtems
#else
    USE messy_main_data_bi,       ONLY: pressi_3d
    USE messy_main_grid_def_mem_bi,ONLY:nlev
    USE messy_main_constants_mem, ONLY: g
#ifdef MESSYTENDENCY
    USE messy_main_grid_def_mem_bi, ONLY: nproma
    USE messy_main_tendency_bi,   ONLY: mtend_add_l, mtend_id_tracer
#else
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte
#endif
#endif
    USE messy_main_tracer,        ONLY: R_MOLARMASS, NUMBERDENSITY
    USE messy_main_constants_mem, ONLY: M_air, pi
    USE messy_made3,              ONLY: i_ss, i_cl, i_du, i_num &
         , acc, cor, sootj, sootc, mfCl, rho

    IMPLICIT NONE
    INTRINSIC EXP, LOG

    ! Local
#ifdef ECHAM5
    REAL(dp), POINTER :: zxtems(:,:)
#else
    REAL(dp)          :: zdp(kproma)
#ifdef MESSYTENDENCY
    REAL(dp)          :: tend(nproma,nlev)
#endif
#endif
    REAL(dp), POINTER :: flux(:,:)
    INTEGER           :: je, jt, idt
! op_cb_20180626+
    INTEGER           :: ii
! op_cb_20180626-
    ! Fraction of DU mass assigned to mode sootj
    REAL(dp), PARAMETER :: massfrac_du_soot = 0.0063_dp
    ! Fraction of DU mass assigned to mode sootc
    REAL(dp), PARAMETER :: massfrac_du_c = 0.9937_dp
    ! Conversion factor from DU mass to number concentration for mode sootj
    REAL(dp) :: DU2SOOT0
    ! Conversion factor from DU mass to number concentration for mode sootc
    REAL(dp) :: DU2CORN
    REAL(dp) :: prefac       ! prefactor for conversion of mass to number conc.
    REAL(dp) :: r            ! radius of log-normal distribution [m]
    REAL(dp) :: zsigma       ! standard deviation of log-normal distribution
! op_cb_20180626+
    REAL(dp), DIMENSION(:), POINTER :: pmea   ! Mean Orography [m]
! op_cb_20180626-

    ! Return if no emissions coupling enabled
    IF (.NOT.ANY(l_emis)) RETURN

! op_cb_20180626+
    pmea => oromea(1:kproma,jrow)
! op_cb_20180626-

#ifdef ECHAM5
    ! units of zxtems:
    ! - mass related: mol(X) mol-1(air) * kg(air) m-2 s-1
    !                 -> emission flux in mol(X) m-2 s-1, multiplied by
    !                    molecular mass of air
    ! - number related: #(X) mol-1(air) * kg(air) m-2 s-1
    !                   -> emission flux in #(X) m-2 s-1, multiplied by
    !                      molecular mass of air
    zxtems => pxtems(:,1,:,jrow)
#else
    ! pressure diff. between top and bottom of currently processed grid cells
    zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) - pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev))
#endif

    ! Assign emissions fluxes to tracers as specified via CPL, depending on the
    ! base model and the use of the TENDENCY submodel
    DO je = 1, n_emis
       DO jt = 1, emis(je)%ntrac

          flux => emis(je)%flux(jt)%ptr
          idt  =  emis(je)%idt(jt)

          IF (ti_gp(idt)%tp%ident%quantity == NUMBERDENSITY) THEN

             ! number tracer
#ifdef ECHAM5
! op_cb_20180626+
          ! Set dust emission flux to zero if orography exceeds threshold
          ! specified in made3.nml. This is to exclude high emission artefacts
          ! in e.g. Himalaya-region when using low horizontal model resolution.
          IF (TRIM(emis(je)%type).EQ.'DU') THEN
             IF (orogr_thr_4DUemis%l) THEN
                DO ii = 1, kproma
                    IF (pmea(ii).GE.orogr_thr_4DUemis%v) THEN
                       flux(ii,jrow) = 0._dp
                    END IF
                END DO
             END IF
          END IF
! op_cb_20180626-
             zxtems(1:kproma,idt) = zxtems(1:kproma,idt) &
                  + emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air * 1.0e-3_dp
#else
#ifdef MESSYTENDENCY
             tend(:,:) = 0._dp
             tend(1:kproma,nlev) = &
                  emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air * g &
                  * 1.0e-3_dp / zdp(1:kproma)
!!$          CALL mtend_add_l(my_handle, mtend_id_tracer, px=tend, idt=idt)
             CALL mtend_add_l(my_handle, idt, px=tend)
#else
             pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) = pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) &
                  + emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air * g &
                  * 1.0e-3_dp / zdp(1:kproma)
#endif
#endif

          ELSE

             ! mass tracer
#ifdef ECHAM5
! op_cb_20180626+
          ! Set dust emission flux to zero if orography exceeds threshold
          ! specified in made3.nml. This is to exclude high emission artefacts
          ! in e.g. Himalaya-region
          IF (TRIM(emis(je)%type).EQ.'DU') THEN
             IF (orogr_thr_4DUemis%l) THEN
                DO ii = 1, kproma
                    IF (pmea(ii).GE.orogr_thr_4DUemis%v) THEN
                       ! print *, 'pmea(ii): ', pmea(ii)
                       flux(ii,jrow) = 0._dp
                       ! print *, 'flux = ', flux(ii,jrow)
                    END IF
                END DO
             END IF
          END IF
! op_cb_20180626-
             zxtems(1:kproma,idt) = zxtems(1:kproma,idt) &
                  + emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air &
                  / ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS)

! op_cb_20180327+
             ! DU_ai mass emission flux in [kg/m2/s]
             IF (TRIM(emis(je)%type).EQ.'DU' .AND. &
                ti_gp(idt)%tp%ident%subname.EQ.'ai') THEN
                p2_DUemfluxai(1:kproma,jrow) = emis(je)%scale(jt) * &
                                               flux(1:kproma,jrow)
             END IF

             ! DU_ci mass emission flux in [kg/m2/s]
             IF (TRIM(emis(je)%type).EQ.'DU' .AND. &
                ti_gp(idt)%tp%ident%subname.EQ.'ci') THEN
                p2_DUemfluxci(1:kproma,jrow) = emis(je)%scale(jt) * &
                                               flux(1:kproma,jrow)
             END IF
! op_cb_20180327-

#else
#ifdef MESSYTENDENCY
             tend(:,:) = 0._dp
             tend(1:kproma,nlev) = &
                  emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air * g &
                  / (ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) * zdp(1:kproma))
!!$          CALL mtend_add_l(my_handle, mtend_id_tracer, px=tend, idt=idt)
             CALL mtend_add_l(my_handle, idt, px=tend)
#else
             pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) = pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) &
                  + emis(je)%scale(jt) * flux(1:kproma,jrow) * M_air * g &
                  / (ti_gp(idt)%tp%meta%cask_r(R_MOLARMASS) * zdp(1:kproma))
#endif
#endif

          END IF  ! number tracer

       END DO
    END DO

  END SUBROUTINE made3_vdiff

!-------------------------------------------------------------------------------

  !> \brief Entry point to MADE3 core routines

  !> \details Calculates current value of relative humidity in cloud-free area
  !>   of grid box, reads or calculates other environmental and tracer data from
  !>   various MESSy variables (and converts their units to ``MADE3 units''),
  !>   calls messy_made3::made3_main, calls the subroutine
  !>   messy_made3_si::update_xte to update MESSy tracer tendencies, and stores
  !>   updated MADE3 variables in the MADE3 channel objects.

  SUBROUTINE made3_physc

  ! External parameters and subroutines
    USE messy_main_constants_mem, ONLY: g, N_A, M_air
    USE messy_main_tools,         ONLY: spec2relhum_q, mass_density
    USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow, nlev, nproma
    USE messy_main_data_bi,       ONLY: press_3d, pressi_3d           &
#ifndef MESSYTENDENCY
         , tm1_3d, tte_3d, qm1_3d, qte_3d                             &
#endif
! op_ck_20120612+
! op_ck_20120612  acdnc, vervel, and tke are no longer required if cloud stuff
!                 is done by CLOUD
!!$                                      aclc, acdnc, &
         , aclc
!!$                                      vervel, tke
! op_ck_20120612-
    USE messy_main_timer,         ONLY: time_step_len
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,   ONLY: mtend_get_start_l, mtend_add_l &
                                      , mtend_id_t, mtend_id_q, mtend_id_tracer
    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, qxtte
#else
    USE messy_main_tracer_mem_bi, ONLY: qxtm1, qxtte
#endif
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_made3,              ONLY: dim1_cblk, dim2_cblk, nmod      &
         , sooti, sootj, sootc, gas          &
         , i_pom, i_bc, i_bctag, i_du        &
         , i_num, i_mom3, i_h2so4                        &
         , MASS2MOM3                             &
! DEBUG+
!!$         , dbg_cblk_1, dbg_cblk_2, dbg_dry_1, dbg_dry_2 &
! DEBUG-
         , made3_main
                                   
    IMPLICIT NONE
    INTRINSIC ASSOCIATED, MAX, MIN, EPSILON, TRIM, SUM

    ! Local
    CHARACTER(LEN=*), PARAMETER :: substr = 'made3_physc'
    INTEGER  :: status                ! error status flag
    ! Interaction between MESSy and MADE3
    REAL(dp) :: zpress(nproma)        ! full level pressure [Pa]
    REAL(dp) :: ztemp(nproma)         ! temperature [K]
    REAL(dp) :: zrelhum(nproma)       ! rel. humidity (0-1) [frac]
    REAL(dp) :: zspechum_m(nproma)    ! specific humidity [kg/kg],
                                    ! grid box mean
    REAL(dp) :: zqs                   ! saturation specific humidity
    REAL(dp) :: zrhoair(nproma)       ! density of air [kg/m3]
    REAL(dp) :: zclcover(nproma)      ! cloud cover [frac]
    REAL(dp) :: zso4rat(nproma)       ! H2SO4(g) production rate [ug/m3/s]
    REAL(dp) :: zsoa(nproma)          ! SOA "emissions" [ug/m3/s]
    REAL(dp) :: zrh_hist(nmod,nproma) ! Deliquescence history
    REAL(dp) :: zcblk(dim1_cblk,dim2_cblk,nproma)    ! tracer conc. array
                                                     ! for MADE3
    REAL(dp) :: zcblk_m1(dim1_cblk,dim2_cblk,nproma) ! tracer conc. array
                                                     ! before MADE3
!  REAL(dp) :: znewsulf(nproma)      ! auxiliary field for H2SO4
    REAL(dp) :: ztmst
    REAL(dp), POINTER :: soa(:,:)
    INTEGER  :: jp, jk, jm, jg, js, je
    INTEGER  :: idt
    ! Additional MADE3 output
    REAL(dp) :: zdg(nmod,nproma)      ! geom. mean diameter (wet) [m]
    REAL(dp) :: zdgdry(nmod,nproma)   ! geom. mean diameter (dry) [m]
    REAL(dp) :: zdens(nmod,nproma)    ! average particle density (wet) [kg m-3]
#ifdef MESSYTENDENCY
    REAL(dp) :: t_curr(nproma,nlev)
    REAL(dp) :: q_curr(nproma,nlev)
    REAL(dp) :: trac_curr(_RI_X_ZN_(nproma,nlev,ntrac_gp))
    REAL(dp) :: tend(_RI_X_ZN_(nproma,nlev,ntrac_gp))
    REAL(dp) :: tend_jk(kproma,ntrac_gp)
#endif

    ! --- 0. Initializations ---------------------------------------------------
    ztmst         = time_step_len
    zsoa(:)       = 0._dp

    ! --- 1. Calculate current values ------------------------------------------

#ifdef MESSYTENDENCY
    CALL mtend_get_start_l(mtend_id_t, v0=t_curr)
    CALL mtend_get_start_l(mtend_id_q, v0=q_curr)
    CALL mtend_get_start_l(mtend_id_tracer, v0t=trac_curr)
    tend(:,:,:) = 0._dp
#endif

    loop_levels: DO jk = 1, nlev

       ! --- 1.1 Meteorological data -------------------------------------------

       ! Current temperature [K]
#ifdef MESSYTENDENCY
       ztemp(1:kproma) = t_curr(1:kproma,jk)
#else
       ztemp(1:kproma) = tm1_3d(_RI_XYZ__(1:kproma,jrow,jk))+tte_3d(_RI_XYZ__(1:kproma,jrow,jk))*ztmst
#endif

       ! Current pressure [Pa]
       zpress(1:kproma) = press_3d(_RI_XYZ__(1:kproma,jrow,jk))

       ! Current specific humidity [kg kg-1], grid box mean
       ! (Based on the usage of qm1_3d in messy/smil/messy_sedi_si.f90, qm1_3d
       ! should be the mass mixing ratio of water vapor to total air, incl.
       ! water vapor).
#ifdef MESSYTENDENCY
       zspechum_m(1:kproma) = MAX(EPSILON(1.0_dp), q_curr(1:kproma,jk))
#else
       zspechum_m(1:kproma) = MAX(EPSILON(1.0_dp), &
            qm1_3d(_RI_XYZ__(1:kproma,jrow,jk)) + qte_3d(_RI_XYZ__(1:kproma,jrow,jk)) * ztmst)
#endif

       ! Current cloud cover [frac]
       zclcover(1:kproma) = MIN(aclc(_RI_XYZ__(1:kproma,jrow,jk)), 0.9999_dp)

       ! Current relative humidity [frac]

       DO jp = 1, kproma


          ! Grid box mean RH
! op_ck_20131218  FIXME: Since EQSAM relies on the assumption that water
!                        activity = RH, the proper definition to use here might
!                        be messy_main_tools::spec2relhum. This should be
!                        checked in the future.
          CALL spec2relhum_q(status, zrelhum(jp), zspechum_m(jp), ztemp(jp) &
               , zpress(jp), liq_only=.TRUE.)
          IF (status == 2) &
               CALL error_bi('Lookup overflow in spec2relhum_q call.', substr)
          IF (status /= 0) &
               CALL error_bi('Could not obtain valid RH value.', substr)

          ! beta test: calculate specific humidity in cloud free area
          !            of grid box instead of average for whole grid box
          !
          ! Q_m = specific humidity, grid box mean
          ! Q_0 = specific humidity, cloud free area of grid box
          ! Q_s = specific humidity, cloudy area of grid box
          !     = saturation specific humidity (liquid water),
          !       this is an assumption for the sake of simplicity
          !
          ! Q_m = (1-cloud_cover)*Q_0 + cloud_cover*Qs
          !
          ! ==> Q_0 = (Q_m - cloud_cover*Q_s) / (1-cloud_cover)
          !     ==> RH_0 = Q_0 / Q_s = (RH_m - cloud cover) / (1 - cloud cover)

          ! Actual relative humidity [frac], cloud free area
          zrelhum(jp) = (zrelhum(jp) - zclcover(jp)) / (1.0_dp - zclcover(jp))

          ! beta test: limit relative humidity to (1%...99%)
          zrelhum(jp) = MAX(0.01_dp,MIN(zrelhum(jp),0.99_dp))

       END DO

       ! current density of air [kg m-3]
       zrhoair(1:kproma) = mass_density(zpress(1:kproma), ztemp(1:kproma) &
            , zspechum_m(1:kproma))

       ! deliquescence history ---> hysteresis
       DO jm = 1, nmod
          zrh_hist(jm,1:kproma)  = p4_rhhist(_RI_XYZN_(1:kproma,jrow,jk,jm))
       END DO

       ! --- 1.2 Chemical data -------------------------------------------------

       soa => NULL()
       DO je = 1, n_emis
          IF (TRIM(emis(je)%type).EQ.'SOA') THEN
             soa => emis(je)%flux(1)%ptr
             EXIT
          END IF
       END DO
       IF ((jk .EQ. nlev) .AND. ASSOCIATED(soa)) THEN
          ! Add SOA "emissions" (and convert from [molecules m-2 s-1] to
          ! [ug m-3 s-1])
          ! NOTE: pressi_3d contains pressure at level interfaces, therefore has
          !       vertical dimension nlev+1
          zsoa(1:kproma) = soa(1:kproma,jrow) * emis(je)%scale(1) * g &
               / (pressi_3d(_RI_XYZ__(1:kproma,jrow,jk+1)) - pressi_3d(_RI_XYZ__(1:kproma,jrow,jk))) &
               * zrhoair(1:kproma) * 1.0e6_dp * aerspec(i_pom)%molm / N_A
       ELSE
          zsoa(1:kproma) = 0.0_dp
       END IF

       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! Production of gas-phase H2SO4: SO2 + OH ---> H2SO4(g)
       ! If no H2SO4(g) production is specified via KPP, the tendency
       ! of the H2SO4(g) tracer will be used as first guess.
       ! In this case, itracprodso4 = itrac(i_h2so4,gas).
       
       ! [mol/mol/s] ---> [ug(H2SO4)/m3/s]
       idt = itrac(i_ph2so4,gas)
! op_ck_20141017  FIXME: Once TENDENCY is able to provide ``total'' tendencies,
!                        the following statement should be replaced.
       zso4rat(1:kproma) = MAX(0.0_dp, qxtte(_RI_X_ZN_(1:kproma,jk,idt)) * &
            gspec(i_h2so4)%molm / M_air * zrhoair(1:kproma) * 1.0e9_dp)
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! Assign ECHAM tracers to array CBLK used in aero. dynamics calculations

       zcblk(:,:,1:kproma)    = 1.0e-30_dp
       zcblk_m1(:,:,1:kproma) = 1.0e-30_dp

       loop_cells: DO jp = 1, kproma

          ! --- 1.2.1 Gases ----------------------------------------------------
          DO jg = 1, ngas
             IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
             idt = itrac(jg,gas)
             IF (idt .NE. -1) THEN   ! tracer exists
#ifdef MESSYTENDENCY
                zcblk(jg,gas,jp) = trac_curr(_RI_X_ZN_(jp,jk,idt))
#else
                zcblk(jg,gas,jp) = qxtm1(_RI_X_ZN_(jp,jk,idt)) + qxtte(_RI_X_ZN_(jp,jk,idt)) * ztmst
#endif
             END IF
          END DO

          loop_modes: DO jm = 1, nmod

             ! --- 1.2.2 Particle mass -----------------------------------------

             DO js = 1, nspec
                idt = itrac(js,jm)
                IF (idt .NE. -1) THEN   ! tracer exists
#ifdef MESSYTENDENCY
                   zcblk(js,jm,jp) = trac_curr(_RI_X_ZN_(jp,jk,idt))
#else
                   zcblk(js,jm,jp) = qxtm1(_RI_X_ZN_(jp,jk,idt)) &
                                     + qxtte(_RI_X_ZN_(jp,jk,idt)) * ztmst
#endif
                END IF
             END DO

             ! --- 1.2.3 Particle numbers --------------------------------------
             idt = itrac(i_num,jm)
             IF (idt .NE. -1) THEN   ! tracer exists
#ifdef MESSYTENDENCY
                zcblk(i_num,jm,jp) = trac_curr(_RI_X_ZN_(jp,jk,idt))
#else
                zcblk(i_num,jm,jp) = qxtm1(_RI_X_ZN_(jp,jk,idt)) &
                                     + qxtte(_RI_X_ZN_(jp,jk,idt)) * ztmst
#endif
             END IF

          END DO loop_modes

       END DO loop_cells

       ! Copy array CBLK for calculating the updated ECHAM tracer tendencies
       ! after aerosol dynamics calculations
       zcblk_m1(:,:,1:kproma) = zcblk(:,:,1:kproma)

! op_va_20090225+
       ! Store externally mixed BC and DU concentrations
       IF (l_bctime) THEN
          p3_BCextIN(1:kproma,jk,jrow) = MAX(0._dp, &
               SUM(zcblk(i_bc,(/sooti,sootj,sootc/),1:kproma), dim=1) &
               + SUM(zcblk(i_bctag,(/sooti,sootj,sootc/),1:kproma), dim=1))
          p3_DUextIN(1:kproma,jk,jrow) = MAX(0._dp, &
               SUM(zcblk(i_du,(/sooti,sootj,sootc/),1:kproma), dim=1))
       END IF
! op_va_20090225-

       ! --- 2. Convert from "ECHAM units" to "MADE3 units" and guarantee
       !        valid concentrations -------------------------------------------
       CALL echam_to_made3(nproma, kproma, zrhoair, zcblk)

       ! --- 3. Calculate aerosol dynamics ---> call MADE3 ---------------------

! DEBUG+
!!$       DO jm = 1, 2
!!$          DO js = 1, 2
!!$             dbg_cblk_1(js,jm)%ptr(1:nproma) => p3arr_dbg_cblk_1(js,jm)%ptr(_RI_XYZ__(1:nproma,jrow,jk))
!!$             dbg_cblk_2(js,jm)%ptr(1:nproma) => p3arr_dbg_cblk_2(js,jm)%ptr(_RI_XYZ__(1:nproma,jrow,jk))
!!$          END DO
!!$       END DO
! DEBUG-
       CALL made3_main(status, nproma, kproma, zpress, ztemp, zrelhum, ztmst &
            , zso4rat, zsoa, zclcover, zcblk, zrh_hist, zdg, zdgdry, zdens   &
! op_va_20090225+
!!$            , BCextin(_RI_XYZ__(:,jrow,jk)), BCextsink(_RI_XYZ__(:,jrow,jk))                       &
! op_va_20090225-
            )
       IF (status /= 0) THEN
          CALL error_bi('Invalid number of modes in ' &
               // 'messy_made3::renaming_fractions', substr)
       END IF
! DEBUG+
!!$       DO jm = 1, 2
!!$          DO js = 1, 2
!!$             IF (ASSOCIATED(dbg_cblk_1(js,jm)%ptr)) NULLIFY(dbg_cblk_1(js,jm)%ptr)
!!$             IF (ASSOCIATED(dbg_cblk_2(js,jm)%ptr)) NULLIFY(dbg_cblk_2(js,jm)%ptr)
!!$          END DO
!!$       END DO
!!$       p3_dbg_dry_1(_RI_XYZ__(1:nproma,jrow,jk)) = dbg_dry_1(:)
!!$       p3_dbg_dry_2(_RI_XYZ__(1:nproma,jrow,jk)) = dbg_dry_2(:)
! DEBUG-

       ! --- 4. Convert from "MADE3 units" to "ECHAM units" --------------------
       CALL made3_to_echam(nproma, kproma, zrhoair, zcblk)

! op_va_20090225+
       IF (l_bctime) THEN
          p3_BCextSINK(1:kproma,jk,jrow) = MAX(0._dp,                   &
               p3_BCextIN(1:kproma,jk,jrow)                             &
               - SUM(zcblk(i_bc,(/sooti,sootj,sootc/),1:kproma), dim=1) &
               - SUM(zcblk(i_bctag,(/sooti,sootj,sootc/),1:kproma), dim=1))
          p3_DUextSINK(1:kproma,jk,jrow) = MAX(0._dp, &
               p3_DUextIN(_RI_XYZ__(1:kproma,jrow,jk)) &
               - SUM(zcblk(i_du,(/sooti,sootj,sootc/),1:kproma), dim=1))
       ENDIF
! op_va_20090225+
     
       ! --- 5. Update ECHAM tracer tendencies for current level ---------------
#ifdef MESSYTENDENCY
       CALL update_xte(nproma, kproma, zcblk, zcblk_m1, jk, ztmst, tend_jk)
       tend(_RI_XYZ__(1:kproma,:,jk)) = tend_jk(:,:)
#else
       CALL update_xte(nproma, kproma, zcblk, zcblk_m1, jk, ztmst)
#endif

       ! --- 6. Store MADE3 quantities in channel objects ----------------------
       DO jp = 1, kproma
          DO jm = 1, nmod
             p4_wetrad(_RI_XYZN_(jp,jrow,jk,jm))  = zdg(jm,jp)    / 2.0_dp
             p4_dryrad(_RI_XYZN_(jp,jrow,jk,jm))  = zdgdry(jm,jp) / 2.0_dp
             p4_densaer(_RI_XYZN_(jp,jrow,jk,jm)) = zdens(jm,jp)
             p4_rhhist(_RI_XYZN_(jp,jrow,jk,jm))  = zrh_hist(jm,jp)
          END DO
       END DO

    END DO loop_levels

#ifdef MESSYTENDENCY
    DO jm = 1, dim2_cblk
       DO js = 1, dim1_cblk
          IF (itrac(js,jm) .NE. -1) THEN
             idt = itrac(js,jm)
             CALL mtend_add_l(my_handle, idt, px=tend(_RI_X_ZN_(:,:,idt)) )
          END IF
       END DO
    END DO
#endif

  END SUBROUTINE made3_physc

!-------------------------------------------------------------------------------

  !> \brief Deallocates memory allocated in messy_made3_si::made3_init_memory
  !> \details Deallocates 5D array pointers that were used for MADE3-internal
  !>   memory management.

  SUBROUTINE made3_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(p5_wetrad))   DEALLOCATE(p5_wetrad)
    IF (ASSOCIATED(p5_dryrad))   DEALLOCATE(p5_dryrad)
    IF (ASSOCIATED(p5_densaer))  DEALLOCATE(p5_densaer)
    IF (ASSOCIATED(p5_rhhist))   DEALLOCATE(p5_rhhist)

  END SUBROUTINE made3_free_memory

!-------------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! PRIVATE ROUTINES
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !> \brief Reads MADE3 CPL namelist
  !> \details Opens, reads, checks, and closes the MADE3 coupling namelist.

  SUBROUTINE made3_read_nml_cpl(status, iou)

    ! External subroutines and parameters
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status flag
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! Local
    CHARACTER(LEN=*), PARAMETER :: substr = 'made3_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    ! Initialize
    status = 1

    ! Open file and write info message (start of submodel initialization) to
    ! standard output
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT. lex) RETURN    ! made3.nml does not exist

    ! Read namelist
    READ(iou, NML=CPL, IOSTAT=fstat)
    ! Write info message to standard output (error/OK)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! Close file and write info message (end of submodel initialization) to
    ! standard output
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE made3_read_nml_cpl

!-------------------------------------------------------------------------------

  !> \brief Converts ``ECHAM units'' to ``MADE3 units''
  !> \details Converts from [mol mol<sup>-1</sup>] and [# mol<sup>-1</sup>] to
  !>   [ug m<sup>-3</sup>] and [# m<sup>-3</sup>], respectively.

  SUBROUTINE echam_to_made3(blksize, numcells, rhoair, cblk)

    ! External parameters
    USE messy_main_constants_mem, ONLY: MWair => M_air
    USE messy_made3,              ONLY : dim1_cblk, dim2_cblk &
         , nspec, ngas, gas, i_h2o, i_num, i_tot              &
         , dgini, nummin, massmin, es36, MASS2MOM3

    IMPLICIT NONE
    INTRINSIC TRIM, MAX, SUM

    ! I/O
    !> Size of input arrays
    INTEGER,  INTENT(in)    :: blksize
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: numcells
    !> Air density [kg m<sup>-3</sup>]
    REAL(dp), INTENT(in)    :: rhoair(blksize)
    !> Tracer array (input in ``ECHAM units'', output in ``MADE3 units'')
    REAL(dp), INTENT(inout) :: cblk(dim1_cblk,dim2_cblk,blksize)
    ! Local
    INTEGER  :: js,jg,jp,jm                     ! loop indices
    REAL(dp) :: mom3                            ! total modal 3rd moment
    REAL(dp) :: conv_a(dim1_cblk), conv_g(ngas) ! conversion factors
    LOGICAL  :: firstime                        ! switch
    DATA firstime / .TRUE. /
    SAVE firstime, conv_a, conv_g


    IF (firstime) THEN
       DO js = 1, nspec
          conv_a(js) = 1.0e9_dp * aerspec(js)%molm / MWair
       END DO
       conv_a(i_num) = 1.0e3_dp / MWair
       DO jg = 1, ngas
          IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
          conv_g(jg) = 1.0e9_dp * gspec(jg)%molm / MWair
       END DO
       firstime = .FALSE.
    END IF

    ! Convert and guarantee minimum concentrations, discard negative values
    loop_cells: DO jp = 1, numcells

       loop_modes: DO jm = 1, nmod

          ! aerosol mass concentrations: [mol mol-1] --> [ug m-3]
          DO js = 1, nspec
             cblk(js,jm,jp) = MAX(1.0e-30_dp, &
                  cblk(js,jm,jp) * conv_a(js) * rhoair(jp))
          END DO

          ! aerosol number concentrations: [# mol-1] --> [# m-3]
          cblk(i_num,jm,jp) = cblk(i_num,jm,jp) * conv_a(i_num) * rhoair(jp)

          ! Ensure minimum values for aerosol number and/or mass concentrations
          IF (cblk(i_num,jm,jp) < nummin(jm)) THEN
             IF (SUM(cblk(1:nspec,jm,jp))-cblk(i_h2o,jm,jp) &
                  < massmin(jm,i_tot)) THEN
                cblk(i_num,jm,jp) = nummin(jm)
                DO js = 1, nspec
                   cblk(js,jm,jp) = massmin(jm,js)
                END DO
             ELSE
                mom3 = 0._dp
                DO js = 1, nspec
                   IF (TRIM(aerspec(js)%bname) == 'H2O') CYCLE
                   mom3 = mom3 + MASS2MOM3(js) * cblk(js,jm,jp)
                END DO
                cblk(i_num,jm,jp) = mom3 &
                     / (dgini(jm) * dgini(jm) * dgini(jm) * es36(jm))
             END IF
          END IF

       END DO loop_modes

       ! gas mass concentrations: [mol mol-1] --> [ug m-3]
       DO jg = 1, ngas
          IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
          cblk(jg,gas,jp) = MAX(1.0e-30_dp, &
               cblk(jg,gas,jp) * conv_g(jg) * rhoair(jp))
       END DO

    END DO loop_cells

  END SUBROUTINE echam_to_made3

!-------------------------------------------------------------------------------

  !> \brief Converts ``MADE3 units'' to ``ECHAM units''
  !> \details Converts from [ug m<sup>-3</sup>] and [# m<sup>-3</sup>] to
  !>   [mol mol<sup>-1</sup>] and [# mol<sup>-1</sup>], respectively.

  SUBROUTINE made3_to_echam(blksize, numcells, rhoair, cblk)

    ! External parameters
    USE messy_main_constants_mem, ONLY: MWair => M_air
    USE messy_made3,              ONLY : dim1_cblk, dim2_cblk                  &
         , nspec, ngas, gas, i_num

    IMPLICIT NONE
    INTRINSIC TRIM, MAX

    ! I/O
    !> Size of input arrays
    INTEGER,  INTENT(in)    :: blksize
    !> Number of grid cells in input arrays
    INTEGER,  INTENT(in)    :: numcells
    !> Air density [kg m<sup>-3</sup>]
    REAL(dp), INTENT(in)    :: rhoair(blksize)
    !> Tracer array (input in ``MADE3 units'', output in ``ECHAM units'')
    REAL(dp), INTENT(inout) :: cblk(dim1_cblk,dim2_cblk,blksize)
    ! Local
    INTEGER :: js,jg,jp,jm                      ! loop indices
    REAL(dp) :: conv_a(dim1_cblk), conv_g(ngas) ! conversion factors
    LOGICAL  :: firstime
    DATA firstime / .TRUE. /
    SAVE firstime, conv_a, conv_g


    IF (firstime) THEN
       DO js = 1, nspec
          conv_a(js) = 1.0e-9_dp * MWair / aerspec(js)%molm
       END DO
       conv_a(i_num) = 1.0e-3_dp * MWair
       DO jg = 1, ngas
          IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
          conv_g(jg) = 1.0e-9_dp * MWair / gspec(jg)%molm
       END DO
       firstime = .FALSE.
    END IF

    DO jp = 1, numcells

       DO jm = 1, nmod
          ! aerosol mass concentrations: [mol mol-1] --> [ug m-3]
          DO js = 1, nspec
             cblk(js,jm,jp) = MAX(1.0e-30_dp, &
                  cblk(js,jm,jp) * conv_a(js) / rhoair(jp))
          END DO
          ! aerosol number concentrations: [# mol-1] --> [# m-3]
          cblk(i_num,jm,jp) = MAX(1.0e-30_dp, &
               cblk(i_num,jm,jp) * conv_a(i_num) / rhoair(jp))
       END DO

       ! gas mass concentrations: [mol mol-1] --> [ug m-3]
       DO jg = 1, ngas
          IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
          cblk(jg,gas,jp) = MAX(1.0e-30_dp, &
               cblk(jg,gas,jp) * conv_g(jg) / rhoair(jp))
       END DO

    END DO

  END SUBROUTINE made3_to_echam

!-------------------------------------------------------------------------------

  !> \brief Updates MADE3 tracer tendencies
  !> \details Updates tracer tendencies (XTE) for all MADE3 tracers for a given
  !>   model level (jk).

#ifdef MESSYTENDENCY
  SUBROUTINE update_xte(blksize, numcells, cblk, cblk_m1, jk, ptmst, tend)
#else
  SUBROUTINE update_xte(blksize, numcells, cblk, cblk_m1, jk, ptmst)
#endif

    ! External parameters
#ifdef MESSYTENDENCY
    USE messy_main_grid_def_mem_bi,ONLY: nproma
    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp
#else
    USE messy_main_tracer_mem_bi, ONLY: qxtte
#endif
    USE messy_made3,              ONLY: dim1_cblk, dim2_cblk, nmod, nspec &
         , gas, i_num

    IMPLICIT NONE
    INTRINSIC TRIM

    ! I/O
    !> Size of input arrays
    INTEGER, INTENT(in) :: blksize
    !> Number of grid cells in input arrays
    INTEGER, INTENT(in) :: numcells
    !> Tracer array (NEW values in ``ECHAM units'')
    REAL(dp), INTENT(in) :: cblk(dim1_cblk,dim2_cblk,blksize)
    !> Tracer array (OLD values in ``ECHAM units'')
    REAL(dp), INTENT(in) :: cblk_m1(dim1_cblk,dim2_cblk,blksize)
    !> Model level index
    INTEGER, INTENT(in) :: jk
    !> Time step [s]
    REAL(dp), INTENT(in) :: ptmst
#ifdef MESSYTENDENCY
    !> Tracer tendencies in given level
    REAL(dp), INTENT(out) :: tend(numcells,ntrac_gp)
#endif
    ! Local
    INTEGER  :: jp,jm,js,jg     ! loop indices
    INTEGER  :: idt             ! tracer index
    REAL(dp) :: pqtmst          ! 1 / ptmst


#ifdef MESSYTENDENCY
    tend(:,:) = 0._dp
#endif
    pqtmst = 1.0_dp / ptmst

    loop_cells: DO jp = 1, numcells

       loop_modes: DO jm = 1, nmod

          ! Aerosol components
          loop_species: DO js = 1, nspec

             idt = itrac(js,jm)

             IF (idt .NE. -1) THEN
#ifdef MESSYTENDENCY
                tend(jp,idt) = (cblk(js,jm,jp) - cblk_m1(js,jm,jp)) * pqtmst
#else
                qxtte(_RI_X_ZN_(jp,jk,idt)) = qxtte(_RI_X_ZN_(jp,jk,idt)) &
                     + (cblk(js,jm,jp) - cblk_m1(js,jm,jp)) * pqtmst
#endif
             END IF

          END DO loop_species

          ! Aerosol numbers
          idt = itrac(i_num,jm)
          IF (idt .NE. -1) THEN   ! tracer exists
#ifdef MESSYTENDENCY
             tend(jp,idt) = (cblk(i_num,jm,jp) - cblk_m1(i_num,jm,jp)) * pqtmst
#else
             qxtte(_RI_X_ZN_(jp,jk,idt)) = qxtte(_RI_X_ZN_(jp,jk,idt)) &
                  + (cblk(i_num,jm,jp) - cblk_m1(i_num,jm,jp)) * pqtmst
#endif
          END IF

       END DO loop_modes

       ! Gas phase species
       DO jg = 1, ngas
          IF (TRIM(gspec(jg)%bname) == 'SOApre') CYCLE
          idt = itrac(jg,gas)
          IF (idt .NE. -1) THEN   ! tracer exists
#ifdef MESSYTENDENCY
             tend(jp,idt) = (cblk(jg,gas,jp) - cblk_m1(jg,gas,jp)) * pqtmst
#else
             qxtte(_RI_X_ZN_(jp,jk,idt)) = qxtte(_RI_X_ZN_(jp,jk,idt)) &
                  + (cblk(jg,gas,jp) - cblk_m1(jg,gas,jp)) * pqtmst
#endif
          END IF
       END DO

    END DO loop_cells

  END SUBROUTINE update_xte

END MODULE messy_made3_si
