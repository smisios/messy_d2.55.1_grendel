{-----------------------------------------------------------------------------}
{-------------------------------- ice phase (i) ------------------------------}
{-----------------------------------------------------------------------------}

{------------------------------- neutral species -----------------------------}

{------------------------------------- O -------------------------------------}

 O2_i         = IGNORE; {@O_2\ice}          {oxygen}
 O3_i         = IGNORE; {@O_3\ice}          {ozone}

{------------------------------------- H -------------------------------------}

 OH_i         = IGNORE; {@OH\ice}           {hydroxyl radical}
 HO2_i        = IGNORE; {@HO_2\ice}         {perhydroxyl radical}
 H2O_i        = IGNORE; {@H_2O\ice}         {water}
 H2O2_i       = IGNORE; {@H_2O_2\ice}       {hydrogen peroxide}

{------------------------------------- N -------------------------------------}

 NH3_i        = IGNORE; {@NH_3\ice}         {ammonia}
 NO_i         = IGNORE; {@NO\ice}           {nitric oxide}
 NO2_i        = IGNORE; {@NO_2\ice}         {nitrogen dioxide}
 NO3_i        = IGNORE; {@NO_3\ice}         {nitrogen trioxide}
 HONO_i       = IGNORE; {@HONO\ice}         {nitrous acid}
 HNO3_i       = IGNORE; {@HNO_3\ice}        {nitric acid}
 HNO4_i       = IGNORE; {@HNO_4\ice}        {pernitric acid}
 N2O5_i       = IGNORE; {@N_2O_5\ice}       {dinitrogen pentoxide}

{------------------------------------- C -------------------------------------}

{1C}
 CH3OH_i      = IGNORE; {@CH_3OH\ice}       {methanol}
 HCOOH_i      = IGNORE; {@HCOOH\ice}        {formic acid}
 HCHO_i       = IGNORE; {@HCHO\ice}         {methanal (formaldehyde)}
 CH3O2_i      = IGNORE; {@CH_3OO\ice}       {methylperoxy radical}
 CH3OOH_i     = IGNORE; {@CH_3OOH\ice}      {}
 CO2_i        = IGNORE; {@CO_2\ice}         {carbon dioxide}

{2C}
 CH3CO2H_i    = IGNORE; {@CH_3COOH\ice}     {acetic acid}
 PAN_i        = IGNORE; {@PAN\ice}          {peroxyacetylnitrate}
 C2H5O2_i     = IGNORE; {@C_2H_5O_2\ice}    {ethylperoxy radical}
 CH3CHO_i     = IGNORE; {@CH_3CHO\ice}      {acetaldehyde}

{3C}
 CH3COCH3_i   = IGNORE; {@CH_3COCH_3\ice}   {acetone}

{------------------------------------- Cl ------------------------------------}

 Cl_i         = IGNORE; {@Cl\ice}           {chlorine atom}
 Cl2_i        = IGNORE; {@Cl_2\ice}         {molecular chlorine}
 HCl_i        = IGNORE; {@HCl\ice}          {hydrogen chloride}
 HOCl_i       = IGNORE; {@HOCl\ice}         {hypochlorous acid}

{------------------------------------- Br ------------------------------------}

 Br_i         = IGNORE; {@Br\ice}           {bromine atom}
 Br2_i        = IGNORE; {@Br_2\ice}         {molecular bromine}
 HBr_i        = IGNORE; {@HBr\ice}          {hydrogen bromide}
 HOBr_i       = IGNORE; {@HOBr\ice}         {hypobromous acid}
 BrCl_i       = IGNORE; {@BrCl\ice}         {bromine chloride}

{------------------------------------- I -------------------------------------}

 I2_i         = IGNORE; {@I_2\ice}          {molecular iodine}
 IO_i         = IGNORE; {@IO\ice}           {iodine monoxide radical}
 HI_i         = IGNORE; {@HI\ice}           {hydrogen iodide}
 HOI_i        = IGNORE; {@HOI\ice}          {hypoiodous acid}
 ICl_i        = IGNORE; {@ICl\ice}          {iodine chloride}
 IBr_i        = IGNORE; {@IBr\ice}          {iodine bromide}
 HIO3_i       = IGNORE; {@HIO_3\ice}        {iodic acid}

{------------------------------------- S -------------------------------------}

 SO2_i        = IGNORE; {@SO_2\ice}         {sulfur dioxide}
 H2SO4_i      = IGNORE; {@H_2SO_4\ice}      {sulfuric acid}
 DMSO_i       = IGNORE; {@DMSO\ice}         {dimethyl sulfoxide: CH3SOCH3}

{----------------------------------- ions ------------------------------------}

{------------------------------------- O -------------------------------------}

 O2m_i        = IGNORE; {@O_2^-\ice}        {}
 OHm_i        = IGNORE; {@OH^-\ice}         {}

{------------------------------------- H -------------------------------------}

 Hp_i         = IGNORE; {@H^+\ice}          {}

{------------------------------------- N -------------------------------------}

 NH4p_i       = IGNORE; {@NH_4^+\ice}       {ammonium}
 NO2m_i       = IGNORE; {@NO_2^-\ice}       {nitrite}
 NO3m_i       = IGNORE; {@NO_3^-\ice}       {nitrate}
 NO4m_i       = IGNORE; {@NO_4^-\ice}       {peroxy nitrate}

{------------------------------------- C -------------------------------------}

{1C}
 CO3m_i       = IGNORE; {@CO_3^-\ice}       {}
 HCOOm_i      = IGNORE; {@HCOO^-\ice}       {formate}
 HCO3m_i      = IGNORE; {@HCO_3^-\ice}      {hydrogen carbonate}

{2C}
 CH3COOm_i    = IGNORE; {@CH_3COO^-\ice}    {acetate}

{------------------------------------- Cl ------------------------------------}

 Clm_i        = IGNORE; {@Cl^-\ice}         {chloride}
 Cl2m_i       = IGNORE; {@Cl_2^-\ice}       {}
 ClOm_i       = IGNORE; {@ClO^-\ice}        {}
 ClOHm_i      = IGNORE; {@ClOH^-\ice}       {}

{------------------------------------- Br ------------------------------------}

 Brm_i        = IGNORE; {@Br^-\ice}         {bromide}
 Br2m_i       = IGNORE; {@Br_2^-\ice}       {}
 BrOm_i       = IGNORE; {@BrO^-\ice}        {}
 BrOHm_i      = IGNORE; {@BrOH^-\ice}       {}
 BrCl2m_i     = IGNORE; {@BrCl_2^-\ice}     {}
 Br2Clm_i     = IGNORE; {@Br_2Cl^-\ice}     {}

{------------------------------------- I -------------------------------------}

 Im_i         = IGNORE; {@I^-\ice}          {iodide}
 IO2m_i       = IGNORE; {@IO_2^-\ice}       {}
 IO3m_i       = IGNORE; {@IO_3^-\ice}       {iodate}
 ICl2m_i      = IGNORE; {@ICl_2^-\ice}      {}
 IClBrm_i     = IGNORE; {@IClBr^-\ice}      {}
 IBr2m_i      = IGNORE; {@IBr_2^-\ice}      {}

{------------------------------------- S -------------------------------------}

 SO3m_i       = IGNORE; {@SO_3^-\ice}       {}
 SO3mm_i      = IGNORE; {@SO_3^<2->\ice}    {sulfite}
 SO4m_i       = IGNORE; {@SO_4^-\ice}       {}
 SO4mm_i      = IGNORE; {@SO_4^<2->\ice}    {sulfate}
 SO5m_i       = IGNORE; {@SO_5^-\ice}       {}
 HSO3m_i      = IGNORE; {@HSO_3^-\ice}      {hydrogen sulfite}
 HSO4m_i      = IGNORE; {@HSO_4^-\ice}      {hydrogen sulfate}
 HSO5m_i      = IGNORE; {@HSO_5^-\ice}      {}
 CH3SO3m_i    = IGNORE; {@CH_3SO_3^-\ice}   {}
 CH2OHSO3m_i  = IGNORE; {@CH_2OHSO_3^-\ice} {}

{-----------------------------------------------------------------------------}
{------------------------------------ dummies --------------------------------}
{-----------------------------------------------------------------------------}

 D1O_i        = IGNORE; {@D_1O\ice}         {}
 D2O_i        = IGNORE; {@D_2O\ice}         {}
 DAHp_i       = IGNORE; {@DAH^+\ice}        {}
 DA_i         = IGNORE; {@DA\ice}           {}
 DAm_i        = IGNORE; {@DA^-\ice}         {}
 DGtAi_i      = IGNORE; {@DGtAi\ice}        {}
 DGtAs_i      = IGNORE; {@DGtAs\ice}        {}
 PROD1_i      = IGNORE; {@PROD1\ice}        {}
 PROD2_i      = IGNORE; {@PROD2\ice}        {}
 Nap_i        = IGNORE; {@Na^+\ice}         {dummy cation}

{ mz_ak_20060516+} 
{ defined in gas.spc to avoid doubling}
{BrSScap       = IGNORE; }{@IGNORE} {yield of Br to SS}
{ mz_ak_20060516-}
