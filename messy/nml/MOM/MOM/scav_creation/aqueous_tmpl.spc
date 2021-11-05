{-----------------------------------------------------------------------------}
{------------------------------ aerosol mode: ## -----------------------------}
{-----------------------------------------------------------------------------}

{------------------------------- neutral species -----------------------------}

{------------------------------------- O -------------------------------------}

 O2_##         = IGNORE; {@O_2\aq}          {oxygen}
 O3_##         = IGNORE; {@O_3\aq}          {ozone}

{------------------------------------- H -------------------------------------}

 OH_##         = IGNORE; {@OH\aq}           {hydroxyl radical}
 HO2_##        = IGNORE; {@HO_2\aq}         {perhydroxyl radical}
 H2O_##        = IGNORE; {@H_2O\aq}         {water}
 H2O2_##       = IGNORE; {@H_2O_2\aq}       {hydrogen peroxide}

{------------------------------------- N -------------------------------------}

 NH3_##        = IGNORE; {@NH_3\aq}         {ammonia}
 NO_##         = IGNORE; {@NO\aq}           {nitric oxide}
 NO2_##        = IGNORE; {@NO_2\aq}         {nitrogen dioxide}
 NO3_##        = IGNORE; {@NO_3\aq}         {nitrogen trioxide}
 HONO_##       = IGNORE; {@HONO\aq}         {nitrous acid}
 HNO3_##       = IGNORE; {@HNO_3\aq}        {nitric acid}
 HNO4_##       = IGNORE; {@HNO_4\aq}        {pernitric acid}
 N2O5_##       = IGNORE; {@N_2O_5\aq}       {dinitrogen pentoxide}

{------------------------------------- C -------------------------------------}

{1C}
 CH3OH_##      = IGNORE; {@CH_3OH\aq}       {methanol}
 HCOOH_##      = IGNORE; {@HCOOH\aq}        {formic acid}
 HCHO_##       = IGNORE; {@HCHO\aq}         {methanal (formaldehyde)}
 CH3O2_##      = IGNORE; {@CH_3OO\aq}       {methylperoxy radical}
 CH3OOH_##     = IGNORE; {@CH_3OOH\aq}      {}
 CO2_##        = IGNORE; {@CO_2\aq}         {carbon dioxide}

{2C}
 CH3CO2H_##    = IGNORE; {@CH_3COOH\aq}     {acetic acid}
 PAN_##        = IGNORE; {@PAN\aq}          {peroxyacetylnitrate}
 C2H5O2_##     = IGNORE; {@C_2H_5O_2\aq}    {ethylperoxy radical}
 CH3CHO_##     = IGNORE; {@CH_3CHO\aq}      {acetaldehyde}
 HOCH2CO2H_##  = IGNORE; {@OHCH_2CO2H\aq}   {glyoxalic acid}
 GLYOX_##      = IGNORE; {@CHOCHO\aq}       {glyoxal}
 OXL_##        = IGNORE; {@HO_2CCO_2H\aq}   {oxalic acid}
 HOCH2CHO_##   = IGNORE; {@HOCH_2CHOaq}     {glycolaldehyde}

{3C}
 CH3COCH3_##   = IGNORE; {@CH_3COCH_3\aq}   {acetone}
 MGLYOX_##     = IGNORE; {@CH_3COCHO\aq}    {methylglyoxal}
 CH3COCO2H_##  = IGNORE; {@CH_3COCO_2H\aq}  {pyruvic acid}
 OLIG_##       = IGNORE; {@oligomeric SOA}

{------------------------------------- Cl ------------------------------------}

 Cl_##         = IGNORE; {@Cl\aq}           {chlorine atom}
 Cl2_##        = IGNORE; {@Cl_2\aq}         {molecular chlorine}
 HCl_##        = IGNORE; {@HCl\aq}          {hydrogen chloride}
 HOCl_##       = IGNORE; {@HOCl\aq}         {hypochlorous acid}

{------------------------------------- Br ------------------------------------}

 Br_##         = IGNORE; {@Br\aq}           {bromine atom}
 Br2_##        = IGNORE; {@Br_2\aq}         {molecular bromine}
 HBr_##        = IGNORE; {@HBr\aq}          {hydrogen bromide}
 HOBr_##       = IGNORE; {@HOBr\aq}         {hypobromous acid}
 BrCl_##       = IGNORE; {@BrCl\aq}         {bromine chloride}

{------------------------------------- I -------------------------------------}

 I2_##         = IGNORE; {@I_2\aq}          {molecular iodine}
 IO_##         = IGNORE; {@IO\aq}           {iodine monoxide radical}
 HI_##         = IGNORE; {@HI\aq}           {hydrogen iodide}
 HOI_##        = IGNORE; {@HOI\aq}          {hypoiodous acid}
 ICl_##        = IGNORE; {@ICl\aq}          {iodine chloride}
 IBr_##        = IGNORE; {@IBr\aq}          {iodine bromide}
 HIO3_##       = IGNORE; {@HIO_3\aq}        {iodic acid}

{------------------------------------- S -------------------------------------}

 SO2_##        = IGNORE; {@SO_2\aq}         {sulfur dioxide}
 H2SO4_##      = IGNORE; {@H_2SO_4\aq}      {sulfuric acid}
 DMSO_##       = IGNORE; {@DMSO\aq}         {dimethyl sulfoxide: CH3SOCH3}

{------------------------------------- Hg ------------------------------------}
 RGM_##        = IGNORE; {@Hg\aq}           {from reactive gaseous Hg}
 Hg_##         = IGNORE; {@Hg\aq}           {from non-reactive gaseous Hg}   

{----------------------------------- ions ------------------------------------}

{------------------------------------- O -------------------------------------}

 O2m_##        = IGNORE; {@O_2^-\aq}        {}
 OHm_##        = IGNORE; {@OH^-\aq}         {}

{------------------------------------- H -------------------------------------}

 Hp_##         = IGNORE; {@H^+\aq}          {}

{------------------------------------- N -------------------------------------}

 NH4p_##       = IGNORE; {@NH_4^+\aq}       {ammonium}
 NO2m_##       = IGNORE; {@NO_2^-\aq}       {nitrite}
 NO3m_##       = IGNORE; {@NO_3^-\aq}       {nitrate}
 NO4m_##       = IGNORE; {@NO_4^-\aq}       {peroxy nitrate}

{------------------------------------- C -------------------------------------}

{1C}
 CO3m_##       = IGNORE; {@CO_3^-\aq}       {}
 HCOOm_##      = IGNORE; {@HCOO^-\aq}       {formate}
 HCO3m_##      = IGNORE; {@HCO_3^-\aq}      {hydrogen carbonate}

{2C}
 CH3COOm_##    = IGNORE; {@CH_3COO^-\aq}    {acetate}
 HOCH2CO2m_##  = IGNORE; {@OHCH_2CO2^-\aq}  {glyoxate}
 OXLm_##       = IGNORE; {@HO_2CCO_2^-\aq}  {hydro-oxalate}
 OXLmm_##      = IGNORE; {@O_2CCO_2^<2->\aq}{oxalate}
 CH3COCO2m_##  = IGNORE; {@CH_3COCO_2^-\aq} {pyruvate}

{------------------------------------- Cl ------------------------------------}

 Clm_##        = IGNORE; {@Cl^-\aq}         {chloride}
 Cl2m_##       = IGNORE; {@Cl_2^-\aq}       {}
 ClOm_##       = IGNORE; {@ClO^-\aq}        {}
 ClOHm_##      = IGNORE; {@ClOH^-\aq}       {}

{------------------------------------- Br ------------------------------------}

 Brm_##        = IGNORE; {@Br^-\aq}         {bromide}
 Br2m_##       = IGNORE; {@Br_2^-\aq}       {}
 BrOm_##       = IGNORE; {@BrO^-\aq}        {}
 BrOHm_##      = IGNORE; {@BrOH^-\aq}       {}
 BrCl2m_##     = IGNORE; {@BrCl_2^-\aq}     {}
 Br2Clm_##     = IGNORE; {@Br_2Cl^-\aq}     {}

{------------------------------------- I -------------------------------------}

 Im_##         = IGNORE; {@I^-\aq}          {iodide}
 IO2m_##       = IGNORE; {@IO_2^-\aq}       {}
 IO3m_##       = IGNORE; {@IO_3^-\aq}       {iodate}
 ICl2m_##      = IGNORE; {@ICl_2^-\aq}      {}
 IClBrm_##     = IGNORE; {@IClBr^-\aq}      {}
 IBr2m_##      = IGNORE; {@IBr_2^-\aq}      {}

{------------------------------------- S -------------------------------------}

 SO3m_##       = IGNORE; {@SO_3^-\aq}       {}
 SO3mm_##      = IGNORE; {@SO_3^<2->\aq}    {sulfite}
 SO4m_##       = IGNORE; {@SO_4^-\aq}       {}
 SO4mm_##      = IGNORE; {@SO_4^<2->\aq}    {sulfate}
 SO5m_##       = IGNORE; {@SO_5^-\aq}       {}
 HSO3m_##      = IGNORE; {@HSO_3^-\aq}      {hydrogen sulfite}
 HSO4m_##      = IGNORE; {@HSO_4^-\aq}      {hydrogen sulfate}
 HSO5m_##      = IGNORE; {@HSO_5^-\aq}      {}
 CH3SO3m_##    = IGNORE; {@CH_3SO_3^-\aq}   {}
 CH2OHSO3m_##  = IGNORE; {@CH_2OHSO_3^-\aq} {}

{-----------------------------------------------------------------------------}
{------------------------------------ dummies --------------------------------}
{-----------------------------------------------------------------------------}
{ mz_ht_20130510+}
 NH50W_##      = IGNORE; {@NH50W}           {aqueous for NH50W tracer}
 SO2t_##       = IGNORE; {@SO2t}            {aqueous for SO2t tracer}
{ mz_ht_20130510-}

 D1O_##        = IGNORE; {@D_1O\aq}         {}
 D2O_##        = IGNORE; {@D_2O\aq}         {}
 DAHp_##       = IGNORE; {@DAH^+\aq}        {}
 DA_##         = IGNORE; {@DA\aq}           {}
 DAm_##        = IGNORE; {@DA^-\aq}         {}
 DGtAi_##      = IGNORE; {@DGtAi\aq}        {}
 DGtAs_##      = IGNORE; {@DGtAs\aq}        {}
 PROD1_##      = IGNORE; {@PROD1\aq}        {}
 PROD2_##      = IGNORE; {@PROD2\aq}        {}
 Nap_##        = IGNORE; {@Na^+\aq}         {dummy cation}

{ mz_ht_20101116+}
 Prod_01_##    = IGNORE; {@ Prod_1\aq }{aqueous phase production rate of sulphate from SO3mm + O3}
 Prod_02_##    = IGNORE; {@ Prod_2\aq }{aqueous phase production rate of sulphate from HSO3m + O3}
 Prod_03_##    = IGNORE; {@ Prod_3\aq }{aqueous phase production rate of sulphate from HSO3m + H2O2}
 Prod_04_##    = IGNORE; {@ Prod_4\aq }{aqueous phase production rate of oxalate  from GLY}
 Prod_05_##    = IGNORE; {@ Prod_5\aq }{aqueous phase production rate of sulphate from GLX + OH}
 Prod_06_##    = IGNORE; {@ Prod_6\aq }{aqueous phase production rate of sulphate from GLX- + OH}
 Prod_07_##    = IGNORE; {@ Prod_7\aq }{aqueous phase production rate of sulphate from GLX + NO3}
 Prod_08_##    = IGNORE; {@ Prod_8\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_09_##    = IGNORE; {@ Prod_9\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_10_##    = IGNORE; {@ Prod_10\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_11_##    = IGNORE; {@ Prod_11\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_12_##    = IGNORE; {@ Prod_12\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_13_##    = IGNORE; {@ Prod_13\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
 Prod_14_##    = IGNORE; {@ Prod_14\aq }{aqueous phase production rate of sulphate from GLX- + NO3}
{ mz_ht_20101116-}

{ mz_ak_20060516+} 
{ defined in gas.spc to avoid doubling}
{BrSScap       = IGNORE; }{@IGNORE} {yield of Br to SS}
{ mz_ak_20060516-}

