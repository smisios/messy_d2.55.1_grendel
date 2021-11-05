#INCLUDE atoms

{ Species are sorted by elements in the following order:                      }
{ O,H,N,C,Cl,Br,I,S                                                           }
{ All peroxides are called ROOH, all peroxy radicals are called RO2           }

{ All species are defined here with #DEFVAR as VARIABLES. Some species        }
{ will be turned into FIXED species with #SETFIX in mecca1.k                  }

#DEFVAR

{-----------------------------------------------------------------------------}
{--------------------------------- gas phase ---------------------------------}
{-----------------------------------------------------------------------------}

{------------------------------------- O -------------------------------------}

O1D           =  O                ; {@O(^1D)}            {O singlet D}
O3P           =  O                ; {@O(^3P)}            {O triplet P}
O2            = 2O                ; {@O_2}               {oxygen}
O3            = 3O                ; {@O_3}               {ozone}

{------------------------------------- H -------------------------------------}

H             =  H                ; {@H}                 {hydrogen atom}
H2            = 2H                ; {@H_2}               {hydrogen}
OH            =  H +  O           ; {@OH}                {hydroxyl radical}
HO2           =  H + 2O           ; {@HO_2}              {hydroperoxy radical}
H2O           = 2H +  O           ; {@H_2O}              {water}
H2O2          = 2H + 2O           ; {@H_2O_2}            {hydrogen peroxide}

{------------------------------------- N -------------------------------------}

N             =            N      ; {@N}                 {nitrogen atom}
N2            =           2N      ; {@N_2}               {nitrogen}
NH3           = 3H      +  N      ; {@NH_3}              {ammonia}
N2O           =       O + 2N      ; {@N_2O}              {nitrous oxide}
NO            =       O +  N      ; {@NO}                {nitric oxide}
NO2           =      2O +  N      ; {@NO_2}              {nitrogen dioxide}
NO3           =      3O +  N      ; {@NO_3}              {nitrogen trioxide}
N2O5          =      5O + 2N      ; {@N_2O_5}            {dinitrogen pentoxide}
HONO          =  H + 2O +  N      ; {@HONO}              {nitrous acid}
HNO3          =  H + 3O +  N      ; {@HNO_3}             {nitric acid}
HNO4          =  H + 4O +  N      ; {@HNO_4}             {peroxynitric acid}
NH2           = 2H      +  N      ; {@NH_2}              {}
HNO           =  H +  O +  N      ; {@HNO}               {}
NHOH          = 2H +  O +  N      ; {@NHOH}              {}
NH2O          = 2H +  O +  N      ; {@NH_2O}             {}
NH2OH         = 3H +  O +  N      ; {@NH_2OH}            {}

{------------------------------------- C -------------------------------------}

{1C}
CH3O2         =  C +  3H + 2O     ; {@CH_3O_2}           {methyl peroxy radical}
CH3OH         =  C +  4H +  O     ; {@CH_3OH}            {methanol}
CH3OOH        =  C +  4H + 2O     ; {@CH_3OOH}           {methyl peroxide}
CH4           =  C +  4H          ; {@CH_4}              {methane}
CO            =  C       +  O     ; {@CO}                {carbon monoxide}
CO2           =  C       + 2O     ; {@CO_2}              {carbon dioxide}
HCHO          =  C +  2H +  O     ; {@HCHO}              {methanal (formaldehyde)}
HCOOH         =  C +  2H + 2O     ; {@HCOOH}             {formic acid}
LCARBON       =  C + IGNORE       ; {@LC_1}              {lumped C1 species}

{2C}
C2H2          = 2C +  2H          ; {@C_2H_2}            {ethyne}
C2H4          = 2C +  4H          ; {@C_2H_4}            {ethene}
C2H5O2        = 2C +  5H + 2O     ; {@C_2H_5O_2}         {ethylperoxy radical}
C2H5OOH       = 2C +  6H + 2O     ; {@C_2H_5OOH}         {ethyl hydro peroxide}
C2H6          = 2C +  6H          ; {@C_2H_6}            {ethane}
CH3CHO        = 2C +  4H +  O     ; {@CH_3CHO}           {acetaldehyde}
CH3CO2H       = 2C +  4H + 2O     ; {@CH_3COOH}          {acetic acid}
CH3CO3        = 2C +  3H + 3O     ; {@CH_3C(O)OO}        {peroxy acetyl radical}
CH3CO3H       = 2C +  4H + 3O     ; {@CH_3C(O)OOH}       {peroxy acetic acid}
ETHGLY        = 2C +  6H + 2O     ; {@ETHGLY}            {HOCH2CH2OH}
ETHOHNO3      = 2C +  5H + 4O + N ; {@ETHOHNO3}          {HOCH2CH2ONO2}
GLYOX         = 2C +  2H + 2O     ; {@GLYOX}             {CHOCHO = glyoxal}
OXL           = 2C +  2H + 4O     ; {@OXL}               {C2H2O4 = oxalic acid}
HCOCO2H       = 2C +  2H + 3O     ; {@HCOCO_2H}          {}
HCOCO3        = 2C +   H + 4O     ; {@HCOCO_3}           {}
HCOCO3H       = 2C +  2H + 4O     ; {@HCOCO_3H}          {}
HOCH2CH2O     = 2C +  5H + 2O     ; {@HOCH_2CH_2O}       {}
HOCH2CH2O2    = 2C +  5H + 3O     ; {@HOCH_2CH_2O_2}     {}
HOCH2CHO      = 2C +  4H + 2O     ; {@HOCH_2CHO}         {glycolaldehyde}
HOCH2CO2H     = 2C +  4H + 3O     ; {@HOCH_2CO_2H}       {}
HOCH2CO3      = 2C +  3H + 4O     ; {@HOCH_2CO_3}        {}
HOCH2CO3H     = 2C +  4H + 4O     ; {@HOCH_2CO_3H}       {}
HYETHO2H      = 2C +  6H + 3O     ; {@HYETHO2H}          {HOCH2CH2OOH}
PAN           = 2C +  3H + 5O + N ; {@PAN}               {CH3C(O)OONO2 = peroxyacetylnitrate}
PHAN          = 2C +  3H + 6O     ; {@PHAN}              {HOCH2C(O)OONO2}
{3C}
ACETOL        = 3C +  6H + 2O     ; {@CH_3COCH_2OH}      {HO-CH2-CO-CH3 = hydroxy acetone}
C3H6          = 3C +  6H          ; {@C_3H_6}            {propene}
C3H8          = 3C +  8H          ; {@C_3H_8}            {propane}
CH3COCH2O2    = 3C +  5H + 3O     ; {@CH_3COCH_2O_2}     {peroxyradical from acetone}
CH3COCH3      = 3C +  6H +  O     ; {@CH_3COCH_3}        {acetone}
HOCH2COCHO    = 3C +  4H + 3O     ; {@HOCH2COCHO}        {}
HOCH2COCO2H   = 3C +  4H + 4O     ; {@HOCH2COCO2H}       {}
HYPERACET     = 3C +  6H + 3O     ; {@CH_3COCH_2O_2H}    {hydroperoxide from CH3COCH2O2}
HYPROPO2      = 3C +  7H + 3O     ; {@HYPROPO2}          {CH3CH(O2)CH2OH}
HYPROPO2H     = 3C +  8H + 3O     ; {@HYPROPO2H}         {CH3CH(OOH)CH2OH}
IC3H7NO3      = 3C +  7H + 3O + N ; {@IC_3H_7ONO_2}      {i-propyl nitrate}
IC3H7O2       = 3C +  7H + 2O     ; {@IC_3H_7O_2}        {n+iso. MCM: NC3H7O2 and IC3H7O2}
IC3H7OOH      = 3C +  8H + 2O     ; {@IC_3H_7OOH}        {n+iso. MCM: NC3H7OOH and IC3H7OOH}
MGLYOX        = 3C +  4H + 2O     ; {@MGLYOX}            {CH3COCHO = methylglyoxal}
NOA           = 3C +  5H + 4O + N ; {@NOA}               {CH3-CO-CH2ONO2 = nitro-oxy-acetone}
PR2O2HNO3     = 3C +  7H + 5O + N ; {@PR2O2HNO3}         {CH3-CH(OOH)-CH2ONO2}
PRONO3BO2     = 3C +  6H + 5O + N ; {@PRONO3BO2}         {CH3-CH(O2)-CH2ONO2}
{4C}
BIACET        = 4C +  6H + 2O     ; {@BIACET}            {CH3-CO-CO-CH3}
BIACETOH      = 4C +  6H + 3O     ; {@BIACETOH}          {CH3-CO-CO-CH2OH}
CO2H3CHO      = 4C +  5H + 3O     ; {@CO2H3CHO}          {CH3-CO-CH(OH)-CHO}
CO2H3CO3      = 4C +  5H + 5O     ; {@CO2H3CO3}          {CH3-CO-CH(OH)-C(O)O2}
CO2H3CO3H     = 4C +  6H + 5O     ; {@CO2H3CO3H}         {CH3-CO-CH(OH)-C(O)OOH}
HO12CO3C4     = 4C +  8H + 3O     ; {@HO12CO3C4}         {CH3-CO-CH(OH)-CH2OH}
LC4H9NO3      = 4C +  9H + 3O + N ; {@LC4H9NO3}          {NC4H9NO3 and SC4H9NO3}
LC4H9O2       = 4C +  9H + 2O     ; {@LC_4H_9O_2}        {CH3-CH2-CH(O2)-CH3 + CH3-CH2-CH2-CH2O2 MCM: NC4H9O2  and SC4H9O2}
LC4H9OOH      = 4C + 10H + 2O     ; {@LC_4H_9OOH}        {CH3-CH2-CH(OOH)-CH3 + CH3-CH2-CH2-CH2OOH MCM: NC4H9OOH and SC4H9OOH}
LHMVKABO2     = 4C +  7H + 4O     ; {@LHMVKABO2}         {HOCH2-CH(O2)-CO-CH3 + CH2(O2)-CH(OH)-CO-CH3}
LHMVKABOOH    = 4C +  8H + 4O     ; {@LHMVKABOOH}        {HOCH2-CH(OOH)-CO-CH3 + CH2(OOH)-CH(OH)-CO-CH3}
LMVKOHABO2    = 4C +  7H + 5O     ; {@LMVKOHABO2}        {HOCH2-CH(O2)-CO-CH2OH + CH2(O2)-CH(OH)-CO-CH2OH}
LMVKOHABOOH   = 4C +  8H + 5O     ; {@LMVKOHABOOH}       {HOCH2-CH(OOH)-CO-CH2OH + CH2(OOH)-CH(OH)-CO-CH2OH}
MACO2H        = 4C +  6H + 2O     ; {@MACO2H}            {CH2=C(CH3)COOH}
MACO3         = 4C +  5H + 3O     ; {@MACO3}             {CH2=C(CH3)C(O)O2}
MACO3H        = 4C +  6H + 2O     ; {@MACO3H}            {CH2=C(CH3)C(O)OOH}
MACR          = 4C +  6H +  O     ; {@MACR}              {CH2=C(CH3)CHO}
MACRO2        = 4C +  7H + 4O     ; {@MACRO2}            {HOCH2C(OO)(CH3)CHO}
MACROH        = 4C +  8H + 3O     ; {@MACROH}            {HOCH2C(OH)(CH3)CHO}
MACROOH       = 4C +  8H + 4O     ; {@MACROOH}           {HOCH2C(OOH)(CH3)CHO}
MEK           = 4C +  8H +  O     ; {@MEK}               {CH3-CO-CH2-CH3 = methyl ethyl ketone}
LMEKO2        = 4C +  7H + 3O     ; {@LMEKO2}            {CH3-CO-CH2-CH2-OO}
LMEKOOH       = 4C +  8H + 3O     ; {@LMEKOOH}           {CH3-CO-CH2-CH2-OOH}
MPAN          = 4C +  5H + 5O + N ; {@MPAN}              {CH2=C(CH3)C(O)OONO2 = peroxymethacryloyl nitrate ; peroxymethacrylic nitric anhydride}
MVK           = 4C +  6H +  O     ; {@MVK}               {CH3-CO-CH=CH2 = methyl vinyl ketone}
MVKOH         = 4C +  6H + 2O     ; {@MVKOH}             {CH2=CHC(=O)CH2OH}
NC4H10        = 4C + 10H          ; {@nC_4H_<10>}        {CH3-CH2-CH2-CH3 = n-butane}
{5C}
C59O2         = 5C +  9H + 5O     ; {@C59O2}             {HOCH2-CO-C(CH3)(O2)-CH2OH}
C59OOH        = 5C + 10H + 5O     ; {@C59OOH}            {HOCH2-CO-C(CH3)(OOH)-CH2OH}
C5H8          = 5C +  8H          ; {@C_5H_8}            {CH2=C(CH3)CH=CH2 = isoprene}
HCOC5         = 5C +  8H + 2O     ; {@HCOC5}             {HOCH2-CO-C(CH3)=CH2}
ISOPAOH       = 5C + 10H + 2O     ; {@ISOPAOH}           {HOCH2-C(CH3)=CH-CH2OH}
ISOPBNO3      = 5C +  9H + 4O + N ; {@ISOPBNO3}          {HOCH2-C(CH3)(ONO2)-CH=CH2}
ISOPBO2       = 5C +  9H + 3O     ; {@ISOPBO2}           {HOCH2-C(CH3)(O2)-CH=CH2}
ISOPBOH       = 5C + 10H + 2O     ; {@ISOPBOH}           {HOCH2-C(CH3)(OH)-CH=CH2}
ISOPBOOH      = 5C + 10H + 3O     ; {@ISOPBOOH}          {HOCH2-C(CH3)(OOH)-CH=CH2}
ISOPDNO3      = 5C +  9H + 4O + N ; {@ISOPDNO3}          {CH2=C(CH3)CH(ONO2)-CH2OH}
ISOPDO2       = 5C +  9H + 3O     ; {@ISOPDO2}           {CH2=C(CH3)CH(O2)-CH2OH}
ISOPDOH       = 5C + 10H + 2O     ; {@ISOPDOH}           {CH2=C(CH3)CH(OH)-CH2OH}
ISOPDOOH      = 5C + 10H + 3O     ; {@ISOPDOOH}          {CH2=C(CH3)CH(OOH)-CH2OH}
LC578O2       = 5C +  9H + 5O     ; {@LC578O2}           {HOCH2-CH(OH)C(CH3)(O2)-CHO + HOCH2-C(CH3)(O2)-CH(OH)-CHO}
LC578OOH      = 5C + 10H + 5O     ; {@LC578OOH}          {HOCH2-CH(OH)C(CH3)(OOH)-CHO + HOCH2-C(CH3)(OOH)-CH(OH)-CHO}
LC5PAN1719    = 5C +  7H + 6O + N ; {@LC5PAN1719}        {HOCH2-C(CH3)=CH-C(O)OONO2 + HOCH2-CH=C(CH3)C(O)OONO2}
LHC4ACCHO     = 5C +  8H + 2O     ; {@LHC4ACCHO}         {HOCH2-C(CH3)=CH-CHO + HOCH2-CH=C(CH3)-CHO}
LHC4ACCO2H    = 5C +  8H + 3O     ; {@LHC4ACCO2H}        {HOCH2-C(CH3)=CH-C(O)OH + HOCH2-CH=C(CH3)-C(O)OH}
LHC4ACCO3     = 5C +  7H + 4O     ; {@LHC4ACCO3}         {HOCH2-C(CH3)=CH-C(O)O2 + HOCH2-CH=C(CH3)-C(O)O2}
LHC4ACCO3H    = 5C +  8H + 4O     ; {@LHC4ACCO3H}        {HOCH2-C(CH3)=CH-C(O)OOH + HOCH2-CH=C(CH3)-C(O)OOH}
LISOPACNO3    = 5C + 10H + 4O + N ; {@LISOPACNO3}        {HOCH2-C(CH3)=CH-CH2ONO2 + HOCH2-CH=C(CH3)-CH2ONO2}
LISOPACO2     = 5C +  9H + 3O     ; {@LISOPACO2}         {HOCH2-C(CH3)=CH-CH2O2 + HOCH2-CH=C(CH3)-CH2O2}
LISOPACOOH    = 5C + 10H + 3O     ; {@LISOPACOOH}        {HOCH2-C(CH3)=CH-CH2OOH + HOCH2-CH=C(CH3)-CH2OOH}
LNISO3        = 5C + IGNORE   + N ; {@LNISO3}            {C510O2+NC4CO3 = CHO-CH(OH)-C(CH3)(O2)-CH2ONO2 + O2NOCH2-C(CH3)=CH-C(O)O2}
LNISOOH       = 5C + IGNORE   + N ; {@LNISOOH}           {CHO-CH(OH)-C(CH3)(OOH)-CH2ONO2 + O2NOCH2-C(CH3)=CH-C(O)OOH}
NC4CHO        = 5C +  7H + 4O + N ; {@nC4CHO}            {O2NOCH2-C(CH3)=CH-CHO}
NISOPO2       = 5C +  8H + 5O + N ; {@nISOPO2}           {O2NOCH2-C(CH3)=CH-CH2O2}
NISOPOOH      = 5C +  9H + 5O + N ; {@nISOPOOH}          {O2NOCH2-C(CH3)=CH-CH2OOH}

{------------------------------------- F -------------------------------------}

{------------------------------------- Cl ------------------------------------}

Cl            = Cl                ; {@Cl}                {chlorine atom}
Cl2           = 2Cl               ; {@Cl_2}              {chlorine}
ClO           = Cl + O            ; {@ClO}               {chlorine oxide}
HCl           = H + Cl            ; {@HCl}               {hydrochloric acid}
HOCl          = H + O + Cl        ; {@HOCl}              {hypochlorous acid}
Cl2O2         = 2Cl + 2O          ; {@Cl_2O_2}           {dichlorine dioxide}
OClO          = Cl + 2O           ; {@OClO}              {chlorine dioxide}
ClNO2         = Cl + 2O + N       ; {@ClNO_2}            {nitryl chloride}
ClNO3         = Cl + N + 3O       ; {@ClNO_3}            {chlorine nitrate}
CCl4          = C + 4Cl           ; {@CCl_4}             {tetrachloro methane}
CH3Cl         = C + 3H + Cl       ; {@CH_3Cl}            {chloromethane}
CH3CCl3       = 2C + 3H + 3Cl     ; {@CH_3CCl_3}         {1,1,1-trichloroethane = methyl chloroform = MCF}
CF2Cl2        = C + 2F + 2Cl      ; {@CF_2Cl_2}          {dichlorodifluoromethane = F12}
CFCl3         = C + F + 3Cl       ; {@CFCl_3}            {trichlorofluoromethane = F11}

{------------------------------------- Br ------------------------------------}

Br            = Br                ; {@Br}                {bromine atom}
Br2           = 2Br               ; {@Br_2}              {bromine}
BrO           = Br + O            ; {@BrO}               {bromine oxide}
HBr           = H + Br            ; {@HBr}               {hydrobromic acid}
HOBr          = H + O + Br        ; {@HOBr}              {hypobromous acid}
BrNO2         = Br + N + 2O       ; {@BrNO_2}            {nitryl bromide}
BrNO3         = Br + N + 3O       ; {@BrNO_3}            {bromine nitrate}
BrCl          = Br + Cl           ; {@BrCl}              {bromine chloride}
CH3Br         = Br + C +3H        ; {@CH_3Br}            {bromomethane}
CF3Br         = Br + 3F + C       ; {@CF_3Br}            {Halon 1301}
CF2ClBr       = Br + 2F + Cl + C  ; {@CF_2ClBr}          {Halon 1211}
CHCl2Br       = C + H + 2Cl + Br  ; {@CHCl_2Br}          {}
CHClBr2       = C + H + Cl + 2Br  ; {@CHClBr_2}          {}
CH2ClBr       = C + 2H + Cl + Br  ; {@CH_2ClBr}          {}
CH2Br2        = C + 2H + 2Br      ; {@CH_2Br_2}          {}
CHBr3         = C + H + 3Br       ; {@CHBr_3}            {}

{------------------------------------- I -------------------------------------}

I             = I                 ; {@I}                 {iodine atomic ground state}
I2            = 2I                ; {@I_2}               {molecular iodine}
IO            = I + O             ; {@IO}                {iodine monoxide radical}
OIO           = I + 2O            ; {@OIO}               {}
I2O2          = 2O + 2I           ; {@I_2O_2}            {}
HI            = H + I             ; {@HI}                {hydrogen iodide}
HOI           = H + O + I         ; {@HOI}               {hypoiodous acid}
HIO3          = H + I + 3O        ; {@HIO_3}             {}
INO2          = I + N + 2O        ; {@INO_2}             {iodine nitrite}
INO3          = I + N + 3O        ; {@INO_3}             {iodine nitrate}
CH3I          = C + 3H + I        ; {@CH_3I}             {iodomethane}
CH2I2         = C + 2H + 2I       ; {@CH_2I_2}           {diiodomethane}
C3H7I         = 3C + 7H + I       ; {@C_3H_7I}           {2-iodopropane}
ICl           = I + Cl            ; {@ICl}               {iodine chloride}
CH2ClI        = C + 2H + Cl + I   ; {@CH_2ClI}           {chloroiodomethane}
IBr           = I + Br            ; {@IBr}               {iodine bromide}

{------------------------------------- S -------------------------------------}

S             = S                 ; {@S}                 {sulfur atomic ground state}
SO            = S + O             ; {@SO}                {sulfur monoxide}
SO2           = S + 2O            ; {@SO_2}              {sulfur dioxide}
SH            = S + H             ; {@SH}                {}
H2SO4         = 2H + S + 4O       ; {@H_2SO_4}           {sulfuric acid}
CH3SO3H       = C + 4H + S + 3O   ; {@CH_3SO_3H}         {MSA: methane sulfonic acid}
DMS           = 2C + 6H + S       ; {@DMS}               {dimethyl sulfide}
DMSO          = 2C + 6H + S + O   ; {@DMSO}              {dimethyl sulfoxide: CH3SOCH3}
CH3SO2        = C + 3H + S + 2O   ; {@CH_3SO_2}          {}
CH3SO3        = C + 3H + S + 3O   ; {@CH_3SO_3}          {}
OCS           = C + S + O         ; {@OCS}               {}
SF6           = S + 6F            ; {@SF_6}              {sulfur hexaflouride}

{--------------------------------- Hg ----------------------------------------}

Hg            = Hg                ; {@Hg}                {}
HgO           = Hg + O            ; {@HgO}               {}
HgCl          = Hg + Cl           ; {@HgCl}              {}
HgCl2         = Hg + 2Cl          ; {@HgCl_2}            {}
HgBr          = Hg + Br           ; {@HgBr}              {}
HgBr2         = Hg + 2Br          ; {@HgBr_2}            {}
ClHgBr        = Hg + Cl + Br      ; {@ClHgBr}            {}
BrHgOBr       = Hg + O + 2Br      ; {@BrHgOBr}           {}
ClHgOBr       = Hg + O + Cl + Br  ; {@ClHgOBr}           {}

{------------------------- Generic Dummy -------------------------------------}

Dummy         = IGNORE            ; {@Dummy}             {just a dummy}

{ mz_pj_20050531+}
{------------------------- Diagnostic Tracers --------------------------------}
O3s           = 3O                ; {@O_3(s)}            {strat. ozone}
LO3s          = IGNORE            ; {@LO_3(s)}           {lost strat. ozone}
Rn_01         = IGNORE            ; {@Rn_01}
Rn_10         = IGNORE            ; {@Rn_10}
Rn_30         = IGNORE            ; {@Rn_30}
CH3I_01       = IGNORE            ; {@CH_3I_01}
CH3I_10       = IGNORE            ; {@CH_3I_10}
CH3I_30       = IGNORE            ; {@CH_3I_30}
{ mz_pj_20050531-}

{ mz_rs_20050831+}
{ozone loss and production reactions}
LossOH   = IGNORE ; {@IGNORE} {O3 + OH = HO2}
LossHO2  = IGNORE ; {@IGNORE} {O3 + HO2 = OH}
LossO1D  = IGNORE ; {@IGNORE} {O1D + H2O = 2 OH}
ProdHO2  = IGNORE ; {@IGNORE} {NO + HO2 = NO2 + OH}
ProdMeO2 = IGNORE ; {@IGNORE} {MeO2 + NO = ...}
ProdRO2  = IGNORE ; {@IGNORE} {RO2 + NO = ...}
{C2H5O2 + NO = ALD + HO2 + NO2  ...}
{CH3CO3 + NO = MeO2 + NO2 + CO2}
{IC3H7O2 + NO = .96 ACET + .96 HO...}
{LHOC3H6O2 + NO = .98 ALD + .98 H...}
{CH3COCH2O2 + NO = NO2 + CH3CO3 + HCHO...}
{LC4H9O2 + NO = .84 NO2 + .56 M...}
{MVKO2 + NO = NO2 + .25 CH3CO3 + ....}
{LMEKO2 + NO = .985 ALD + .985 ...}
{ISO2 + NO = .88 NO2 + .88 MVK...}
UNITY = IGNORE ; {UNITY=1}
{ mz_rs_20050831-}

{ op_vg_20060411+}
ProdO3   = IGNORE  ; {@IGNORE} {...+ XX = O3 }
LossO3   = IGNORE  ; {@IGNORE} {O3 + XX = YYY}
{ op_vg_20060411-}
{ op_cf_20110405+}
{ o3 prod and loss, additional diagnostic for airtrac }
LossO3Y = IGNORE ; {@IGNORE}  {O3 Loss via Y}
ProdHNO3 = IGNORE ; {@IGNORE}  {HNO3 Prod}
LossHNO3 = IGNORE ; {@IGNORE}  {HNO3 Loss}
{ op_cf_20110405-}
{ op_rd_20120124+}
DiagO3Or  = IGNORE  ; {@IGNORE} {yield from photolysis of methyl chloride}
{ op_rd_20120124-}
{ mz_ak_20060516+}
Clozone = IGNORE  ; {@IGNORE} {ozone removed by chloride}
Cldiag  = IGNORE  ; {@IGNORE} {yield from photolysis of methyl chloride}
Brhcs   = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbons}
Brcoh   = IGNORE  ; {@IGNORE} {yield from reaction with OH of bromocarbons}
BrSSrel = IGNORE  ; {@IGNORE} {Br released from SS into gasphase}
BrSScap = IGNORE  ; {@IGNORE} {yield of Br to SS}
CH3Brhcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CF3Brhcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CF2ClBrhcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CHCl2Brhcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CHClBr2hcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CH2ClBrhcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CH2Br2hcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CHBr3hcs = IGNORE  ; {@IGNORE} {yield from photolysis of bromocarbon}
CH3Brcoh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
CHCl2Brcoh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
CHClBr2coh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
CH2ClBrcoh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
CH2Br2coh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
CHBr3coh = IGNORE  ; {@IGNORE} {yield from OH reaction of bromocarbon}
{ mz_ak_20060516-}

{ op_kg_20091208+}
{ application of L41 Upper Boundary Parameterisation}
RH2O          = 2H + O           ; {@RH2O}              {H2O reservoir}
RNOy          = N                ; {@RNO_y}             {NOy reservoir}
RCly          = Cl               ; {@RCl_y}             {Cly reservoir}
RBr           = Br               ; {@RBr}               {Br reservoir}
{ for calculation of L41 UBP coefficients}
LossN2O    = IGNORE ; {@IGNORE} {N2O + O1D = RNOy}
LossN2Oa   = IGNORE ; {@IGNORE} {N2O + O1D = 2NO}
LossN2Ob   = IGNORE ; {@IGNORE} {N2O + O1D = N2 + O2}
LossCH4    = IGNORE ; {@IGNORE} {CH4 + O1D = YYY}
LossCH4OH = IGNORE ; {@IGNORE} {CH4 + OH = YYY}
LossCH4Cl = IGNORE ; {@IGNORE} {CH4 + Cl = YYY}
LossCFC    = IGNORE ; {@IGNORE} {CFC + O1D = YYY}
LossCFCOH = IGNORE ; {@IGNORE} {CFC + OH = YYY}
LossBrG    = IGNORE ; {@IGNORE} {CH3Br + OH = YYY}
LossBrJ    = IGNORE ; {@IGNORE} {CFBr  + hv = YYY}
LossN2OJ   = IGNORE ; {@IGNORE} {N2O   + hv = YYY}
LossCH4J   = IGNORE ; {@IGNORE} {CH4   + hv = YYY}
LossCFCJ   = IGNORE ; {@IGNORE} {CFC   + hv = YYY}
LossRH2O   = IGNORE ; {@IGNORE} {RH2O + M = H2O}
LossRNOy   = IGNORE ; {@IGNORE} {RNOy + M = .7 HNO3 + .3 NO2}
LossRCly   = IGNORE ; {@IGNORE} {RCly + M = HCl}
LossRBr    = IGNORE ; {@IGNORE} {RBr + M = HBr}
{ op_kg_20091208-}

{ mz_rs_20100227+}
{only for MIM1, not used in MIM2:}
LHOC3H6O2  = 3C + 7H + 3O     ; {@CH_3CH(O_2)CH_2OH} {hydroxyperoxyradical from propene+OH}
LHOC3H6OOH = 3C + 8H + 3O     ; {@CH_3CH(OOH)CH_2OH} {C3H6OHOOH = hydroxyhydroperoxides from C3H6}
ISO2       = 5C + 9H + 3O     ; {@ISO2}              {isoprene (hydroxy) peroxy radicals}
ISON       = IGNORE + N       ; {@ISON}              {organic nitrates from ISO2 and C5H8+NO3}
ISOOH      = 5C + 10H + 3O    ; {@ISOOH}             {isoprene (hydro) peroxides}
MVKO2      = 4C + 7H + 4O     ; {@MVKO2}             {MVK/MACR peroxy radicals}
MVKOOH     = 4C + 8H + 4O     ; {@MVKOOH}            {MVK hydroperoxides}
NACA       = 2C + 3H + 4O + N ; {@NACA}              {nitro-oxy acetaldehyde}
{ mz_rs_20100227-}

{ mz_ht_20130510+}
NH50W      = IGNORE; {@NH50W}           {gaseous for NH50W tracer}
SO2t       = IGNORE; {@SO2t}            {gaseous for SO2t tracer}
{ mz_ht_20130510-}

{ op_vg_20100608+}
{ ClOx cycles }
{ Cl  + O3  -> ClO +O2 }
{ ClO + O   -> Cl + O2 }
{ ClO + NO  -> Cl + NO2  (Null cycle)  }
LossO3Cl = IGNORE ; {@IGNORE}  {Loss by Cl}
LossO3H = IGNORE ; {@IGNORE}  {Loss by HOx}
LossO3R  = IGNORE ; {@IGNORE}  {Loss by RHs}
LossO3Br = IGNORE ; {@IGNORE}  {Loss by BrO}
LossO3N = IGNORE ; {@IGNORE}  {Loss by NOx}
LossO3N2 = IGNORE ; {@IGNORE}  {Loss by NOx}
LossO3O = IGNORE ; {@IGNORE}  {Loss by Ox}
{ op_vg_20100608-}

