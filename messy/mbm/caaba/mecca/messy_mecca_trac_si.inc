! -*- f90 -*- created automatically by xmecca, DO NOT EDIT!
CALL new_tracer(status, setname, 'BrNO2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_BrNO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_BrNO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF2ClBr',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF2ClBr)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF2ClBr, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClBr, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClBr, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClBr, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF3Br',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF3Br)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF3Br, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3Br, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3Br, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3Br, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3I',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3I)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3I, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3I, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3I, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3I, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'O3s',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'O3', &
  medium               = AIR, &
  idx                  = idt_O3s)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_O3s, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3s, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3s, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3s, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3s, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LCARBON',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'C', &
  medium               = AIR, &
  idx                  = idt_LCARBON)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LCARBON, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LCARBON, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LCARBON, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3SO3H',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3SO3H)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3SO3H, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3H, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3H, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3H, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3H, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NO3m',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  subname              = 'cs', &
  refspec              = 'NO3', &
  medium               = AEROSOL, &
  idx                  = idt_NO3m_cs)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_sedi              , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_aerosol_method    , MODAL)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, S_aerosol_model     , c_pa_asm)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3m_cs, I_aerosol_mode      , i_pa_amode)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Hp',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  subname              = 'cs', &
  refspec              = 'H', &
  medium               = AEROSOL, &
  idx                  = idt_Hp_cs)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Hp_cs, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_sedi              , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_aerosol_method    , MODAL)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, S_aerosol_model     , c_pa_asm)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Hp_cs, I_aerosol_mode      , i_pa_amode)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LO3s',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'O3', &
  medium               = AIR, &
  idx                  = idt_LO3s)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LO3s, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LO3s, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LO3s, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LO3s, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LO3s, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LossCH4',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'ONE', &
  medium               = AIR, &
  idx                  = idt_LossCH4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LossCH4, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCH4, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCH4, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCH4, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCH4, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LossCO',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'ONE', &
  medium               = AIR, &
  idx                  = idt_LossCO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LossCO, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCO, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCO, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCO, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossCO, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LossN2O',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'ONE', &
  medium               = AIR, &
  idx                  = idt_LossN2O)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LossN2O, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossN2O, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossN2O, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossN2O, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossN2O, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LossO3',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'ONE', &
  medium               = AIR, &
  idx                  = idt_LossO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LossO3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossO3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossO3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossO3, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LossO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ProdO3',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'ONE', &
  medium               = AIR, &
  idx                  = idt_ProdO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ProdO3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ProdO3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ProdO3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ProdO3, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ProdO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'BIACET',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_BIACET)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_BIACET, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BIACET, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BIACET, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BIACET, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BIACET, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Cl2O2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_Cl2O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Cl2O2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2O2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2O2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'S',              modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_S)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_S, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_S, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_S, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_S, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C2H6',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C2H6)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C2H6, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H6, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H6, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H6, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NC4H10',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NC4H10)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NC4H10, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NC4H10, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NC4H10, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NC4H10, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NC4H10, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CCl4',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CCl4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CCl4, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CCl4, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CCl4, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CCl4, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF2Cl2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF2Cl2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF2Cl2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2Cl2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2Cl2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2Cl2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF2ClCF2Cl',     modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF2ClCF2Cl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF2ClCF2Cl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCF2Cl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCF2Cl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCF2Cl, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF3CF2Cl',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF3CF2Cl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF3CF2Cl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3CF2Cl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3CF2Cl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF3CF2Cl, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CFCl3',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CFCl3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CFCl3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CFCl3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CFCl3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CFCl3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH2Br2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH2Br2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH2Br2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Br2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Br2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Br2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Br2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHBr3',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHBr3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHBr3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHBr3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHBr3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHBr3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHBr3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3SO3',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3SO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3SO3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NH3',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NH3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NH3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C3H8',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C3H8)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C3H8, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H8, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H8, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H8, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CF2ClCFCl2',     modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CF2ClCFCl2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CF2ClCFCl2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCFCl2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCFCl2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CF2ClCFCl2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ClNO2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ClNO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ClNO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'OClO',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_OClO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_OClO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OClO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OClO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OClO, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH2ClBr',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH2ClBr)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH2ClBr, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2ClBr, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2ClBr, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2ClBr, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2ClBr, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3Br',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3Br)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3Br, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Br, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Br, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Br, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Br, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHCl2Br',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHCl2Br)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHCl2Br, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl2Br, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl2Br, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl2Br, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl2Br, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHClBr2',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHClBr2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHClBr2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHClBr2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHClBr2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHClBr2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHClBr2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'H2SO4',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_H2SO4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_H2SO4, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2SO4, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2SO4, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2SO4, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2SO4, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'SH',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_SH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_SH, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SH, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SH, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SH, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHF3',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHF3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHF3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH2Cl2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH2Cl2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH2Cl2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Cl2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Cl2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2Cl2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CCl3',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CCl3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CCl3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CCl3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CCl3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CCl3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHCl3',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHCl3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHCl3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHCl3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'OCS',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_OCS)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_OCS, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OCS, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OCS, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OCS, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NACA',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NACA)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NACA, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NACA, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NACA, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NACA, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'N',              modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_N)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_N, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'N2O',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_N2O)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_N2O, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NH2OH',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NH2OH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NH2OH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2OH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2OH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2OH, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2OH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'SO3',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_SO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_SO3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3Cl',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3Cl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3Cl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Cl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Cl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Cl, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3Cl, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NHOH',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NHOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NHOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NHOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NHOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NHOH, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NHOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MPAN',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MPAN)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MPAN, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MPAN, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MPAN, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MPAN, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MPAN, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LMEKOOH',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LMEKOOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LMEKOOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKOOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKOOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKOOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKOOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CHF2',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CHF2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CHF2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHF2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHF2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHF2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C2H5OOH',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C2H5OOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C2H5OOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5OOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5OOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5OOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5OOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CO3H',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CO3H)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CO3H, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3H, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3H, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3H, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3H, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'SO',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_SO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_SO, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'IC3H7NO3',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_IC3H7NO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_IC3H7NO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7NO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7NO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7NO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7NO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ISOOH',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ISOOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ISOOH, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISOOH, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISOOH, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISOOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISOOH, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LHOC3H6OOH',     modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LHOC3H6OOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LHOC3H6OOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6OOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6OOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6OOH, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'IC3H7OOH',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_IC3H7OOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_IC3H7OOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7OOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7OOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7OOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7OOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HNO4',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HNO4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HNO4, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO4, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO4, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO4, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO4, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CO2H',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CO2H)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CO2H, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO2H, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO2H, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO2H, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO2H, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'DMSO',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_DMSO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_DMSO, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMSO, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMSO, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMSO, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMSO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HYPERACET',      modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HYPERACET)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HYPERACET, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HYPERACET, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HYPERACET, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HYPERACET, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HYPERACET, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HONO',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HONO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HONO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HONO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HONO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HONO, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HONO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH2FCF3',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH2FCF3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH2FCF3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2FCF3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2FCF3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2FCF3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C2H2',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C2H2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C2H2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CFCl2',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CFCl2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CFCl2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CFCl2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CFCl2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CFCl2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHF2CF3',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHF2CF3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHF2CF3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2CF3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2CF3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2CF3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CF3',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CF3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CF3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CF3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CF3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CF3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH2F2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH2F2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH2F2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2F2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2F2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH2F2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'PAN',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_PAN)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_PAN, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_PAN, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_PAN, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_PAN, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_PAN, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HCOOH',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HCOOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HCOOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCOOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCOOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCOOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCOOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CHF2Cl',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CHF2Cl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CHF2Cl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2Cl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2Cl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CHF2Cl, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Cl2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_Cl2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Cl2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LC4H9OOH',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LC4H9OOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LC4H9OOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9OOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9OOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9OOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9OOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3SO2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3SO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3SO2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3SO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MEK',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MEK)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MEK, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MEK, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MEK, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MEK, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MEK, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MVKOOH',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MVKOOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MVKOOH, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKOOH, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKOOH, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKOOH, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'N2O5',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_N2O5)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_N2O5, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O5, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O5, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O5, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2O5, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NH2O',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NH2O)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NH2O, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2O, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2O, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2O, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2O, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HNO',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HNO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HNO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3COCH3',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3COCH3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3COCH3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'H2O2',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_H2O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_H2O2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2O2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2O2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2O2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2O2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3OH',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3OH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3OH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'BrCl',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_BrCl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_BrCl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrCl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrCl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrCl, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrCl, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NH2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NH2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NH2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2, I_drydep            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NH2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ISON',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ISON)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ISON, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISON, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISON, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISON, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISON, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'IC3H7O2',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_IC3H7O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_IC3H7O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_IC3H7O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3COCH2O2',     modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3COCH2O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3COCH2O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH2O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH2O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3COCH2O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MGLYOX',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MGLYOX)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MGLYOX, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MGLYOX, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MGLYOX, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MGLYOX, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MGLYOX, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CO',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C3H6',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C3H6)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C3H6, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H6, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H6, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C3H6, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LMEKO2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LMEKO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LMEKO2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKO2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKO2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LMEKO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'H2',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_H2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_H2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH4',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH4, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH4, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH4, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH4, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C2H4',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C2H4)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C2H4, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H4, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H4, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H4, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Br2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_Br2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Br2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HNO3',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HNO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HNO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HNO3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LC4H9O2',        modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LC4H9O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LC4H9O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'BrNO3',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_BrNO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_BrNO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrNO3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'SO2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_SO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_SO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_SO2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3OOH',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3OOH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3OOH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OOH, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OOH, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OOH, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3OOH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C5H8',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C5H8)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C5H8, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C5H8, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C5H8, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C5H8, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C5H8, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ACETOL',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ACETOL)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ACETOL, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ACETOL, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ACETOL, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ACETOL, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ACETOL, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ISO2',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ISO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ISO2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISO2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISO2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ISO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MVK',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MVK)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MVK, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVK, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVK, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVK, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVK, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LC4H9NO3',       modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LC4H9NO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LC4H9NO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9NO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9NO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9NO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LC4H9NO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'MVKO2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_MVKO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_MVKO2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKO2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKO2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_MVKO2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HOCl',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HOCl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HOCl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOCl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOCl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOCl, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOCl, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'DMS',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_DMS)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_DMS, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMS, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMS, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_DMS, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'LHOC3H6O2',      modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_LHOC3H6O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_LHOC3H6O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_LHOC3H6O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ClNO3',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ClNO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ClNO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClNO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'C2H5O2',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_C2H5O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_C2H5O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_C2H5O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HOBr',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HOBr)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HOBr, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOBr, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOBr, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOBr, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HOBr, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CHO',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CHO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CHO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHO, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CHO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'O1D',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'O', &
  medium               = AIR, &
  idx                  = idt_O1D)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_O1D, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O1D, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O1D, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O1D, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3CO3',         modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3CO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3CO3, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3CO3, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'H',              modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_H)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_H, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_H, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HBr',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HBr)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HBr, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HBr, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HBr, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HBr, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HBr, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'ClO',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_ClO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_ClO, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClO, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClO, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_ClO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'O3',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_O3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_O3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'O3P',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  refspec              = 'O', &
  medium               = AIR, &
  idx                  = idt_O3P)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_O3P, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3P, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3P, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O3P, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'OH',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_OH)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_OH, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OH, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OH, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_OH, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HO2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HO2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HO2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HO2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HO2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HCl',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HCl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HCl, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCl, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCl, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCl, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCl, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'BrO',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_BrO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_BrO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_BrO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CH3O2',          modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CH3O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CH3O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3O2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CH3O2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'HCHO',           modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_HCHO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_HCHO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCHO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCHO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCHO, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_HCHO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NO',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NO)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NO, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NO2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO2, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO2, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'NO3',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_NO3)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_NO3, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3, I_drydep            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_NO3, I_scav              , ON)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Br',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_Br)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Br, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Br, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'Cl',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_Cl)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_Cl, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_Cl, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'O2',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_O2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_O2, R_vini              , 0.21_dp)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_O2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'N2',             modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_N2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_N2, R_vini              , 0.78_dp)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2, I_advect            , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2, I_convect           , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2, I_vdiff             , OFF)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_N2, I_scav              , OFF)
CALL tracer_halt(substr, status)
CALL new_tracer(status, setname, 'CO2',            modstr, &
  quantity             = AMOUNTFRACTION, &
  unit                 = 'mol/mol', &
  medium               = AIR, &
  idx                  = idt_CO2)
CALL tracer_halt(substr,status)
CALL set_tracer(status, setname, idt_CO2, R_vini              , 3.5E-04_dp)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO2, I_advect            , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO2, I_convect           , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO2, I_vdiff             , ON)
CALL tracer_halt(substr, status)
CALL set_tracer(status, setname, idt_CO2, I_scav              , ON)
CALL tracer_halt(substr, status)
