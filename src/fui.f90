!**********************************************************************!
FUNCTION fui(arg)                                                      !i like this name, but it's not informative
!----------------------------------------------------------------------!
!dN_lep/(dE_lep*dE_nu)=2*m_ini*P_lep/E_lep*deltax*fui                  !is it right and clear?
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

real:: &
    Flux_Get_dF,dsQESCC_dQ2_SM,MA_QES_eff                              !do we have NC QES processes? don't you want to insert "Get" into section program name?

common/NuAnu/NuAnu                                                     !Switch: neutrino/antineutrino
common/Flavor/Flavor                                                   !Switch fot lepton flavor
common/Target/Target                                                   !Switch for target nucleus
common/x_lim/x_ini,deltax                                              !Integration limits
common/P_lep/P_lep,E_lep                                               !Charged lepton momentum
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon

integer &
    NuAnu,Flavor,Target
real &
    x_ini,deltax,P_lep,E_lep,m_ini,mm_ini,&
    arg,x,E_nu,Q2,&
    flux,section,Jacobian_v,fui

!change of variables---------------------------------------------------!
    x=deltax*arg+x_ini
    E_nu=P_lep/x
    Q2=2.0*m_ini*(E_nu-E_lep)
!----------------------------------------------------------------------!
    flux=Flux_Get_dF(NuAnu,Flavor,E_nu)
    section=dsQESCC_dQ2_SM(Flavor,NuAnu,Target,E_nu,Q2,MA_QES_eff(E_nu))
    Jacobian_v=E_nu*E_nu/P_lep                                         !variable part of the Jacobian
    fui=flux*section*Jacobian_v
    return
!----------------------------------------------------------------------!
endFUNCTION fui
!**********************************************************************!
