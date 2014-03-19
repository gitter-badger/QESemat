!**********************************************************************!
FUNCTION fui(arg)                                                      !i like this name, but it's not informative
!----------------------------------------------------------------------!
!QES d2N_lep/(dE_lep*dE_nu)=2*m_ini*P_lep/E_lep*delta_x*fui            !is it right and clear?
!(on one nucleon)
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

real:: &
    fui,Flux_Get_dF,QESNuc_dQ2,MA_QES_eff

common/NuAnu/NuAnu                                                     !Switch: neutrino/antineutrino
common/Flavor/Flavor                                                   !Switch fot lepton flavor
common/Target/Target                                                   !Switch for target nucleus
common/x_lim/x_ini,delta_x                                             !Integration limits
common/P_lep/P_lep,E_lep                                               !Charged lepton momentum
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon

integer &
    NuAnu,Flavor,Target
real &
    x_ini,delta_x,P_lep,E_lep,m_ini,mm_ini,&
    arg,x,E_nu,Q2,&
    flux,section,Jacob_var

!change of variables---------------------------------------------------!
    x=delta_x*arg+x_ini
    E_nu=P_lep/x
    Q2=2.*m_ini*(E_nu-E_lep)
!calculation-----------------------------------------------------------!
    flux=Flux_Get_dF(NuAnu,Flavor,E_nu)
    section=QESNuc_dQ2(NuAnu,Flavor,Target,E_nu,Q2,MA_QES_eff(Target,E_nu))
    Jacob_var=E_nu*E_nu/P_lep                                          !variable part of the Jacobian
    fui=flux*section*Jacob_var                                         !*(-1)
    return
!----------------------------------------------------------------------!
endFUNCTION fui
!**********************************************************************!

!**********************************************************************!
FUNCTION fuisect(arg)
!----------------------------------------------------------------------!
!QES d2N_lep/(dE_lep*dE_nu)=delta_x*fui (on one nucleon)
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

real:: &
    fuisect,QESNuc_dQ2,MA_QES_eff

common/NuAnu/NuAnu                                                     !Switch: neutrino/antineutrino
common/Flavor/Flavor                                                   !Switch fot lepton flavor
common/Target/Target                                                   !Switch for target nucleus
common/x_lim/x_ini,delta_x                                             !Integration limits
common/E_nu/E_nu,E_lep                                               !Charged lepton momentum and energy
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon
common/m_lep/m_lep,mm_lep                                              !Mass and squared mass of final charged lepton

integer &
    NuAnu,Flavor,Target
real &
    x_ini,delta_x,P_lep,E_lep,m_ini,mm_ini,m_lep,mm_lep,&
    arg,x,E_nu,Q2

!change of variables---------------------------------------------------!
    x=delta_x*arg+x_ini
    Q2=x
    E_lep=E_nu-0.5*Q2/m_ini
    P_lep=sqrt(E_lep*E_lep-mm_lep)
!calculation-----------------------------------------------------------!
    fuisect=QESNuc_dQ2(NuAnu,Flavor,Target,E_nu,Q2,MA_QES_eff(Target,E_nu))   !*(-1)
    return
!----------------------------------------------------------------------!
endFUNCTION fuisect
!**********************************************************************!
