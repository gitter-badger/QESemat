!**********************************************************************!
FUNCTION GeM_FV_SM(arg)
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                               ???, O.Petrova!
!**********************************************************************!

implicit none

real:: &
    GeM_FV_SM,rho_SM

common/n_SM_TMD/n_SM_TMD                                               !Switch for target nucleon momentum distribution
                                                                       !name!
common/P_FeMAX/P_FeMAX                                                 !Maximum value of Fermi momentum of target nucleon
                                                                       !what's the difference between P_FeMAX and P_Fermi?! name!

integer &
    n_SM_TMD
real &
    arg,p,P_FeMAX

!change of variables---------------------------------------------------!
    p=arg*P_FeMAX
!----------------------------------------------------------------------!
    GeM_FV_SM=p*p*rho_SM(n_SM_TMD,p)
    return
!----------------------------------------------------------------------!
endFUNCTION GeM_FV_SM
!**********************************************************************!
