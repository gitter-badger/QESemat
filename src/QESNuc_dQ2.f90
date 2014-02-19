!**********************************************************************!
FUNCTION QESNuc_dQ2_Init(n_FF_QES,i_SM_TMD,i_SM_PEF)!n_FF_QES - Switch for nucleon form factor model in QES reactions
!----------------------------------------------------------------------!
!ds/dQ^2 of QES (on one nucleon)
!----------------------------------------------------------------------!
!edited by                                  ???, O.Petrova, A.Sheshukov!
!**********************************************************************!
use PhysMathConstants

implicit none

logical:: &
    QESFree_dQ2_Init,&
    QESNuc_Print_TgtParam,QESNuc_Get_TgtName,QESNuc_dQ2_Init
integer:: &
    QESNuc_Get_TgtNumb
real:: &
    MassNucleus,QESFree_dQ2,QESDeut_dQ2,&
    QESNuc_Get_TgtAWght,QESNuc_Get_TgtMass,QESNuc_Get_TgtNcln,QESNuc_dQ2
real,external:: &
    MuL_QESNuc_dQ2,GeM_FV_SM

common/MulLim/Xlow(3),Xupp(3)                                          !MuL integration limits
                                                                       !stupid MuL interface! unexterminatable common-block?!
common/Target/Target                                                   !Switch for target nucleus
common/n_SM_TMD/n_SM_TMD                                               !Switch for momentum distribution type for target nucleon
common/n_SM_PEF/n_SM_PEF                                               !Switch for momentum distribution type for outgoing nucleon (Pauli exclusion factor)
                                                                       !why this name? what does mean 'Pauli exclusion factor'?
common/Q2/Q2                                                           !Square of momentum transfer (Q^2=-q^2)
common/E_nu/E_nu                                                       !Neutrino energy
common/P_FeMax/P_FeMax                                                 !Maximum value of Fermi momentum of target nucleon
common/P_Fermi/P_Fermi                                                 !Fermi momentum of target nucleon
common/T_Fermi/T_Fermi                                                 !Effective Fermi temperature of target nucleon
common/E_nuBin/E_nuBin                                                 !Neutrino binding energy
common/FV_SM/FV_SM                                                     !Normalization factor for nuclear volume !phase space integral?
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon
common/m_lep/m_lep,mm_lep                                              !Mass and squared mass of final charged lepton
common/m_fin/m_fin,mm_fin                                              !Mass and squared mass of final hadron (system)
common/m_tar/m_tar,mm_tar                                              !Mass and squared mass of target nucleus
common/m_rnu/m_rnu,mm_rnu                                              !Mass and squared mass of residual nucleus
common/MA_ELS/MA_ELS                                                   !Mass of axial-vector in ELS reactions
common/MA_QES/MA_QES                                                   !Mass of axial-vector in QES reactions
common/MM_QES/MM_QES                                                   !Mass of monopole axial-vector in QES reactions
common/MS_QES/MS_QES                                                   !Mass of scalar in QES reactions
common/MT_QES/MT_QES                                                   !Mass of tensor in QES reactions
common/MV_QES/MV_QES                                                   !Mass of isovector in QES reactions
common/n_AG_QES/n_AG                                                   !Switch for model of QES reactions
common/n_AP_QES/n_AP                                                   !Switch for model of axial form factor in QES reactions
common/n_GE_QES/n_GE                                                   !Switch for parametrization of Sachs electric form factor of neutron
common/n_MC_QES/n_MC                                                   !Switch for correction of masses of initial and final nucleons
common/n_MS_QES/n_MS                                                   !Switch for value of axial mass in Sehgal's model
common/n_PT/n_PT                                                       !Switch for lepton polarization type
common/n_NT/n_NT                                                       !Switch for neutrino type
common/xi_A/xi_A                                                       !Normalization of axial form factor
common/xi_M/xi_M                                                       !Normalization of monopole form factor
common/xi_P/xi_P                                                       !Normalization of pseudoscalar form factor
common/xi_S/xi_S                                                       !Normalization of scalar form factor
common/xi_T/xi_T                                                       !Normalization of tensor form factor
common/xi_V/xi_V                                                       !Normalization of vector form factor
common/phi_S/phi_S                                                     !Phase of scalar form factor
common/phi_T/phi_T                                                     !Phase of tensor form factor
                                                                       !are they all indeed needed here? it seems no

integer,parameter:: &
    NNuAnu=2, NFlavor=3,&
    MinCal_GeM=100, MinCal_MuL=100
real,parameter:: &
    RelErr_GeM=1.0d-05, RelErr_MuL=5.0d-05,&
    factor=G_Fermi**2*hc2*1.0d+38/pi
character(*),parameter:: &
    hf='****************************************',vf='*'
logical &
    bufL
integer &
    NuAnu,Flavor,Target,&
    n_AG,n_AP,n_GE,n_MC,n_MS,n_PT,&
    n_NT,n_FF_QES,n_SM_PEF,n_SM_TMD,&
    iNuAnu,iFlavor,iTarget,i_SM_PEF,i_SM_TMD
real &
    E_nu,Q2,S,&
    iE_nu,iQ2,iMA_QES,&
    m_lep,mm_lep,m_ini,mm_ini,m_fin,mm_fin,m_tar,mm_tar,m_rnu,mm_rnu,&
    Z,A,E_nuBin,P_FeMax,P_Fermi,T_Fermi,&
    FV_SM,NVF_inv,NVF_var,&                                            !Normalization factor for nuclear volume !phase space integral?
    E_nu_thr,Q2_min,Q2_max
character*2 &
    TgtName

save
real &
    MA_ELS,MA_QES,MM_QES,MS_QES,MT_QES,MV_QES,&                        !save MA_QES? nonsense!
    phi_S,phi_T,xi_A,xi_M,xi_P,xi_S,xi_T,xi_V,&
    mass_lep(NFlavor)/m_e,m_mu,m_tau/,&                                !maybe into PMC? or special module? procedure-module?
    mass_fin(NNuAnu)/m_p,m_n/,&
    mass_ini(NNuAnu)/m_n,m_p/,&
    NucMass(0:NNuAnu,18),&                                             !Mass of nuclei: (A,Z),(A,Z-1),(A-1,Z-1)
    NucVolFact(NNuAnu,NFlavor,18),&                                    !Normalization factor for nuclear volume !phase space integral?
    E_nu_FrThr(NNuAnu,NFlavor),&                                       !Neutrino energy threshold for scattering on free nucleon
    E_nu_NucThr(NNuAnu,NFlavor,18),&                                   !Neutrino energy threshold for scattering on nucleus
    Xlow/3*0./,Xupp/3*1./
character*2 &
    NucName(-1:18)/'D2','H','Li','C','O','Mg','Si','Ca','Fe','Ni','Ar','Sn','Pb','N','F','Ne','Al','Cl','Cu','Br'/

!free nucleon cross section initialization-----------------------------!
    bufL=QESFree_dQ2_Init()                                            !can give E_nu_FrThr!
!settings for d3sQES_dQ2dnudkF_SM--------------------------------------!maybe, common-blocks -> arguments?
    n_SM_TMD=i_SM_TMD
    n_SM_PEF=i_SM_PEF
!settings: nucleon form factors----------------------------------------!
    call NucQESFF_Init(n_FF_QES)
!settings: integrators-------------------------------------------------!
    call GeMSet(RelErr_GeM,MinCal_GeM,*101)                            !stupid GeM interface! unused arguments!
    Xlow=0.; Xupp=1.                                                   !stupid MuL interface! unexterminatable common-block?!
    call MuLSet(RelErr_MuL,MinCal_MuL,2,*102)                          !stupid MuL interface! unused arguments!
!table filling and neutrino energy threshold calculation---------------!
    do Target=1,18
        Z=NucParam(1,Target); A=NucParam(2,Target)
        NucMass(0,Target)=MassNucleus(Z,A)                             !Initial nucleus mass
        NucMass(1,Target)=MassNucleus(Z,A-1)                           !Residual nucleus mass for neutrino
        NucMass(2,Target)=MassNucleus(Z-1,A-1)                         !Residual nucleus mass for antineutrino
!settings: common-blocks m_* for E_nu_th_SM----------------------------!
        m_tar=NucMass(0,Target); mm_tar=m_tar**2                       !this squaring is everywhere!!
        do NuAnu=1,2
            m_ini=mass_ini(NuAnu);       mm_ini=m_ini**2
            m_fin=mass_fin(NuAnu);       mm_fin=m_fin**2
            m_rnu=NucMass(NuAnu,Target); mm_rnu=m_rnu**2
!----------------------------------------------------------------------!
            P_Fermi=NucParam(4+NuAnu,Target); P_FeMax=P_Fermi          !why P_FeMax=P_Fermi?!
!nucleon volume factor calculation-------------------------------------!
            NVF_inv=4.*pi*P_FeMax                                      !invariant part of NucVolFact
            do Flavor=1,3
                call GeMInt(GeM_FV_SM,NVF_var,0.,1.,*103)              !variable part of NucVolFact
                NucVolFact(NuAnu,Flavor,Target)=NVF_inv*NVF_var
!settings: common-blocks m_* for E_nu_th_SM----------------------------!
                m_lep=mass_lep(Flavor); mm_lep=m_lep**2
!----------------------------------------------------------------------!
                call E_nu_th_SM(E_nu_NucThr(NuAnu,Flavor,Target))
            enddo
        enddo
    enddo
!free nucleon neutrino energy threshold calculation--------------------!
    do NuAnu=1,2
        do Flavor=1,3
            E_nu_FrThr(NuAnu,Flavor)=&
            0.5*((mass_fin(NuAnu)+mass_lep(Flavor))**2-mass_ini(NuAnu)**2)/mass_ini(NuAnu)
        enddo
    enddo
!----------------------------------------------------------------------!
    QESNuc_dQ2_Init=.true.
    return

!**********************************************************************!
INCLUDE 'QESNuc_Entries.f90'
!----------------------------------------------------------------------!

!**********************************************************************!
ENTRY QESNuc_dQ2(iNuAnu,iFlavor,iTarget,iE_nu,iQ2,iMA_QES)
!----------------------------------------------------------------------!
    E_nu=iE_nu
    Q2=iQ2
    MA_QES=iMA_QES
!initialization of common-block variables------------------------------!
    m_ini=mass_ini(iNuAnu);  mm_ini=m_ini**2
    m_fin=mass_fin(iNuAnu);  mm_fin=m_fin**2
    m_lep=mass_lep(iFlavor); mm_lep=m_lep**2
    if(iNuAnu==1)then
        n_NT=+1
    else
        n_NT=-1
    endif
!special case for neutrino scattering on H-----------------------------!
    if((iTarget==0).and.(iNuAnu==1))then
        QESNuc_dQ2=0.                                                  !no neutrons in Hydrogen
        return
    endif
!kinematic limit initialization----------------------------------------!
    if((iTarget==0).or.(iTarget==-1))then
        E_nu_thr=E_nu_FrThr(iNuAnu,iFlavor)
        call Q2QES_lim(E_nu,Q2_min,Q2_max)
    else
        E_nu_thr=E_nu_NucThr(iNuAnu,iFlavor,iTarget)
        call Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
    endif
!kinematic limit applying----------------------------------------------!
    if(E_nu<=E_nu_thr)then
        QESNuc_dQ2=0.
        return
    else
        if((Q2<=Q2_min).or.(Q2>=Q2_max))then
            QESNuc_dQ2=0.
            return
        else
!QES cross section calculation-----------------------------------------!
            selectcase(iTarget)
                case( 0)                                               !Hydrogen
                    QESNuc_dQ2=QESFree_dQ2(iNuAnu,iFlavor,E_nu,Q2,MA_QES)
                case(-1)                                               !Deuterium
                    QESNuc_dQ2=QESDeut_dQ2(iNuAnu,iFlavor,E_nu,Q2,MA_QES)
                case default
!initialization of common-block variables depending on iTarget---------!for every E_nu, Q2 - the same actions!!
                    m_tar=NucMass(0,iTarget);      mm_tar=m_tar**2
                    m_rnu=NucMass(iNuAnu,iTarget); mm_rnu=m_rnu**2
                    E_nuBin=NucParam(2+iNuAnu,iTarget)
                    P_Fermi=NucParam(4+iNuAnu,iTarget)
                    P_FeMax=P_Fermi                                    !why P_FeMax=P_Fermi?!
                    T_Fermi=NucParam(6+iNuAnu,iTarget)
                    FV_SM=NucVolFact(iNuAnu,iFlavor,iTarget)
!----------------------------------------------------------------------!
                    call MuLInt(MuL_QESNuc_dQ2,S,*104)
                    QESNuc_dQ2=factor*S
            endselect
!----------------------------------------------------------------------!
        endif
    endif
    return
!emergency exits-------------------------------------------------------!
101 stop 'QESNuc_dQ2_Init ERROR: GeMSet failed!'
102 stop 'QESNuc_dQ2_Init ERROR: MuLSet failed!'
103 stop 'QESNuc_dQ2_Init ERROR: GeMInt failed!'
104 stop 'QESNuc_dQ2 ERROR: MuLInt failed!'
!----------------------------------------------------------------------!
endFUNCTION QESNuc_dQ2_Init
!**********************************************************************!

!**********************************************************************!
SUBROUTINE QESNucPrintAll()
!----------------------------------------------------------------------!
!Prints properties of all nuclei
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!
logical:: &
    QESNuc_Print_TgtParam

logical &
    bufL
integer &
    Target

    do Target=-1,18
        bufL=QESNuc_Print_TgtParam(Target)
    enddo
!----------------------------------------------------------------------!
endSUBROUTINE QESNucPrintAll
!**********************************************************************!
