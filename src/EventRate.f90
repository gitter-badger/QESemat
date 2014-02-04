!**********************************************************************!
FUNCTION EventRate_Init_Flux(fluxfile,iNuAnu,iFlavor)
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!
use PhysMathConstants

implicit none

logical:: &
    EventRate_Init_Flux,EventRate_Init_Section,EventRate_Set_Tgt,&
    Flux_Init,Flux_Open_File,Flux_Read_Hdr,Flux_Read_Table,Flux_Close_File,&
    Flux_Has_Table,Flux_Print_Table,Flux_Calc_Spline,&
    MA_QES_Init,QES_print
real:: &
    EventRate_Calc,&
    QES_Get_TgtNcln,&
    Flux_Get_Emin,Flux_Get_Emax,&
    QES_init,dsqescc_dq2_fn_init,dsqescc_dq2_fp_set                    !rename dsqescc_dq2_*!
real,external:: &
    fui                                                                !Integrand

common/NuAnu/NuAnu                                                     !Switch: neutrino/antineutrino
common/Flavor/Flavor                                                   !Switch fot lepton flavor
common/Target/Target                                                   !Switch for target nucleus
common/P_lep/P_lep,E_lep                                               !Charged lepton momentum and energy
common/x_lim/x_ini,deltax                                              !Integration limits
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon
common/MA_QES/MA_QES                                                   !Mass of axial-vector in QES reactions
common/MA_ELS/MA_ELS                                                   !Mass of axial-vector in ELS reactions
common/n_PT/n_PT                                                       !Switch for lepton polarization type
common/n_AG_QES/n_AG                                                   !Switch for model of QES reactions
common/n_MC_QES/n_MC                                                   !Switch for correction of masses of initial and final nucleons
common/n_AP_QES/n_AP                                                   !Switch for model of axial form factor in QES reactions
common/n_MS_QES/n_MS                                                   !Switch for value of axial mass in Sehgal's model
common/n_GE_QES/n_GE                                                   !Switch for parametrization of Sachs electric form factor of neutron
common/MV_QES/MV_QES                                                   !Mass of isovector in QES reactions
common/MM_QES/MM_QES                                                   !Mass of monopole axial-vector in QES reactions
common/MS_QES/MS_QES                                                   !Mass of scalar in QES reactions
common/MT_QES/MT_QES                                                   !Mass of tensor in QES reactions
common/xi_V/xi_V                                                       !Normalization of vector form factor
common/xi_M/xi_M                                                       !Normalization of monopole form factor
common/xi_S/xi_S                                                       !Normalization of scalar form factor
common/xi_A/xi_A                                                       !Normalization of axial form factor
common/xi_P/xi_P                                                       !Normalization of pseudoscalar form factor
common/xi_T/xi_T                                                       !Normalization of tensor form factor
common/phi_T/phi_T                                                     !Phase of tensor form factor
common/phi_S/phi_S                                                     !Phase of scalar form factor

integer,parameter:: &
    n_FF_QES=8                                                         !NucQESFF setting: Bodek,Avvakumov,Bradford&Budd form factor
logical &
    bufL
integer &
    NuAnu,Flavor,Target,CorV,&
    n_AG,n_AP,n_GE,n_MC,n_MS,n_PT,&
    iNuAnu,iFlavor,iTarget
real &
    bufR,Jacobian_c,MA_QES,Res,&
    x_ini,x_fin,deltax,P_lep,E_lep,m_ini,mm_ini,&
    EmO_lep,EpP_lep,E_nu_ini,E_nu_fin,E_nu_low,E_nu_upp,&
    MA_ELS,MM_QES,MS_QES,MT_QES,MV_QES,phi_S,phi_T,xi_A,xi_M,xi_P,xi_S,xi_T,xi_V,&
    iMA_QES,iP_lep
character*80 &
    fluxfile

save
real &
    E_nu_min,E_nu_max,m_lep,mm_lep,&
    O_lep,P_kth

    NuAnu=iNuAnu
    Flavor=iFlavor
!settings: neutrino flux-----------------------------------------------!settings FROM and INTO Flux
    bufL=Flux_init()
    bufL=Flux_open_file(fluxfile)
    do while(Flux_Read_Hdr())
        bufL=Flux_read_table()
    enddo
    bufL=Flux_close_file()
    if(Flux_has_table(NuAnu,Flavor).eqv..FALSE.)then
        write(*,*) 'QESemat ERROR: Flux table doesn''t exist!'
        stop
    endif
    bufL=Flux_print_table(NuAnu,Flavor)                                !do we really need this?
    bufL=Flux_calc_spline(NuAnu,Flavor)
         
    E_nu_min=Flux_Get_Emin(NuAnu,Flavor)
    E_nu_max=Flux_Get_Emax(NuAnu,Flavor)
    write(*,'(2(A9,F10.3,1X))') 'E_nu_min=',E_nu_min,'E_nu_max=',E_nu_max
    EventRate_Init_Flux=.true.
    return
!**********************************************************************!
ENTRY EventRate_Init_Section(iNuAnu,iFlavor,CorV,iMA_QES)
!----------------------------------------------------------------------!
    NuAnu=iNuAnu
    Flavor=iFlavor
    MA_QES=iMA_QES
!settings: MA_QES_eff--------------------------------------------------!setting INTO MA_QES_eff
    bufL=MA_QES_Init(CorV,0,MA_QES)
!settings: *QES*-------------------------------------------------------!setting INTO *QES*
    bufR=QES_init(n_FF_QES,1,2)
!    CALL QES_PRINT_ALL()
    bufR=dsQESCC_dQ2_fN_init()
    bufR=dsQESCC_dQ2_FP_set(0.,0.,MA_QES)                              !rewrite like dsQESCC_dQ2_fN_init!
!settings: other-------------------------------------------------------!
    MA_ELS=MA_QES
    n_PT=0; n_AG=1; n_MC=1; n_AP=1; n_MS=1; n_GE=5
    MV_QES=0.84; MM_QES=0.8; MS_QES=1.; MT_QES=1.5
    xi_V=1.; xi_M=1.; xi_S=0.; xi_A=1.; xi_P=1.; xi_T=0.
    phi_T=0.; phi_T=phi_T*dtr
    phi_S=0.; phi_S=phi_S*dtr
    if(NuAnu==1)then
        m_ini=m_n; mm_ini=mm_n
    else
        m_ini=m_p; mm_ini=mm_p
    endif
    selectcase(Flavor)
        case(1)
            m_lep=m_e; mm_lep=mm_e
        case(2)
            m_lep=m_mu; mm_lep=mm_mu
        case(3)
            m_lep=m_tau; mm_lep=mm_tau
    endselect
    O_lep=0.5*mm_lep/m_I
    P_kth=0.5*m_I-O_lep
    EventRate_Init_Section=.true.
    return
!**********************************************************************!
ENTRY EventRate_Set_Tgt(iTarget)
!----------------------------------------------------------------------!
    Target=iTarget
    bufL=QES_print(Target)
    EventRate_Set_Tgt=.true.
    return
!**********************************************************************!
ENTRY EventRate_Calc(iP_lep)
!----------------------------------------------------------------------!
    P_lep=iP_lep
    E_lep=sqrt(P_lep*P_lep+mm_lep)
    EpP_lep=E_lep+P_lep
    EmO_lep=E_lep-O_lep
    E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2.*O_lep)
    E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
    E_nu_ini=max(E_nu_low,E_nu_min)
    if(P_lep<=abs(P_kth))then
        if(P_kth<0.)then
            E_nu_fin=E_nu_ini
        else
            E_nu_fin=min(E_nu_upp,E_nu_max)
        endif
    else
        E_nu_fin=E_nu_max
    endif
    x_ini=P_lep/E_nu_fin
    x_fin=P_lep/E_nu_ini
    deltax=x_fin-x_ini                                                 !-
    call GeMInt(fui,Res,0.,1.,*100)
    Jacobian_c=2.*m_ini*P_lep/E_lep
    EventRate_Calc=QES_Get_TgtNcln(NuAnu,Target)*deltax*Res*Jacobian_c
    return
!emergency exits-------------------------------------------------------!
100 stop 'EventRate ERROR: GeMInt failed!'
!----------------------------------------------------------------------!
endFUNCTION EventRate_Init_Flux
!**********************************************************************!
