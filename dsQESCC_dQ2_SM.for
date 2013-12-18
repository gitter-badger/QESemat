************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
*                                                                      *
*                                RedASh edition 2013/12                *
************************************************************************
      FUNCTION dsQESCC_dQ2_SM_init(n_FF,n_TMD,n_PEF)
************************************************************************
* Arguments:  n_FF - FormFactor model (passed in NucQESFF_init)
*             n_TMD - Switch for momentum distribution type for target nucleon
*             n_PEF - Switch for momentum distribution type for outgoing nucleon
************************************************************************
         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER(2) R,T

         SAVE

         LOGICAL(2),PARAMETER::
     #              OLDSCHEME =.FALSE.
            INTEGER,PARAMETER::
     #              MinCal_GeM= 100,
     #              MinCal_MuL= 100
               REAL,PARAMETER::
     #              RelErr_GeM= 1.0d-05,
     #              RelErr_MuL= 5.0d-05,
     #              factor    = G_Fermi**2*hc2*1.0d+38/pi
     
! ===============  setup arrays for target nuclei description =============
        CHARACTER(2) TGT_NAME(0:18)/
     #         'H','Li','_C','_0','Mg','Si','Ca','Fe','Ni',
     #         'Ar','Sn','Pb','_N','_F','Ne','Al','Cl','Cu','Br'/

        REAL TGT_PARAM(8,18)/ !Z,A,E_BIN[nu/an],P_Fermi[nu/an],T_Fermi[nu/an]
     #   3,   7, 0.0151, 0.0151, 0.1690, 0.1690, 0.01, 0.01,             !Li
     #   6,  12, 0.0256, 0.0257, 0.2210, 0.2210, 0.01, 0.01,             !C
     #   8,  16, 0.0266, 0.0266, 0.2250, 0.2250, 0.01, 0.01,             !O
     #   12, 24, 0.0289, 0.0290, 0.2350, 0.2350, 0.01, 0.01,             !Mg
     #   14, 28, 0.0299, 0.0300, 0.2390, 0.2390, 0.01, 0.01,             !Si
     #   20, 40, 0.0329, 0.0330, 0.2510, 0.2510, 0.01, 0.01,             !Ca
     #   26, 56, 0.0361, 0.0330, 0.2630, 0.2510, 0.01, 0.01,             !Fe
     #   28, 59, 0.0361, 0.0346, 0.2630, 0.2570, 0.01, 0.01,             !Ni
     #   18, 40, 0.0350, 0.0307, 0.2590, 0.2420, 0.01, 0.01,             !Ar
     #   50,119, 0.0391, 0.0300, 0.2740, 0.2390, 0.01, 0.01,             !Sn
     #   82,207, 0.0417, 0.0315, 0.2830, 0.2450, 0.01, 0.01,             !Pb
     #    7, 14, 0.0261, 0.0262, 0.2230, 0.2230, 0.01, 0.01,             !N
     #    9, 19, 0.0283, 0.0284, 0.2325, 0.2325, 0.01, 0.01,             !F
     #   10, 20, 0.0278, 0.0278, 0.2300, 0.2300, 0.01, 0.01,             !Ne
     #   13, 27, 0.0294, 0.0295, 0.2390, 0.2390, 0.01, 0.01,             !Al
     #   17, 36, 0.0314, 0.0315, 0.2450, 0.2450, 0.01, 0.01,             !Cl
     #   29, 64, 0.0361, 0.0346, 0.2630, 0.2570, 0.01, 0.01,             !Cu
     #   35, 80, 0.0381, 0.0315, 0.2703, 0.2450, 0.01, 0.01/             !Br
     
        REAL TGT_MASS(0:3,18)      !(A,Z)_mass, (A,Z-1)_mass, (A-1,Z-1)_mass
        REAL TGT_E_THR(3,2,18)   !E_nu threshold for Nu and ANu
        REAL TGT_FV_SM(3,2,18)   !Phase space integral for Nu and ANu
        REAL MAS_LEP(3) /m_e, m_mu,m_tau/
        REAL MAS_FIN(2) /m_p, m_n/
        REAL MAS_INI(2) /m_n, m_p/
        REAL NUC_E_THR(3,2)     !neutrino energy threshold for nucleon

         COMMON       /n_NT/n_NT                                         Switch for neutrino type
         COMMON       /n_PT/n_PT                                         Switch for lepton polarization type
         COMMON       /n_TT/n_TT                                         Switch for nuclear target type
         COMMON   /n_AP_QES/n_AP                                         Switch for model of axial form factor in QES reactions
         COMMON   /n_MS_QES/n_MS                                         Switch for value of axial mass in Sehgal's model
         COMMON   /n_GE_QES/n_GE                                         Switch for parametrization of Sachs electric form factor of neutron
         COMMON   /n_AG_QES/n_AG                                         Switch for model of QES reactions
         COMMON   /n_MC_QES/n_MC                                         Switch for correction of masses of initial and final nucleons
         COMMON   /n_SM_TMD/n_SM_TMD                                     Switch for target nucleon momentum distribution
         COMMON   /n_SM_PEF/n_SM_PEF                                     Switch for Pauli exclusion factor for outgoing nucleon
         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON      /m_tar/m_tar,mm_tar                                 Mass of target nucleus
         COMMON      /m_rnu/m_rnu,mm_rnu                                 Mass of residual nucleus
         COMMON         /Q2/Q2_tmp                                       Square of mometum transfer (Q^2=-q^2)
         COMMON       /E_nu/E_nu_tmp                                     Neutrino energy
         COMMON     /MA_QES/MA_QES_tmp                                   Mass of axial-vector in QES CC reactions
         COMMON     /MV_QES/MV_QES                                       Mass of isovector in QES reactions
         COMMON     /MM_QES/MM_QES                                       Mass of monopole axial-vector in QES reactions
         COMMON     /MS_QES/MS_QES                                       Mass of scalar in QES reactions
         COMMON     /MT_QES/MT_QES                                       Mass of tensor in QES reactions
         COMMON    /P_FeMAX/P_FeMAX                                      Maximum value of Fermi momentum of target nucleon
         COMMON    /P_Fermi/P_Fermi                                      Fermi momentum of target nucleon
         COMMON    /T_Fermi/T_Fermi                                      Effective Fermi temperature of target nucleon
         COMMON    /E_nuBIN/E_nuBIN                                      Neutrino binding energy
         COMMON      /FV_SM/FV_SM                                        Normalization factor for nuclear volume
         COMMON       /xi_V/xi_V                                         Normalization of vector form factor
         COMMON       /xi_M/xi_M                                         Normalization of monopole form factor
         COMMON       /xi_S/xi_S                                         Normalization of scalar form factor
         COMMON       /xi_A/xi_A                                         Normalization of axial form factor
         COMMON       /xi_P/xi_P                                         Normalization of pseudoscalar form factor
         COMMON       /xi_T/xi_T                                         Normalization of tensor form factor
         COMMON      /phi_T/phi_T                                        Phase of tensor form factor
         COMMON      /phi_S/phi_S                                        Phase of scalar form factor
         COMMON     /MulLim/Xlow(3),Xupp(3)                              MuL integration limits

         
*************** initialization *************************
        EXTERNAL MuL_dsQESCC_dQ2_SM,GeM_FV_SM

        n_SM_TMD= n_TMD  !1                                                Switch for target nucleon momentum distribution
        n_SM_PEF= n_PEF  !2                                                Switch for Pauli exclusion factor for outgoing nucleon

        CALL NucQESFF(n_FF)

        Xlow= zero
        Xupp= one
************** initialize integrators (GeM and MuL) ************
        CALL MuLSet(MuL_dsQESCC_dQ2_SM,R1,RelErr_MuL,MinCal_MuL,2,*101)
        CALL GeMSet(GeM_FV_SM,R1,zero,one,RelErr_GeM,MinCal_GeM,  *102)
**************  Fill tables (TGT_FV_SM, TGT_E_THR)  ***********
        DO n_TGT=1,18
            TGT_MASS(0,n_TGT) = MassNucleus(Z,A)     !Initial nucleus mass
            TGT_MASS(1,n_TGT) = MassNucleus(Z,A-1)   !Residual nucl mass (Nu)
            TGT_MASS(2,n_TGT) = MassNucleus(Z-1,A-1) !Residual nucl mass (ANu)
            m_tar=TGT_MASS(0,n_TGT)
            mm_tar= m_tar**2
            
            !T_Fermi= 1.0d-02
            DO n_NA=1,2
                m_ini=MAS_INI(n_NA)  !initial nucleon mass
                mm_ini= m_ini**2
                m_fin=MAS_FIN(n_NA)  !final nucleon mass
                mm_fin= m_fin**2
                m_rnu=TGT_MASS(n_NA,n_TGT) !residual nucleus mass
                mm_rnu= m_rnu**2
                P_Fermi=TGT_PARAM(4+n_NA,n_TGT)
                P_FeMAX=P_Fermi
                DO n_Flv=1,3
                    m_lep=MAS_LEP(n_Flv)
                    mm_lep= m_lep**2
                    !calculate FV_SM
                    CALL GeMInt(GeM_FV_SM,FV_SM,zero,one,*103)
                    TGT_FV_SM(n_Flv,n_NA,n_TGT)=FV_SM
                    !calculate Nu energy threshold
                    CALL E_nu_th_SM(E_nu_thr)
                    TGT_E_THR(n_Flv,n_NA,n_TGT)=E_nu_thr
                    !TGT_PARAM(6+n_NA,n_TGT)=T_Fermi
                endDO
            endDO
!             dsQESCC_PRINT(n_TGT)
        endDO
        DO n_NA=1,2
            DO n_Flv=1,3
              NUC_E_THR(n_Flv,n_NA)= 
     #         ((MAS_FIN(n_NA)+MAS_LEP(n_Flv))**2-
     #          MAS_INI(n_NA)**2)/(2*MAS_INI(n_NA))
            endDO
        endDO
        dsQESCC_dQ2_SM_init= one
        
      RETURN
  101    STOP 'ERROR IN MuLSet. FUNCTION dsQESCC_dQ2_SM_init'
  102    STOP 'ERROR IN GeMSet. FUNCTION dsQESCC_dQ2_SM_init'
  103    STOP 'ERROR IN MuLInt. FUNCTION dsQESCC_dQ2_SM_init'
  104    STOP 'ERROR IN GeMInt. FUNCTION dsQESCC_dQ2_SM_init'
     
*     ==================================================================
      ENTRY dsQESCC_PRINT(n_TARGET)
*     ==================================================================
        WRITE(*,'(A4,A7,I3,A1,A2)')
     #    "*** ","Target#",n_TARGET,"=",TGT_NAME(n_TARGET)
        WRITE(*,'(A4,A2,F3.0,A3,F3.0)')
     #    "*** ","Z=",TGT_PARAM(1,n_TARGET),"A=",TGT_PARAM(2,n_TARGET)
        WRITE(*,'(A4,A8,F5.3,A1,F5.3)')
     #    "*** ","Enu_BIN=",TGT_PARAM(3,n_TARGET),
     #    "/",TGT_PARAM(4,n_TARGET)
        WRITE(*,'(A4,A8,F5.3,A1,F5.3)')
     #    "*** ","P_Fermi=",TGT_PARAM(5,n_TARGET),
     #    "/",TGT_PARAM(6,n_TARGET)
        WRITE(*,'(A4,A8,F5.3,A1,F5.3)')
     #    "*** ","T_Fermi=",TGT_PARAM(7,n_TARGET),
     #    "/",TGT_PARAM(8,n_TARGET)
      RETURN


*     ==================================================================
      ENTRY dsQESCC_dQ2_SM(n_Fl,n_NuAnu,n_TARGET,E_nu,Q2,MA_QES)
*     ==================================================================
*      ONE function for cross-section calculation, 
*      given neutrino type and flavour

*     Init kinematic variables
         m_ini     = MAS_INI(n_NuAnu); mm_ini    = m_ini**2
         m_fin     = MAS_FIN(n_NuAnu); mm_fin    = m_fin**2
         m_lep     = MAS_LEP(n_Lep);   mm_lep    = m_lep**2
*     Fill common block vars 
         E_nu_tmp  = E_nu
         Q2_tmp    = Q2
         MA_QES_tmp= MA_QES
         IF(n_NuAnu.EQ.1)THEN
            n_NT=+1
         ELSE
            n_NT=-1
         endIF
         
         IF(n_TARGET.EQ.0)THEN
            ! Hydrogen:
            IF(n_NuAnu.EQ.1)THEN
            ! no neutrons in Hydrogen
                dsQESCC_dQ2_SM=Precision
            ELSE
                SELECTCASE(n_Fl)
                  CASE(1) !nu_e 
                    dsQESCC_dQ2_SM= dsQESCC_dQ2_fN_ea(E_nu,Q2,MA_QES)
                  CASE(2) !nu_e 
                    dsQESCC_dQ2_SM= dsQESCC_dQ2_fN_ma(E_nu,Q2,MA_QES)
                  CASE(3) !nu_e 
                    dsQESCC_dQ2_SM= dsQESCC_dQ2_fN_ta(E_nu,Q2,MA_QES)
                endSELECT
            endIF
         ELSE
             IF (E_nu.le.TGT_E_THR(n_Fl,n_nuAnu,n_TARGET)) THEN                             O  (OXYGEN),    Z= 8
                 dsQESCC_dQ2_SM=Precision
                                        ELSE
!                 CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
!                 IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
!                   dsQESCC_dQ2_SM=Precision
!                                                     ELSE
            m_tar=TGT_MASS(0,n_TARGET); mm_tar =m_tar**2
            m_rnu=TGT_MASS(n_NuAnu,n_TARGET); mm_rnu =m_rnu**2
            
            n_Z=TGT_PARAM(1,n_TARGET)
            n_A=TGT_PARAM(2,n_TARGET)
            IF(n_nuanu.eq.1) THEN
              n_N=n_A-n_Z
            ELSE
              n_N=n_Z
            endIF
            
            E_nuBIN=TGT_PARAM(2+n_NuAnu,n_TARGET)
            P_Fermi=TGT_PARAM(4+n_NuAnu,n_TARGET)
            P_FeMAX=TGT_PARAM(4+n_NuAnu,n_TARGET)
            T_Fermi=TGT_PARAM(6+n_NuAnu,n_TARGET)
            FV_SM=TGT_FV_SM(n_Fl,n_NuAnu,n_TARGET)
            CALL MuLInt(MuL_dsQESCC_dQ2_SM,S,*104)
            dsQESCC_dQ2_SM=factor*n_N*S
!              endIF
            endIF
         endIF
        RETURN
      END FUNCTION dsQESCC_dQ2_SM_init
****************************************************************
      FUNCTION dsQESCC_PRINT_ALL()
        DO N=1,18
            X=dsQESCC_PRINT(N)
        endDO
      END FUNCTION dsQESCC_PRINT_ALL

