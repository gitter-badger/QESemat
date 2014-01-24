************************************************************************
      FUNCTION EventNumbers()
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER*7 MatName
         CHARACTER*4
     #                MAn
        CHARACTER*3 tmp_string
        
              INTEGER,PARAMETER::
     #                Nel=10,
     #                MinCal   = 100,
                 REAL,PARAMETER::
     #                Xlow     = zero,
     #                Xupp     = one,
     #                RelErr   = 1.0d-13,
     #                cff_flux = 4*pi,                                   Coefficient for neutrino flux
     #                cff_sctn = 1.0d-38,                                Coefficient for section (section is multiplied by 1.00d+38)     
     #                cff_mass = N_Avogadro*1.0d+03,                     Nuclei in mol, gramms in kg - nuclei in kg multiplyed by molar mass
     #                cff_time = 60*60*24*365.25,                        Seconds in year
     #                factor   = cff_flux*cff_sctn*cff_mass*cff_time,    Coefficient for number of events per kg of detector per second multiplied by molar mass
     #                E_nu_min = 1.0d-01,                                Minimal energy given by AN spectrum
     #                E_nu_max = 1.0d+03,                                Maximal energy given by AN spectrum
     #                P_lep_min= 9.0d-02,
     #                P_lep_max= 5.0d+01

     REAL
     #                ValP(NP_lep),m_frac(Nel)/Nel*0/,
     #                Intel(Nel,NP_lep)/NPtotal*0/,
     #                     R(NP_lep)/NP_lep*0/
              INTEGER
     #                Flavor,Target,nm_TT(Nel)/Nel*0/
         CHARACTER*2 name_TT(Nel)
         COMMON   /Target/Target !Target type
         COMMON   /NuAnu/NuAnu                                           Switch for neutrino type
         COMMON     /n_NT/n_NT                                           Switch for neutrino type
         COMMON     /Flavor/Flavor                                           Switch fot lepton flavor
         COMMON    /P_lep/P_lep,E_lep                                    Charged lepton momentum
         COMMON    /x_lim/x_ini,deltax                                   Limits (for neutrino energy)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON /E_nu_thr/E_nu_thr,P_cher,O_lep,P_kth                    Neutrino energy threshold, Cherenkov threshold of lepton momentum
         COMMON      /n_l/n_l
         COMMON      /n_b/n_b
         COMMON     /MA_cen/MA_cen
         
         EXTERNAL fui

!*************************************************************************
      entry EventNumbers_Init()
!*************************************************************************
         
         CALL GeMSet(fui,one,Xlow,Xupp,RelErr,MinCal,*99)

         n_l    = 3
         N       = 1
!  ***  init neutino spectrum
         set=dFANom_dE_init()

!  ***  init QES crossections:
         n_b=2
         n_FF_QES= 8                                                     (Bodek,Avvakumov,Bradford&Budd form factor)
         buSM=dsQESCC_dQ2_SM_init(n_FF_QES,1,2)
         bufN=dsQESCC_dQ2_fN_init()
         buFP=dsQESCC_dQ2_FP_set(zero,zero,MA_QES)
         CALL setEds
! *** done
         return 
 
!*************************************************************************
      entry EventNumbers_SetTarget(TargetType)
!*************************************************************************        
         ! *** this part is for QES only:
         Target=TargetType
         X=dsQESCC_PRINT(Target)
         return
!*************************************************************************
      entry EventNumbers_CalcQES()
!*************************************************************************        
         E_lep= sqrt(P_lep**2+mm_lep)
         ! *** set neutrino energy limits
         EpP_lep= E_lep+P_lep
         EmO_lep= E_lep-O_lep
         E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2*O_lep)
         E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
         E_nu_ini=max(E_nu_low,E_nu_thr,E_nu_min)
         IF (P_lep.LE.abs(P_kth)) THEN
         IF (P_kth.LT.zero) THEN
             E_nu_fin=E_nu_ini
                            ELSE
             E_nu_fin=min(E_nu_upp,E_nu_max)
         endIF
                            ELSE
             E_nu_fin=E_nu_max
        endIF
        !**** setup integration over Enu
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-
        !**** integration over Enu
           CALL GeMInt(fui,Res,Xlow,Xupp,*100)
           Jacobianc   = 2*m_ini*P_lep/E_lep
           Intel(n_el,n_NP_lep)=m_frac(n_el)*deltax*Res*Jacobianc
           !WRITE(Nfilof,102) P_lep,Intel(n_el,n_NP_lep)
      endDO
         R=R+Intel(n_el,:)
      endIF
      endDO
        R=fact*R
        Intel=fact*Intel
         
      END FUNCTION


