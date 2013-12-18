************************************************************************
      PROGRAM qesemat
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         CHARACTER*80
     #                arg
         CHARACTER*4
     #                MAn
 
              INTEGER,PARAMETER::
     #                Nfilof   = 100,
     #                Nfilobn  = 101,
     #                Nfiloba  = 102,
     #                Nfiloe   = 103,
     #                NP_lep   = 100,
     #                MinCal   = 100
                 REAL,PARAMETER::
     #                Xlow     = zero,
     #                Xupp     = one,
     #                RelErr   = 1.0d-13,
     #                f1       = c2C/(8*pi)*mm_I*G_Fermi**2*hc2,         Coefficient for section
     #                f2       = 4*pi,                                   Coefficient for neutrino flux     
     #                f3       = 1.055d+39,                              Coefficient for number of events
     #                factorf  = 2*f1*f2*f3*mm_W**2,                     (m_W is from W-boson propagator)
     #                factorb  = 8*f2*f3*1.00d-38,                       (section is multiplied by 1.00d+38)
     #                E_nu_min = 1.0d-01,                                Minimal energy given by AN spectrum
     #                E_nu_max = 1.0d+03,                                Maximal energy given by AN spectrum
     #                P_lep_min= 9.0d-02,
     #                P_lep_max= 5.0d+01
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
         CHARACTER*128
     #                namfof,namfobn,namfoba,namfoe

         CHARACTER*3
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
         CHARACTER*1
     #                CorV(2)/'c','v'/,
     #                nb(3)/'l','c','r'/,
     #                nl(5)/'0','n','m','x','1'/,
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/
                 REAL
     #                ValP(NP_lep),
     #                Intf(NP_lep),Intbn(NP_lep),Intba(NP_lep)

         COMMON     /n_MA/n_MA                                           Switch for MA_QES
         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON   /N_CorV/N_CorV
         COMMON   /N_TT/N_TT
         COMMON   /NuAnu/NuAnu                                           Switch for neutrino type
         COMMON     /n_DM/n_DM                                           Name of the Earth density model
         COMMON     /n_hi/n_hi                                           Switch for neutrino mass hierarchy
         !COMMON     /n_NT/n_NT                                           Switch for neutrino type
         COMMON     /n_fl/n_fl                                           Switch fot lepton flavor
         !COMMON /n_FF_QES/n_FF_QES                                       Switch for model of nucleon form factors in QES reactions
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

! Read arguments:
         IF (IARGC().LT.5) THEN
           WRITE(*,*) 'ERROR: Missing arguments!'
           WRITE(*,*) 'Usage: ./qesemat NuAnu[1,2] N_Fl[1,2,3] CorV[1,2]
     # MA Target[1,10,23...]'
           STOP
      endIF
         CALL GETARG(1,arg)
         READ(arg,*) NuAnu
         CALL GETARG(2,arg)
         READ(arg,*) N_Fl
         CALL GETARG(3,arg)
         READ(arg,*) N_CorV
         CALL GETARG(4,arg)
         READ(arg,*) MA_cen
         WRITE(MAn,'(F4.2)') MA_cen
         CALL GETARG(5,arg)
         READ(arg,*) N_TT
         WRITE(*,*) 'Set to ',NTn(NuAnu),' ',fln(N_Fl),' ',CorV(N_CorV),
     #' MA_cen=',MAn,' (',N_TT,'. events)'
         
         CALL GeMSet(fui,one,Xlow,Xupp,RelErr,MinCal,*99)

         n_l    = 3
!         DO n_DM=1,4
         n_DM    = 4
         SELECTCASE(n_DM)
               CASE(0)
                     WRITE(*,*) ' with no oscillation '
               CASE(1)
                     WRITE(*,*) ' vacuum case '
               CASE(2)
                     WRITE(*,*) ' 2LEM '
               CASE(3)
                     WRITE(*,*) ' 7LEM '
               CASE(4)
                     WRITE(*,*) ' PREM '
      endSELECT
         DO n_hi=1,2
!         n_hi    = 1
         SELECTCASE(n_hi)
               CASE(1)
                     WRITE(*,*) ' normal hierarchy '
               CASE(2)
                     WRITE(*,*) ' inverse hierarchy '
      endSELECT
!         DO N=1,2
         N       = 1
         set=dFANom_dE(one)

         n_FF_QES= 8                                                     (Bodek,Avvakumov,Bradford&Budd form factor)

         n_b=2
         
         buSM=dsQESCC_dQ2_SM_init(n_FF_QES,1,2)
!          CALL dsQESCC_PRINT_ALL()
         bufN=dsQESCC_dQ2_fN_init()
         buFP=dsQESCC_dQ2_FP_set(zero,zero,MA_QES)

         P_lep_ini=P_lep_min
         P_lep_fin=P_lep_max
         lgP_lep_ini= log10(P_lep_ini)
         lgP_lep_fin= log10(P_lep_fin)
         steplgP_lep= (lgP_lep_fin-lgP_lep_ini)/(NP_lep-1)

         IF (NuAnu.EQ.1) THEN
             n_NT= 1
                          ELSE
             n_NT=-1
        endIF

         CALL setEds
         X=dsQESCC_PRINT(n_TT)
         namfof=Out//'QESnewP/'//DM(n_DM)//hin(n_hi)//
     #               fln(N_Fl)//NTn(NuAnu)//MAn//'_'//nl(n_l)//nb(n_b)//
     #                                                 CorV(N_CorV)//ext
         WRITE(*,*)"Output to file ",namfof
         OPEN(Nfilof,FILE=namfof)

         IF (N_ForB.EQ.1) THEN
           factor=factorf
                          ELSE
           factor=factorb
      endIF
      
      
         DO n_NP_lep=1,NP_lep
           P_lep= 10**(lgP_lep_ini+(n_NP_lep-1)*steplgP_lep)
           ValP(n_NP_lep)=P_lep

           E_lep= sqrt(P_lep**2+mm_lep)
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
           x_ini =P_lep/E_nu_fin
           x_fin =P_lep/E_nu_ini
           deltax=x_fin-x_ini                                           !-

           CALL GeMInt(fui,Res,Xlow,Xupp,*100)
           Jacobianc   = 2*m_ini*P_lep/E_lep
           Intf(n_NP_lep)=factor*deltax*Res*Jacobianc
           WRITE(Nfilof,102) P_lep,Intf(n_NP_lep)
      endDO
         CLOSE(Nfilof)
         CALL GeMInf

      
!      endDO
      endDO
!      endDO
         
         STOP 'THE END OF PROGRAM qesemat'

   99    STOP 'ERROR WITH GeMSet. PROGRAM qesemat'
  100    STOP 'ERROR WITH GeMInt. PROGRAM qesemat'

  101 FORMAT(A,$)
  102 FORMAT(2(1PE16.8))
  103 FORMAT(I4,$)
      END PROGRAM qesemat
