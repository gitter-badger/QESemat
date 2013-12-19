************************************************************************
      FUNCTION dFANom_dE_init()
************************************************************************
*                                                                      *
*     This FUNCTION returns the values of atmospheric electron and     *
*     muon (anti)neutrino fluxes dF/dE [1/(cm^2 s sr GeV)] modulated   *
*     neutrino oscillations and averaged over the cosine of            *
*     zenith angle for any E_\nu from range 50 MeV to 10^9 GeV.        *
*                                                                      *
*            Data neutrino energy range is 50 MeV to 1 EeV:            *
*                                                                      *
*    0     1   10    10^2  10^3  10^4  10^5  10^6  10^7  10^8  10^9    *
*    |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|     *
*     |----|-----|---|-------|-----|-----|-----|-----|-----------|     *
*    0.05  1   10   70     10^3  10^4  10^5  10^6  10^7        10^9    *
*                            |--------------ISU_HE---------------|     *
*     |---------------------------CORT---------------------------|     *
*      |----------Honda11----------|                                   *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: one
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE
         
           LOGICAL(2),PARAMETER::
     #                Tie  =.TRUE.,
     #                Quiz =.TRUE.
              INTEGER,PARAMETER::
     #                Mult =   0,
     #                Mode =   1,
     #                NP   = 100,
     #                L    =   1
                 REAL,PARAMETER::
     #                P_min = 1.0d-01,
     #                P_max = 1.0d+04,
     #                Eps   = 1.0d-12
         CHARACTER(*),PARAMETER::
     #                ext='.dat'
        
        CHARACTER*59 Fn_cur
        CHARACTER*3
     #                Sp,
     #                DM(0:4)/'wno','vac','2lm','7lm','mat'/
        CHARACTER*1
     #                SA,
     #                nl(5)/'0','n','m','x','1'/,
     #                fln(3)/'e','m','t'/,
     #                NTn(2)/'n','a'/,
     #                hin(2)/'n','i'/


         DIMENSION FF(3,2,NP),
     #             CC(3,2,NP+2),
     #             F_cur(NP),
     #             C_cur(NP+2)

        INTEGER Iss(3,2)/1,2,3,4,5,6/
         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON     /n_DM/n_DM                                           Name of the Earth density model
         COMMON     /n_hi/n_hi                                           Switch for neutrino mass hierarchy
         COMMON      /n_l/n_l

         SELECTCASE(N)
               CASE(1)
                     WRITE(*,*) ' AN_Honda11 + AN_ISU_HE + AN_SHE '
                     Sp='H11'
                     WRITE(*,*) ' Maximal solar activity '
                     SA='x'
               CASE(2)
                     WRITE(*,*) ' CORTout '
                     Sp='CRT'
                     WRITE(*,*) ' Minimal solar activity '
                     SA='n'
      endSELECT
      
        lgP_min = log10(P_min)
        lgP_max = log10(P_max)
        DO n_Fl=1,3
          DO n_NuAnu=1,2
            Fn_cur=Out//'spectra/'//DM(n_DM)//hin(n_hi)//Sp//SA//
     #             fln(n_Fl)//NTn(n_NuAnu)//'_'//nl(n_l)//ext
*            read data from file to ->F_cur
            Nf_cur=301
            OPEN(Nf_cur,FILE=Fn_cur)
            WRITE(*,*)"Opened file ",Fn_cur
            DO n_NP=1,NP
                READ(Nf_cur,101) E,FF(n_Fl,n_NuAnu,n_NP)
            endDO
            CLOSE(Nf_cur)
*            calculate coefficients ->C_cur
            CALL Coeff1(Mult,Mode,Tie,Eps,Iss(n_Fl,n_NuAnu),
     #               NP,lgP_min,lgP_max,
     #               FF(n_Fl,n_NuAnu,:),CC(n_Fl,n_NuAnu,:),Quiz,L)
          endDO
          
        endDO
      
         !FIXME - pass F_cur to F() and C_cur to C()
        
        dFANom_dE_init=one
        RETURN
*     ==================================================================
      ENTRY dFANom_dE(nn_Fl,nn_NuAnu,P)
*     ==================================================================
        dFANom_dE=Sp1(Iss(nn_Fl,nn_NuAnu),CC(nn_Fl,nn_NuAnu,:),log10(P))
        RETURN
         
  101 FORMAT(2(1PE16.8))
      END FUNCTION dFANom_dE_init
