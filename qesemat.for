************************************************************************
      PROGRAM qesemat
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants
         USE InpOutUnits

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

      LOGICAL :: Flux_init,Flux_open_file,Flux_set_sbit,
     #  Flux_read_hdr,Flux_read_table,Flux_print_table,
     #       Flux_calc_spline,Flux_close_file

         logical buf,hdr
         CHARACTER*80
     #                arg
        CHARACTER*80
     #             formula
         CHARACTER*7 MatName
         CHARACTER*4
     #                MAn
        CHARACTER*3 tmp_string
        
              INTEGER,PARAMETER::
     #                Nel=10,
     #                Nfilof   = 100,
     #                Nfilobn  = 101,
     #                Nfiloba  = 102,
     #                Nfiloe   = 103,
     #                NP_lep   = 100,
     #                MinCal   = 100,
     #                NPtotal  = Nel*NP_lep
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
     #                ValP(NP_lep),m_frac(Nel)/Nel*0/,
     #                Intel(Nel,NP_lep)/NPtotal*0/,
     #                     R(NP_lep)/NP_lep*0/
              INTEGER
     #                nm_TT(Nel)/Nel*0/
         CHARACTER*2 name_TT(Nel)
             character(50) FileName
             
         COMMON     /n_MA/n_MA                                           Switch for MA_QES
         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON   /N_CorV/N_CorV
         COMMON   /N_TT/N_TT
         COMMON   /NuAnu/NuAnu                                           Switch for neutrino type
         COMMON     /n_DM/n_DM                                           Name of the Earth density model
         COMMON     /n_hi/n_hi                                           Switch for neutrino mass hierarchy
         COMMON     /n_NT/n_NT                                           Switch for neutrino type
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
         
         
         FileName="2lmnH11xe.dat"

! Read arguments:
         IF (IARGC().LT.5) THEN
           WRITE(*,*) 'ERROR: Missing arguments!'
           WRITE(*,*) 'Usage: ./qesemat NuAnu[1,2] N_Fl[1,2,3] CorV[1,2]
     #   MA "mixture"  
     #   "formula[element1 fraction1 element2 fraction2...]"
     #    e.g. ./qesemat 2 3 1 1.0 "water" "H 2 O 2"'
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
         CALL GETARG(5,MatName)
         CALL GETARG(6,arg)
         !READ(arg,*) N_TT
         READ(arg,'(A80)') formula
********** number-based formula ****************************
!         READ(formula,*,END=990)(nm_TT(n_el),m_frac(n_el),n_el=1,Nel)
********** name-based formula ****************************
        READ(formula,*,END=990)
     #        (name_TT(n_el),m_frac(n_el),n_el=1,Nel)
990      WRITE(*,*) 'Set to ',NTn(NuAnu),' ',fln(N_Fl),' ',CorV(N_CorV),
     #   ' MA_cen=',MAn,', formula: ',formula!', N_TT=',N_TT
!*         convert names to numbers:
         
       NelLast=1
         DO n_el=1,Nel
           nm_TT(n_el)=GET_TGT_NUMBER(name_TT(n_el))
           WRITE(*,*)"element=",name_TT(n_el),nm_TT(n_el),m_frac(n_el)
           if(m_frac(n_el).gt.0)NelLast=n_el
         endDO

c~ * write the chemical formula!
c~          DO n_el=1,Nel
c~         
c~            frac=m_frac(n_el)
c~            IF(frac.GT.0)THEN
c~              X=GET_TGT_NAME(nm_TT(n_el),tmp_string)
c~              WRITE(*,'(A2$)')adjustr(tmp_string)
c~              IF(frac.NE.1)THEN 
c~                WRITE(tmp_string,'(F3.1)')frac
c~                WRITE(*,'(A2$)')adjustl(tmp_string)
c~              endIf
c~            endIf
c~          endDO
c~          WRITE(*,*)"."
c~ *********** done *****************************    
         mu=0 
         DO n_el=1,Nel
           IF(m_frac(n_el).GT.0)
     #         mu=mu+GET_TGT_A(nm_TT(n_el))*m_frac(n_el)                 molar mass, numerically equals to atomic mass, num.app.eq. nucleon number
         endDO
         fact=factor/mu
    
         CALL GeMSet(fui,one,Xlow,Xupp,RelErr,MinCal,*99)

         n_l    = 3
         
         buf=Flux_init()
         buf=Flux_open_file(FileName)
         hdr=.true.
         do while(Flux_read_hdr())
           buf=Flux_read_table()
           buf=Flux_print_table()
           buf=Flux_calc_spline()
         end do
         buf=Flux_close_file()

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
      
      DO n_el=1,Nel
         IF(m_frac(n_el).GT.0)THEN
         n_TT=nm_TT(n_el)
         X=dsQESCC_PRINT(n_TT)

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
           Intel(n_el,n_NP_lep)=m_frac(n_el)*deltax*Res*Jacobianc
           !WRITE(Nfilof,102) P_lep,Intel(n_el,n_NP_lep)
           WRITE(*,*)n_NP_lep,"/",NP_lep,"E_lep=",
     #      E_lep,Intel(n_el,n_NP_lep)
      endDO
         R=R+Intel(n_el,:)
      endIF
      endDO
        R=fact*R
        Intel=fact*Intel
         namfof=Out//'QESnewP/'//TRIM(MatName)//
     #               fln(N_Fl)//NTn(NuAnu)//MAn//'_'//nl(n_l)//nb(n_b)//
     #                                                 CorV(N_CorV)//ext
         WRITE(*,*)"Output to file ",namfof
         OPEN(Nfilof,FILE=namfof)
         WRITE(Nfilof,'(11A25)')adjustl("#Plep, GeV"),
     #    (adjustl(name_TT(n_el)),n_el=1,NelLast),adjustl(formula)
          DO n_NP_lep=1,NP_lep
          WRITE(*,*)n_NP_lep,"/",NP_lep,"E_lep=",
     #      ValP(n_NP_lep),R(n_NP_lep)
           WRITE(Nfilof,*) ValP(n_NP_lep),R(n_NP_lep),
     #     (Intel(n_el,n_NP_lep),n_el=1,NelLast)
           WRITE(Nfilof,*)
      endDO
         CLOSE(Nfilof)
         CALL GeMInf
         
         STOP 'THE END OF PROGRAM qesemat'

   99    STOP 'ERROR WITH GeMSet. PROGRAM qesemat'
  100    STOP 'ERROR WITH GeMInt. PROGRAM qesemat'

  101 FORMAT(A,$)
  102 FORMAT(1PE16.8,$)
  103 FORMAT(I4,$)
      END PROGRAM qesemat
