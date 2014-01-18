************************************************************************
      FUNCTION Flux_init()
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2014/01/10  *
************************************************************************
*                                                                      *
*     Neutrino fluxes dF/dE [1/(cm^2 s sr GeV)]                        *
************************************************************************
       implicit none

                 logical Flux_init               
                 logical Flux_open_file
                 logical Flux_read_hdr
                 logical Flux_read_table
                 logical Flux_close_file
                 logical Flux_calc_spline
                 logical Flux_print_table
                 logical Flux_has_table
                 real Flux_get_dF
                 real Flux_GetEmin,Flux_GetEmax
                 real Flux_GetZmin,Flux_GetZmax
                 
                 integer ioer,ne_cur,Issue
                 real Energy
                 
                 real Sp1
                 
      CHARACTER(*),PARAMETER::
     #     FName="<<<FluxReader>>>:"
              INTEGER,PARAMETER::
     #                str  = 300,                                         Tsiferka dlya failika, chtoby chitat
     #                MaxNE=1000,
     #                MaxNcos=1000,
     #                NFlv=3,
     #                NNuType=2,
     #                NNeutrinos=NFlv*NNuType,
     #                Size_dF=NNeutrinos*MaxNE,
     #                Size_d2F=NNeutrinos*MaxNE*MaxNcos
      SAVE      
                 logical htable(NNuType,NFlv)/NNeutrinos*.FALSE./
              INTEGER length,I,
     #                n_NE,NE(NNuType,NFlv)/NNeutrinos*0/,
     #                n_NCos,NCos(NNuType,NFlv)/NNeutrinos*0/,
     #                n_NuAnu,n_Fl,
     #                NuAnu,Flavor

                 REAL
     #                 E(NNuType,NFlv,MaxNE),
     #                 E_max(NNuType,NFlv),E_min(NNuType,NFlv),
     #                 Z(NNuType,NFlv,MaxNcos),
     #                 Z_max(NNuType,NFlv),Z_min(NNuType,NFlv),
     #                 dF(NNuType,NFlv,MaxNE)/Size_dF*0/,
     #                 CdF(NNuType,NFlv,MaxNE+2),
     #                 d2F(NNuType,NFlv,MaxNE,MaxNcos)/Size_d2F*0/,
     #                 Cd2F(NNuType,NFlv,MaxNE+2,MaxNcos+2)

         CHARACTER(20) filename
         CHARACTER(*)  file_name                                           Imya faila. Dolzhno by chitat'sya glavnoy programoy
         CHARACTER(100) lline
         CHARACTER(2) Anutype
!     initialization:
         WRITE(*,*)FName,"Init flux routine"
         Flux_init=.TRUE.
         RETURN
!*************************************************************************
         ENTRY Flux_GetEmin(NuAnu,Flavor)
           Flux_GetEmin=E_min(NuAnu,Flavor)
         RETURN
!*************************************************************************
         ENTRY Flux_GetEmax(NuAnu,Flavor)
           Flux_GetEmax=E_max(NuAnu,Flavor)
         RETURN
!*************************************************************************
        ENTRY Flux_GetZmin(NuAnu,Flavor)
           Flux_GetZmin=Z_min(NuAnu,Flavor)
         RETURN
!*************************************************************************
         ENTRY Flux_GetZmax(NuAnu,Flavor)
           Flux_GetZmax=Z_max(NuAnu,Flavor)
         RETURN
!*************************************************************************
      ENTRY Flux_open_file(file_name)
!************ open & read file with dF/dE table **************************
         filename=trim(file_name)
         WRITE(*,*) FName,'Open file ',filename
         !length = len_trim(filename) - 1
         Flux_open_file=.FALSE.
         OPEN(str,STATUS='OLD',FILE=filename,ERR=2001)
         Flux_open_file=.TRUE.
         RETURN
!*************************************************************************
      ENTRY Flux_read_hdr()
!************ read file with dF/dE table *********************************
         WRITE(*,*) FName,'Read table header'
         !************** read header ***************
         do 
             read(str,*,iostat=ioer),lline
!             write(*,*),"line=",lline
             if(ioer.ne.0)then
                 Flux_read_hdr=.FALSE.
                 goto 2003
             end if
             if(lline(1:1).ne."#")exit
             
             if(lline(1:7).eq.'#Ntype=')then
                 read(lline,'(7x,A2)')Anutype
                 write(*,'("#Ntype=",A2)')Anutype
                 selectcase(Anutype(1:1))
                     case("n");n_NuAnu=1;
                     case("a");n_NuAnu=2;
                     case default
                         Flux_read_hdr=.FALSE.
                         goto 2002
                 end select
                 selectcase(Anutype(2:2))
                     case("e");n_Fl=1;
                     case("m");n_Fl=2;
                     case("t");n_Fl=3;
                     case default
                         Flux_read_hdr=.FALSE.
                         goto 2002
                 end select
             elseif(lline(1:4).eq."#NE=")then
                 read(lline,'(4x,I4)')NE(n_NuAnu,n_Fl)
                 write(*,'("#NE=",I4)')NE(n_NuAnu,n_Fl)
             elseif(lline(1:6).eq."#Ncos=")then
                read(lline,'(6x,I4)')Ncos(n_NuAnu,n_Fl)
                write(*,'("#Ncos=",I4)')Ncos(n_NuAnu,n_Fl)         
             else
                write(*,*)'Unknown line ',lline
             endif
         end do
         
         write(*,*)Fname,"Read header complete"
         Flux_read_hdr=.TRUE.
         htable(n_NuAnu,n_Fl)=.TRUE.
         RETURN
         
!*************************************************************************
      ENTRY Flux_read_table()
!************ read file with dF/dE table *********************************
        write(*,*),Fname,"Read flux table from file"
!         RETURN
         DO n_NE=1,NE(n_NuAnu,n_Fl)
            READ(str,*,iostat=ioer)
     #            E(n_NuAnu,n_FL,n_NE),dF(n_NuAnu,n_Fl,n_NE)
      endDO
      ! get limits:
      E_min(n_NuAnu,n_Fl)=E(n_NuAnu,n_Fl,1)
      E_max(n_NuAnu,n_Fl)=E(n_NuAnu,n_Fl,NE(n_NuAnu,n_Fl))
      if(E_max(n_NuAnu,n_Fl).lt.E_min(n_NuAnu,n_Fl))then
          E_min(n_NuAnu,n_Fl)=E_max(n_NuAnu,n_Fl)
          E_max(n_NuAnu,n_Fl)=E(n_NuAnu,n_Fl,1)
      end if
      Flux_read_table=.TRUE.
      RETURN
!*************************************************************************
      ENTRY Flux_close_file()
!************ read file with dF/dE table *********************************     
        write(*,*)Fname,"Close file"
        CLOSE(str)
        RETURN
!*************************************************************************
      ENTRY Flux_has_table(NuAnu,Flavor)
!************ read file with dF/dE table *********************************     
        Flux_has_table=htable(NuAnu,Flavor)
        RETURN
!*************************************************************************
      ENTRY Flux_print_table(NuAnu,Flavor)
!************ read file with dF/dE table *********************************     
        write(*,*)"#Neutrino type = [",NuAnu,Flavor,"]"
        write(*,*)"#Table size ",NE(NuAnu,Flavor),"x",Ncos(NuAnu,Flavor)
        write(*,'(2E12.5)')(E(NuAnu,Flavor,I),dF(NuAnu,Flavor,I),
     #  I=1,NE(NuAnu,Flavor))
        RETURN
!*************************************************************************
      ENTRY Flux_calc_spline(NuAnu,Flavor)
!************ read file with dF/dE table *********************************
        write(*,*)Fname,"Calc spline"
        
        Issue=Flavor+NFlv*(NuAnu-1)
        CALL Coeff1(0,1,.TRUE.,1.0d-12,Issue,NE(NuAnu,Flavor),
     #  log10(E_min(NuAnu,Flavor)),log10(E_max(NuAnu,Flavor)),
     #  dF(NuAnu,Flavor,:),CdF(NuAnu,Flavor,:),.TRUE.,1)
        Flux_calc_spline=.TRUE.
        RETURN
!*************************************************************************
      ENTRY Flux_get_dF(NuAnu,Flavor,Energy)
!************ dF/dE spline *********************************    
         Issue=Flavor+NFlv*(NuAnu-1)
        if((NE(NuAnu,Flavor).eq.0).or.
     #   (Energy.lt.E_min(NuAnu,Flavor)).or.
     #   (Energy.gt.E_max(NuAnu,Flavor))) then
            Flux_get_dF=0;
        else
          ne_cur=NE(NuAnu,Flavor)+2
          Flux_get_dF=Sp1(Issue,CdF(NuAnu,Flavor,1:ne_cur),
     #        log10(Energy))
        end if
        RETURN
!*************************************************************************
        

2001    write(*,*)FName,'ERROR: ','cannot open file ',filename
        RETURN
2002    write(*,*)FName,'ERROR: ','Unknown neutrino type=[',
     #                   Anutype(1:2),'] !=[n,a](e,m,t)'
        RETURN
2003    write(*,*)FName,'INFO: ','EOF reached'
        RETURN

      END FUNCTION Flux_init
