************************************************************************
      FUNCTION Flux_Init()
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2014/01/10  *
************************************************************************
*                                                                      *
*     Neutrino fluxes dF/dE [1/(cm^2 s sr GeV)]                        *
************************************************************************
       implicit none

                 logical:: Flux_Init               
                 logical:: Flux_Open_File
                 logical:: Flux_Read_Hdr
                 logical:: Flux_Read_Table
                 logical:: Flux_Close_File
                 logical:: Flux_Calc_Spline
                 logical:: Flux_Print_Table
                 logical:: Flux_Has_Table
                 real:: Flux_get_dF
                 real:: Flux_Get_Emin,Flux_Get_Emax
                 real:: Flux_Get_Zmin,Flux_Get_Zmax
                 logical:: Flux_Get_Last_Nu
                 integer:: ioer,ne_cur,Issue
                 
                 real Energy
                 
                 real:: Sp1
                 
      CHARACTER(*),PARAMETER::
     #     FName="<<<FluxReader>>>:"
              INTEGER,PARAMETER::
     #                str  = 300,                                         !Tsiferka dlya failika, chtoby chitat
     #                MaxNE=1000,
     #                MaxNcos=1000,
     #                NFlv=3,
     #                NNuType=2,
     #                NNeutrinos=NFlv*NNuType,
     #                Size_dF=NNeutrinos*MaxNE,
     #                Size_d2F=NNeutrinos*MaxNE*MaxNcos
        REAL,PARAMETER::
     #  Eps=1d-12
      SAVE      
                 logical htable(NNuType,NFlv)/NNeutrinos*.FALSE./
                 logical LogE,LogE_table(NNuType,NFlv)
              INTEGER I,
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
         REAL SpLim(4),SpVar(2)
         CHARACTER(80) filename
         CHARACTER(*)  file_name                                           !Imya faila. Dolzhno by chitat'sya glavnoy programoy
         CHARACTER(100) lline
         CHARACTER(2) Anutype
!     initialization:
         WRITE(*,*)FName,"Init flux routine"
         Flux_Init=.TRUE.
         RETURN
!*************************************************************************
         ENTRY Flux_Get_Emin(NuAnu,Flavor)
           Flux_Get_Emin=E_min(NuAnu,Flavor)
         RETURN
!*************************************************************************
         ENTRY Flux_Get_Emax(NuAnu,Flavor)
           Flux_Get_Emax=E_max(NuAnu,Flavor)
         RETURN
!*************************************************************************
        ENTRY Flux_Get_Zmin(NuAnu,Flavor)
           Flux_Get_Zmin=Z_min(NuAnu,Flavor)
         RETURN
!*************************************************************************
         ENTRY Flux_Get_Zmax(NuAnu,Flavor)
           Flux_Get_Zmax=Z_max(NuAnu,Flavor)
         RETURN
!*************************************************************************
         ENTRY Flux_Get_Last_Nu(NuAnu,Flavor)
         write(*,*)n_NuAnu,n_Fl
           NuAnu=n_NuAnu; Flavor=n_Fl ! output 
         RETURN
!*************************************************************************
      ENTRY Flux_Open_File(file_name)
!************ open & read file with dF/dE table **************************
         filename=trim(file_name)
         WRITE(*,*) FName,'Open file ',filename
         LogE=.TRUE.
         Flux_Open_File=.FALSE.
         OPEN(str,STATUS='OLD',FILE=filename,ERR=2001)
         Flux_Open_File=.TRUE.
         RETURN
!*************************************************************************
      ENTRY Flux_Read_Hdr()
!************ read file with dF/dE table *********************************
         WRITE(*,*) FName,'Read table header'
         !************** read header ***************
         do 
             read(str,*,iostat=ioer),lline
             write(*,*),"line=",lline
             if(ioer.ne.0)then
                 Flux_Read_Hdr=.FALSE.
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
                         Flux_Read_Hdr=.FALSE.
                         goto 2002
                 end select
                 selectcase(Anutype(2:2))
                     case("e");n_Fl=1;
                     case("m");n_Fl=2;
                     case("t");n_Fl=3;
                     case default
                         Flux_Read_Hdr=.FALSE.
                         goto 2002
                 end select
             elseif(lline(1:4).eq."#NE=")then
                 read(lline,'(4x,I4)')NE(n_NuAnu,n_Fl)
                 write(*,'("#NE=",I4)')NE(n_NuAnu,n_Fl)
             elseif(lline(1:6).eq."#Ncos=")then
                read(lline,'(6x,I4)')Ncos(n_NuAnu,n_Fl)
                write(*,'("#Ncos=",I4)')Ncos(n_NuAnu,n_Fl)         
             elseif(lline(1:5).eq.'#RegE')then
                 write(*,*)"#RegE"
                 LogE=.FALSE.
             elseif(lline(1:5).eq.'#LogE')then
                 write(*,*)"#LogE"
                 LogE=.TRUE.
             else
                write(*,*)'Unknown line ',lline
             endif
         end do
         
         write(*,*)Fname,"Read header complete"
         htable(n_NuAnu,n_Fl)=.TRUE.
         LogE_table(n_NuAnu,n_Fl)=LogE
         Flux_Read_Hdr=.TRUE.
         RETURN
         
!*************************************************************************
      ENTRY Flux_Read_Table()
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
      Flux_Read_Table=.TRUE.
      RETURN
!*************************************************************************
      ENTRY Flux_Close_File()
!************ read file with dF/dE table *********************************     
        write(*,*)Fname,"Close file"
        CLOSE(str)
        RETURN
!*************************************************************************
      ENTRY Flux_Has_Table(NuAnu,Flavor)
!************ read file with dF/dE table *********************************     
        Flux_Has_Table=htable(NuAnu,Flavor)
        RETURN
!*************************************************************************
      ENTRY Flux_Print_Table(NuAnu,Flavor)
!************ read file with dF/dE table *********************************     
        write(*,*)"#Neutrino type = [",NuAnu,Flavor,"]"
        write(*,*)"#Table size ",NE(NuAnu,Flavor),"x",Ncos(NuAnu,Flavor)
        if(LogE_table(NuAnu,Flavor))then 
            write(*,*)"#E step is LOG"
        else 
            write(*,*)"#E step is REGULAR"
        end if
        write(*,'(2E12.5)')(E(NuAnu,Flavor,I),dF(NuAnu,Flavor,I),
     #  I=1,NE(NuAnu,Flavor))
        RETURN
!*************************************************************************
      ENTRY Flux_Calc_Spline(NuAnu,Flavor)
!************ read file with dF/dE table *********************************
        write(*,*)Fname,"Calc spline"
        
        Issue=Flavor+NFlv*(NuAnu-1)
        !set spline limits
        SpLim(1)=E_min(NuAnu,Flavor)
        SpLim(2)=E_max(NuAnu,Flavor)
        if(LogE_table(NuAnu,Flavor))then
          SpLim(1)=log10(SpLim(1))
          SpLim(2)=log10(SpLim(2))
        endif
        
        CALL Coeff1(0,1,.TRUE.,Eps,Issue,NE(NuAnu,Flavor),
     #  SpLim(1),SpLim(2),
     #  dF(NuAnu,Flavor,:),CdF(NuAnu,Flavor,:),.TRUE.,1)
        Flux_Calc_Spline=.TRUE.
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
          SpVar(1)=Energy
          if(LogE_table(NuAnu,Flavor))SpVar(1)=log10(SpVar(1))
          Flux_get_dF=Sp1(Issue,CdF(NuAnu,Flavor,1:ne_cur),
     #        SpVar(1))
          if(Flux_get_dF.le.0)Flux_get_dF=0
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

      END FUNCTION Flux_Init
