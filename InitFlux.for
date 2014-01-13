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
                 logical Flux_set_sbit
                 logical Flux_print_table
                 
                 logical Stop_bits(5)/4*.FALSE.,.TRUE./
                 logical Sbit
                 integer ioer, Nbit
                 real lgE_max,lgE_min
                 
      CHARACTER(*),PARAMETER::
     #     FName="<<<FluxReader>>>:"
              INTEGER,PARAMETER::
     #                str  = 300,                                         Tsiferka dlya failika, chtoby chitat
     #                MaxNE=1000,
     #                MaxNcos=1000,
     #                NFlv=3
      SAVE
              INTEGER length,I,
     #                n_NE,NE(3),
     #                n_NCos,NCos(3),
     #                n_NuAnu,n_Fl

                 REAL
     #                 E(NFlv,MaxNE),Z(NFlv,MaxNcos),
     #                 dF(NFlv,MaxNE),CdF(NFlv,MaxNE+2),
     #                 d2F(NFlv,MaxNE,MaxNcos),
     #                 Cd2F(NFlv,MaxNE+2,MaxNcos+2)

         CHARACTER(20) filename
         CHARACTER(*)  file_name                                           Imya faila. Dolzhno by chitat'sya glavnoy programoy
         CHARACTER(100) lline
         CHARACTER(2) Anutype
!     initialization:
         WRITE(*,*)FName,"Init flux routine"
         Flux_init=.TRUE.
         RETURN

!*************************************************************************
      ENTRY Flux_set_sbit(Nbit,Sbit)
!************ open & read file with dF/dE table **************************
         Stop_Bits(Nbit)=Sbit
         return

!*************************************************************************
      ENTRY Flux_open_file(file_name)
!************ open & read file with dF/dE table **************************
         filename=trim(file_name)
         WRITE(*,*) FName,'Open file ',filename
         !length = len_trim(filename) - 1
         Flux_open_file=.FALSE.
         OPEN(str,STATUS='OLD',FILE=filename,ERR=2001)
         Flux_open_file=.TRUE.
         if(Stop_Bits(1))RETURN
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
                 selectcase(Anutype(2:2))
                     case("e");n_Fl=1;
                     case("m");n_Fl=2;
                     case("t");n_Fl=3;
                     case default
                         Flux_read_hdr=.FALSE.
                         goto 2002
                 end select
             elseif(lline(1:4).eq."#NE=")then
                 read(lline,'(4x,I4)')NE(n_Fl)
                 write(*,'("#NE=",I4)')NE(n_Fl)
             elseif(lline(1:6).eq."#Ncos=")then
                read(lline,'(6x,I4)')Ncos(n_Fl)
                write(*,'("#Ncos=",I4)')Ncos(n_Fl)         
             else
                write(*,*)'Unknown line ',lline
             endif
         end do
         
         write(*,*)Fname,"Read header complete"
         Flux_read_hdr=.TRUE.
         if(Stop_Bits(2))RETURN
         
!*************************************************************************
      ENTRY Flux_read_table()
!************ read file with dF/dE table *********************************
        write(*,*),Fname,"Read flux table from file"
!         RETURN
         DO n_NE=1,NE(n_Fl)
            READ(str,*)E(n_FL,n_NE),dF(n_Fl,n_NE)
      endDO
      Flux_read_table=.TRUE.
      if(Stop_Bits(3))RETURN
!*************************************************************************
      ENTRY Flux_calc_spline()
!************ read file with dF/dE table *********************************
        write(*,*)Fname,"Calc spline"
        lgE_min=log10(E(n_Fl,1))
        lgE_max=log10(E(n_Fl,NE(n_Fl)))
         
        CALL Coeff1(0,1,.TRUE.,1.0d-12,1,NE(n_Fl),lgE_min,lgE_max,
     #                       dF(n_Fl,:),CdF(n_Fl,:),.TRUE.,1)
        
        if(Stop_Bits(4))RETURN
!*************************************************************************
      ENTRY Flux_close_file()
!************ read file with dF/dE table *********************************     
        write(*,*)Fname,"Close file"
        CLOSE(str)
        if(Stop_Bits(5))RETURN

!*************************************************************************
      ENTRY Flux_print_table()
!************ read file with dF/dE table *********************************     
        write(*,*)"#Neutrino type = ",n_Fl
        write(*,*)"#Table size ",NE(n_Fl),"x",Ncos(n_Fl)
        write(*,'(2E12.5)')(E(n_Fl,I),dF(n_Fl,I),I=1,NE(n_Fl))
        RETURN
!*************************************************************************

2001    write(*,*)FName,'ERROR: ','cannot open file ',filename
        RETURN
2002    write(*,*)FName,'ERROR: ','Unknown neutrino flavour=[',
     #                   Anutype(2:2),'] !=(e,m,t)'
        RETURN
2003    write(*,*)FName,'INFO: ','EOF reached'
        RETURN

      END FUNCTION Flux_init
