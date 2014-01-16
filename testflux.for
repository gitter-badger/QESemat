C
C File:   testflux.for
C Author: ash
C
C Created on 13 Январь 2014 г., 14:20
C
      LOGICAL :: Flux_init,Flux_open_file,Flux_set_sbit
      LOGICAL :: Flux_read_hdr,Flux_print_table
      LOGICAL :: Flux_close_file
      REAL :: Flux_get_dF
      logical res,tmp
      integer,parameter:: Npts=100
      real p
      real Emin(2,3), Emax(2,3), Emin1,Emax1
      character(50) arg
           
      real flx(2,3)
      character(50) FileName
      FileName="test.flux"
      
      IF (IARGC().GT.0) THEN
        ! read filename from argument
        CALL GETARG(1,arg)
        FileName=arg
      END IF
      
      
      PRINT*, 'Hello World'
      res=Flux_init()
      res=Flux_set_sbit(1,.TRUE.)
      res=Flux_set_sbit(2,.FALSE.)
      res=Flux_set_sbit(3,.FALSE.)
      res=Flux_set_sbit(4,.TRUE.)
      res=Flux_open_file(FileName)
      do while(res)
          res=Flux_read_hdr()
          if(res)tmp=Flux_print_table()
      end do
      res=Flux_close_file()
      write(*,'(A16," |",3A16," |",3A16)')"P, GeV",
     #    "NuE flux","NuMU flux","NuTau flux",
     #    "AntiNuE flux","AntiNuMU flux","AntiNuTau flux"
      write(*,'(116A1)')("=",i=1,116)
      ! **** find spectrum limits ****
      
      do nf=1,3
              do nanu=1,2
                Emin(nanu,nf)=Flux_GetEmin(nanu,nf)
                Emax(nanu,nf)=Flux_GetEmax(nanu,nf)
              end do
      end do
      Emin1=minval(Emin)
      Emax1=maxval(Emax)
!      WRITE(*,*)"Spectrum Emin=[",Emin,"] => Emin1=",Emin1
 !     WRITE(*,*)"Spectrum Emax=[",Emax,"] => Emax1=",Emax1
      p=Emin1
      dp=(Emax1-Emin1)/20
      do 
          p=p+dp
          if(p.GT.Emax1)EXIT
          do nf=1,3
              do nanu=1,2
                flx(nanu,nf)=Flux_get_dF(nanu,nf,p)
              end do
          end do
          write(*,'(e16.8," |",3e16.8," |",3e16.8)'),p,flx(1,:),flx(2,:)
         
      end do
!      flx(1)=Flux_get_dF(1,p)
!      flx(2)=Flux_get_dF(2,p)
!      flx(3)=Flux_get_dF(3,p)
!      write(*,'("Flux E(p=",e12.3," GeV)=",e12.6)')p,flx(1)
!      write(*,'("Flux M(p=",e12.3," GeV)=",e12.6)')p,flx(2)
!      write(*,'("Flux T(p=",e12.3," GeV)=",e12.6)')p,flx(3)
      END

