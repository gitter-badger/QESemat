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
      character(50) arg
           
      real flx(3)
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
      res=Flux_set_sbit(4,.TRUE.)
      res=Flux_open_file(FileName)
      do while(res)
          res=Flux_read_hdr()
!          if(res)tmp=Flux_print_table()
      end do
      res=Flux_close_file()
      write(*,'(A16," |",3A16)')"P, GeV",
     #    "NuE flux","NuMU flux","NuTau flux"
      write(*,'(66A1)')("=",i=1,66)
      do i=1,Npts
          p=2e-1+i*0.01
          do nf=1,3
            flx(nf)=Flux_get_dF(nf,p)
          end do
          write(*,'(e16.8," |",3e16.8)'),p,flx
      end do
!      flx(1)=Flux_get_dF(1,p)
!      flx(2)=Flux_get_dF(2,p)
!      flx(3)=Flux_get_dF(3,p)
!      write(*,'("Flux E(p=",e12.3," GeV)=",e12.6)')p,flx(1)
!      write(*,'("Flux M(p=",e12.3," GeV)=",e12.6)')p,flx(2)
!      write(*,'("Flux T(p=",e12.3," GeV)=",e12.6)')p,flx(3)
      END

