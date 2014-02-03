C
C File:   test_fluxtable.for
C Author: ash
C
C Created on 20 Январь 2014 г., 22:50
C
      LOGICAL :: Flux_init,Flux_open_file
      LOGICAL :: Flux_read_head,Flux_read_table,Flux_print_table
      LOGICAL :: Flux_close_file,Flux_calc_spline,Flux_get_last_nu
      REAL :: Flux_get_dF
      logical res
      integer,parameter:: Npts=100
      real p
      real Emin1,Emax1
      character(80) arg
      real flx
      character(80) FluxFile
      FluxFile="test.flux"
      
      IF (IARGC().GT.0) THEN
        ! read FLUX filename from argument
        CALL GETARG(1,arg)
        FluxFile=arg
      END IF
      res=Flux_init()
      res=Flux_open_file(FluxFile)
      res=Flux_read_head()
      res=Flux_read_table()
      res=Flux_close_file()
      res=Flux_get_last_nu(nanu,nf)
      write(*,*)nanu,nf
      ! **** find spectrum limits ****
      res=Flux_calc_spline(nanu,nf)
      res=Flux_print_table(nanu,nf)
      Emin1=Flux_Get_Emin(nanu,nf)
      Emax1=Flux_Get_Emax(nanu,nf)
      
      !WRITE(*,*)"Spectrum Emin=[",Emin,"] => Emin1=",Emin1
      !WRITE(*,*)"Spectrum Emax=[",Emax,"] => Emax1=",Emax1
      p=Emin1
      dp=(Emax1-Emin1)/Npts
       write(*,'(A16," |",A16)')"P, GeV","flux"
      write(*,'(34A1)')("=",i=1,34)
      do
          p=p+dp
          if(p.GE.Emax1)EXIT
          flx=Flux_get_dF(nanu,nf,p)
          write(*,'(e16.8,"  ",e16.8)'),p,flx
          write(901,'(e16.8,"  ",e16.8)'),p,flx
      end do
      END



