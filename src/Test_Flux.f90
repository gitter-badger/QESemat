!**********************************************************************!
PROGRAM Test_Flux
!----------------------------------------------------------------------!
!Test for function 'Flux' and file contained flux                      !am i right?
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

logical:: &
    Flux_Init,Flux_Open_File,Flux_Read_Hdr,Flux_Read_Table,Flux_Close_File,&
    Flux_Print_Table,Flux_Calc_Spline
real:: &
    Flux_Get_Emin,Flux_Get_Emax,Flux_Get_dF

integer,parameter:: &
    Npts=100
logical &
    bufL
integer &
    NuAnu,Flavor,i
real &
    E,dE,Eminmin,Emaxmax,Emin(2,3),Emax(2,3),flx(2,3)
character*80 &
    FluxFile

!file name setting-----------------------------------------------------!
    if(iargc()==0)then
        FluxFile='../input/flux.sng'                                   !Default
    else
        call GetArg(1,FluxFile)                                        !reading
    endif
!flux file reading-----------------------------------------------------!
    bufL=Flux_Init()
    bufL=Flux_Open_File(FluxFile)
    do while(bufL)
        bufL=Flux_Read_Hdr()
        if(bufL)bufL=Flux_Read_Table()
    enddo
    bufL=Flux_Close_File()
!spectrum limit finding------------------------------------------------!
    do NuAnu=1,2
        do Flavor=1,3
            bufL=Flux_Calc_Spline(NuAnu,Flavor)
            bufL=Flux_Print_Table(NuAnu,Flavor)
            Emin(NuAnu,Flavor)=Flux_Get_Emin(NuAnu,Flavor)
            Emax(NuAnu,Flavor)=Flux_Get_Emax(NuAnu,Flavor)
        enddo
    enddo
    Eminmin=minval(Emin)
    Emaxmax=maxval(Emax)
    write(*,'("Spectrum Emin=[",6(1PE9.2),"] => Eminmin=",1PE9.2)') Emin,Eminmin
    write(*,'("Spectrum Emax=[",6(1PE9.2),"] => Emaxmax=",1PE9.2)') Emax,Emaxmax
!summary flux table printing-------------------------------------------!
    E=Eminmin
    dE=(Emaxmax-Eminmin)/Npts
    write(*,'(A16,1X,"||",3(A13,3X),1X,"|",3(A15,1X))')'E [GeV]',&
    '\nu_e','\nu_\mu','\nu_\tau','\bar\nu_e','\bar\nu_\mu','\bar\nu_\tau'
    write(*,'(118A1)')('=',i=1,118)
    do while(E<Emaxmax)
        do NuAnu=1,2
            do Flavor=1,3
                flx(NuAnu,Flavor)=Flux_Get_dF(NuAnu,Flavor,E)
            enddo
        enddo
        write(*,'(1PE16.8,1X,"||",3E16.8,1X,"|",3E16.8)'),E,flx(1,:),flx(2,:)
        E=E+dE
    enddo
    write(*,'(118A1)')('=',i=1,118)
    write(*,'(A16,1X,"||",3(A13,3X),1X,"|",3(A15,1X))')'E [GeV]',&
    '\nu_e','\nu_\mu','\nu_\tau','\bar\nu_e','\bar\nu_\mu','\bar\nu_\tau'
!----------------------------------------------------------------------!
    stop
endPROGRAM Test_Flux
