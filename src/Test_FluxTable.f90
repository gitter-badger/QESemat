!**********************************************************************!
PROGRAM Test_FluxTable
!----------------------------------------------------------------------!
!Test for first flux table in file with fluxes                         !don't understand it!
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

logical:: &
    Flux_Init,Flux_Open_File,Flux_Read_Hdr,Flux_Read_Table,Flux_Close_File,&
    Flux_Print_Table,Flux_Calc_Spline,Flux_Get_LastNu
real:: &
    Flux_Get_Emin,Flux_Get_Emax,Flux_Get_dF
    
integer,parameter:: &
    Npts=20
logical &
    bufL
integer &
    NuAnu,Flavor,i
real &
    E,dE,E_min,E_max,flux
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
    bufL=Flux_Read_Hdr()
    bufL=Flux_Read_Table()
    bufL=Flux_Close_File()
    bufL=Flux_Get_LastNu(NuAnu,Flavor)
    write(*,*) NuAnu,Flavor
!spectrum limit finding------------------------------------------------!
    bufL=Flux_Calc_Spline(NuAnu,Flavor)
    bufL=Flux_Print_Table(NuAnu,Flavor)
    E_min=Flux_Get_Emin(NuAnu,Flavor)
    E_max=Flux_Get_Emax(NuAnu,Flavor)
    write(*,'("Spectrum Emin=",1PE9.2)') E_min
    write(*,'("Spectrum Emax=",1PE9.2)') E_max
!flux table printing---------------------------------------------------!and what we need *first* table from the file for?!
    E=E_min
    dE=(E_max-E_min)/Npts
    write(*,'(8X,"E [GeV]",1X,"|",10X,"flux")')
    write(*,'(34A1)')('=',i=1,34)
    do while(E<E_max)
        flux=Flux_Get_dF(NuAnu,Flavor,E)
        write(*,'(2(E16.8))'),E,flux
        E=E+dE
    enddo
!----------------------------------------------------------------------!
    stop
endPROGRAM Test_FluxTable
