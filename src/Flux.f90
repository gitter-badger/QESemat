!**********************************************************************!
FUNCTION Flux_Init() !bind(C,name='Flux_Init')
!----------------------------------------------------------------------!
!Neutrino flux dF/dE [1/(cm^2 s sr GeV)] spline interpolation
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!

implicit none

logical:: &
    Flux_Init,Flux_Open_File,Flux_Read_Hdr,Flux_Read_Table,Flux_Close_File,&
    Flux_Has_Table,Flux_Has_Spline,Flux_Print_Table,Flux_Calc_Spline,Flux_Get_LastNu
real:: &
    Sp1,&
    Flux_Get_Emin,Flux_Get_Emax,Flux_Get_Zmin,Flux_Get_Zmax,&
    Flux_Get_dF

integer,parameter::&
    flxf=300,&                                                         !Output device number
    MaxNE=1000, MaxNZ=1000,&
    NNuAnu=2, NFlavor=3, NNuType=NFlavor*NNuAnu,&                      !why not just 2 and 3 like all over the program?
    Size_dF=NNuType*MaxNE, Size_d2F=NNuType*MaxNE*MaxNZ
real,parameter::&
    Eps=1.0d-12
character(*),parameter:: &
    hf='-----------------------------',vf='|'
logical &
    LogE
integer &
    ioer,NCdF,&
    iNuAnu,iFlavor
real &
    Energy

save 
logical &
    HasTable (NNuAnu,NFlavor)/NNuType*.false./,&
    HasSpline(NNuAnu,NFlavor)/NNuType*.false./,&
    LogTable (NNuAnu,NFlavor)/NNuType*.true./
integer &
    NuAnu,Flavor,&
    n_E,NE(NNuAnu,NFlavor)/NNuType*0/,&
    n_Z,NZ(NNuAnu,NFlavor)/NNuType*0/,&
    Issue(2,3)                                                         !"numerates different sets of interpolation ranges" - no idea what does it mean
real &
    E(NNuAnu,NFlavor,MaxNE),&
    E_max(NNuAnu,NFlavor),E_min(NNuAnu,NFlavor),&
    Z(NNuAnu,NFlavor,MaxNZ),&
    Z_max(NNuAnu,NFlavor),Z_min(NNuAnu,NFlavor),&
    dF(NNuAnu,NFlavor,MaxNE)/Size_dF*0./,&
    CdF(NNuAnu,NFlavor,MaxNE+2),&
    d2F(NNuAnu,NFlavor,MaxNE,MaxNZ)/Size_d2F*0./,&
    Cd2F(NNuAnu,NFlavor,MaxNE+2,MaxNZ+2),&
    SpLim(4),SpVar(2)                                                  !why array? i suspect some thought was here
character*132 &
    fline
character*80 &
    file_name,filename
character*2 &
    NuType

    Flux_Init=.true.
    return

!**********************************************************************!
INCLUDE 'Flux_Entries.f90'
!----------------------------------------------------------------------!

!**********************************************************************!
ENTRY Flux_Calc_Spline(iNuAnu,iFlavor) !bind(C,name='Flux_Calc_Spline')
!----------------------------------------------------------------------!
    if(HasSpline(iNuAnu,iFlavor))then
        Flux_Calc_Spline=.true.
        return
    endif
    write(*,*) 'Flux spline coefficient calculation...'
    Issue(iNuAnu,iFlavor)=iFlavor+NFlavor*(iNuAnu-1)
!spline limit setting--------------------------------------------------!
    SpLim(1)=E_min(iNuAnu,iFlavor)                                     !why array? i suspect some thought was here
    SpLim(2)=E_max(iNuAnu,iFlavor)
    if(LogTable(iNuAnu,iFlavor))SpLim(1:2)=log10(SpLim(1:2))           !using arrays for limits is not very clear for me
!spline coefficient calculation----------------------------------------!
    call Coeff1(0,1,.true.,Eps,Issue(iNuAnu,iFlavor),NE(iNuAnu,iFlavor),&!what about Mode in 'spline1'? what does mean 'logarithmic'?
    SpLim(1),SpLim(2),dF(iNuAnu,iFlavor,:),CdF(iNuAnu,iFlavor,:),.true.,1)
!----------------------------------------------------------------------!
    write(*,*) 'Flux spline coefficients have been calculated'
    Flux_Calc_Spline=.true.
    return

!**********************************************************************!
ENTRY Flux_Get_dF(iNuAnu,iFlavor,Energy) !bind(C,name='Flux_Get_dF')
!----------------------------------------------------------------------!
!default for flux which doesn't have a table in file-------------------!we have Flux_Has_Table for it!
    if(NE(iNuAnu,iFlavor)==0)then
        Flux_Get_dF=0.
        return
!flux beyond its energy limits-----------------------------------------!
    elseif((Energy<E_min(iNuAnu,iFlavor)).or.(Energy>E_max(iNuAnu,iFlavor)))then
        Flux_Get_dF=0.
        return
!----------------------------------------------------------------------!
    else
!settings--------------------------------------------------------------!
        NCdF=NE(iNuAnu,iFlavor)+2                                      !each time!
        SpVar(1)=Energy                                                !why array? i suspect some thought was here
        if(LogTable(iNuAnu,iFlavor))SpVar(1)=log10(SpVar(1))
!flux interpolation----------------------------------------------------!
        Flux_Get_dF=Sp1(Issue(iNuAnu,iFlavor),CdF(iNuAnu,iFlavor,1:NCdF),SpVar(1))
!patch for interpolation errors----------------------------------------!
        if(Flux_Get_dF<0.)Flux_Get_dF=0.                               !rough!
!----------------------------------------------------------------------!
    endif
    return
!emergency exits-------------------------------------------------------!
 98 write(*,*) 'Flux_Open_File ERROR: Cannot open file ',filename,'!'
    stop
!----------------------------------------------------------------------!
endFUNCTION Flux_Init
!**********************************************************************!
