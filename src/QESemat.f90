!**********************************************************************!
PROGRAM QESemat
!----------------------------------------------------------------------!
!QES neutrino event rate vs. charged lepton momentum
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!
use PhysMathConstants, only: pi,N_Avogadro

implicit none

logical:: &
    EventRate_Init_GeM, &
    EventRate_Init_Flux, &
    EventRate_Init_Section, &
    EventRate_Set_Tgt
integer:: &
    QESNuc_Get_TgtNumb
real:: &
    QESNuc_Get_TgtAWght,EventRate
real,external:: &
    fui                                                                !Integrand

integer,parameter:: &
    outf=100,&                                                         !Output device number
    NElmax=20,&!Different nucleus types that we have                   !into PMC?
    NP_lep=100, Npoint=NElmax*NP_lep, &
    MinCal=100                                                         !GeM setting: minimal number of integrand calls
real,parameter:: &
    RelErr = 1.0d-4,&                                                 !GeM setting: maximal relative error
    cff_flux = 4*pi,&                                                  !Coefficient for neutrino flux averaged over whole sphere
    cff_sect = 1.0d-38,&                                               !Coefficient for section (section is multiplied by 1.00d+38)     
    cff_mass = N_Avogadro*1.0d+03,&                                    !Nuclei in mol, gramms in kg - nuclei in kg, multiplyed by molar mass
    cff_time = 60*60*24*365.25,&                                       !Seconds in year
    cff = cff_sect*cff_mass*cff_time                                   !Coefficient for number of events per kg of detector per second multiplied by molar mass
character(*),parameter:: &
    usage='Usage: ./qesemat "outputfile" "fluxfile" NuAnu[1,2] Flavor[1,2,3] CorV[1,2] M_A[GeV] "mixture" &
&"formula[el1 frac1 el2 frac2...]" P_lep_min[GeV] P_lep_max[GeV]',&
    example='e.g. ./qesemat "output.dat" "flux.sng" 2 3 1 1.0 "water" "H 2 O 1" 0.09 50'
logical &
    bufL
integer &
    NuAnu,Flavor,CorV,&!CorV need more clear name!
    numb_TN(NElmax)/NElmax*0/,&
    NEl,n_El,n_P_lep,&
    ieof
real &
    P_lep,&
    P_lep_min,P_lep_max,lgP_lep_min,lgP_lep_max,steplgP_lep,&
    factor,MA_QES,mu,&
    ValP(NP_lep),frac(NElmax)/NElmax*0/,IntEl(NElmax,NP_lep)/Npoint*0/,R(NP_lep)/NP_lep*0/
character*80 &
    arg,formula,outfile,fluxfile,mixture
character*4 MAn
character*2 name_TN(NElmax)
character*1 &
    NAn(2)/'n','a'/,Fln(3)/'e','m','t'/,CorVn(2)/'c','v'/

!reading arguments-----------------------------------------------------!
    if(iargc()<10)then
        write(*,*) 'QESemat ERROR: Missing arguments!'
        write(*,*) usage
        write(*,*) example
        stop
    endif
    call GetArg( 1,outfile)
    call GetArg( 2,fluxfile)
    call GetArg( 3,arg); read(arg,*) NuAnu
    call GetArg( 4,arg); read(arg,*) Flavor
    call GetArg( 5,arg); read(arg,*) CorV
    call GetArg( 6,arg); read(arg,*) MA_QES; write(MAn,'(F4.2)') MA_QES
    call GetArg( 7,mixture)
    call GetArg( 8,arg); read(arg,'(A80)') formula
    call GetArg( 9,arg); read(arg,*) P_lep_min
    call GetArg(10,arg); read(arg,*) P_lep_max
!echo------------------------------------------------------------------!
    write(*,*) 'Output to file: ',outfile
    write(*,'("Set to:",1X,2A1,1X,"on",1X,A80)') NAn(NuAnu),Fln(Flavor),formula
    write(*,'(A1,1X,"MA_QES=",F6.3)') CorVn(CorV),MA_QES
    write(*,*) 'Flux file: ',fluxfile
    write(*,'(2(A10,F8.3,1X))') 'P_lep_min=',P_lep_min,'P_lep_max=',P_lep_max
!formula processing----------------------------------------------------!
    read(formula,*,iostat=ieof) (name_TN(n_El),frac(n_El),n_El=1,NElmax)!error report?.. but how?
    NEl=1
    do n_El=1,NElmax
        if(frac(n_El)>0.)then
            numb_TN(n_El)=QESNuc_Get_TgtNumb(name_TN(n_El))
            if(numb_TN(n_El)==-999)stop                                !don't like this 999! why not stopping in QESNuc_Get_TgtNumb?
            write(*,'("element #",I2,":",1X,A2,1X,"(",I2,")",1X,F7.3)')&
            n_El,name_TN(n_El),numb_TN(n_El),frac(n_El)
            NEl=n_El
        endif
    enddo
!mixture molar mass calculation----------------------------------------!
    mu=0.
    do n_El=1,NEl
        if(frac(n_El)>0.)mu=mu+QESNuc_Get_TgtAWght(numb_TN(n_El))*frac(n_El)!molar mass [g/mol], numerically equals to atomic mass, num.app.eq. nucleon number
                                                                       !not actually A! why though?..
    enddo
    factor=cff/mu                                                      !Coefficient for number of events per kg of detector per second
!settings--------------------------------------------------------------!
    bufL=EventRate_Init_Flux(fluxfile,NuAnu,Flavor)
    bufL=EventRate_Init_Section(NuAnu,Flavor,CorV,MA_QES)
    bufL=EventRate_Init_GeM(RelErr,MinCal)
!calculation of calculation points-------------------------------------!English?..
    lgP_lep_min=log10(P_lep_min)
    lgP_lep_max=log10(P_lep_max)
    steplgP_lep=(lgP_lep_max-lgP_lep_min)/(NP_lep-1)
    do n_P_lep=1,NP_lep
        ValP(n_P_lep)=10**(lgP_lep_min+(n_P_lep-1)*steplgP_lep)
    enddo
!----------------------------------------------------------------------!
    do n_El=1,NEl
        if(frac(n_El)>0.)then
            bufL=EventRate_Set_Tgt(numb_TN(n_El))
            do n_P_lep=1,NP_lep
                P_lep=ValP(n_P_lep)
                IntEl(n_El,n_P_lep)=frac(n_El)*EventRate(P_lep)
                write(*,'(I3,"/",I3,3X,"P_lep=",F8.3,2X,"IntEl=",1PE15.8)')&
                n_P_lep,NP_lep,P_lep,IntEl(n_El,n_P_lep)
            enddo
            R=R+IntEl(n_El,:)
        endif
    enddo
    IntEl=factor*IntEl
    R=factor*R
!output----------------------------------------------------------------!
    write(*,*) 'Output to file ',outfile
    open(outf,file=outfile)
    write(outf,'("# P_lep [GeV]",20(21X,A2))',advance='no') (name_TN(n_El),n_El=1,NEl)
    write(outf,*) adjustl(formula)
    do n_P_lep=1,NP_lep
        write(*,'(I3,"/",I3,3X,"P_lep=",F8.3,2X,"R=",1PE15.8)')&
        n_P_lep,NP_lep,ValP(n_P_lep),R(n_P_lep)
        write(outf,'(22E23.15)') ValP(n_P_lep),(IntEl(n_El,n_P_lep),n_El=1,NEl),R(n_P_lep)
    enddo
    call GeMInf
    stop 'QESemat finished'
!----------------------------------------------------------------------!
endPROGRAM QESemat
!**********************************************************************!
