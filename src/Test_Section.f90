!**********************************************************************!
PROGRAM Test_Section
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!
use PhysMathConstants, only: pi,N_Avogadro

implicit none

logical:: &
    CrossSection_Init,CrossSection_Set_Tgt
integer:: &
    QESNuc_Get_TgtNumb
real:: &
    QESNuc_Get_TgtAWght,CrossSection
real,external:: &
    fui                                                                !Integrand

integer,parameter:: &
    outf=100,&                                                         !Output device number
    NElmax=20,&!Different nucleus types that we have                   !into PMC?
    NE_nu=100, Npoint=NElmax*NE_nu,&
    MinCal=100                                                         !GeM setting: minimal number of integrand calls
real,parameter:: &
    RelErr = 1.0d-13,&                                                 !GeM setting: maximal relative error
    cff_sect = 1.0d-38,&                                               !Coefficient for section (section is multiplied by 1.00d+38)     
    cff_mass = N_Avogadro*1.0d+03,&                                    !Nuclei in mol, gramms in kg - nuclei in kg, multiplyed by molar mass
    cff = cff_sect*cff_mass                                            !Coefficient for number of events per kg of detector per second multiplied by molar mass
character(*),parameter:: &
    usage='Usage: ./Test_Section "outputfile" NuAnu[1,2] Flavor[1,2,3] CorV[1,2] M_A[GeV] "mixture" &
&Mode[0,1] "formula[el1 frac1 el2 frac2...]" E_nu_min[GeV] E_nu_max[GeV]',&
    example='e.g. ./Test_Section "output.dat" 2 3 1 1.0 0 "hydrogen" "H 1" 3 20'
logical &
    bufL
integer &
    NuAnu,Flavor,CorV,Mode,&!CorV,Mode need more clear name!
    numb_TN(NElmax)/NElmax*0/,&
    NEl,n_El,n_E_nu,&
    ieof
real &
    E_nu,&
    E_nu_min,E_nu_max,lgE_nu_min,lgE_nu_max,steplgE_nu,&
    factor,MA_QES,mu,&
    ValE(NE_nu),frac(NElmax)/NElmax*0/,IntEl(NElmax,NE_nu)/Npoint*0/,R(NE_nu)/NE_nu*0/
character*80 &
    arg,formula,outfile,mixture
character*4 MAn
character*2 name_TN(NElmax)
character*1 &
    NAn(2)/'n','a'/,Fln(3)/'e','m','t'/,CorVn(2)/'c','v'/

!reading arguments-----------------------------------------------------!
    if(iargc()<10)then
        write(*,*) 'Test_Section ERROR: Missing arguments!'
        write(*,*) usage
        write(*,*) example
        stop
    endif
    call GetArg( 1,outfile)
    call GetArg( 2,arg); read(arg,*) NuAnu
    call GetArg( 3,arg); read(arg,*) Flavor
    call GetArg( 4,arg); read(arg,*) CorV
    call GetArg( 5,arg); read(arg,*) MA_QES; write(MAn,'(F4.2)') MA_QES
    call GetArg( 6,mixture)
    call GetArg( 7,arg); read(arg,*) Mode
    call GetArg( 8,arg); read(arg,'(A80)') formula
    call GetArg( 9,arg); read(arg,*) E_nu_min
    call GetArg(10,arg); read(arg,*) E_nu_max
!echo------------------------------------------------------------------!
    write(*,*) 'Output to file: ',outfile
    write(*,'("Set to:",1X,2A1,1X,"on",1X,A80)') NAn(NuAnu),Fln(Flavor),formula
    write(*,'(A1,1X,"MA_QES=",F6.3)') CorVn(CorV),MA_QES
    write(*,'(2(A10,F8.3,1X))') 'E_nu_min=',E_nu_min,'E_nu_max=',E_nu_max
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
        if(frac(n_El)>0.)then
            if(Mode==0)then
                mu=mu+QESNuc_Get_TgtAWght(numb_TN(n_El))*frac(n_El)    !molar mass [g/mol], numerically equals to atomic mass, num.app.eq. nucleon number
                                                                       !not actually A! why though?..
            else
                mu=mu+frac(n_El)
                write(*,'(A2,1X,F7.3)',advance='no') name_TN(n_El),frac(n_El)
                frac(n_El)=frac(n_El)/QESNuc_Get_TgtAWght(numb_TN(n_El))
                write(*,*) '->',frac(n_El)
            endif
        endif
    enddo
    factor=cff/mu                                                      !Coefficient for number of events per kg of detector per second
!settings--------------------------------------------------------------!
    call GeMSet(RelErr,MinCal,*99)
    bufL=CrossSection_Init(NuAnu,Flavor,CorV,MA_QES)
!calculation of calculation points-------------------------------------!English?..
    lgE_nu_min=log10(E_nu_min)
    lgE_nu_max=log10(E_nu_max)
    steplgE_nu=(lgE_nu_max-lgE_nu_min)/(NE_nu-1)
    do n_E_nu=1,NE_nu
        ValE(n_E_nu)=10**(lgE_nu_min+(n_E_nu-1)*steplgE_nu)
    enddo
!----------------------------------------------------------------------!
    do n_El=1,NEl
        if(frac(n_El)>0.)then
            bufL=CrossSection_Set_Tgt(numb_TN(n_El))
            do n_E_nu=1,NE_nu
                E_nu=ValE(n_E_nu)
                IntEl(n_El,n_E_nu)=frac(n_El)*CrossSection(E_nu)
                write(*,'(I3,"/",I3,3X,"E_nu=",F8.3,2X,"IntEl=",1PE15.8)')&
                n_E_nu,NE_nu,E_nu,IntEl(n_El,n_E_nu)
            enddo
            R=R+IntEl(n_El,:)
        endif
    enddo
    IntEl=factor*IntEl
    R=factor*R
!output----------------------------------------------------------------!
    write(*,*) 'Output to file ',outfile
    open(outf,file=outfile)
    write(outf,'("#  E_nu [GeV]",20(21X,A2))',advance='no') (name_TN(n_El),n_El=1,NEl)
    write(outf,*) adjustl(formula)
    do n_E_nu=1,NE_nu
        write(*,'(I3,"/",I3,3X,"E_nu=",F8.3,2X,"R=",1PE15.8)')&
        n_E_nu,NE_nu,ValE(n_E_nu),R(n_E_nu)
        write(outf,'(22E23.15)') ValE(n_E_nu),(IntEl(n_El,n_E_nu),n_El=1,NEl),R(n_E_nu)
    enddo
    close(outf)
    call GeMInf
    stop 'Test_Section finished'
!emergency exits-------------------------------------------------------!
 99 stop 'Test_Section ERROR: GeMSet failed!'
!----------------------------------------------------------------------!
endPROGRAM Test_Section
!**********************************************************************!
