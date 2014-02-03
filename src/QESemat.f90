!**********************************************************************!
PROGRAM QESemat
!----------------------------------------------------------------------!
!
!----------------------------------------------------------------------!
!edited by                                       O.Petrova, A.Sheshukov!
!**********************************************************************!
use PhysMathConstants

implicit none

logical:: &
    Flux_Init,Flux_Open_File,Flux_Read_Head,Flux_Read_Table,Flux_Close_File,&
    Flux_Has_Table,Flux_Print_Table,Flux_Calc_Spline,&
    MA_QES_init,QES_print
integer:: &
    QES_Get_tgtNumb                                                    !until we invent a better name
real:: &
    QES_Get_tgtA,QES_Get_tgtNcln,&
    Flux_Get_Emin,Flux_Get_Emax,&
    QES_init,dsqescc_dq2_fn_init,dsqescc_dq2_fp_set                    !rename dsqescc_dq2_*!
real,external:: &
    fui                                                                !Integrand

common/NuAnu/NuAnu                                                     !Switch: neutrino/antineutrino
common/Flavor/Flavor                                                   !Switch fot lepton flavor
common/Target/Target                                                   !Switch for target nucleus
common/P_lep/P_lep,E_lep                                               !Charged lepton momentum and energy
common/x_lim/x_ini,deltax                                              !Integration limits
common/m_ini/m_ini,mm_ini                                              !Mass and squared mass of initial nucleon
common/MA_QES/MA_QES                                                   !Mass of axial-vector in QES reactions
common/MA_ELS/MA_ELS                                                   !Mass of axial-vector in ELS reactions
common/n_PT/n_PT                                                       !Switch for lepton polarization type
common/n_AG_QES/n_AG                                                   !Switch for model of QES reactions
common/n_MC_QES/n_MC                                                   !Switch for correction of masses of initial and final nucleons
common/n_AP_QES/n_AP                                                   !Switch for model of axial form factor in QES reactions
common/n_MS_QES/n_MS                                                   !Switch for value of axial mass in Sehgal's model
common/n_GE_QES/n_GE                                                   !Switch for parametrization of Sachs electric form factor of neutron
common/MV_QES/MV_QES                                                   !Mass of isovector in QES reactions
common/MM_QES/MM_QES                                                   !Mass of monopole axial-vector in QES reactions
common/MS_QES/MS_QES                                                   !Mass of scalar in QES reactions
common/MT_QES/MT_QES                                                   !Mass of tensor in QES reactions
common/xi_V/xi_V                                                       !Normalization of vector form factor
common/xi_M/xi_M                                                       !Normalization of monopole form factor
common/xi_S/xi_S                                                       !Normalization of scalar form factor
common/xi_A/xi_A                                                       !Normalization of axial form factor
common/xi_P/xi_P                                                       !Normalization of pseudoscalar form factor
common/xi_T/xi_T                                                       !Normalization of tensor form factor
common/phi_T/phi_T                                                     !Phase of tensor form factor
common/phi_S/phi_S                                                     !Phase of scalar form factor

logical &
    bufL
integer,parameter:: &
    outf=100,&                                                         !Output device number
    NElmax=19,&                                                        !Different nucleus types that we have
    NP_lep=100,NPtotal=NElmax*NP_lep,&!don't like name 'NPtotal'
    MinCal=100,&                                                       !GeM setting: minimal number of integrand calls
    n_FF_QES=8                                                         !NucQESFF setting: Bodek,Avvakumov,Bradford&Budd form factor
real,parameter:: &
    RelErr = 1.0d-13,&                                                 !GeM setting: maximal relative error
    cff_flux = 4*pi,&                                                  !Coefficient for neutrino flux averaged over whole sphere
    cff_sect = 1.0d-38,&                                               !Coefficient for section (section is multiplied by 1.00d+38)     
    cff_mass = N_Avogadro*1.0d+03,&                                    !Nuclei in mol, gramms in kg - nuclei in kg multiplyed by molar mass
    cff_time = 60*60*24*365.25,&                                       !Seconds in year
    factor   = cff_flux*cff_sect*cff_mass*cff_time                     !Coefficient for number of events per kg of detector per second multiplied by molar mass
character(*),parameter:: &
    usage='Usage: ./qesemat "outputfile" "fluxfile" NuAnu[1,2] Flavor[1,2,3] CorV[1,2] M_A[GeV] "mixture" &
&"formula[el1 frac1 el2 frac2...]" P_lep_min[GeV] P_lep_max[GeV]',&
    example='e.g. ./qesemat "output.dat" "flux.sng" 2 3 1 1.0 "water" "H 2 O 1" 0.09 50'
integer &
    NuAnu,Flavor,Target,CorV,&
    numb_TN(NElmax)/NElmax*0/,&
    NEl,n_El,n_P_lep,&
    ieof,&
    n_AG,n_AP,n_GE,n_MC,n_MS,n_PT
real &
    bufR,fact,Jacobian_c,MA_QES,mu,Res,&
    x_ini,x_fin,deltax,P_lep,E_lep,m_ini,mm_ini,m_lep,mm_lep,&
    P_lep_min,P_lep_max,lgP_lep_min,lgP_lep_max,steplgP_lep,&
    E_nu_min,E_nu_max,E_nu_ini,E_nu_fin,E_nu_low,E_nu_upp,&
    EmO_lep,EpP_lep,O_lep,P_kth,&
    ValP(NP_lep),frac(NElmax)/NElmax*0/,Intel(NElmax,NP_lep)/NPtotal*0/,R(NP_lep)/NP_lep*0/,&
    MA_ELS,MM_QES,MS_QES,MT_QES,MV_QES,phi_S,phi_T,xi_A,xi_M,xi_P,xi_S,xi_T,xi_V
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
    write(*,'(A8,2A1,1X,A3,A80)') 'Set to: ',NAn(NuAnu),Fln(Flavor),'on ',formula
    write(*,'(A1,1X,A7,F6.3)') CorVn(CorV),'MA_QES=',MA_QES
    write(*,*) 'Flux file: ',fluxfile
    write(*,'(2(A10,F8.3,1X))') 'P_lep_min=',P_lep_min,'P_lep_max=',P_lep_max
    write(*,*) 'Output to file: ',outfile
!formula processing----------------------------------------------------!
    read(formula,*,iostat=ieof) (name_TN(n_El),frac(n_El),n_El=1,NElmax)
    NEl=1
    do n_El=1,NElmax
        if(frac(n_El)>0.)then
            numb_TN(n_El)=QES_Get_tgtNumb(name_TN(n_El))
            if(numb_TN(n_El)==-999)stop                                !don't like this 999
            write(*,'(A9,I2,A2,A2,1X,A1,I2,A1,1X,F7.3)')'element #',n_El,': ',name_TN(n_El),'(',numb_TN(n_El),')',frac(n_El)
            NEl=n_El
        endif
    enddo
!mixture molar mass calculation----------------------------------------!
    mu=0.
    do n_El=1,NEl
        if(frac(n_El)>0.)mu=mu+QES_Get_tgtA(numb_TN(n_El))*frac(n_El)  !molar mass, numerically equals to atomic mass, num.app.eq. nucleon number
                                                                       !not actually A!
    enddo
    fact=factor/mu                                                     !comment it!
!settings: GeM---------------------------------------------------------!setting INTO GeM
    call GeMSet(fui,1.,0.,1.,RelErr,MinCal,*99)
!settings: neutrino flux-----------------------------------------------!settings FROM and INTO InitFlux
    bufL=Flux_init()
    bufL=Flux_open_file(fluxfile)
    do while(Flux_read_head())
        bufL=Flux_read_table()
    enddo
    bufL=Flux_close_file()
    if(Flux_has_table(NuAnu,Flavor).eqv..FALSE.)then
        write(*,*) 'QESemat ERROR: Flux table doesn''t exist!'
        stop
    endif
    bufL=Flux_print_table(NuAnu,Flavor)                                !do we really need this?
    bufL=Flux_calc_spline(NuAnu,Flavor)
         
    E_nu_min=Flux_Get_Emin(NuAnu,Flavor)
    E_nu_max=Flux_Get_Emax(NuAnu,Flavor)
    write(*,'(2(A9,F10.3,1X))') 'E_nu_min=',E_nu_min,'E_nu_max=',E_nu_max
!settings: MA_QES_eff--------------------------------------------------!setting INTO MA_QES_eff
    bufL=MA_QES_init(CorV,0,MA_QES)
!settings: *QES*-------------------------------------------------------!setting INTO *QES*
    bufR=QES_init(n_FF_QES,1,2)
!    CALL QES_PRINT_ALL()
    bufR=dsQESCC_dQ2_fN_init()
    bufR=dsQESCC_dQ2_FP_set(0.,0.,MA_QES)                              !rewrite like dsQESCC_dQ2_fN_init!
!settings: other-------------------------------------------------------!
    MA_ELS=MA_QES
    n_PT=0; n_AG=1; n_MC=1; n_AP=1; n_MS=1; n_GE=5
    MV_QES=0.84; MM_QES=0.8; MS_QES=1.; MT_QES=1.5
    xi_V=1.; xi_M=1.; xi_S=0.; xi_A=1.; xi_P=1.; xi_T=0.
    phi_T=0.; phi_T=phi_T*dtr
    phi_S=0.; phi_S=phi_S*dtr
    if(NuAnu==1)then
        m_ini=m_n; mm_ini=mm_n
    else
        m_ini=m_p; mm_ini=mm_p
    endif
    selectcase(Flavor)
        case(1)
            m_lep=m_e; mm_lep=mm_e
        case(2)
            m_lep=m_mu; mm_lep=mm_mu
        case(3)
            m_lep=m_tau; mm_lep=mm_tau
    endselect
    O_lep=half*mm_lep/m_I
    P_kth=half*m_I-O_lep
!----------------------------------------------------------------------!
    lgP_lep_min=log10(P_lep_min)
    lgP_lep_max=log10(P_lep_max)
    steplgP_lep=(lgP_lep_max-lgP_lep_min)/(NP_lep-1)
!----------------------------------------------------------------------!
    do n_El=1,NEl
        if(frac(n_El)>0.)then
            Target=numb_TN(n_El)
            bufL=QES_print(Target)
            do n_P_lep=1,NP_lep
                P_lep=10**(lgP_lep_min+(n_P_lep-1)*steplgP_lep)
                ValP(n_P_lep)=P_lep
                E_lep=sqrt(P_lep*P_lep+mm_lep)
                EpP_lep=E_lep+P_lep
                EmO_lep=E_lep-O_lep
                E_nu_low=EmO_lep*EpP_lep/(EpP_lep-2.*O_lep)
                E_nu_upp=EmO_lep*m_I/(m_I-EpP_lep)
                E_nu_ini=max(E_nu_low,E_nu_min)
                if(P_lep<=abs(P_kth))then
                    if(P_kth<0.)then
                        E_nu_fin=E_nu_ini
                    else
                        E_nu_fin=min(E_nu_upp,E_nu_max)
                    endif
                else
                    E_nu_fin=E_nu_max
                endif
                x_ini=P_lep/E_nu_fin
                x_fin=P_lep/E_nu_ini
                deltax=x_fin-x_ini                                     !-
                call GeMInt(fui,Res,0.,1.,*100)
                Jacobian_c=2.*m_ini*P_lep/E_lep
                IntEl(n_El,n_P_lep)=frac(n_El)*QES_Get_tgtNcln(NuAnu,Target)*deltax*Res*Jacobian_c
                WRITE(*,*)n_P_lep,'/',NP_lep,'P_lep=',P_lep,Intel(n_El,n_P_lep)
            enddo
            R=R+Intel(n_El,:)
        endif
    enddo
    R=fact*R
    Intel=fact*Intel
    write(*,*) 'Output to file ',outfile
    open(outf,file=outfile)
    write(outf,'(11A25)')adjustl('#Plep, GeV'),(adjustl(name_TN(n_El)),n_El=1,NEl),adjustl(formula)
    do n_P_lep=1,NP_lep
        write(*,*) n_P_lep,'/',NP_lep,'P_lep=',ValP(n_P_lep),R(n_P_lep)
        write(outf,*) ValP(n_P_lep),R(n_P_lep),(Intel(n_El,n_P_lep),n_El=1,NEl)
        write(outf,*)
    enddo
    call GeMInf
stop 'QESemat finished'
!emergency exits-------------------------------------------------------!
 99 stop 'QESemat ERROR: GeMSet failed!'
100 stop 'QESemat ERROR: GeMInt failed!'
!----------------------------------------------------------------------!
endPROGRAM QESemat
!**********************************************************************!
