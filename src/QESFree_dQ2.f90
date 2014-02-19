!**********************************************************************!
FUNCTION QESFree_dQ2_Init()
!----------------------------------------------------------------------!
!ds/dQ^2 of QES on free nucleon
!----------------------------------------------------------------------!
!edited by                                               ???, O.Petrova!
!**********************************************************************!
use PhysMathConstants

implicit none

logical:: &
    QESFree_dQ2_Init
real:: &
    dsQESCC_dQ2,&
    QESFree_dQ2

common/MA_QES/MA_QES!Mass of axial-vector in QES reactions             !dsQESCC_dQ2 and NucQESFF need the common-block

integer,parameter:: &
    NNuAnu=2, NFlavor=3
real,parameter:: &
    factor=G_Fermi**2*c2C*mm_I/(8*pi)*hc2*1.0d+38                      !multiply each time! and in *Nuc* is a part of this coefficient
integer &
    NuAnu,Flavor,&
    iNuAnu,iFlavor
real &
    E_nu,Q2,Q2_min,Q2_max,MA_QES,&
    iMA_QES

save
real &
    mass_lep(NFlavor)/m_e,m_mu,m_tau/,&                                !maybe into PMC? or special module? procedure-module?
    mass_fin(NNuAnu)/m_p,m_n/,&
    mass_ini(NNuAnu)/m_n,m_p/,&
    E_nu_FrThr(NNuAnu,NFlavor)                                         !Neutrino energy threshold for scattering on free nucleon

!neutrino energy threshold calculation---------------------------------!now in *Nuc*, but if we want to use this function separately...
!    do NuAnu=1,2
!        do Flavor=1,3
!            E_nu_FrThr(NuAnu,Flavor)=&
!            0.5*((mass_fin(NuAnu)+mass_lep(Flavor))**2-mass_ini(NuAnu)**2)/mass_ini(NuAnu)
!        enddo
!    enddo
!----------------------------------------------------------------------!
    QESFree_dQ2_Init=.true.
    return

!**********************************************************************!
ENTRY QESFree_dQ2(iNuAnu,iFlavor,E_nu,Q2,iMA_QES)
!----------------------------------------------------------------------!
    MA_QES=iMA_QES
!    if(E_nu<=E_nu_FrThr(iNuAnu,iFlavor))then
!        QESFree_dQ2=0.
!    else
!        call Q2QES_lim(E_nu,Q2_min,Q2_max)
!        if((Q2<=Q2_min).or.(Q2>=Q2_max))then
!            QESFree_dQ2=0.
!        else
            QESFree_dQ2=factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2            !what's 'dsQESCC_dQ2'? why it's not multiple by factor and E_nu^2 itself?
!        endif
!    endif
    return
!----------------------------------------------------------------------!
endFUNCTION QESFree_dQ2_Init
!**********************************************************************!

!**********************************************************************!
FUNCTION QESDeut_dQ2_Init()
!----------------------------------------------------------------------!
!ds/dQ^2 of QES on one nucleon of Deuterium
!----------------------------------------------------------------------!
!edited by                                                    O.Petrova!
!**********************************************************************!
implicit none

logical:: &
    QESFree_dQ2_Init,&
    QESDeut_dQ2_Init
real:: &
    FactorPauli_D2,QESFree_dQ2,&
    QESDeut_dQ2

logical &
    bufL
integer &
    NuAnu,Flavor,&
    iNuAnu,iFlavor
real &
    E_nu,Q2,Q2_min,Q2_max,MA_QES
    
!free nucleon cross section initialization-----------------------------!
    bufL=QESFree_dQ2_Init()
!----------------------------------------------------------------------!
    QESDeut_dQ2_Init=.true.
    return

!**********************************************************************!
ENTRY QESDeut_dQ2(iNuAnu,iFlavor,E_nu,Q2,MA_QES)
!----------------------------------------------------------------------!
    QESDeut_dQ2=FactorPauli_D2(Q2)*QESFree_dQ2(iNuAnu,iFlavor,E_nu,Q2,MA_QES)
    return
!----------------------------------------------------------------------!
endFUNCTION QESDeut_dQ2_Init
!**********************************************************************!
