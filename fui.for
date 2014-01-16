************************************************************************
      FUNCTION fui(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         EXTERNAL MA_QES_EFF

         COMMON        /N/N                                              Atmospheric neutino spectrum
         COMMON     /n_TT/n_TT
         COMMON    /NuAnu/NuAnu                                           Switch for neutrino type
         COMMON     /N_Fl/N_Fl                                           Switch fot lepton flavor
         COMMON   /N_CorV/N_CorV
         COMMON    /P_lep/P_lep,E_lep                                    Charged lepton momentum
         COMMON    /x_lim/x_ini,deltax                                   Limits (for neutrino energy)
         COMMON    /m_ini/m_ini,mm_ini                                   Mass and square of the mass of initial nuclon
         COMMON    /m_lep/m_lep,mm_lep                                   Mass and square of the mass of charged lepton
         COMMON    /m_fin/m_fin,mm_fin                                   Mass of final hadron or hadron system
         COMMON    /m_tar/m_tar,mm_tar                                   Mass of target nucleus
         COMMON   /MA_QES/MA_QES                                         Mass of axial-vector in QES CC reactions

************************************************************************
         x=deltax*var+x_ini
         E_nu=P_lep/x
         Q2=2*m_ini*(E_nu-E_lep)
         den=mm_W+Q2                                                     (is from W-boson propagator)
!****************************************** select spectrum ******************************
        spectrum=Flux_get_dF(NuAnu,n_Fl,E_nu)
!****************************************** select section ******************************

        section=dsQESCC_dQ2_SM(n_Fl,NuAnu,n_TT,E_nu,Q2,MA_QES_EFF(E_nu))
        
        Jacobianv=E_nu**2/P_lep
        fui=spectrum*section*Jacobianv
       ! WRITE(*,*)"MA=",MA_QES,"Q2=",Q2
       ! WRITE(*,*)"conf=",NuAnu,N_Fl," Enu=",E_nu
       ! WRITE(*,*)"spectrum=",spectrum
       ! WRITE(*,*)"section=",section
       ! WRITE(*,*)"Jacobian=",Jacobianv
       ! WRITE(*,*)" Enu=",E_nu,"fui=",fui
        RETURN

      END FUNCTION fui
