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
         COMMON   /NuAnu/NuAnu                                           Switch for neutrino type
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

         SELECTCASE(N_Fl)
               CASE(1)
                     IF (NuAnu.EQ.1) THEN
                       spectrum=dFANomen_dE(E_nu)
                                     ELSE
                       spectrum=dFANomea_dE(E_nu)
                  endIF
               CASE(2)
                     IF (NuAnu.EQ.1) THEN
                       spectrum=dFANommn_dE(E_nu)
                                     ELSE
                       spectrum=dFANomma_dE(E_nu)
                  endIF
               CASE(3)
                     IF (NuAnu.EQ.1) THEN
                       spectrum=dFANomtn_dE(E_nu)
                                     ELSE
                       spectrum=dFANomta_dE(E_nu)
                  endIF
      endSELECT
!****************************************** select section ******************************

           SELECTCASE(N_Fl)
                 CASE(1)
                       IF (NuAnu.EQ.1) THEN
                         section=dsQESCC_dQ2_SM_en(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                                       ELSE
                         section=dsQESCC_dQ2_SM_ea(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                    endIF
                 CASE(2)
                       IF (NuAnu.EQ.1) THEN
                         section=dsQESCC_dQ2_SM_mn(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                                       ELSE
                         section=dsQESCC_dQ2_SM_ma(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                    endIF
                 CASE(3)
                       IF (NuAnu.EQ.1) THEN
                         section=dsQESCC_dQ2_SM_tn(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                                       ELSE
                         section=dsQESCC_dQ2_SM_ta(E_nu,Q2,
     #                                                 MA_QES_EFF(E_nu))
                    endIF
        endSELECT
           Jacobianv=E_nu**2/P_lep
         fui=spectrum*section*Jacobianv
         WRITE(*,*)"MA=",MA_QES,"Q2=",Q2
         WRITE(*,*)"conf=",N_ForB,NuAnu,N_Fl," Enu=",E_nu
         WRITE(*,*)"fui=",fui
         RETURN

      END FUNCTION fui