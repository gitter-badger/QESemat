************************************************************************
      FUNCTION dsQESCC_dQ2_fN_init()
************************************************************************
*                                                                      *
*                                BLTP JINR, Dubna, Russia, 2013/06/01  *
************************************************************************

         USE PhysMathConstants

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         SAVE

         REAL,PARAMETER:: factor= G_Fermi**2*c2C*mm_I/(8*pi)*
     #                            hc2*1.0d+38
         REAL E_thr(3,2)

         COMMON      /m_ini/m_ini,mm_ini                                 Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 Mass of final hadron or hadron system
         COMMON     /MA_QES/MA_QES_tmp                                   Mass of axial-vector in QES CC reactions

         E_thr(1,1)= ((m_p+m_e  )**2-mm_n)/(2*m_n)
         E_thr(2,1)= ((m_p+m_mu )**2-mm_n)/(2*m_n)
         E_thr(3,1)= ((m_p+m_tau)**2-mm_n)/(2*m_n)

         E_thr(1,2)= ((m_n+m_e  )**2-mm_p)/(2*m_p)
         E_thr(2,2)= ((m_n+m_mu )**2-mm_p)/(2*m_p)
         E_thr(3,2)= ((m_n+m_tau)**2-mm_p)/(2*m_p)

         dsQESCC_dQ2_fN_init= one
         RETURN

*     ==================================================================
      ENTRY dsQESCC_dQ2_fN(n_Fl,n_NuAnu,E_nu,Q2,MA_QES)                            ea + p --> e^+   + n
*     ==================================================================
         MA_QES_tmp= MA_QES
         !WRITE(*,*)"Enu=",E_nu,"EnuMIN",E_ea_thr
         IF (E_nu.le.E_thr(n_Fl,n_NuAnu)) THEN
           dsQESCC_dQ2_fN=Precision
                                 ELSE
!            CALL Q2QES_lim(E_nu,Q2_min,Q2_max)
!            WRITE(*,*)"Q2=",Q2,"Q2_lims=",Q2_min,Q2_max
!            IF (Q2.le.Q2_min .or. Q2.ge.Q2_max) THEN
!              dsQESCC_dQ2_fN= Precision
!                                                ELSE
            !WRITE(*,*)"dsQESCC_dQ2=",dsQESCC_dQ2(E_nu,Q2)
             dsQESCC_dQ2_fN= factor*dsQESCC_dQ2(E_nu,Q2)/E_nu**2
!         endIF
      endIF
         RETURN

      END FUNCTION dsQESCC_dQ2_fN_init
