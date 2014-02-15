************************************************************************
      FUNCTION MuL_d3sQES_dQ2dnudkF_SM(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         COMMON       /E_nu/E_nu                                         !Neutrino energy

         DIMENSION var(*)

         CALL Q2QES_SM_lim(E_nu,Q2_min,Q2_max)
         Q2=(Q2_max-Q2_min)*var(1)+Q2_min
         CALL nuQES_SM_lim(E_nu,Q2,nu_min,nu_max)
         nu=(nu_max-nu_min)*var(2)+nu_min
         CALL kFQES_SM_lim(Q2,nu,kF_min,kF_max)
         kF=(kF_max-kF_min)*var(3)+kF_min

         MuL_d3sQES_dQ2dnudkF_SM=d3sQES_dQ2dnudkF_SM(E_nu,Q2,nu,kF)*
     #   (Q2_max-Q2_min)*(nu_max-nu_min)*(kF_max-kF_min)

         RETURN
      END FUNCTION MuL_d3sQES_dQ2dnudkF_SM

************************************************************************
      FUNCTION MuL_QESNuc_dQ2(var)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero

         IMPLICIT REAL (A-Z)

         COMMON         /Q2/Q2                                           !Square of mometum transfer (Q^2=-q^2)
         COMMON       /E_nu/E_nu                                         !Neutrino energy

         DIMENSION var(*)

         CALL nuQES_SM_lim(E_nu,Q2,nu_min,nu_max)
         IF (nu_min.ge.nu_max) THEN
           MuL_QESNuc_dQ2=zero
                               ELSE
           nu=(nu_max-nu_min)*var(1)+nu_min
           CALL kFQES_SM_lim(Q2,nu,kF_min,kF_max)
           IF (kF_min.ge.kF_max) THEN
             MuL_QESNuc_dQ2=zero
                                 ELSE
             kF=(kF_max-kF_min)*var(2)+kF_min
             MuL_QESNuc_dQ2=d3sQES_dQ2dnudkF_SM(E_nu,Q2,nu,kF)*
     #                          (nu_max-nu_min)*(kF_max-kF_min)
        endIF
      endIF
         RETURN
      END FUNCTION MuL_QESNuc_dQ2

************************************************************************
      FUNCTION MuL_d3sQES_dQ2dnudkF_SM_TEST(var)               
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-Z)

         COMMON       /E_nu/E_nu                                         !Neutrino energy
         COMMON     /Q2_lim/Q2_min,Q2_max                                !Q^2 kinematic limits

         DIMENSION var(*)

         Q2=(Q2_max-Q2_min)*var(1)+Q2_min
         CALL nuQES_SM_lim(E_nu,Q2,nu_min,nu_max)
         nu=(nu_max-nu_min)*var(2)+nu_min
         CALL kFQES_SM_lim(Q2,nu,kF_min,kF_max)
         kF=(kF_max-kF_min)*var(3)+kF_min

         MuL_d3sQES_dQ2dnudkF_SM_TEST=
     #   d3sQES_dQ2dnudkF_SM(E_nu,Q2,nu,kF)*
     #   (Q2_max-Q2_min)*(nu_max-nu_min)*(kF_max-kF_min)

         RETURN
      END FUNCTION MuL_d3sQES_dQ2dnudkF_SM_TEST
