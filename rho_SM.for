************************************************************************
      FUNCTION rho_SM(N,p)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: zero,one

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON    /P_Fermi/P_Fermi                                      !Fermi momentum of target nucleon
         COMMON    /T_Fermi/T_Fermi                                      !Effective Fermi temperature of target nucleon

*           ---------------------------------------------------------- *
          IF (N.eq.0) THEN                                               !No Pauli blocking
*           ---------------------------------------------------------- *
            rho_SM= zero
*           ---------------------------------------------------------- *
      ELSEIF (N.eq.1) THEN                                               !Pure Fermi gaz with T_Fermi=0
*           ---------------------------------------------------------- *
            IF (p.le.P_Fermi) THEN; rho_SM= one
                              ELSE; rho_SM= zero
         endIF
*           ---------------------------------------------------------- *
      ELSEIF (N.eq.2) THEN                                               !Fermi-Dirac distribution
*           ---------------------------------------------------------- *
            var=-(P_Fermi-p)/T_Fermi
            if(var.gt.700)then
                rho_SM=0
            else
                rho_SM= one/(one+exp(var))
            end if
*           ---------------------------------------------------------- *
                      ELSE
*           ---------------------------------------------------------- *
            rho_SM= one
*           ---------------------------------------------------------- *
       endIF

         RETURN
      END FUNCTION rho_SM