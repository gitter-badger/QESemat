************************************************************************
      FUNCTION rho_SM(N,p)
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON    /P_Fermi/P_Fermi                                      !Fermi momentum of target nucleon
         COMMON    /T_Fermi/T_Fermi                                      !Effective Fermi temperature of target nucleon

         real,parameter:: MAXEXP=log(huge(real))  !maximum exponent
*           ---------------------------------------------------------- *
          IF (N.eq.0) THEN                                               !No Pauli blocking
*           ---------------------------------------------------------- *
            rho_SM= 0.0
*           ---------------------------------------------------------- *
      ELSEIF (N.eq.1) THEN                                               !Pure Fermi gaz with T_Fermi=0
*           ---------------------------------------------------------- *
            IF (p.le.P_Fermi) THEN; rho_SM= 1.0
                              ELSE; rho_SM= 0.0
         endIF
*           ---------------------------------------------------------- *
      ELSEIF (N.eq.2) THEN                                               !Fermi-Dirac distribution
*           ---------------------------------------------------------- *
            var=-(P_Fermi-p)/T_Fermi
            !write(*,*)"MAXEXP=",MAXEXP; stop
            if(var.ge.MAXEXP)then
                rho_SM=0
            else
                rho_SM= 1.0/(1.0+exp(var))
            end if
*           ---------------------------------------------------------- *
                      ELSE
*           ---------------------------------------------------------- *
            rho_SM= 1.0
*           ---------------------------------------------------------- *
       endIF

         RETURN
      END FUNCTION rho_SM
