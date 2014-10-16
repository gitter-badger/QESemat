************************************************************************
      FUNCTION dsQESCC_dQ2(E_nu,Q2)
************************************************************************
*                                                                      *
*                                                                      *
*     REFERENCES                                                       *
*                                                                      *
*     [1] C.H. Llewellyn Smith, "Neutrino reactions at accelerator     *
*         energies,"  Phys. Rept. 3 (1972) 261-379 [Preprint SLAC-     *
*         PUB-958 (1971)].                                             *
*     [2] E.A. Paschos and J.Y. Yu,  "Neutrino interactions in os-     *
*         cillation  experiments,"  Phys. Rev. D 65 (2002)  033002     *
*         [arXiv: hep-ph/0107261].                                     *
*                                                                      *
*     [3] K.S. Kuzmin,  "Neutrino scattering  off nucleons and po-     *
*         larization of charged  leptons  in quasielastic  reacti-     *
*         ons,"  Ph.D. Thesis, JINR, Dubna, 2009/04/01 (Ph.D. The-     *
*         sis advisor V.A. Naumov, BLTP JINR).                         *
*                                                                      *
************************************************************************

         USE PhysMathConstants, ONLY: m_I,mm_I

         IMPLICIT REAL (A-M,O-Z), INTEGER (N)

         COMMON       /n_NT/n_NT                                         !Switch for neutrino type
         COMMON       /n_PT/n_PT                                         !Switch for lepton polarization type
         COMMON   /n_AG_QES/n_AG                                         !Switch for model of QES reactions
         COMMON   /n_MC_QES/n_MC                                         !Switch for correction of masses of initial and final nucleons
         COMMON      /m_ini/m_ini,mm_ini                                 !Mass of target nucleon
         COMMON      /m_lep/m_lep,mm_lep                                 !Mass of final charged lepton
         COMMON      /m_fin/m_fin,mm_fin                                 !Mass of final hadron or hadron system

         SELECTCASE(n_AG)
*              ------------------------------------------------------- *
               CASE(   1)                                                !V.A. Naumov et al. [3]
*              ------------------------------------------------------- *
               x= Q2/(4*mm_I)
               s= mm_ini+2*m_ini*E_nu
               u= mm_fin-2*m_ini*E_nu+mm_lep+Q2
               z= mm_lep/(4*mm_I)
               IF (n_MC.eq.0) THEN
                 r= 0.0
                              ELSE
                 r= (m_fin-m_ini)/(2*m_I)
            endIF
               CALL FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                      ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)

               a0=+(1.0+x)*(ReF_A**2+ImF_A**2-4*x*(ReF_T**2+ImF_T**2))
     #            -(1.0-x)*(ReF_V**2+ImF_V**2-  x*(ReF_M**2+ImF_M**2))
     #            +4*x*(ReF_V*ReF_M+ImF_V*ImF_M)
     #            -z*(ReF_V**2+ImF_V**2+ReF_M**2+ImF_M**2
     #            +2*(ReF_V*ReF_M+ImF_V*ImF_M)
     #            +ReF_A**2+ImF_A**2+4*(ReF_P**2+ImF_P**2)
     #            +4*(ReF_A*ReF_P+ImF_A*ImF_P)
     #            -4*(1.0+x)*(ReF_P**2+ImF_P**2+ReF_S**2+ImF_S**2))

               a1=+(1.0+x)*x*(ReF_T*ReF_A+ImF_T*ImF_A)
     #            +z*(-n_NT*(ReF_A*ReF_V+ReF_A*ReF_M
     #                   +ImF_A*ImF_V+ImF_A*ImF_M)
     #            +    x*(ReF_T*ReF_A+ImF_T*ImF_A
     #                   +2*(ReF_T*ReF_P+ImF_T*ImF_P))
     #            +(1.0+x)*(ReF_S*ReF_V+ImF_S*ImF_V))
     #            +z**2*(ReF_S*ReF_V+ImF_S*ImF_V
     #                  +ReF_S*ReF_M+ImF_S*ImF_M)

               a2=+(1.0+x)*(ReF_A**2+ImF_A**2+8*x*(ReF_T**2+ImF_T**2))
     #            +(1.0-x)*(ReF_V**2+ImF_V**2-  x*(ReF_M**2+ImF_M**2))
     #            -4*x*(ReF_V*ReF_M+ImF_V*ImF_M)
     #            -z*(ReF_V**2+ImF_V**2+ReF_M**2+ImF_M**2
     #               +2*(ReF_V*ReF_M+ImF_V*ImF_M)
     #               -(ReF_A**2+ImF_A**2)-4*(ReF_P**2+ImF_P**2)
     #               -4*(ReF_A*ReF_P+ImF_A*ImF_P)
     #               -  (1.0+  x)*(ReF_M**2+ImF_M**2)
     #               +4*(1.0+  x)*(ReF_P**2+ImF_P**2)
     #               -4*(1.0+2*x)*(ReF_T**2+ImF_T**2))
     #            -4*z**2*(ReF_P**2+ImF_P**2)

               a3=+(1.0+x)*(ReF_T*ReF_A+ImF_T*ImF_A)
     #            +z*(ReF_T*ReF_A+ImF_T*ImF_A
     #            +2*(ReF_T*ReF_P+ImF_T*ImF_P))

               a4=+(1.0+x+z)*(ReF_T**2+ImF_T**2)

               b0=-n_NT*x*(ReF_A*ReF_V+ImF_A*ImF_V
     #                    +ReF_A*ReF_M+ImF_A*ImF_M)
     #            +z*(ReF_T*ReF_A+ImF_T*ImF_A
     #               -2*x*(ReF_T*ReF_P+ImF_T*ImF_P)
     #            -ReF_S*ReF_V-ImF_S*ImF_V
     #            +x*(ReF_S*ReF_M+ImF_S*ImF_M))

               b1=+z*(ReF_M**2+ImF_M**2+ReF_V*ReF_M+ImF_V*ImF_M
     #               +2*(ReF_A*ReF_P+ImF_A*ImF_P))

               b2=+z*(ReF_T*ReF_P+ImF_T*ImF_P)

               c0=+0.25*(ReF_V**2+ImF_V**2
     #                  +x*(ReF_M**2+ImF_M**2)+ReF_A**2+ImF_A**2)
     #            +x*(ReF_T**2+ImF_T**2)

               c1=+ReF_T*ReF_A+ImF_T*ImF_A

               c2=+ReF_T**2+ImF_T**2

               A =+4*((x+z)*a0-r*(4*a1+  r*(a2+4*r*(a3+r*a4))))
               B =    +2*(2*b0-r*(  b1+4*r*b2))
               C =         +c0+r*(  c1+  r*c2)

               dsQESCC_dQ2= A+(s-u)*(B+(s-u)*C/mm_I)/mm_I
*              ------------------------------------------------------- *
               CASE(   2)                                                !C.H. Llewellyn Smith [1]
*              ------------------------------------------------------- *
               CALL FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                      ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
               t =-Q2/(4*mm_I)
               su= (4*m_ini*E_nu-Q2-mm_lep)/mm_ini
               A = ((mm_lep+Q2)/mm_ini)*
     #             ((ReF_A**2+ImF_A**2)*(1.0-t)
     #            -(ReF_V**2+ImF_V**2)*(1.0+t)
     #            -(ReF_M**2+ImF_M**2)*(1.0+t)*t
     #            +(ReF_T**2+ImF_T**2)*(1.0-t)*4*t
     #            -(ReF_V*ReF_M+ImF_V*ImF_M)*(4*t+mm_lep/(2*mm_ini))
     #            -mm_lep/(4*mm_ini)*
     #             (ReF_M**2+ImF_M**2+ReF_V**2
     #             +ImF_V**2+ReF_A**2+ImF_A**2
     #             +4*(ReF_P**2+ImF_P**2)
     #             +4*(ReF_A*ReF_P+ImF_A*ImF_P)
     #             -4*(1.0-t)*(ReF_S**2+ImF_S**2+ReF_P**2+ImF_P**2)))
               B =-4*t*(ReF_V*ReF_A+ImF_V*ImF_A+ReF_M*ReF_A+ImF_M*ImF_A)
     #            -(mm_lep/mm_ini)*(ReF_S*(ReF_V+t*ReF_M)
     #             +ImF_S*(ImF_V+  t*ImF_M)
     #             -ReF_T*(ReF_A+2*t*ReF_P)
     #             -ImF_T*(ImF_A+2*t*ImF_P))
               C = 0.25*(ReF_V**2+ImF_V**2+ReF_A**2+ImF_A**2
     #                 +0.25*Q2*(ReF_M**2+ImF_M**2)/mm_ini
     #                 +     Q2*(ReF_T**2+ImF_T**2)/mm_ini)
               dsQESCC_dQ2= mm_ini*(A-su*(n_NT*B-C*su))/(4*E_nu**2)
*              ------------------------------------------------------- *
               CASE(   3)                                                !E.A. Paschos and J.Y. Yu [2]
*              ------------------------------------------------------- *
               CALL FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                      ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
               dsQESCC_dQ2=
     #         (ReF_V**2*
     #          (Q2**2-4*mm_ini*(mm_lep+Q2)-mm_lep**2)/(4*mm_ini)
     #         +ReF_M**2*(4*mm_ini*(Q2**2-mm_lep**2)-Q2**2*(mm_lep+Q2))/
     #          (16*mm_ini**2)
     #         +ReF_A**2*
     #          (Q2**2+4*mm_ini*(mm_lep+Q2)-mm_lep**2)/(4*mm_ini)
     #         +ReF_P**2*mm_lep*Q2*(mm_lep+Q2)/(4*mm_ini**2)
     #         +ReF_V*ReF_M*(2*Q2**2-Q2*mm_lep+mm_lep**2)/(2*mm_ini)
     #         -ReF_A*ReF_P*mm_lep*(mm_lep+Q2)/(2*mm_ini)
     #         -n_NT*ReF_A*
     #          (ReF_V+ReF_M)*Q2*(4*m_ini*E_nu-mm_lep-Q2)/mm_ini
     #         +(ReF_V**2+ReF_M**2*Q2/(4*mm_ini)+ReF_A**2)*
     #          (4*m_ini*E_nu-mm_lep-Q2)**2/(4*mm_ini))/(4*E_nu**2)
*              ------------------------------------------------------- *
               CASE(   4)                                                !V.A. Naumov et al., OBSOLETE
*              ------------------------------------------------------- *
               y= Q2/(2*m_ini*E_nu)
               a= 0.5*mm_lep/(E_nu*m_ini)
               b= 0.5*y
               c= b*E_nu/m_ini
               d= y+a
               CALL FFCC(Q2,ReF_V,ReF_M,ReF_A,ReF_P,ReF_T,ReF_S,
     #                      ImF_V,ImF_M,ImF_A,ImF_P,ImF_T,ImF_S)
               IF (n_PT.eq.0) THEN
                 GOTO 1
                              ELSE
                 E_lep= E_nu-y*E_nu
                 P_lep= sqrt(E_lep**2-mm_lep)
                 h    = E_lep+P_lep
                 p1   = y+a*(h-2*E_nu)/h
                 p2   = 2.0-(y+a*(h+2*E_nu)/h)
                 m1   = y+a*(mm_lep-2*E_nu*h)/mm_lep
                 m2   = 2.0-(y+a*(mm_lep+2*E_nu*h)/mm_lep)
                 IF (n_NT.eq.1) THEN
                   IF (n_PT.eq.1) THEN
                     GOTO 2
                                  ELSE
                     GOTO 3
                endIF
                              ELSE
                   IF (n_PT.eq.1) THEN
                     GOTO 3
                                  ELSE
                     GOTO 2
                endIF
              endIF
            endIF
*              ------------------------------------------------------- *
*              NO POLARIZATION                                         *
*              ------------------------------------------------------- *
    1          f = 1.0
               A1= 1.0-(1.0-b+0.5*m_ini/E_nu)*d
               A2= 1.0-(1.0-b-0.5*m_ini/E_nu)*d
               A3= b*(b+(E_nu-y*E_nu)/m_ini)-
     #             a*(c-0.25*c*y+0.25*a*(1.0-c))
               A4= 0.25*mm_lep/mm_ini*y*d
               A5= d*(y-0.5*a)
               A6=-a*d
               A7=-y*(1.0-b-0.5*a)
               GOTO 4
*              ------------------------------------------------------- *
*              "CORRECT" POLARIZATION                                  *
*              ------------------------------------------------------- *
    2          f = 0.5*h/P_lep
               A1= 0.5*((E_nu-m_ini)*p1/E_nu+P_lep*p2/E_nu)
               A2= 0.5*((E_nu+m_ini)*p1/E_nu+P_lep*p2/E_nu)
               A3= b*((1.0-b)*P_lep/m_ini+b*(1.0+E_lep/m_ini))
     #               +0.5*a*(y*(1.0-0.5*c-E_nu*
     #                       (1.0+h/m_ini+P_lep*(1.0-c)/E_nu-a)/h))
               A4=-0.25*y*m1*(mm_lep/(m_ini*h))**2
               A5= (y+0.5*a)*p1-a*d*P_lep/h
               A6= a*mm_lep*m1/h**2
               A7= (E_nu+P_lep)*p1/(2*E_nu)
               GOTO 4
*              ------------------------------------------------------- *
*              "UNCORRECT" POLARIZATION                                *
*              ------------------------------------------------------- *
    3          f = 0.5*mm_lep/(P_lep*h)
               A1= 0.5*(-(E_nu-m_ini)*m1/E_nu+P_lep*m2/E_nu)
               A2= 0.5*(-(E_nu+m_ini)*m1/E_nu+P_lep*m2/E_nu)
               A3= b*((1.0-b)*P_lep/m_ini-b*(m_ini+E_lep)/m_ini)
     #               -0.5*a*(y*(1.0-0.5*c-E_nu*(1.0+mm_lep/(m_ini*h)
     #                      -P_lep*(1.0-c)/E_nu-a)*h/mm_lep))
               A4= 0.25*y*p1*(h/m_ini)**2
               A5=-(y+0.5*a)*m1-a*d*P_lep*h/mm_lep
               A6=-0.5*p1*h**2/(E_nu*m_ini)
               A7=-0.5*(E_nu-P_lep)*m1/E_nu
*              ------------------------------------------------------- *
    4          dsQESCC_dQ2= f*(A1*ReF_V**2+A2*ReF_A**2+A3*ReF_M**2
     #                        +A4*ReF_P**2+A5*ReF_V*ReF_M+A6*ReF_A*ReF_P
     #                        +n_NT*A7*2*ReF_A*(ReF_V+ReF_M))
      endSELECT

         RETURN
      END FUNCTION dsQESCC_dQ2
