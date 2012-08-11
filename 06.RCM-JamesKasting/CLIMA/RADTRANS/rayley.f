
      SUBROUTINE RAYLEY(SIGR,SIGRUV)
C
* Subroutine modified to calculate Rayleigh scattering in the near UV
* wavelenghts (Antigona Segura, 08/2005)
      INCLUDE '../INCLUDECLIM/parND.inc'
      INCLUDE '../INCLUDECLIM/parNSOL_NSOLUV.inc'
C
       INCLUDE '../INCLUDECLIM/comSOLARBLK.inc'
       INCLUDE '../INCLUDECLIM/comCBLOK.inc'
       INCLUDE '../INCLUDECLIM/comCONS.inc'

      DIMENSION A(3),B(3),DEL(3),SIG(3),SIG1(3)
      dimension SIGR(NSOL),SIGRUV(NSOLUV)

C
      DATA A/29.06,26.63,43.9/
      DATA B/7.7,5.07,6.4/
      DATA DEL/.0305,.054,.0805/

**** NEAR UV ******
      DO 1120 I=1,NSOLUV
         ALUV2=XLAMBUV(I)*XLAMBUV(I)
         ALUV4=ALUV2*ALUV2
C
         DO 1125 J=1,3
            PAREN1=(1.E-5*A(J)*(1.+1.E-3*B(J)/ALUV2))**2
            SMDEL1=(6.+3.*DEL(J))/(6.-7.*DEL(J))
            SIG1(J)=4.577E-21*SMDEL1*PAREN1/ALUV4
 1125    CONTINUE
C
         SIGRUV(I) = FCO2*SIG1(3) + (1.-FCO2) * (FO2*SIG1(2) +
     &     (1.-FO2)*SIG1(1))
 1120 CONTINUE 

***** VISIBLE AND NEAR IR *******
      DO 1135 I=1,NSOL
         AL2=XLAMBDA(I)*XLAMBDA(I)
         AL4=AL2*AL2
C
         DO 1140 J=1,3
            PAREN=(1.E-5*A(J)*(1.+1.E-3*B(J)/AL2))**2
            SMDEL=(6.+3.*DEL(J))/(6.-7.*DEL(J))
            SIG(J)=4.577E-21*SMDEL*PAREN/AL4
 1140    CONTINUE
C
         SIGR(I) = FCO2*SIG(3) + (1.-FCO2) * (FO2*SIG(2) +
     &     (1.-FO2)*SIG(1))
 1135 CONTINUE 
C
      RETURN
      END  
