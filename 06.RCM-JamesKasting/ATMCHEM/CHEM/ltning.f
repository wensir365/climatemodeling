
      SUBROUTINE LTNING(FO2,FCO2,P0CGS)
      
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      
      INCLUDE '../INCLUDECHEM/comBBLOK.inc'
      INCLUDE '../INCLUDECHEM/comLTBLOK.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'

C     THIS SUBROUTINE CALCULATES LIGHTNING PRODUCTION RATES FOR O2
C     AND NO IN AN N2-O2-CO2-H2O ATMOSPHERE USING THE EQUATIONS IN
C     THESIS APPENDIX C.
C
C     EQUILIBRIUM CONSTANTS AT 3500 K
      AK1 = .103
      AK2 = .619
      AK3 = 5.3
      AK4 = .22
C
      P0 = P0CGS/1.013E6
      FN2 = 1. - FO2 - FCO2
      PN2 = FN2*P0
      PO2 = FO2*P0
      PCO2 = FCO2*P0
      PH2O = USOL(LH2O,1)*P0
      PH2 = USOL(LH2,1)*P0
      PCO = USOL(LCO,1)*P0
C
      O2T = PO2 + PCO2 + 0.5*(PH2O + PCO)
      H2T = PH2 + PH2O
      CT = PCO2 + PCO
      ALPHA = AK2*SQRT(O2T)
      BETA = AK3*SQRT(O2T)
      A = (AK1*SQRT(PN2) + AK4)/(2.*SQRT(O2T))
      B = 0.5*CT/O2T
      C = 0.5*H2T/O2T
C
C     INITIAL GUESS FOR XO2 AT 3500 K
      X = 0.1 + 0.9*PO2/O2T + 0.2*PCO2/O2T
C
C     NEWTON STEP
      DO 1 N=1,20
      NS = N
      XS = X
      X2 = SQRT(X)
      FX = X + A*X2 - B/(1.+ALPHA*X2) - C/(1.+BETA*X2) + 2.*B + C - 1.
      FPX = 1. + (A + ALPHA*B/(1.+ALPHA*X2)**2 + BETA*C/(1.+BETA*X2)
     2  **2)/(2.*X2)
      X = X - FX/FPX
      ERR = ABS((X-XS)/X)
      IF(ERR.LT.1.E-5) GO TO 2
   1  CONTINUE
   2  PO2 = X*O2T
      PNO = AK1*SQRT(PN2*PO2)
C
C     SCALE AGAINST ESTIMATED COLUMN PRODUCTION OF NO IN THE PRESENT-
C     DAY TROPOSPHERE.  DISTRIBUTE PRODUCTION OVER LOWEST 6 KM.
C
      PRNOX = 3.E9/7.942E5
      PNONOW = 3.574E-2
      ZAPNO = PRNOX*PNO/PNONOW

c ZAPO2 is only needed in high CO2 atmospheres
c      ZAPO2 = ZAPNO*PO2/PNO
c      write(90,100)PO2,PNO,ZAPO2,ZAPNO
c 100  FORMAT(/1X,'PO2=',1PE10.3,2X,'PNO=',1PE10.3,2X,'ZAPO2=',
c     2  1PE10.3,2X,'ZAPNO=',1PE10.3)
      write(90,100)PO2,PNO,ZAPNO
 100  FORMAT(/1X,'PO2=',1PE10.3,2X,'PNO=',1PE10.3,2X,'ZAPNO=',1PE10.3)
      write(90,101)NS,X,ERR
 101  FORMAT(1X,'N=',I2,5X,'X=',1PE12.5,2X,'ERR=',1PE12.5)
      RETURN
      END
      FUNCTION SIGRAY(W)
      W1 = 1.E-4 * W
      W2 = W1 * W1
      W4 = W2 * W2
      SIGRAY = 4.006E-28*(1. + .0113/W2 + .00013/W4)/W4
      RETURN
      END

      FUNCTION E1(X)
      COMMON/EBLOK1/A(6),B(4),C(4)
      E1 = 0.
      IF(X.EQ.0.) RETURN
C
      X2 = X*X
      X3 = X2*X
      X4 = X3*X
      IF(X.GT.1) GO TO 1
      X5 = X4*X
      E1 = -ALOG(X) + A(1)*X + A(2)*X2 + A(3)*X3 + A(4)*X4 + A(5)*X5
     2  + A(6)
      RETURN
C
   1  SUM1 = X4 + B(1)*X3 + B(2)*X2 + B(3)*X + B(4)
      SUM2 = X4 + C(1)*X3 + C(2)*X2 + C(3)*X + C(4)
      E1 = EXP(-X)*SUM1/SUM2/X
      RETURN
      END
