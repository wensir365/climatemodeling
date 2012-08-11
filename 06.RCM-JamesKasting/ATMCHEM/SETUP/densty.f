
      SUBROUTINE DENSTY(FO2,FCO2,P0)
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comPRESS1.inc'
C
      G0 = G 
      RGAS = 8.3143E7
      BK = 1.38054E-16
      R0 = 6.371E8
      FT = FO2 + FCO2 + 0.01
      WT = FO2*32. + FCO2*44. + (1.-FT)*28. + 0.4
      ROVERM = RGAS/WT
      PG = 1.013E6
      P0 = PG
C     P0 = PG GIVES YOU A ONE BAR ATMOSPHERE
      DZ = Z(2) - Z(1)
      T0 = T(1) + (T(1)-T(2))/2.
      HA = ROVERM*0.5*(T0 + T(1))/G0
      P1 = P0 * EXP(-0.5*DZ/HA)
      DEN(1) = P1/(BK*T(1))
C
C ***** FIND DENSITY FROM HYDROSTATIC EQUILIBRIUM *****
      DO 1 I=2,NZ
      DZ = Z(I) - Z(I-1)
      R = R0 + Z(I)
      G = G0 * (R0/R)*(R0/R)
      TAV = 0.5*(T(I) + T(I-1))
      HA = ROVERM*TAV/G
   1  DEN(I) = DEN(I-1)*EXP(-DZ/HA)*T(I-1)/T(I)
C****** FIND PRESSURE FROM THIS DENSITY ****************
      DO 334 I=1,NZ
        PRESS(I)=DEN(I)*BK*T(I)
 334  CONTINUE
C
      RETURN
      END
