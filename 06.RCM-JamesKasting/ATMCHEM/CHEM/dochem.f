
      SUBROUTINE DOCHEM(FVAL,N)
      
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/parNF.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'

      DIMENSION FVAL(NQ,NZ),XP(NZ),XL(NZ),D(NSP2,NZ)

      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comBBLOK.inc'
      INCLUDE '../INCLUDECHEM/comCBLOK.inc'
      INCLUDE '../INCLUDECHEM/comDBLOK.inc'
      INCLUDE '../INCLUDECHEM/comGBLOK.inc'
      INCLUDE '../INCLUDECHEM/comLTBLOK.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'
      INCLUDE '../INCLUDECHEM/comZBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRRATS.inc'
      INCLUDE '../INCLUDECHEM/comSATBLK.inc'
      INCLUDE '../INCLUDECHEM/comSULBLK.inc'

      SRAIN = 1.E-4
      DO 1 I=1,NQ
      DO 1 J=1,NZ
   1  D(I,J) = USOL(I,J) * DEN(J)
C
      DO 2 J=1,NZ
      D(LSO4AER,J) = SO4AER(J) * DEN(J)
      D(NSP-2,J) = O2(J) * DEN(J)
      D(NSP-1,J) = CO2(J) * DEN(J)
      D(NSP,J) = (1. - O2(J) - CO2(J)) * DEN(J)
      D(NSP,J) = AMAX1(D(NSP,J),1.E-60)
      D(NSP1,J) = 1.
   2  D(NSP2,J) = DEN(J) 
C-AP   2  D(NSP1,J) = 1.
C-AP  Added M and SO4AER 
C
C ***** SOLVE FOR THE PHOTOCHEMICAL EQUILIBRIUM SPECIES *****
C
      NQT1 = NQT + 1
      NSP3 = NSP - 4 
      DO 3 I=NQT1,NSP3
      CALL CHEMPL(D,XP,XL,I)
      DO 3 J=1,NZ
   3  D(I,J) = XP(J)/XL(J)

C ***** LONG-LIVED SPECIES CHEMISTRY *****
      DO 4 I=1,NQ
      CALL CHEMPL(D,XP,XL,I)
      DO 4 J=1,NZ
      XLJ = XL(J) + RAINGC(I,J)
      FVAL(I,J) = XP(J)/DEN(J) - XLJ*USOL(I,J)
      YP(I,J) = XP(J)
   4  YL(I,J) = XLJ
C-AP ***** TRIDIAGONAL SPECIES (H2SO4 aerosol) *****
      CALL CHEMPL(D,XP,XL,LSO4AER)
      DO 16 J=1,NZ
      YL(LSO4AER,J) = XL(J) + RAINGC(LH2SO4,J)
  16  YP(LSO4AER,J) = XP(J)
C-AP ***********************************************
C   ZERO OUT H2O TERMS IN THE TROPOSPHERE AND INCLUDE LIGHTNING
C   PRODUCTION OF NO, O2, AND CO.  (MUST INCLUDE CO IN ORDER TO
C   BALANCE THE HYDROGEN BUDGET)
C-AS This is a leftover of the low O2, high CO2 version of the 
c    photochemical code, all the commented lines on this do loop
c    are not nedeed for this highO2 version of the code
      DO J=1,JTROP
      FVAL(LH2O,J) = 0.
      SCALE = RAIN(J)/RAIN(1)
      ZAP = ZAPNO * SCALE
      FVAL(LNO,J) = FVAL(LNO,J) + ZAP/DEN(J)
      YP(LNO,J) = YP(LNO,J) + ZAP
c     ZAP = ZAPO2 * SCALE
c     FVAL(LO2,J) = FVAL(LO2,J) + ZAP/DEN(J)
c     YP(LO2,J) = YP(LO2,J) + ZAP
c      ZAP = (ZAPNO + 2.*ZAPO2)*SCALE
c      FVAL(LCO,J) = FVAL(LCO,J) + ZAP/DEN(J)
c      YP(LCO,J) = YP(LCO,J) + ZAP
      enddo
C
C  VOLCANIC OUTGASSING OF SO2
C    Note: PRSO2 = 3.5E9, distributed over lower 10 km.
C        (JTEN=5 covers five 2-km levels on Mars)
C-AP Value from Archean we divide by 3      PRSO2 = 3.5E9/1.E6
      PRSO2 = 3.5E9/1.E6
      JTEN = 1.E6/DELZ + 0.01
      DO 351 J=1,JTEN
      FVAL(LSO2,J) = FVAL(LSO2,J) + PRSO2/DEN(J)
 351  YP(LSO2,J) = YP(LSO2,J) + PRSO2
C
C  VOLCANIC OUTGASSING OF H2S 
C    Note: PRH2S = 3.5E8, distributed over lower 10 km.
C        (JTEN=5 covers five 2-km levels on Mars)
C-AP Value from Archean we divide by 3   PRH2S = 3.5E8/1.E6
      PRH2S = 3.5E8/1.E6
      JTEN = 1.E6/DELZ + 0.01
      DO 352 J=1,JTEN
      FVAL(LH2S,J) = FVAL(LH2S,J) + PRH2S/DEN(J)
 352  YP(LH2S,J) = YP(LH2S,J) + PRH2S
C
C   H2O CONDENSATION IN THE STRATOSPHERE
C   (RHCOLD IS THE ASSUMED RELATIVE HUMIDITY AT THE COLD TRAP
      RHCOLD = 0.1
      JT1 = JTROP + 1
      CONFAC = 1.6E-5
      DO 13 J=JT1,NZ
      H2OCRT = RHCOLD * H2OSAT(J)
      IF (USOL(LH2O,J) .LT. H2OCRT) GO TO 13
      CONDEN(J) = CONFAC * (USOL(LH2O,J) - H2OCRT)
      FVAL(LH2O,J) = FVAL(LH2O,J) - CONDEN(J)
  13  CONTINUE
C
C   H2SO4 CONDENSATION
      DO 14 J=1,NZ
      CONSO4(J) = CONFAC * (USOL(LH2SO4,J) - H2SO4S(J))
C-AP
      CONSO4(J) = AMAX1(CONSO4(J),0.)
C-AP
      FVAL(LH2SO4,J) = FVAL(LH2SO4,J) - CONSO4(J)
      YL(LH2SO4,J) = YL(LH2SO4,J) + CONFAC
C-AP For now comment out ficticious production 
      YP(LH2SO4,J) = YP(LH2SO4,J) + CONFAC*H2SO4S(J)*DEN(J)
      YP(LSO4AER,J) = YP(LSO4AER,J) + CONSO4(J)*DEN(J)
  14  CONTINUE
C
      DO  J=1,NZ
      O3(J) = D(LO3,J)/DEN(J)
      ENDDO                                        
C-AP
C   SPECIAL TREATMENT FOR H2O AND NO
      DO 5 J=1,NZ
      IF(Z(J).GT.ZTROP) GO TO 7
      FVAL(LH2O,J) = 0.
      YP(LH2O,J) = 0.
      YL(LH2O,J) = 1.E-99
   7  CONTINUE
      IF(Z(J).GT.10.E5) GO TO 5
      ZAP = ZAPNO*RAIN(J)/RAIN(1)
      FVAL(LNO,J) = FVAL(LNO,J) + ZAP/DEN(J)
      YP(LNO,J) = YP(LNO,J) + ZAP
   5  CONTINUE
      IF(N.LT.1) RETURN
C
C ***** CALCULATE COLUMN-INTEGRATED PRODUCTION AND LOSS *****
      O3COL = 0.
      H2SCOL = 0.
      SO2COL = 0.
      DO 10 L=1,NR
  10  RAT(L) = 0.
      DO 11 K=1,NQT
      SR(K) = 0.
      TP(K) = 0.
      TL(K) = 0.
  11  CONTINUE
C
      DO 6 J=1,NZ
      RELH(J) = USOL(LH2O,J)/H2OSAT(J) 
      H2SCOL = H2SCOL + D(LH2S,J)*DZ
      SO2COL = SO2COL + D(LSO2,J)*DZ
   6  O3COL = O3COL + USOL(7,J)*DEN(J)*DZ
C
      DO 12 L=1,NR
      M = JCHEM(1,L)
      K = JCHEM(2,L)
      DO 12 J=1,NZ
  12  RAT(L) = RAT(L) + AR(L,J)*D(M,J)*D(K,J)*DZ
C
      DO 8 I=1,NQT
      XLG(I) = YL(I,1)
      DEQ(I) = YP(I,1)/(YL(I,1) + 1.E-99)
      DO 8 J=1,NZ
      TP(I) = TP(I) + YP(I,J)*DZ
      TL(I) = TL(I) + YL(I,J)*D(I,J)*DZ
   8  SR(I) = SR(I) + RAINGC(I,J)*D(I,J)*DZ
C
C ***** SAVE THESE DENSITIES FOR PRINTOUT *****
      NQT1 = NQT + 1
C-AP
      DO 9 I=NQT1,NSP
      DO 9 J=1,NZ
   9  SL(I,J) = D(I,J)
      RETURN
      END
