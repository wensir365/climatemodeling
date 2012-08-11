      SUBROUTINE CONVEC(T1,T2,P1,P2,FH1,FH2,FC1,FC2,DZ,ITROP,cflag,Idry)

C   THIS SUBROUTINE FINDS THE CONVECTIVE LAPSE RATE BETWEEN GRID
C   POINTS J1 AND J. IT CONSIDERS THE CASE FOR HIGH CO2.

C JK  Idry is a parameter that has been added in order to compute the dry adiabatic lapse rate
c    
c-as This subroutine has been modified to take into account CO2 condensation (sept-2004)
c  Tag code: 1-Main loop
c            Units-Options IN the main loop   
c            1x-Moist adiabatic lapse rate
c            2x-CO2 condensation   
     
      INCLUDE '../INCLUDECLIM/parNT.inc'
      INCLUDE '../INCLUDECLIM/parMT.inc'
      
      INCLUDE '../INCLUDECLIM/comCBLOK.inc'
      INCLUDE '../INCLUDECLIM/comEBLOK.inc'
      INCLUDE '../INCLUDECLIM/comFBLOK.inc'
      INCLUDE '../INCLUDECLIM/comGBLOK1.inc'
      INCLUDE '../INCLUDECLIM/comCONS.inc'
      INCLUDE '../INCLUDECLIM/comSBLOK.inc'
      INCLUDE '../INCLUDECLIM/comCO2BLOK.inc'

      AMN = DM
      NDIV = 10                  !number of sublevels
      PLT = DZ
      IF(P2.LT.P1) PLT = - DZ
      DLP = PLT/FLOAT(NDIV)
      AC = 0.
      CALL SATRAT(T1,PSAT)  
      IF (IMOIST.EQ.1) FH1 = PSAT/P1
      FH1 = amin1(fh1,0.99)
      AV = 18./AMN * FH1/(1.-FH1)   !THIS VALUE OF AV IS ONLY USED BELOW 273 K
      T = T1
      P = P1
      TL = ALOG(T)
      PL1 = ALOG(P1)

      cflag=0.      !convection flag
c  Flag to indicate that CO2 has reach the saturation pressure
c (0= not saturated, 1=saturated)
      imco2 = 0     
C
C-KK	Begin Sublevel Integration
      DO 1 L=1,NDIV
c Calculation of heat capacities (cal/mol/K)
      CPCO2 = 7.7 + 5.3E-3*T - 8.3E-7*T*T
      CPN2 = 6.76 + 6.06E-4*T + 1.3E-7*T*T
      CPO2 = 8.27 + 2.58E-4*T - 1.877E5/T/T
      CPO2 = AMAX1(CPO2,CPN2)
      CPCH4 = 8.3
      CPN = FC1*CPCO2 + FN2*CPN2 + FO2*CPO2 + FAR*4.97 + FCH4*CPCH4
      CALL SATCO2(T,PSCO2)
      FCSAT = PSCO2/P     
c Once CO2 reach the saturation pressure it stays saturated up to JTROP
      if(FCSAT.lt.FC1) imco2=1   
      
      PL = PL1 + (L-1)*DLP
      P = EXP(PL)
      TC = T - 273.15
      If (Idry.eq.1) go to 40
      if(imco2.eq.0) go to 3

c   Moist CO2 adiabat
      cflag = 3. 
      N = (TC + 130.)/5. + 1
      N = MAX0(N,1)
      IF (TC .GT. -56.595) N = (TC + 130.)/5. + 3

c-mm  Next line added to allow for surface temps above 303K, since this is
c-mm  the limit of the table.
      IF (TC.GT.29.) N = MT - 1
      N1 = N + 1
      FRC = (TC - TCTAB(N))/(TCTAB(N1) - TCTAB(N))
      BET = FRC*BETASC(N1) + (1.-FRC)*BETASC(N)
c-mm  Treat CO2 as ideal if above 303K
      IF (TC.GT.29.) BET = 1.
      AVC = 44./AMN * BET * FCSAT/(1. - FCSAT)
      IF (L.EQ.1) BETAC1 = BET
      TSQ = T*T
      IF (T .LT. 216.56) GO TO 21
C   Lapse rate in pure CO2 over liquid CO2
      if (FCO2.gt.0.9) then
      DLPVLT = 2.303*T*(867.2124/TSQ + 18.65612E-3 - 2.*72.4882E-6*T
     &  + 3.*93.E-9*TSQ)
      endif
      GO TO 22
C   Lapse rate in pure CO2 over solid CO2
   21 if(FCO2.gt.0.9) then
       T47 = T - 4.718
       DLPVLT = 2.303*T*(1284.07/T47/T47 + 1.256E-4)
      endif
C   Note that, entropy for CO2 is expressed in units of cal/g-K,
C   so that the conversion factor 4.184 disappears from
C   the expression for DLADLT
   22 DSVDLT = FRC*DSVC(N1) + (1.-FRC)*DSVC(N)
      SVCM = FRC*SVSC(N1) + (1.-FRC)*SVSC(N)
      BET = FRC*BETASC(N1) + (1.-FRC)*BETASC(N)
      DLRVLT = FRC*DRCVAP(N1) + (1.-FRC)*DRCVAP(N)
      CVN = CPN - R
      SUM1 = (R*DLRVLT - CVN)/AMN - AVC*DSVDLT
      SUM2 = AVC*SVCM + R/AMN
      DLADLT = SUM1/SUM2
      DLPDLT = DLPVLT- DLADLT/(1. + AVC*AMN/44./BET)
C
      GO TO 2
c   WATER MOIST ADIABAT

   3   CALL SATRAT(T,PSAT)
       cflag =1.
      IF (ITROP .EQ. 0) THEN 
	DLPDLT = CPN/R
        cflag = 2.
        GO TO 2
      END IF
      IF (IMOIST.EQ.1) GO TO 11
C
      IF (L.EQ.1) BET = BETA2
      FP = 1./(1. + BET*(1.-FVDRY)/FVDRY)
      PH2O = FP * P
      IF (PH2O.LT.PSAT) GO TO 12
  11  IMOIST = 1
      IF(T.LT.273.16) GO TO 13
C
C   INGERSOLL'S FORMULATION BETWEEN 0 AND 374 C
      N = TC/5. + 1
      N1 = N + 1
      FR = (TC - TTAB(N))/(TTAB(N1) - TTAB(N))
      DLPVLT = FR*DPVAP(N1) + (1.-FR)*DPVAP(N)
      DSVDLT = FR*DSV(N1) + (1.-FR)*DSV(N)
      SVCM = FR*SVC(N1) + (1.-FR)*SVC(N)
      PV = FR*PVAP(N1) + (1.-FR)*PVAP(N)
      RV = FR*RHOV(N1) + (1.-FR)*RHOV(N)
      DLRVLT = FR*DRHOV(N1) + (1.-FR)*DRHOV(N)
      BET = 41.84*RV*R*T/(18.*PV)
      CVN = CPN - R
      IF (L.NE.1) GO TO 14
      AV = 18./AMN * BET * FH1/(1.-FH1)
      BETA1 = BET
  14  CONTINUE
C
      SUM1 = 4.184*(R*DLRVLT - CVN)/AMN - AV*DSVDLT
      SUM2 = AV*SVCM + 4.184*R/AMN
      DLADLT = SUM1/SUM2
      DADLT = AV*DLADLT
      SUM3 = 1. + DLRVLT - DLADLT
      FAC = BET*18./AV/AMN
      DLPDLT = PV/P * (DLPVLT + FAC*SUM3)
      GO TO 2
C
C   DRY ADIABAT ABOVE 374 C
  12  DADLT = 0.
      DLADLT = 0.
      AV = 18./AMN * FVDRY/(1.-FVDRY)
      TC = T - 273.15
      IMOIST = 0
      PDRY = P
      N = TC/10. + 2
      IF (TC.GT.600.) N = (TC-600.)/100. + 62
      PV = FH1 * P
      M = PV/5. + 1
      M1 = M - 1
      IF (M.EQ.1) M1 = M
C
      FP = 0.
      IF (M.GT.1) FP = (PV - PCP(M))/(PCP(M) - PCP(M1))
      FT = (TC - TCP(N))/(TCP(N+1) - TCP(N))
      DLPVLT = (1.+FP-FT)*DPDTL(M,N) - FP*DPDTL(M1,N) + FT*
     &  DPDTL(M,N+1)
      BET = (1.+FP-FT)*BETAM(M,N) - FP*BETAM(M1,N) + FT*BETAM(M,N+1)
      IF (L.EQ.1) BETA1 = BET
      DLPDLT = FH1*DLPVLT + (1.-FH1)*CPN/R
      GO TO 2
C
C   CLAUSIUS-CLAPEYRON EQUATION BELOW 0 C
  13  HL = SUBL
      HLT = HL/T
      DDLT = 0.
      CC = 0.5
      SUM1 = (HLT*18. - CPN)/AMN - AC*CC
      SUM2 = DDLT + CC - HLT
      SUM3 = HLT + R/AV/AMN
      DADLT = (SUM1 - AV*SUM2)/SUM3
      DLPDLT = 18.*HLT/R - DADLT/AV/(1. + AV*AMN/18.)
      BET = 1.
      GO TO 2
C
C Dry adiabat
  40  DLPDLT = CPN*AMN/4.184*R
C
   2  DLT = DLP/DLPDLT     
      TL = TL + DLT
      T = EXP(TL)

      if(imco2.eq.0)then
      DAV = DADLT*DLT
      AV = AV + DAV
      endif

   1  CONTINUE       !End Sublevel integration
      T2 = T
      
      FH2 = 1./(1. + BET*18./AMN/AV)
      BETA2 = BET
      CALL SATRAT(T2,PSAT)
      IF (PSAT.GT.P2) GOTO 30
      FH2 = RELHUM(P2) * PSAT/P2
         
  30  IF(IMCO2.EQ.1) THEN
       BETA2 = 1.
       BETAC2 = BET
      ELSE
       BETA2=BET
      ENDIF  
   
      FC2 = FCO2
      CALL SATCO2(T2,PSCO2)
      FCSAT = PSCO2/P2
      IF (IMCO2.EQ.1) FC2 = FCSAT
      
      RETURN
      END
