
      SUBROUTINE RATES

       INCLUDE '../INCLUDECHEM/parNZ.inc'
       INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE '../INCLUDECHEM/parNR.inc'
       INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE '../INCLUDECHEM/parNMAX.inc'

      DIMENSION A3(NZ),A4(NZ),A5(NZ),A6(NZ),A10(NZ),A11(NZ),A12(NZ),
     2  A13(NZ),A14(NZ),A15(NZ),A16(NZ),A17(NZ),A18(NZ),A19(NZ),
     3  A20(NZ),A21(NZ),A22(NZ),A30(NZ),A31(NZ),A32(NZ),A41(NZ),
     4  A43(NZ),A46(NZ),A47(NZ),A48(NZ)
      
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comCBLOK.inc'     
      INCLUDE '../INCLUDECHEM/comGBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRBLOK.inc'

C
C ***** TEMPERATURE-DEPENDENT RATE COEFFICIENTS *****
C
      DO 1 I=1,NZ
      PATM = DEN(I)*1.38E-16*T(I)/1.013E6
      A3(I) = 3.E-14*T(I)*EXP(-4480./T(I))
      A4(I) = 5.5E-12*EXP(-2000./T(I))
      A5(I) = 1.4E-10*EXP(-470./T(I))
      A10(I) = 2.2E-11*EXP(120./T(I))
      A12(I) = 1.6E-12*EXP(-940./T(I))
      A13(I) = 3.0E-11*EXP(200./T(I))
      A14(I) = 1.1E-14*EXP(-500./T(I))
      A15(I) = 2.3E-13*EXP(600./T(I)) + 1.7E-33*EXP(1000./T(I))
     2  *DEN(I)
      A16(I) = 2.9E-12*EXP(-160./T(I))
      A17(I) = 2.76E-34*EXP(710./T(I))*DEN(I)
      A19(I) = 8.0E-12 * EXP(-2060./T(I))
      A20(I) = 4.2E-12*EXP(-240./T(I))
      A21(I) = 1.8E-11*EXP(110./T(I))
      A22(I) = 3.2E-11*EXP(70./T(I))
      A30(I) = 1.5E-13 * (1. + 0.6*PATM)
      A31(I) = 6.5E-33*EXP(-2180./T(I))*DEN(I)
      A32(I) = 2.0E-33*EXP(-850./T(I))*DEN(I)
      A41(I) = 2.8E-11*EXP(-1540./T(I))
      A43(I) = 2.6E-33*EXP(375./T(I))*DEN(I)
      A46(I) = 6.1E-26/T(I)/T(I) *DEN(I)
      A48(I) = 3.4E-11*EXP(-1600./T(I))
   1  CONTINUE
C
C ***** THREE-BODY COEFFICIENTS *****
      DO 3 I=1,NZ
      TT = T(I)
      DN = DEN(I)
C     AR(J,I) = TBDY(K0,KI,N,M,T,DEN)
C
C  HO2
      AR(6,I) = TBDY(5.7E-32,7.5E-11,1.6,0.,TT,DN)
C
C  O3
      AR(18,I) = TBDY(6.E-34,1.E-10,2.3,0.,TT,DN)
C
C  H2O2
      AR(47,I) = TBDY(6.9E-31,1.5E-11,0.8,0.,TT,DN)
C
C   CH3O2
      AR(67,I) = TBDY(4.5E-31,1.8E-12,3.,1.7,TT,DN)
C
C   NO2
      AR(84,I) = TBDY(9.E-32,3.E-11,1.5,0.,TT,DN)
C
C   HNO2
      AR(86,I) = TBDY(7.E-31,1.5E-11,2.6,0.5,TT,DN)
C
C   HNO3
      AR(88,I) = TBDY(2.6E-30,2.4E-11,3.2,1.3,TT,DN)
C
C   HO2NO2
      AR(91,I) = TBDY(1.8E-31,4.7E-12,3.2,1.4,TT,DN)
C
C   NOCL
      AR(114,I) = TBDY(9.E-32,1.E-10,1.6,0.,TT,DN)
C
C   CLONO
      AR(115,I) = TBDY(1.3E-30,1.E-10,2.,1.,TT,DN)
C
C   CLO2
      AR(117,I) = TBDY(2.7E-33,1.E-10,1.5,0.,TT,DN)
c
C   CLONO2
      AR(122,I) = TBDY(1.8E-31,1.5E-11,3.4,1.9,TT,DN)
C
C   CL2O2
      AR(141,I) = TBDY(1.9E-32,7.E-12,3.9,0.,TT,DN)
C
C   NO3
      AR(150,I) = TBDY(9.E-32,2.2E-11,2.,0.,TT,DN)
C
C   N2O5
      AR(152,I) = TBDY(2.2E-30,1.5E-12,3.9,0.7,TT,DN)
C
C   HSO3
      AR(163,I) = TBDY(3.E-31,1.5E-12,3.3,0.,TT,DN)   
c   N2O
      AR(218,I) = TBDY(3.5E-37,1.E-10,0.6,0.,TT,DN)                          

   3  CONTINUE
C
C ***** FILL UP RATE MATRIX *****
      DO 4 I=1,NZ
      AR(1,I) = 2.2E-10
      AR(2,I) = 1.0E-10
      AR(3,I) = A3(I)
      AR(4,I) = A4(I)
      AR(5,I) = A5(I)
      AR(7,I) = 8.1E-11 * 0.08
      AR(8,I) = 8.1E-11 * 0.02
      AR(9,I) = 8.1E-11 * 0.90
      AR(10,I) = A10(I)
      AR(11,I) = 4.8E-11*EXP(250./T(I))
      AR(12,I) = A12(I)
      AR(13,I) = A13(I)
      AR(14,I) = A14(I)
      AR(15,I) = A15(I)
      AR(16,I) = A16(I)
      AR(17,I) = A17(I)
      AR(19,I) = A19(I)
      AR(20,I) = A20(I)
      AR(21,I) = A21(I)
      AR(22,I) = A22(I)
      AR(30,I) = A30(I)
      AR(31,I) = A31(I)
      AR(32,I) = A32(I)
      AR(33,I) = 1.2E-10
      AR(34,I) = 0.
      AR(35,I) = 5.0E-11
      AR(36,I) = 1.0E-10
      AR(37,I) = 1.0E-10
      AR(40,I) = 1.0E-2
      AR(41,I) = A41(I)
      AR(43,I) = A43(I)
      AR(44,I) = 3.5E-12 * EXP(140./T(I))
      AR(45,I) = 1.E-11
      AR(46,I) = A46(I)
      AR(48,I) = A48(I)
      AR(49,I) = 1.4E-12 * EXP(-2000./T(I))
      AR(58,I) = 2.9E-12 * EXP(-1820./T(I))
      AR(59,I) = 1.4E-10
      AR(60,I) = 1.4E-11
      AR(61,I) = 1.9E-12
      AR(62,I) = 3.E-11
      AR(63,I) = 5.E-13
      AR(64,I) = 5.E-14
      AR(65,I) = 5.E-14
      AR(66,I) = 1.5E-12
      AR(68,I) = 1.E-10
      AR(69,I) = 1.1E-10
      AR(70,I) = 5.4E-12 * EXP(-220./T(I))
      AR(71,I) = 3.8E-13 * EXP(800./T(I))
      AR(72,I) = 2.5E-13 * EXP(190./T(I))
      AR(73,I) = 4.2E-12 * EXP(180./T(I))
      AR(74,I) = 3.9E-14*EXP(-900./T(I))
      AR(75,I) = 1.E-14
   4  CONTINUE
C
      DO 5 I=1,NZ
      AR(76,I) = 5.E-11
      AR(77,I) = 6.7E-11
      AR(78,I) = 4.9E-11
      AR(79,I) = 1.5E-11*EXP(-3600./T(I))
      AR(80,I) = 1.E-16
      AR(81,I) = 5.3E-11
      AR(82,I) = 3.4E-11
      AR(83,I) = 2.0E-12 * EXP(-1400./T(I))
      AR(85,I) = 3.7E-12 * EXP(250./T(I))
      AR(87,I) = 6.5E-12 * EXP(120./T(I))
      AR(89,I) = 4.8E-10 * EXP(-340./T(I))
      AK0 = 7.2E-15*EXP(785./T(I))
      AK2 = 4.1E-16*EXP(1440./T(I))
      AK3M = 1.9E-33*EXP(725./T(I))*DEN(I)
      AR(90,I) = AK0 + AK3M/(1. + AK3M/AK2)
      AR(92,I) = 1.3E-12 * EXP(380./T(I))
      AR(93,I) = 7.8E-11 * EXP(-3400./T(I))
      AR(94,I) = AR(91,I)/(2.33E-27*EXP(10870./T(I)))
      AR(96,I) = 3.8E-12 * EXP(200./T(I))
      AR(97,I) = AR(11,I)
      AR(98,I) = 1.2E-13*EXP(-2450./T(I))
      AR(99,I) = 4.5E-14*EXP(-1260./T(I))
      AR(100,I) = 1.E-11
      AR(102,I) = 1.5E-11*EXP(170./T(I))
      AR(103,I) = 2.3E-11
      AR(104,I) = 2.1E-12*EXP(-1150./T(I))
      AR(105,I) = 2.9E-11*EXP(-260./T(I))
      AR(106,I) = 3.7E-11*EXP(-2300./T(I))
      AR(107,I) = 1.1E-11*EXP(-1400./T(I))
      AR(108,I) = 3.3E-11*EXP(-1250./T(I))
      AR(109,I) = 8.1E-11*EXP(-30./T(I))
      AR(110,I) = 1.1E-11*EXP(-980./T(I))
      AR(111,I) = 1.8E-11*EXP(170./T(I))
      AR(112,I) = 4.1E-11*EXP(-450./T(I))
      AR(113,I) = 6.8E-12*EXP(160./T(I))
      AR(116,I) = 5.8E-11*EXP(100./T(I))
      AR(118,I) = 2.3E-10
      AR(119,I) = 1.2E-11
      AR(120,I) = 3.0E-11*EXP(70./T(I))
      AR(121,I) = 6.4E-12*EXP(290./T(I))
      AR(123,I) = 4.8E-13*EXP(700./T(I))
      AR(124,I) = 1.1E-11*EXP(120./T(I))
      AR(125,I) = 2.6E-12*EXP(-350./T(I))
      AR(126,I) = 3.E-12*EXP(-500./T(I))
      AR(127,I) = 1.2E-12*EXP(-330./T(I))
      AR(128,I) = 1.E-11*EXP(-3300./T(I))
      AR(129,I) = 1.E-11*EXP(-2200./T(I))
      AR(130,I) = 2.9E-12*EXP(-800./T(I))
      AR(131,I) = 1.4E-12*EXP(-900./T(I))
      AR(139,I) = AR(117,I)/(2.43E-25*EXP(2979./T(I)))
      AR(140,I) = 4.1E-12
   5  CONTINUE
C
      DO 6 I=1,NZ
      AR(143,I) = AR(141,I)/(3.E-27*EXP(8450./T(I)))
      AR(144,I) = AR(117,I)/(5.7E-25*EXP(2500./T(I)))
      AR(145,I) = 2.6E-11
      AR(146,I) = 3.E-12*EXP(-130./T(I))
      AR(147,I) = 4.E-13
      AR(148,I) = 3.5E-14
      AR(149,I) = 2.5E-12*EXP(-950./T(I))
      AR(154,I) = AR(152,I)/(4.E-27*EXP(10930./T(I)))
      AR(155,I) = 2.E-19
C  AR(155,I) should be viewed as a tuning parameter. The gas phase upper 
C  limit is 2.e-21. A value 100 times higher than this yields an effective
C  first-order rate of ~2.e-5 in the lower stratosphere.
 6    CONTINUE
C-AP
C-AP Adding sulfur rate constants
C
C-AP
      DO 7 I=1,NZ
      AR(159,I) = 2.4E-13 * EXP(-2370./T(I))
      AR(160,I) = 2.8E-11
      AR(161,I) = 6.0E-31 * DEN(I)
      AR(162,I) = 8.6E-11
      AR(164,I) = 3.4E-32 * EXP(-1130./T(I)) * DEN(I)
      AR(165,I) = 6.0E-15
      AR(166,I) = 1.3E-12 *EXP(-330./T(I))
      AR(167,I) = 1.0E-11
      AR(168,I) = 1.0E-11
      AR(169,I) = 1.0E-11
      AR(170,I) = 6.0E-12 * EXP(-75./T(I))
      AR(171,I) = 1.3E-11 * EXP(-860./T(I))
      AR(172,I) = 9.2E-12 * EXP(-1800./T(I))
      AR(173,I) = 1.6E-10
      AR(174,I) = 4.0E-19
      AR(175,I) = 3.0E-11
      AR(176,I) = 1.2E-11
      AR(177,I) = 5.0E-11
      AR(178,I) = 1.0E-11
      AR(179,I) = 2.2E-11 * EXP(120./T(I))
      AR(180,I) = 2.3E-12
      AR(181,I) = 6.6E-11
C-AP We need to zero out this reaction S +HCO ->HS +CO because both S and HCO are short-lived
C-AP      AR(182,I) = 5.0E-11
C-KK	This reaction has been removed from the primo3s rxn list.
C-KK      AR(182,I) = 1.E-99
      AR(183,I) = 1.5E-11
      AR(184,I) = 1.5E-11
      AR(185,I) = 1.7E-11 * EXP(-800./T(I))
   7  CONTINUE
C-AP 
C
      DO 8 I=1,NZ
      AR(190,I) = 1.0E-12
      AR(191,I) = 1.0E-11
      AR(192,I) = 1.5E+3
      AR(193,I) = 2.2E+4
      AR(194,I) = 1.0E-16
      AR(195,I) = 4.0E-12
      AR(196,I) = 1.5E-13
      AR(197,I) = 1.13E+3
      AR(198,I) = 7.0E-14
      AR(199,I) = 1.4E-11
      AR(200,I) = 3.6E-12 * EXP(-1100./T(I))
      AR(201,I) = 0.
      AR(202,I) = 9.0E-12 * EXP(-280./T(I))
      AR(203,I) = 2.9E-11 * EXP(240./T(I))
      AR(204,I) = 1.2E-11
      AR(205,I) = 8.3E-15
      AR(206,I) = 2.0E-15
      AR(207,I) = 1.0E-20
      AR(208,I) = 0.
      AR(209,I) = AR(44,I)
      AR(210,I) = AR(6,I)
      AR(212,I) = AR(11,I)
      AR(213,I) = AR(9,I)
      AR(214,I) = AR(7,I)
      AR(215,I) = 1.E-12
      AR(216,I) = AR(13,I)
      AR(217,I) = 1.E-11
      AR(219,I) = 5.03e-7 * (298./T(I))**2.16 * exp(-18701./T(I))
      AR(220,I) = 2.92e-13 * (T(I)/300.)**2.23 * exp(-23292./T(I))
   8  CONTINUE


C ***** GIORGI AND CHAMEIDES RAINOUT RATE *****
      GAM15 = 8.64E+05/2.0
      GAM8 = 7.0E+06/2.0
      AV = 6.02E+23
      WL = 1.0
      R = 1.36E-22
c     print*,ZTROP, DZ
      NH = ZTROP/DZ + 0.01
c      print*,NH
C  Loop over altitude
      DO 10 I=1,NH
      ZKM = Z(I)/1.E5
      TEMP = T(I)
C
C  Find appropriate GAMMA
      IF (ZKM.LE.1.51) THEN
         GAMMA = GAM15
      ELSE IF (ZKM.LT.8.) THEN
         GAMMA = GAM15 + (GAM8-GAM15)*((ZKM-1.5)/6.5)
      ELSE
         GAMMA = GAM8
      END IF
C
C  Find WH2O
      IF (ZKM.LE.1.) THEN
         X = 11.35 + 0.1*ZKM
      ELSE
         X = 11.5444 - 0.085333*ZKM - 9.1111E-03*ZKM*ZKM
      END IF
      WH2O = 10.0**X
C
C  Find F(Z)
      IF (ZKM.LE.1.51) THEN
         F = 0.1
      ELSE
         F = 0.16615 - 0.04916*ZKM + 3.37451E-03*ZKM*ZKM
      END IF
C
C  Loop over species
      DO 10 J=1,NQ
      RKJ = (WH2O/55.)/(AV*WL*1.0E-9 + 1./(H(J)*R*TEMP))
      QJ = 1. - F + F/(GAMMA*RKJ) * (1.0 - EXP(-RKJ*GAMMA))
  10  RAINGC(J,I) = (1. - EXP(-RKJ*GAMMA))/(GAMMA*QJ)
C
      NH1 = NH + 1
      DO 11 I=NH1,NZ
      DO 11 J=1,NQ
  11  RAINGC(J,I) = 0.
C
C ***** RAINOUT RATE *****
      DO 2 I=1,NZ
      ZKM = Z(I)/1.E5
      RAIN(I) = 0.
      IF(ZKM.LT.10.) RAIN(I) = 2.4E-6*EXP((6. - ZKM)/2.42)
      IF(ZKM.LT.6.) RAIN(I) = 2.4E-6
   2  CONTINUE

      RETURN
      END
