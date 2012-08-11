        
         program atm_chem
 
C     This program was created by James Kasting and modified mainly by
c     Alex Pavlov (AP) and Kara Kerelove (KK). The user-friendly version 
c     of this code was created by Antigona Segura (AS) (2005).
c     Some modifications are identified with the author's initials.
c
c     The code is mostly written in f77 but is compiled in f90 and it 
c     contains some f90 features.
c
C     This program was created from PRESO3.F in June, 1995, for the
C     purpose of calculating ozone levels on planets around stars of
C     different spectral types. It differs from PRESO3 by including 
C     NO3, N2O5, and CL2O2, and by not including the anthropogenic
C     CFC's. Sulfur chemistry was added by Alex Pavlov.
c     
c     The UV fluxes for stars that are not the Sun were provided by
c     Martin Cohen from new IUE (+model) data. 
c
c     Check the notes along the main program before use it. Look for 
c     the word 'NOTE:'
C
C     This program has been updated to the JPL '92 rate constant
C     recommendations. I have not checked all the photolysis cross
C     sections to see if they are up to date.
C
C         THIS PROGRAM IS DESIGNED SPECIFICALLY FOR AN O2 LEVEL
C     OF 1 PAL.  IT DIFFERS FROM THE OLD PRIMO3 BY INCLUDING AN
C     UPDATED ALGORITHM FOR CALCULATING PHOTOLYSIS OF O2 IN THE
C     SCHUMANN-RUNGE BANDS.  THE PHOTOLYSIS ALGORITHMS FOR O2 AND
C     NO MAY NOT BE APPLICABLE TO LOWER O2 LEVELS.  (THIS HAS NOT
C     BEEN CHECKED, BUT IT WOULD BE SURPRISING IF THEY WERE.)
C
C     THIS VERSION OF PRESO3 ALSO CONTAINS THE MODIFICATIONS USED IN 
C     PRIMO2 TO CALCULATE O2 SCHUMANN-RUNGE PHOTOLYSIS BY EXPONENTIAL
C     SUMS AND NO PHOTOLYSIS BY THE MODIFIED CIESLIK AND NICOLET METHOD.
C     THIS EDITION ALSO HAS GIORGI AND CHAMEIDES RAINOUT RATES AND A 
C     MANABE-WETHERALD RELATIVE HUMIDITY DISTRIBUTION. 
C
C       THIS PROGRAM IS A ONE-DIMENSIONAL MODEL OF THE PRIMORDIAL
C     ATMOSPHERE.  THE MIXING RATIOS OF THE LONG-LIVED SPECIES
C     ARE CALCULATED FROM THE EQUATION
C
C     DF/DT  =  (1/N)*D/DZ(KN*DF/DZ) + P/N - LF
C
C     WHERE
C     F = MIXING RATIO (USOL)
C     K = EDDY DIFFUSION COEFFICIENT (EDD)
C     N = TOTAL NUMBER DENSITY (DEN)
C     L = CHEMICAL LOSS FREQUENCY (XL)
C     P = CHEMICAL PRODUCTION RATE (XP)
C
C          THE SYSTEM OF PDES IS SOLVED USING THE REVERSE EULER
C     METHOD.  LETTING THE TIME STEP GO TO INFINITY GIVES YOU NEWTONS
C     METHOD, E.G. IT REVERTS TO AN EFFICIENT STEADY-STATE SOLVER.
C
C          THE LIST OF SUBROUTINES IS AS FOLLOWS:
C     (1) GRID   -  SETS UP THE ALTITUDE GRID
C     (2) RATES  -  DEFINES CHEMICAL REACTION RATES AND RAINOUT RATE
C     (3) PHOTO  -  COMPUTES PHOTOLYSIS RATES (CALLS TWOSTR)
C     (4) DENSTY -  COMPUTES ATMOSPHERIC DENSITIES FROM HYDROSTATIC
C                   EQUILIBRIUM AND INITIALIZES ABSORBER PROFILES
C     (5) DIFCO  -  COMPUTES DK = K*N BETWEEN GRID POINTS
C     (6) OUTPUTP -  PRINTS OUT RESULTS
C     (7) DOCHEM - DOES CHEMISTRY FOR ALL SPECIES AT ALL GRID
C                  POINTS BY CALLING CHEMPL
C     (8) CHEMPL - COMPUTES CHEMICAL PRODUCTION AND LOSS RATES
C                  FOR ONE SPECIES AT ALL GRID POINTS
C     (9) LTNING -  COMPUTES LIGHTNING PRODUCTION RATES FOR O2 AND
C                   N2 BASED ON CHAMEIDES' RESULTS
C    (10) SCALE  -  CAN BE USED TO SCALE THE EQUATIONS BY ARBITRARY
C                   FACTORS TO IMPROVE CONDITIONING OF MATRIX
C
C          OTHER DEFINED FUNCTIONS INCLUDE:
C     (1) TBDY   -  COMPUTES 3-BODY REACTION RATES
C     (2) E1     - EXPONENTIAL INTEGRAL OF ORDER ONE
C
C ***** REACTION LIST *****
C     1)  H2O + O(1D) = 2OH
C     2)  H2 + O(1D) = OH + H
C     3)  H2 + O = OH + H
C     4)  H2 + OH = H2O + H
C     5)  H + O3 = OH + O2
C     6)  H + O2 + M = HO2 + M
C     7)  H + HO2 = H2 + O2
C     8)  H + HO2 = H2O + O
C     9)  H + HO2 = OH + OH
C    10)  OH + O = H + O2
C    11)  OH + HO2 = H2O + O2
C    12)  OH + O3 = HO2 + O2
C    13)  HO2 + O = OH + O2
C    14)  HO2 + O3 = OH + 2O2
C    15)  HO2 + HO2 = H2O2 + O2
C    16)  H2O2 + OH = HO2 + H2O
C    17)  O + O + M = O2 + M
C    18)  O + O2 + M = O3 + M
C    19)  O + O3 = 2O2
C    20)  OH + OH = H2O + O
C    21)  O(1D) + N2 = O(3P) + N2
C    22)  O(1D) + O2 = O(3P) + O2
C    23)  O2 + HV = O(3P) + O(1D)
C    24)  O2 + HV = O(3P) + O(3P)
C    25)  H2O + HV = H + OH
C    26)  O3 + HV = O2 + O(1D)
C    27)  O3 + HV = O2 + O(3P)
C    28)  H2O2 + HV = 2OH
C    29)  CO2 + HV = CO + O(3P)
C    30)  CO + OH = CO2 + H
C    31)  CO + O + M = CO2 + M
C    32)  H + CO + M = HCO + M
C    33)  H + HCO = H2 + CO
C    34)  HCO + HCO = H2CO + CO
C    35)  OH + HCO = H2O + CO
C    36)  O + HCO = H + CO2
C    37)  O + HCO = OH + CO
C    38)  H2CO + HV = H2 + CO
C    39)  H2CO + HV = HCO + H
C    40)  HCO + HV = H + CO
C    41)  H2CO + H = H2 + HCO
C    42)  CO2 + HV = CO + O(1D)
C    43)  H + H + M = H2 + M
C    44)  HCO + O2 = HO2 + CO
C    45)  H2CO + OH = H2O + HCO
C    46)  H + OH + M = H2O + M
C    47)  OH + OH + M = H2O2 + M
C    48)  H2CO + O = HCO + OH
C    49)  H2O2 + O = OH + HO2
C    50)  HO2 + HV = OH + O
C    51)  CH4 + HV  =  1CH2 + H2
C    52)  CH3OOH + HV  =  H3CO + OH
C    53)  N2O + HV  =  N2 + O
C    54)  HNO2 + HV  = NO + OH
C    55)  HNO3 + HV  = NO2 + OH
C    56)  NO + HV  =  N + O
C    57)  NO2 + HV  =  NO + O
C    58)  CH4 + OH  =  CH3 + H2O
C    59)  CH4 + O(1D)  =  CH3 + OH
C    60)  CH4 + O(1D)  =  H2CO + H2
C    61)  1CH2 + CH4  =  2 CH3
C    62)  1CH2 + O2  =  H2CO + O
C    63)  1CH2 + N2  =  3CH2 + N2
C    64)  3CH2 + H2  =  CH3 + H
C    65)  3CH2 + CH4  =  2 CH3
C    66)  3CH2 + O2  =  H2CO + O
C    67)  CH3 + O2 + M  =  CH3O2 + M
C    68)  CH3 + OH  =  H2CO + H2
C    69)  CH3 + O  =  H2CO + H
C    70)  CH3 + O3  =  H2CO + HO2
C    71)  CH3O2 + HO2  =  CH3OOH + O2
C    72)  CH3O2 + CH3O2  =  2 H3CO + O2
C    73)  CH3O2 + NO  =  H3CO + NO2
C    74)  H3CO + O2  =  H2CO + HO2
C    75)  H3CO + O  =  H2CO + OH
C    76)  H3CO + OH  =  H2CO + H2O
C    77)  N2O + O(1D)  =  NO + NO
C    78)  N2O + O(1D)  =  N2 + O2
C    79)  N + O2  =  NO + O
C    80)  N + O3  =  NO + O2
C    81)  N + OH  =  NO + H
C    82)  N + NO  =  N2 + O
C    83)  NO + O3  =  NO2 + O2
C    84)  NO + O + M  =  NO2 + M
C    85)  NO + HO2  =  NO2 + OH
C    86)  NO + OH + M  =  HNO2 + M
C    87)  NO2 + O  =  NO + O2
C    88)  NO2 + OH + M  =  HNO3 + M
C    89)  NO2 + H  =  NO + OH
C    90)  HNO3 + OH  =  H2O + NO3
C    91)  HO2 + NO2 + M  =  HO2NO2 + M
C    92)  HO2NO2 + OH  =  NO2 + H2O + O2
C    93)  HO2NO2 + O  =  NO2 + OH + O2
C    94)  HO2NO2 + M  =  HO2 + NO2 + M
C    95)  HO2NO2 + HV  =  HO2 + NO2
C    96)  CH3OOH + OH  =  CH3O2 + H2O
C    97)  CH3O2 + OH  = H3CO + HO2
C    98)  O3 + NO2  =  O2 + NO3
C    99)  NO2 + NO3  =  NO + NO2 + O2
C   100)  O + NO3  =  O2 + NO2
C   101)  CH3Cl + hv  =  CH3 + Cl
C   102)  NO + NO3  =  NO2 + NO2
C   103)  OH + NO3  =  HO2 + NO2
C   104)  CH3Cl + OH  =  Cl + H2O
C   105)  CL + O3  =  ClO + O2
C   106)  Cl + H2  = HCl + H
C   107)  Cl + CH4  =  HCl + CH3
C   108)  Cl + CH3Cl  =  Cl + HCl
C   109)  Cl + H2CO  =  HCl + HCO
C   110)  Cl + H2O2  =  HCl + HO2
C   111)  Cl + HO2  =  HCl + O2
C   112)  Cl + HO2  =  ClO + OH
C   113)  Cl + ClONO2  =  Cl + Cl + NO2 (+O)
C   114)  Cl + NO + M  = NOCl + M
C   115)  Cl + NO2 + M  =  ClONO + M
C   116)  Cl + NOCl  =  NO + Cl2
C   117)  Cl + O2 + M  =  ClO2 + M
C   118)  Cl + ClO2  =  Cl2 + O2
C   119)  Cl + ClO2  =  ClO + ClO
C   120)  ClO + O  =  Cl + O2
C   121)  ClO + NO  =  Cl + NO2
C   122)  ClO + NO2 + M  =  ClONO2 + M
C   123)  ClO + HO2  =  HOCl + O2
C   124)  ClO + OH  =  CL + HO2
C   125)  HCl + OH  =  Cl + H2O
C   126)  HOCl + OH  =  ClO + H2O
C   127)  ClONO2 + OH  =  Cl + HO2 + NO2
C   128)  HCl + O  =  Cl + OH
C   129)  HOCl + O  =  ClO + OH
C   130)  ClONO2 + O  =  Cl + O2 + NO2
C   131)  Cl2 + OH  =  HOCl + Cl
C   132)  Cl2 + hv  =  Cl + Cl
C   133)  ClO2 + hv  =  ClO + O
C   134)  HCl + hv  =  H + Cl
C   135)  HOCl + hv  =  OH + Cl
C   136)  NOCl + hv  =  Cl + NO
C   137)  ClONO + hv  =  Cl + NO2
C   138)  CLONO2 + hv  =  Cl + NO3
C   139)  ClO2 + hv  =  Cl + O2
C   140)  HO2 + NO3  =  HNO3 + O2
C   141)  ClO + ClO + M  = Cl2O2 + M
C   142)  Cl2O2 + hv  =  ClO2 + Cl
C   143)  Cl2O2 + M  =  ClO + ClO
C   144)  ClO2 + M  =  Cl + O2
C   145)  Cl + NO3  =  ClO + NO2
C   146)  Cl + HOCl  =  Cl2 + OH
C   147)  ClO + NO3  =  ClONO + O2
C   148)  ClONO + OH  =  HOCl + NO2
C   149)  ClO2 + O  =  ClO + O2
C   150)  NO2 + O + M  =  NO3 + M
C   151)  NO3 + hv  =  NO2 + O
C   152)  NO3 + NO2 + M  =  N2O5 + M
C   153)  N2O5 + hv  =  NO2 + NO3
C   154)  N2O5 + M  =  NO2 + NO3 + M
C   155)  N2O5 + H2O  =  2 HNO3
C************** Added sulfur chemistry ********
C   156)  SO   + HV   =     S    +     O
C   157)  SO2  + HV   =     SO   +     O
C   158)  H2S  + HV   =     HS   +     H
C   159)  SO   + O2   =     O    +     SO2
C   160)  SO   + HO2  =     SO2  +     OH
C   161)  SO   + O    =     SO2
C   162)  SO   + OH   =     SO2  +     H
C   163)  SO2  + OH   =     HSO3
C   164)  SO2  + O    =     SO3
C   165)  SO3  + H2O  =     H2SO4
C   166)  HSO3 + O2   =     HO2  +     SO3
C   167)  HSO3 + OH   =     H2O  +     SO3
C   168)  HSO3 + H    =     H2   +     SO3
C   169)  HSO3 + O    =     OH   +     SO3
C   170)  H2S  + OH   =     H2O  +     HS
C   171)  H2S  + H    =     H2   +     HS
C   172)  H2S  + O    =     OH   +     HS
C   173)  HS   + O    =     H    +     SO
C   174)  HS   + O2   =     OH   +     SO
C   175)  HS   + HO2  =     H2S  +     O2
C   176)  HS   + HS   =     H2S  +     S
C   177)  HS   + HCO  =     H2S  +     CO
C   178)  HS   + H    =     H2   +     S
C   179)  HS   + S    =     H    +     S2
C   180)  S    + O2   =     SO   +     O
C   181)  S    + OH   =     SO   +     H
c-as This reaction was deleted and subtituted by 182a)
C   182)  S    + HCO  =     HS   +     CO
c   182a) SO2  + HV   =     S    +     O2
C   183)  S    + HO2  =     HS   +     O2
C   184)  S    + HO2  =     SO   +     OH
C   185)  HS   + H2CO =     H2S  +     HCO
C   186)  SO2  +     HV  =      SO21
C   187)  SO2  +     HV  =      SO23
C   188)  H2SO4 +     HV =       SO2   +   OH  +    OH
C   189)  SO3   +    HV  =      SO2  +     O
C   190)  SO21  +    M   =      SO23 +     M
C   191)  SO21  +    M   =      SO2  +     M
C   192)  SO21  +    HV  =      SO23 +     HV
C   193)  SO21  +    HV  =      SO2  +     HV
C   194)  SO21  +    O2  =      SO3  +     O
C   195)  SO21  +    SO2 =      SO3  +     SO
C   196)  SO23  +    M   =      SO2  +     M
C   197)  SO23  +    HV  =      SO2  +     HV
C   198)  SO23  +    SO2 =      SO3  +     SO
C   199)  SO    +    NO2 =      SO2  +     NO
C   200)  SO    +    O3  =      SO2  +     O2
C   201)  SO2   +    HO2 =      SO3  +     OH
C   202)  HS    +    O3  =      HSO  +     O2
C   203)  HS    +    NO2 =      HSO  +     NO
C   204)  S     +    O3  =      SO   +     O2
C   205)  SO    +    SO  =      SO2  +     S
C   206)  SO3   +    SO  =      SO2  +     SO2
C   207)  S     +    CO2 =      SO   +     CO
C   208)  SO    +    HO2 =      HSO  +     O2
C   209)  SO    +    HCO =      HSO  +     CO
C   210)  H     +    SO  =      HSO
C   211)  HSO   +    HV  =      HS   +     O
C   212)  HSO   +    OH  =      H2O  +     SO
C   213)  HSO   +    H   =      HS   +     OH
C   214)  HSO   +    H   =      H2   +     SO
C   215)  HSO   +    HS  =      H2S  +     SO
C   216)  HSO   +    O   =      OH   +     SO
C   217)  HSO   +    S   =      HS   +     SO
c-as Reactions added for N2O 
c   218)  N2 + O1D = N2O
c   219)  N2O + H = NO + NO + OH
c   220)  N2O + NO = NO2 + N2
C***********************************************************
C
C        THIS PROGRAM DOES THE CHEMISTRY AUTOMATICALLY.  THE CHEMICAL
C     REACTIONS ARE ENTERED ON DATA CARDS IN FIVE 10-DIGIT COLUMNS
C     STARTING IN COLUMN 11, I.E.
C
C         REAC1     REAC2     PROD1     PROD2     PROD3
C
C     THE IMPORTANT PARAMETERS DESCRIBING THE CHEMISTRY ARE
C        NR   = NUMBER OF REACTIONS
C        NSP  = NUMBER OF CHEMICAL SPECIES
C        NSP1 = NSP + 1 (INCLUDES HV)
C        NQ   = NUMBER OF SPECIES FOR WHICH A DIFFUSION EQUATION
C               IS SOLVED
C        NMAX = MAXIMUM NUMBER OF REACTIONS IN WHICH AN INDIVIDUAL
C               SPECIES PARTICIPATES
C
C        PHOTOLYSIS REACTIONS ARE IDENTIFIED BY THE SYMBOL HV (NOT
C     COUNTED IN EVALUATING NSP).  THREE-BODY REACTIONS ARE WRITTEN
C     IN TWO-BODY FORM, SO THE DENSITY FACTOR MUST BE INCLUDED IN
C     THE RATE CONSTANT.
C        THE CHEMICAL REACTION SCHEME IS STORED IN THE FOLLOWING MATRICE
C
C     ISPEC(NSP2) = VECTOR CONTAINING THE HOLLERITH NAMES OF THE
C                  CHEMICAL SPECIES.  THE LAST ENTRY MUST BE HV.
C     JCHEM(5,NR) = MATRIX OF CHEMICAL REACTIONS.  THE FIRST TWO ARE
C                   REACTANTS, THE LAST THREE ARE PRODUCTS.
C     ILOSS(2,NSP,NMAX) = MATRIX OF LOSS PROCESSES.  ILOSS(1,I,L)
C                         HOLDS REACTION NUMBER J, ILOSS(2,I,L) HOLDS
C                         REACTANT NUMBER.
C     IPROD(NSP,NMAX) = MATRIX OF PRODUCTION PROCESSES.  IPROD(I,L)
C                       HOLDS REACTION NUMBER J.
C     NUML(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF ILOSS
C     NUMP(NSP) = NUMBER OF NON-ZERO ELEMENTS FOR EACH ROW OF IPROD
C
c      PARAMETER(NZ=64, NQ=34, NQT=NQ+1)
c      PARAMETER(NEQ=NQ*NZ,LDA=3*NQ+1)
c      PARAMETER(NR=220, NF=34)
c      PARAMETER(NSP=55, NSP1=NSP+1, NSP2=NSP+2, NMAX=70)

       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNEQ_LDA.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'

      DIMENSION FVAL(NQ,NZ),FV(NQ,NZ),DJAC(LDA,NEQ),RHS(NEQ),IPVT(NEQ)
     2  ,SGFLUX(NQ),SMFLUX(NQ),VDEP(NQ),DD(NZ),DL(NZ),DU(NZ),
     3  USAVE(NQ,NZ),R(NZ),U(NQ)
      DIMENSION DPU(NZ,3),DPL(NZ,3)
      DIMENSION TA(NZ),TB(NZ),TC(NZ),TY(NZ)                                    
      DIMENSION TSAV(NZ)
      DIMENSION water(NZ),FLOW(NQT),fluxsave(108),sfxsave(10)

      CHARACTER :: STARR*3,DIRDATA*4, AA*11,DIRIO*2
     
c Name of the star
      INCLUDE 'INCLUDECHEM/comSTR.inc'
      INCLUDE 'INCLUDECHEM/comFLUXPHOTO.inc'
c DIRP contains DIRCOUP, this variable is defined as a character 
c in comDIRP.inc
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comEBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comGBLOK.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comPRESS1.inc'
      INCLUDE 'INCLUDECHEM/comQBLOK.inc'
      INCLUDE 'INCLUDECHEM/comRBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'

        DATA LH2CO,LO,LH2O,LOH,LHO2,LH2O2,LO3,LH,LH2,LCH4,LCO,
     2  LCH3OOH,LCH3O2,LN2O,LNO,LNO2,LHNO2,LHNO3,LHO2NO2,LNO3,LN2O5,
     3  LCL2O2,LCH3CL,LHOCL,LCL,LCLO,LHCL,LCLONO2,LH2S,LHS,LSO,LSO2,
     4  LH2SO4,LHSO,LSO4AER,LCH21,LCH23,LO1D,LCH3,LH3CO,LHCO,LN,LNOCL,
     5  LCLONO,LCLO2,LCL2,LS,LSO21,LSO23,LHSO3,LSO3,LS2,LO2,LCO2,LN2/
     6  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     7  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
     8  44,45,46,47,48,49,50,51,52,53,54,55/
C-AP
C
C   NO PREDISSOCIATION COEFFICIENTS (ALLEN AND FREDERICK, 1982)
C
      DATA ANO/-1.790868E+1, -1.924701E-1, -7.217717E-2, 5.648282E-2,
     2  4.569175E-2, 8.353572E-3, 3*0.,
     3  -1.654245E+1, 5.836899E-1, 3.449436E-1, 1.700653E-1,
     4  -3.324717E-2, -4.952424E-2, 1.579306E-2, 1.835462E-2,
     5  3.368125E-3/
C
      DATA BNO/7.836832E+3, -1.549880E+3, 1.148342E+2, -3.777754E+0,
     2  4.655696E-2, 1.297581E+4, -2.582981E+3, 1.927709E+2,
     3  -6.393008E+0, 7.949835E-2/
C
      DATA LLNO/3*0, 2*2, 3*0, 2*1, 25*0/
      DATA RNO2/60*0., .79, .83, .66, .15, 4*0./
C
C   CONSTANTS FOR 1ST EXPONENTIAL INTEGRAL
      DATA AI/.99999193, -.24991055, .05519968, -.00976004,
     2  .00107857, -.57721566/
      DATA BI/8.5733287401, 18.0590169730, 8.6347608925,
     2  .2677737343/
      DATA CI/9.5733223454, 25.6329561486, 21.0996530827,
     2  3.9584969228/
      DATA NUML,NUMP/NSP*0,NSP*0/
C
C ***** SOLUBILITY (GIORGI AND CHAMEIDES) *****
      DATA H/1.3E+04, 1.0E-99, 1.0E-99, 1.0E+05, 3.3E+04,
     2       2.0E+05, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99,
     3       1.0E-99, 2.0E+05, 3.3E+04, 1.0E-99, 1.9E-03,
     4       7.0E-03, 7.0E+11, 7.0E+11, 7.0E+11, 7.0E-03,
     5       7.0E+11, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99,
     6       1.0E-99, 7.0E+11, 1.0E-99, 
     7       0.14, 1.E+5, 1.9E-3, 1.E+4, 7E+11,
     8       9E+3/
C-AP  I have changed the Henry constant similar to the Archean code
C-AP  only for the sulfur species  from Archean. Note that Henry(SO2) 
C-AP  in the table is higher because we use effective Henry constant
C-AP  from Archean.

c  Temperature from the US Standard Atmosphere 1976. Used when the 
c  code is not coupled to the climate model 
      DATA T/288.15, 281.65, 275.15, 268.66, 262.17,
     &       255.68, 249.19, 242.70, 236.21, 229.73,
     &       223.25, 216.77, 216.65, 216.65, 216.65,
     &       216.65, 216.65, 216.65, 216.65, 216.65,
     &       216.65, 217.58, 218.57, 219.57, 220.56, 
     &       221.55, 222.54, 223.54, 224.53, 225.52,
     &       226.51, 227.50, 228.49, 230.97, 233.74, 
     &       236.51, 239.28, 242.05, 244.82, 247.58,
     &       250.35, 253.14, 255.88, 258.64, 261.40,
     &       264.16, 266.96, 269.68, 270.65, 270.65,
     &       270.65, 270.65, 269.03, 266.27, 263.52, 
     &       260.77, 258.02, 255.27, 252.52, 249.77, 
     &       247.02, 244.27, 241.53, 230.78/

c NOTE: Boundary conditions are defined here
C ***** UPPER BOUNDARY FLUXES *****
      DATA SMFLUX/NQ*0./
C
C ***** EFFUSION VELOCITIES *****
      DATA VEFF/NQ*0./
C
C ***** UPPER BOUNDARY CONDITIONS *****
      DATA MBOUND/0, 1, 8*0, 1, 24*0/
C
C   0 = CONSTANT EFFUSION VELOCITY (VEFF)
C   1 = CONSTANT FLUX (SMFLUX)
C
C ***** LOWER BOUNDARY FLUXES *****
c NOTE: SGFLUX is redefined later in the main code when the O2 
c mixing ratio is less than 0.21 or when the star is not the Sun
      DATA SGFLUX/NQ*0./
C
C ***** DEPOSITION VELOCITIES (NQ)*****
      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1., 3*0., 0.2, 1.,
     2  0., 2*0.2, 6*0.5, 0., 0.5, 1., 3*0.5, 0.02, 1, 3.E-4, 
     3  3*1./
c Deposition velocities used for planets around quiescent M dwarfs
c      DATA VDEP/0.2, 1., 0., 1., 1., 0.2, 0.1, 1.,2.4e-4,0.,1.2e-4, 0.2,
c     2  1.,0., 2*0.2, 6*0.5, 0., 0.5, 1., 3*0.5, 0.02, 1, 3.E-4, 
c     3  3*1./
C
C ***** LOWER BOUNDARY CONDITIONS (NQT)*****
c NOTE: Mixing ratios for present Earth for H2,CO,CH4, N2O and CH3Cl
c are defined after the model parameters.
c fixed mixing ratios
        DATA LBOUND/2*0, 1, 5*0, 1, 1, 1, 2*0, 1, 5*0, 3*0, 1, 12*0/
c fixed fluxes
c       DATA LBOUND/2*0, 1, 5*0, 2, 2, 2, 2*0, 2, 5*0, 3*0, 2, 12*0/
c Boundary conditions for planets around quiescent M dwarfs
c      DATA LBOUND/2*0, 1, 5*0, 0, 1, 0, 2*0, 2, 5*0, 3*0, 2, 12*0/
c Boundary conditions for AD Leo
c      DATA LBOUND/2*0, 1, 5*0, 0, 2, 0, 2*0, 2, 5*0, 3*0, 2, 12*0/

C
C   0 = CONSTANT DEPOSITION VELOCITY (VDEP)
C   1 = CONSTANT MIXING RATIO
C   2 = CONSTANT UPWARD FLUX (SGFLUX)
C
C Modify Henry constants for sensitivity run comparison to Archean atmosphere
      H(LHO2) = 9.E3
      H(LH2O2) = 6.2E5
      H(LH2CO) =  4.25E4

C============== FILE SECTION =============================

      DIRDATA = 'DATA'
      DIRIO = 'IO'

c **** INPUT FILES
c
c This file contains the chemical reactions used in the code
      OPEN(unit=61, file= DIRDATA//'/primo3s.chm', status='old')
c-as Next file contains solar flux and atmospheric data
c-as it MUST be read for all the cases  
      OPEN(unit=62, file= DIRDATA//'/photos.pdat', status='old') 
      OPEN(unit=63,file= DIRDATA//'/h2so4.pdat',status='old')
      OPEN(unit=64,file= DIRDATA//'/eddy.pdat',status='old')

c Files with the input parameters
c NOTE: Check this two before runnig the program
       OPEN(unit=65,file= DIRIO//'/input_atmchem.dat')
       OPEN(unit=66,file= DIRIO//'/planet.dat')  

c-as Next file was formely named atm_chem_coefs.dat.
c-as I has the same format as atm_composition.out   
      OPEN(unit=67, file= DIRIO//'/atm_composition.dat') 

c Next file is generated by the climate model contains altitude,
c temperature and water. Formerly called photo_input.dat.
c Only used when ICOUPLE=1
      OPEN(unit=71,file= DIRIO//'/fromClima2Photo.dat') 
      
c-as  Unit 72 is an input and output file that is shared with the climate code
c-as  when ICOUPLE= 1. It is OPEN later in the program to WRITE on it.
c Only used when ICOUPLE=1    
      OPEN(unit=72,file= DIRIO//'/mixing_ratios.dat')

      OPEN(unit=73, file= DIRDATA//'/faruvs.pdat')
      
c  NOTE: IMPORTANT files to read the far UV of all the stars (including 
c  the Sun) and the UV fluxes for stars others than the Sun.
c        74     fluxesKGF_photo.pdat
c        75     M star flux (name it as you like)
c        76     far UV flux (name depends on the star)

C **** OUTPUT FILES
c Main output
      OPEN(unit=90,file= DIRIO//'/outchem.dat')
c Files commented would contain information that is not
c needed for now but the write commands for them still in 
c program. To activate them just remove all c here and in 
c the write commands
c      OPEN(unit=10,file='primo3.plt')
c      OPEN(unit=14,file='primotemp.dat')
c      OPEN(unit=15,file='o3graph.dat')  

c-as MAIN output file MOVED to couple.f
C       OPEN(unit=90,file='outchem.dat')  

c Written on subroutine TWOSTR
c      OPEN(unit=82,file= DIRIO//'/wave_means.dat')

c This file contains the altitude in cm vs the ozone number density (cm^-3)
c      OPEN(unit=83,file= DIRIO//'/O3numdens.out')

c  This is an output file containing altitude, H2O, O3.
c  To be used as input of the climate model, formerly called Pass2SurfMP.dat
      OPEN(unit=84,file= DIRIO//'/fromPhoto2Clima.dat')

C These OTPUT files are open along the program 
c	UNIT   	NAME  
c        19     mixing_ratios.dat (main program)
c This file contains the altitude (cm) vs. water mixing ratio
c	 86	 DIRIO/h20mixing.out   (subroutine OUTPUTP)
c=================================================================

********* SET MODEL PARAMETERS *****

C     ZY = SOLAR ZENITH ANGLE (IN DEGREES)
C     AGL = DIURNAL AVERAGING FACTOR FOR PHOTORATES
C     ISEASON = TELLS WHETHER P AND T VARY WITH TIME (THEY DON'T FOR
C               ISEASON < 3)
C     IZYO2 = TELLS WHETHER SOLAR ZENITH ANGLE VARIES WITH TIME (0 SAYS
C             IT DOESN'T; 1 SAYS IT DOES)
C     IO2 = 0 FOR ALLEN AND FREDERICK O2 SCHUMANN-RUNGE COEFFICIENTS
C         = 1 FOR EXPONENTIAL SUM FITS (FOR LOW-O2 ATMOSPHERES)
C     INO = 0 FOR ALLEN AND FREDERICK NO PREDISSOCIATION COEFFICIENTS
C         = 1 FOR MODIFIED CIESLIK AND NICOLET FORMULATION
C     EPSJ = AMOUNT BY WHICH TO PERTURB THINGS FOR JACOBIAN CALCULATION
C     ZTROP = TROPOPAUSE HEIGHT (ABOVE WHICH H2O BEHAVES AS A NONCONDENS
C             ABLE GAS
C     STARR - Character variable to choose a star, it can be:
c             Sun, F2V, K2V,dMV 
c             DO NOT FORGET quotation marks
c     ICOUPLE - 1 = Coupled to the climate model              
c               0 = Not coupled
C     FCO2 = CO2 mixing ratio when ICOUPLE = 0
C     FO2 = O2 mixing ratio when ICOUPLE = 0
C     DT = INITIAL TIME STEP
C     TSTOP = TIME AT WHICH CALCULATION IS TO STOP
C     NSTEPS = NUMBER OF TIME STEPS TO RUN (IF TSTOP IS NOT REACHED)
      read(65,555)
      read(65,*)AA,STARR
      read(65,*)AA,FLUXFAC
      read(65,*)AA,INIT
      read(65,*)AA,TSTOP
      read(65,*)AA,DT
      read(65,*)AA,NSTEPS
      read(65,*)AA,ZY
      read(65,*)AA,AGL
      read(65,*)AA,ISEASON
      read(65,*)AA,IZYO2
      read(65,*)AA,IO2
      read(65,*)AA,INO
      read(65,*)AA,EPSJ
      read(65,*)AA,ZTROP
      read(65,*)AA,FCO2
      read(65,*)AA,FO2
  555  format(3/)
      close(65)
      ICOUPLE = 0

*** Read the flux from a star
      call readstar(FLUXFAC)
C-AP
***** READ THE PLANET PARAMETER DATAFILE *****
      READ(66,502) G,FSCALE,ALB,DELZ,ZTROP,JTROP
 502  FORMAT(F5.1/,F4.2/,F5.3/,E5.1/,E5.1/,I2)
      close(66)
C
C-AP
***** READ THE ATMOSPHERIC COMPOSITION
      READ(67,400) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 400  FORMAT(1P8E12.5)
      close(67)

      if (ICOUPLE.eq.1) then
	DO J=1, NZ
c Read the temperature and water profiles from the climate code
	  READ(71,*) Z(J),T(J),water(J)
	END DO
      close(71)
      endif
 351  FORMAT(I3,1PE10.3, 1PE12.3, 1PE12.3)

C-KK  Surface mixing ratios to share with the climate code
c  FAR is not needed on this code but it must be transfered to 
c  the climate model.
        READ(72,*) FAR                  !Argon
	READ(72,*) FCH4			!Methane
	READ(72,*) FCO2			!Carbon dioxide	
	READ(72,*) FO2			!Oxygen	
	READ(72,*) JTROP		!Tropopause layer
        READ(72,*) O3OLD                !Former O3 column depth
      close(72)

C   Initial constant mixing ratios used for present Earth
c       if (STARR=="Sun".and.FO2.gt.0.20) then
     
       if(INIT.eq.1) goto 77
c      Surface mixing ratios
         if(LBOUND(9).eq.1)FH2 = 5.5E-7
         if(LBOUND(11).eq.1)FCO = 9.0E-8
         if(LBOUND(14).eq.1)FN2O = 3.0E-7
         if(LBOUND(23).eq.1)FCH3CL = 5.0E-10

       DO i = 1, NZ
        if(LBOUND(9).eq.1) then
           USOL(LH2, i) = USOL(LH2,i)*(FH2/USOL(LH2,1))
        endif
       if(LBOUND(10).eq.1) then
           USOL(LCH4,i) = USOL(LCH4,i)*(FCH4/USOL(LCH4,1))
        endif
        if(LBOUND(11).eq.1)then
            USOL(LCO,i) = USOL(LCO,i)*(FCO/USOL(LCO,1))
        endif
        if(LBOUND(14).eq.1) then
            USOL(LN2O,i) = USOL(LN2O,i)*(FN2O/USOL(LN2O,1))
        endif
        if(LBOUND(23).eq.1)then
            USOL(LCH3CL,i) = USOL(LCH3CL,i)*(FCH3CL/USOL(LCH3CL,1))
        endif
        END DO
   77   continue       
c        endif
       
C-KK	Constant flux for long-lived gases. For details see the model 
c       description on: Segura et al. (2003) Astrobiology, 3(4), 689-708.
C-KK	Only used at O2 levels lower than 1 PAL and for planets around 
C-KK    other stars
C-KK	9 - H2   10 - CH4   11 - CO  14 - N2O   23 - CH3Cl
         if(LBOUND(9).eq.2)SGFLUX(9) = -7.821E9
         if(LBOUND(10).eq.2)SGFLUX(10) = 1.709E11
c        if(LBOUND(10).eq.2)SGFLUX(10) = 1.5E11  !value for AD Leo
         if(LBOUND(11).eq.2)SGFLUX(11) = 2.944e11
         if(LBOUND(14).eq.2)SGFLUX(14) = 1.134E9
         if(LBOUND(23).eq.2)SGFLUX(23) = 4.276E8


       write(90,*)"SOLAR ZENITH ANGLE (deg) = ", ZY
	 write(90,*)"FIXED SPECIES MIXING RATIOS:"
	 write(90,*)"   O2 = ",FO2," CO2 = ",FCO2," CH4 = ",FCH4
C	 print*,"   H2 = ",FH2," CO = ",FCO," N2O = ",FN2O
	
C ***** SPECIES DEFINITIONS *****
      ISPEC(1) = 4HH2CO
      ISPEC(2) = 1HO
      ISPEC(3) = 3HH2O
      ISPEC(4) = 2HOH
      ISPEC(5) = 3HHO2
      ISPEC(6) = 4HH2O2
      ISPEC(7) = 2HO3
      ISPEC(8) = 1HH
      ISPEC(9) = 2HH2
      ISPEC(10) = 3HCH4
      ISPEC(11) = 2HCO
      ISPEC(12) = 6HCH3OOH
      ISPEC(13) = 5HCH3O2
      ISPEC(14) = 3HN2O
      ISPEC(15) = 2HNO
      ISPEC(16) = 3HNO2
      ISPEC(17) = 4HHNO2
      ISPEC(18) = 4HHNO3
      ISPEC(19) = 6HHO2NO2
      ISPEC(20) = 3HNO3
      ISPEC(21) = 4HN2O5
      ISPEC(22) = 5HCL2O2
      ISPEC(23) = 5HCH3CL
      ISPEC(24) = 4HHOCL
      ISPEC(25) = 2HCL
      ISPEC(26) = 3HCLO
      ISPEC(27) = 3HHCL
      ISPEC(28) = 6HCLONO2
      ISPEC(29) = 3HH2S
      ISPEC(30) = 2HHS
      ISPEC(31) = 2HSO
      ISPEC(32) = 3HSO2
      ISPEC(33) = 5HH2SO4
      ISPEC(34) = 3HHSO
C
C   TRIDIAGONAL SOLVER
      ISPEC(35) = 6HSO4AER
C
C   Short-lived species
      ISPEC(36) = 4HCH21
      ISPEC(37) = 4HCH23
      ISPEC(38) = 3HO1D
      ISPEC(39) = 3HCH3
      ISPEC(40) = 4HH3CO
      ISPEC(41) = 3HHCO
      ISPEC(42) = 1HN
      ISPEC(43) = 4HNOCL
      ISPEC(44) = 5HCLONO
      ISPEC(45) = 4HCLO2
      ISPEC(46) = 3HCL2
      ISPEC(47) = 1HS
      ISPEC(48) = 4HSO21
      ISPEC(49) = 4HSO23
      ISPEC(50) = 4HHSO3
      ISPEC(51) = 3HSO3
C
C   Inert species
      ISPEC(52) = 2HS2
      ISPEC(53) = 2HO2
      ISPEC(54) = 3HCO2
      ISPEC(55) = 2HN2
      ISPEC(56) = 2HHV
      ISPEC(57) = 1HM

      do i=1,NSP
       NUML(i) = 0
       NUMP(i) = 0
      enddo

C
C ***** READ THE CHEMISTRY DATA CARDS *****
      READ (61,200) JCHEM
 200  FORMAT(10X,A8,2X,A8,2X,A8,2X,A8,2X,A8)
      close(61)
     
C NOTE: This command will print the list of chemical reactions in the 
c       output file 
C     PRINT 201,(J,(JCHEM(M,J),M=1,5),J=1,NR)
c 201  FORMAT(1X,I3,1H),5X,A8,4H +  ,A8,7H  =    ,A8,4H +  ,A8,4X,A8)
      KJAC = LDA*NEQ
C      PRINT 202,NQ,NZ,KJAC
c 202  FORMAT(//1X,'NQ=',I2,5X,'NZ=',I3,5X,'KJAC=',I7)
C     

C ***** REPLACE HOLLERITH LABELS WITH SPECIES NUMBERS IN JCHEM *****
      DO 5 J=1,NR
      DO 5 M=1,5
      IF(JCHEM(M,J).EQ.' ' ) GO TO 5
C-AP Change NSP1 to NSP2 because we added M to the species list
      DO 6 I=1,NSP2
      IF(JCHEM(M,J).NE.ISPEC(I)) GO TO 6
      JCHEM(M,J) = I
      GO TO 5
   6  CONTINUE
      IERR = J
      GO TO 25
   5  CONTINUE
C 
C ***** FILL UP CHEMICAL PRODUCTION AND LOSS MATRICES *****
      DO 7 M=1,2
      N = 3-M
      DO 7 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 7
      NUML(I) = NUML(I) + 1
      IF(NUML(I).GT.NMAX) GO TO 20
      K = NUML(I)
      ILOSS(1,I,K) = J
      ILOSS(2,I,K) = JCHEM(N,J)
   7  CONTINUE
C
      DO 8 M=3,5
      DO 8 J=1,NR
      I = JCHEM(M,J)
      IF(I.LT.1.OR.I.GT.NSP) GO TO 8
      NUMP(I) = NUMP(I) + 1
      IF(NUMP(I).GT.NMAX) GO TO 20
      K = NUMP(I)
      IPROD(I,K) = J
   8  CONTINUE
C-AP
C
** Read the eddy diffusion profile  
	do J=1,NZ
	 READ(64,245) Z(J),EDD(J)	
	end do 
 245  FORMAT(2(1PE11.3, 1X))
      close(64)

C-AP variable methane flux
C-AP      SGFLUX(LCH4) = SGFLUX(LCH4)*8
C-AP
C-AP      DO I=1,NZ
C-AP      USOL(LCH4,I) =  USOL(LCH4,I)*(1.E-3/USOL(LCH4,1))
C-AP      ENDDO 
C-AP      USOL(LH2,I) =  USOL(LH2,I)*(5.5E-7/USOL(LH2,1)) 
C-AP      USOL(LCO,I) =  USOL(LCO,I)*(9.0E-8/USOL(LCO,1)) 
C-AP      USOL(LN2O,I) =  USOL(LN2O,I)*(3.0E-7/USOL(LN2O,1)) 
C-AP  40  USOL(LCH3CL,I) =  USOL(LCH3CL,I)*(6.0E-10/USOL(LCH3CL,1)) 
C-AP
C-AP      NQOLD1 = NQOLD + 1
C-AP      DO 41 K=NQOLD1,NQ
C-AP      DO 41 I=1,NZ
C-AP  41  USOL(K,I) = 1.E-10
C-AP
C-AP      DO K=1,3
C-AP      DO I=1,NZ
C-AP        RPAR(I,K) = 1.E-9
C-AP        AERSOL(I,K) = 1.E-12
C-AP        WFALL(I,K) = 1.
C-AP        SO4AER(I) = 1.E-7
C-AP      ENDDO
C-AP      ENDDO 
C


      LTIMES = 0       !Counter for the photorate subroutine
      ISULF = 0
      VOLFLX = 3.E9    !Volcanic flux
      CALL GRID
      DZ = Z(2) - Z(1)
      CALL DENSTY(FO2,FCO2,P0) 
      CALL RATES
      CALL DIFCO

c  If the star is other than the Sun or there is a flare
c  the fluxes are saved here
      do i=1,108
       fluxsave(i)=FLUX(i)
      enddo 
     
c  Read the data for the photochemistry and the UV flux of the Sun
      CALL READPHOTO
      if(STARR.ne.'Sun') then
        do i=1,108
         FLUX(i)= fluxsave(i)
        enddo       
      endif
     
** Water calculation

      JTROP = ZTROP/DZ + 0.01

      CALL PSATRAT(H2O)

       DO 23 J=1,JTROP 
  23    USOL(LH2O,J) = H2O(J)
      
       if(ICOUPLE.eq.1) then
C-KK	This is added to make sure that tropospheric water is being
C-KK	handled consistently. The #s are imported from the climate model.
        DO J = 1, NZ
         USOL(LH2O,J) = water(J)
        END DO
       endif       !end of water calculation

      CALL LTNING(FO2,FCO2,P0)
      CALL AERTAB
      NZ1 = NZ - 1
      HA = 1.38E-16*T(NZ)/(1.67E-24*28.*980.)

      DTINV = 1./DT
      TIME = 0.
C
C ***** PRINT OUT INITIAL DATA *****
      CALL OUTPUTP(0,NSTEPS,0.,FLOW)	
C     
C ***** STORE CONSTANT JACOBIAN COEFFICIENTS *****
      DZ2 = DZ*DZ
      DU(1) = DK(1)/DEN(1)/DZ2
      DL(NZ) = DK(NZ1)/DEN(NZ)/DZ2
      DO 2 J=2,NZ1
      DU(J) = DK(J)/DEN(J)/DZ2
      DL(J) = DK(J-1)/DEN(J)/DZ2
   2  DD(J) = DU(J) + DL(J)
C
      KD = 2*NQ + 1
      KU = KD - NQ
      KL = KD + NQ
C
C   PRINT OUT RESULTS EVERY NPR TIME STEPS
      NPR = NSTEPS
      PRN = NPR
C
C   DO PHOTORATES EVERY MP TIME STEPS
      NPHOT = 0
      MP = 3
      PM = MP
      NN = 0
C
C ***** START THE TIME-STEPPING LOOP *****
C	STARSHIP
      DO 1 N=1,NSTEPS
      TIME = TIME + DT
      NN = NN + 1
      MS = (N-1)/MP
      SM = (N-1)/PM
      IF(NN.EQ.NSTEPS) SM = MS
      IF(SM-MS.GT.0.01) GO TO 18
      IF(N.GT.1 .AND. TIME.LT.1.E4) GO TO 18
C
      DO 35 I=1,NZ
      H2O(I) = ABS(USOL(LH2O,I))
      O3(I) = ABS(USOL(LO3,I))
      O2(I) = FO2
      CO2(I) = FCO2
C-AP
      FSO2(I) = ABS(USOL(LSO2,I))
      H2S(I) = ABS(USOL(LH2S,I))
C-AP
      CH4(I) = ABS(USOL(LCH4,I))
  35  CONTINUE
       
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1      
      CALL PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,IDO)
      CALL AERCON(H2O)
    
C
      VCO2 = (PCO2(NZ) + PCO2D(NZ)) * HA
      SMFLUX(11) = - VCO2*CO2(NZ)*DEN(NZ)
      SMFLUX(2) = SMFLUX(11)
      NMP = NSTEPS - MP
      IF (NN.gt.1 .and. nn.LT.NSTEPS) GO TO 18
      write(90,97)
  97  FORMAT(//1X,'PHOTOLYSIS RATES')
      write(90,98) 
  98  FORMAT(/5X,'Z',7X,'PO2',6X,'PO2D',5X,'PCO2',5X,'PCO2D',4X,
     2  'PH2O',5X,'PO3',6X,'PO3D',5X,'PH2O2',4X,'PHCO',5X,'PH2',
     3  6X,'PHO2')
      write(90,99)(Z(I),PO2(I),PO2D(I),PCO2(I),PCO2D(I),PH2O(I),
     2  PO3(I),PO3D(I),PH2O2(I),PHCO(I),PH2(I),PHO2(I),I=1,NZ,
     3  3)
  99  FORMAT(2X,1P12E9.2)
      write(90,198)
 198  FORMAT(/5X,'Z',6X,'PCH4',5X,'PCH3OOH',2X,'PN2O',5X,'PHNO3',4X,
     2  'PNO',6X,'PNO2',5X,'PHNO4',4X,'PCCL3F',3X,'PCCL2F2',2X,
     3  'PCCL4',4X,'PCH3CL')
      write(90,99) (Z(I),PCH4(I),PMOOH(I),PN2O(I),PHNO3(I),PNO(I),
     2  PNO2(I),PHNO4(I),PCCL3F(I),PCCL2F2(I),PCCL4(I),PCH3CL(I),
     3  I=1,NZ,3)
      write(90,197)
 197  FORMAT(/5X,'Z',6X,'PMCCL3',3X,'PCL2',5X,'PHOCL',4X,'PNOCL',4X,
     2  'PCLONO',3X,'PCLONO2',2X,'PCLO2',4X,'PHCL')
      write(90,199) (Z(I),PMCCL3(I),PCL2(I),PHOCL(I),PNOCL(I),PCLONO(I),
     2  PCLONO2(I),PCLO2(I),PHCL(I),I=1,NZ,3)
 199  FORMAT(2X,1P9E9.2)
      write(90,298)
 298  format(/5x,'Z',6x,'PNO3',5x,'PN2O5',4x,'PCL2O2')
      write(90,299) (z(i),pno3(i),pn2o5(i),pcl2o2(i),i=1,nz,3)
 299  format(2x,1p4e9.2)
  18  CONTINUE

      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL SEDMNT(FSULF,IDO)

      DO J=1,NZ
        AERSOL(J,1) = SO4AER(J)*DEN(J)/CONVER(J,1)
      ENDDO                
                                                   
C
C ***** SET UP THE JACOBIAN MATRIX AND RIGHT-HAND SIDE *****
      DO 17 J=1,LDA
      DO 17 K=1,NEQ
  17  DJAC(J,K) = 0.
      DO 19 K=1,NEQ
  19  RHS(K) = 0.
C
C     (DJAC IS EQUAL TO (1/DT)*I - J, WHERE J IS THE JACOBIAN MATRIX)
C
C   COMPUTE CHEMISTRY TERMS AT ALL GRID POINTS
      IDO = 0
      IF (NN.EQ.NSTEPS) IDO = 1
      CALL DOCHEM(FVAL,IDO)

      DO 9 I=1,NQ
      DO 9 J=1,NZ
      K = I + (J-1)*NQ
      RHS(K) = FVAL(I,J)
   9  USAVE(I,J) = USOL(I,J)
C
      DO 3 I=1,NQ
      DO 11 J=1,NZ
      R(J) = EPSJ * ABS(USOL(I,J))
  11  USOL(I,J) = USAVE(I,J) + R(J)
      CALL DOCHEM(FV,0)
C
      DO 12 M=1,NQ
      MM = M - I + KD
      DO 12 J=1,NZ
      K = I + (J-1)*NQ
  12  DJAC(MM,K) = (FVAL(M,J) - FV(M,J))/R(J)
C
      DO 10 J=1,NZ
  10  USOL(I,J) = USAVE(I,J)
   3  CONTINUE
C
C   COMPUTE TRANSPORT TERMS AT INTERIOR GRID POINTS
      DO 13 I = 1,NQ
      DO 14 J=2,NZ1
      K = I + (J-1)*NQ
      RHS(K) = RHS(K) - DD(J)*USOL(I,J) + DU(J)*USOL(I,J+1) +
     2  DL(J)*USOL(I,J-1)
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DD(J)
      DJAC(KU,K+NQ) = - DU(J)
  14  DJAC(KL,K-NQ) = - DL(J)
  13  CONTINUE
C
C ***** LOWER BOUNDARY CONDITIONS *****
      DO 15 K=1,NQ
      U(K) = USOL(K,1)
      LB = LBOUND(K)
C
C   CONSTANT DEPOSITION VELOCITY
      IF(LB.NE.0) GO TO 16
      RHS(K) = RHS(K) + DU(1)*(USOL(K,2) - U(K)) - VDEP(K)*U(K)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(1) + VDEP(K)/DZ
      DJAC(KU,K+NQ) = - DU(1)
      GO TO 15
C
C   CONSTANT MIXING RATIO
  16  IF(LB.NE.1) GO TO 31
      RHS(K) = 0.
      DO 36 M=1,NQ
      MM = KD + K - M
  36  DJAC(MM,M) = 0.
      DJAC(KU,K+NQ) = 0.
      DJAC(KD,K) = DTINV + DU(1)
      GO TO 15
C
C   CONSTANT UPWARD FLUX
  31  CONTINUE
      RHS(K) = RHS(K) + DU(1)*(USOL(K,2) - U(K)) + SGFLUX(K)/DEN(1)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DU(1)
      DJAC(KU,K+NQ) = - DU(1)
  15  CONTINUE
C
C ***** UPPER BOUNDARY CONDITIONS *****
      DO 30 I=1,NQ
      U(I) = USOL(I,NZ)
      K = I + NZ1*NQ
      MB = MBOUND(I)
C
C   CONSTANT EFFUSION VELOCITY
      IF(MB.NE.0) GO TO 29
      RHS(K) = RHS(K) + DL(NZ)*(USOL(I,NZ1) - U(I)) - VEFF(I)*U(I)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(NZ) + VEFF(I)/DZ
      DJAC(KL,K-NQ) = - DL(NZ)
      GO TO 30
C
C   CONSTANT DOWNWARD FLUX
  29  RHS(K) = RHS(K) + DL(NZ)*(USOL(I,NZ1) - U(I)) - SMFLUX(I)/
     2  DEN(NZ)/DZ
      DJAC(KD,K) = DJAC(KD,K) + DTINV + DL(NZ)
      DJAC(KL,K-NQ) = - DL(NZ)
  30  CONTINUE
C
      DO 33 J=1,NZ
      IF(Z(J).GT.ZTROP) GO TO 34
      K = 3 + (J-1)*NQ
      RHS(K) = 0.
      DO 32 M=1,NQ
      MM = M - 3 + KD
  32  DJAC(MM,K) = 0.
      DJAC(KD,K) = DTINV
      DJAC(KU,K+NQ) = 0.
      IF(J.EQ.1) GO TO 33
      DJAC(KL,K-NQ) = 0.
  33  CONTINUE
  34  CONTINUE
C
C ***** FACTOR THE JACOBIAN AND SOLVE THE LINEAR SYSTEM *****
      CALL SGBFA(DJAC,LDA,NEQ,NQ,NQ,IPVT,INDEX)
      IF(INDEX.NE.0.) write(90,103)N,INDEX
 103  FORMAT(/1X,'N =',I3,5X,'INDEX =',I3)
      CALL SGBSL(DJAC,LDA,NEQ,NQ,NQ,IPVT,RHS,0)
C
C   COMPUTE NEW CONCENTRATIONS
      EMAX = 0.
      DO 26 I=1,NQ
      DO 26 J=1,NZ
      K = I + (J-1)*NQ
C-KK	For all runs, inc. standard O2
      IF((I.EQ.LH2S).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHS).AND.(Z(J).GT.1.2E6)) GOTO 26
      IF((I.EQ.LHSO).AND.(Z(J).GT.2.E6)) GOTO 26
      IF((I.EQ.LCH3CL).AND.(Z(J).GT.2.8E6)) GOTO 26
      IF((I.EQ.LH).AND.(Z(J).LT.3.E6)) GOTO 26
      IF (I.EQ.LSO) GOTO 26
      IF (I.EQ.LCL2O2) GOTO 26
c For AD Leo
c      IF(I.EQ.LH2) GOTO 26
C-KK	Added for low-O2 environments
c      IF((I.EQ.LCLONO2).AND.(Z(J).GT.1.4E6)) GOTO 26
c       IF((I.EQ.LH2SO4).AND.(Z(J).GT.3.4E6)) GOTO 26
c       IF(I.EQ.LH2) GOTO 26
c       IF((I.EQ.LHSO).AND.(Z(J).GT.2.3E6)) GOTO 26
c      IF((I.EQ.LSO2).AND.(Z(J).GT.4.0E6)) GOTO 26
c       IF(I.EQ.LN2O5) GOTO 26
C      IF((I.EQ.LHNO2).AND.(Z(J).GT.4.4E6)) GOTO 26
c        IF((I.EQ.LCH3OOH).AND.(Z(J).GT.5.0E6)) GOTO 26
c       IF((I.EQ.LHO2NO2).AND.(Z(J).GT.5E6)) GOTO 26
C       IF((I.EQ.LHNO3).AND.(Z(J).GT.3.4E6)) GOTO 26
C       IF((I.EQ.LH2O).AND.(Z(J).GT.5.9E6)) GOTO 26
C       IF((I.EQ.LN2O).AND.(Z(J).GT.1.4E6)) GOTO 26
C       IF((I.EQ.LH2O2).AND.(Z(J).GT.2.3E6)) GOTO 26
C      IF((I.EQ.LNO3).AND.(Z(J).GT.4.9E6)) GOTO 26
c       IF((I.EQ.LHOCL).AND.(Z(J).GT.4.8E6)) GOTO 26
c        IF((I.EQ.LH2CO).AND.(Z(J).GT.5.E6)) GOTO 26
c        IF((I.EQ.LCH3O2).AND.(Z(J).GT.6.3E6)) GOTO 26
C       IF((I.EQ.LHO2).AND.(Z(J).GT.5.2E6)) GOTO 26
C	IF((I.EQ.LNO2).AND.(Z(J).GT.6.0E6)) GOTO 26
C	IF((I.EQ.LCLO).AND.(Z(J).GT.5.8E6)) GOTO 26
C	IF((I.EQ.LO3).AND.(Z(J).GT.6.1E6)) GOTO 26
c	IF(I.EQ.LCH3CL) GOTO 26
C	IF((I.EQ.LH2S).AND.(Z(J).GT.8.E5)) GOTO 26
 

      REL(I,J) = RHS(K)/USOL(I,J)
      EREL = ABS(REL(I,J))
      EMAX = AMAX1(EMAX,EREL)
      IF(EREL.LT.EMAX) GO TO 26
      IS = I
      JS = J
      UMAX = USOL(I,J)
      RMAX = RHS(K)
  26  USOL(I,J) = USOL(I,J) + RHS(K)
C
      DO 4 J=1,NZ
      IF(Z(J).LT.ZTROP) USOL(3,J) = H2O(J)
   4  CONTINUE

C-AP Adding tridiagonal solver
C       TRIDIAGONAL INVERSION *****
      L=1
      I = NQ + L
      IF(I.EQ.LSO4AER) MZ = 50
C-AP      IF(I.EQ.LS8TESTAER) MZ = 40
C-AP      IF(I.EQ.LHCAER) MZ = NZ
      MZ1 = MZ - 1
      MZP1 = MZ + 1
C
C   COMPUTE ADVECTION TERMS FOR PARTICLES
      DPU(1,L) = WFALL(2,L)*DEN(2)/DEN(1)/(2.*DZ)
      DPL(NZ,L) = WFALL(NZ1,L)*DEN(NZ1)/DEN(NZ)/(2.*DZ)
      DO 38 J=2,NZ1
      DPU(J,L) = WFALL(J+1,L)*DEN(J+1)/DEN(J)/(2.*DZ)
  38  DPL(J,L) = WFALL(J-1,L)*DEN(J-1)/DEN(J)/(2.*DZ)
C                                                                              
C
C   TA = LOWER DIAGONAL, TB = DIAGONAL, TC = UPPER DIAGONAL, TY =
C   RIGHT-HAND SIDE
      DO 70 J=1,NZ
      TA(J) = 0.
      TB(J) = 0.
      TC(J) = 0.
  70  TY(J) = 0.
C
      DO 44 J=1,MZ
      TB(J) = YL(I,J)
  44  TY(J) = YP(I,J)/DEN(J)
C
      DO 45 J=2,MZ1
      TA(J) = - DL(J) + DPL(J,L)
      TB(J) = TB(J) + DD(J)
  45  TC(J) = - DU(J) - DPU(J,L)
C                                                                              
C   BOUNDARY CONDITIONS
      TA(MZ) = - DL(MZ) + DPL(MZ,L)
      TB(MZ) = TB(MZ) + DL(MZ) + 0.5*WFALL(MZ,L)/DZ
      TB(1) = TB(1) + DU(1) + (.01 - 0.5*WFALL(1,L))/DZ
      TC(1) = - DU(1) - DPU(1,L)
C
      CALL SGTSL(MZ,TA,TB,TC,TY,NFLAG)
C-AP      STOP
      IF (NFLAG.NE.0) write(90,401) N,NFLAG,I
 401  FORMAT(//1X,'TRIDIAGONAL SOLVER FAILED AT N =',I3,2X,
     2  'NFLAG =',I2,2X,'SPECIES #',I2)
C
      IF(I.EQ.LSO4AER) THEN
        DO 59 J=1,MZ
   59     SO4AER(J) = TY(J)
C-AP      ELSEIF(I.EQ.LS8AER) THEN
C-AP        DO 46 J=1,MZ
C-AP   46     S8(J) = TY(J)
C-AP      ELSEIF(I.EQ.LHCAER) THEN
C-AP        DO 60 J=1,MZ
C-AP          HCAER(J) = TY(J)
C-AP   60   CONTINUE
      ENDIF
C
C   FILL UP UPPER PORTION WITH APPROXIMATE ANALYTIC SOLUTION
      IF(I.EQ.LSO4AER .AND. MZ.NE.NZ) THEN
        DO 61 J=MZP1,NZ
          SO4AER(J) = SO4AER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
   61     SO4AER(J) = AMAX1(SO4AER(J),1E-100)
C-AP      ELSEIF(I.EQ.LS8AER .AND. MZ.NE.NZ) THEN
C-AP        DO 47 J=MZP1,NZ
C-AP          S8(J) = S8(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
C-AP  47      S8(J) = AMAX1(S8(J),1E-100)
C-AP      ELSEIF(I.EQ.LHCAER .AND. MZ.NE.NZ) THEN
C-AP        DO 62 J=MZP1,NZ
C-AP          HCAER(J) = HCAER(J-1) * EXP(-WFALL(J,L)*DZ/EDD(J))
C-AP          HCAER(J) = AMAX1(HCAER(J),1E-100)
C-AP   62   CONTINUE
      ENDIF
C-AP   58 CONTINUE

      if(NN.eq.NSTEPS) GOTO 2010
C   AUTOMATIC TIME STEP CONTROL
      DTSAVE = DT
      IF(EMAX.GT.0.20)  DT = 0.7*DTSAVE
      IF(EMAX.GT.0.15)  DT = 0.9*DTSAVE
      IF(EMAX.LT.0.10)  DT = 1.1*DTSAVE
      IF(EMAX.LT.0.05)  DT = 1.3*DTSAVE
      IF(EMAX.LT.0.03)  DT = 1.5*DTSAVE
      IF(EMAX.LT.0.01)  DT = 2.0*DTSAVE
      IF(EMAX.LT.0.003) DT = 5.0*DTSAVE
      IF(EMAX.LT.0.001) DT = 10.*DTSAVE

c Adjusting time step to assure that TIME<=TSTOP
      TIME1 = TIME + DT
      if(TIME1.ge.TSTOP) then
        DT = ABS(TIME-TSTOP)
        NN = NSTEPS -1
      endif

      DTINV = 1./DT
C
2010  ISP = ISPEC(IS)
      ZMAX = Z(JS)
      IF(SM-MS.GT.0.01) GO TO 317
      write(90,100)N,EMAX,ISP,ZMAX,
     & UMAX,RMAX,DT,TIME
 100  FORMAT(1X,'N =',I4,2X,'EMAX =',1PE9.2,' FOR ',A8,
     2  'AT Z =',E9.2,1X,'U =',E9.2,1X,'RHS =',E9.2,
     3  2X,'DT =',E9.2,2X,'TIME =',E9.2)
C-AP
C   COMPUTE ATMOSPHERIC OXIDATION STATE
      DO 42 I=1,NQ
      SR(I) = 0.
      DO 43 J=1,JTROP
  43  SR(I) = SR(I) + RAINGC(I,J)*USOL(I,J)*DEN(J)*DZ
      PHIDEP(I) = VDEP(I)*USOL(I,1)*DEN(1)
  42  TLOSS(I) = SR(I) + PHIDEP(I)
C
      SR(LSO4AER) = 0.
      DO 48 J=1,JTROP
      SR(LSO4AER) = SR(LSO4AER) + RAINGC(LH2SO4,J)*SO4AER(J)*DEN(J)*DZ
  48  CONTINUE
      PHIDEP(LSO4AER) = (WFALL(1,1) + .01) * SO4AER(1) * DEN(1)
      TLOSS(LSO4AER) = SR(LSO4AER) + PHIDEP(LSO4AER)
C
C   COMPUTE SULFUR BUDGET AND READJUST SO2 (H2S) OUTGASSING RATE IF SO
C   DESIRED (PROGRAM IS SET UP FOR PURE SO2 OUTGASSING)
      SLOSS = TLOSS(LH2S) + TLOSS(LHS) + TLOSS(LSO) +
     2  TLOSS(LSO2) + TLOSS(LH2SO4) + TLOSS(LHSO) + 
     3  TLOSS(LSO4AER)
      SLOSSP = SLOSS - TLOSS(LSO2)
      IF (ISULF.EQ.0 .OR. TIME.LT.1.E6) GO TO 316
      IF (LBOUND(LSO2).EQ.2) SGFLUX(LSO2) = SGFLUX(LSO2) * VOLFLX/SLOSS
      IF (LBOUND(LH2S).EQ.2) SGFLUX(LH2S) = SGFLUX(LH2S) * VOLFLX/SLOSS
 316  CONTINUE
C
C-KK      PRINT 101, SLOSS,SLOSSP
C-KK 101  FORMAT(10X,'SLOSS =',E10.3,2X,'SLOSSP =',
C-KK     2  E10.3/)
 317  CONTINUE
C
C   RETRY TIME STEP IF EMAX EXCEEDS 30 PERCENT
      IF(EMAX.LT.0.3) GO TO 28
      write(90,*)'RETRY TIME STEP BECAUSE EMAX=',EMAX
      DT = 0.5*DTSAVE
      TIME = TIME - DTSAVE
      if(TIME.lt.0) TIME = 0.
      DO 27 I=1,NQ
      DO 27 J=1,NZ
  27  USOL(I,J) = USAVE(I,J)
  28  CONTINUE
C
      NS = N/NPR
      SN = N/PRN
      IF(NN.EQ.NSTEPS) SN = NS
      IF(SN-NS.GT.1.E-3) GO TO 37
C
      write(90,*)'BEFOREOUT'
      CALL OUTPUTP(NN,NSTEPS,TIME,FLOW)
  37  CONTINUE
      IF(INDEX.NE.0) STOP
      IF(NN.EQ.NSTEPS) GO TO 22
      IF(TIME.GE.TSTOP) NN = NSTEPS - 1
      if(TIME.lt.TSTOP.and.N.eq.NSTEPS) then
       write(*,'(A,I4,A)')'Time not reached after',NSTEPS,
     &  ' of the photochemical model. The run is stopping now.'
       STOP
      endif
          
   1  CONTINUE      
     
C ***** END THE TIME-STEPPING LOOP *****

  22  CONTINUE
      if(ICOUPLe.eq.1) NSTEPS =N

      OPEN(unit=81,file= DIRIO//'/atm_composition.out')
      WRITE(81,399) USOL,T,EDD,DEN,SO4AER,AERSOL,WFALL,RPAR
 399  FORMAT(1P8E12.5)
      close(81) 
        

c Transfer results to the climate model (ICOUPLE=1)
      if (ICOUPLE.eq.1) then
       DO 255 I=1,NZ
c Transfer O3 and H2O
        WRITE(84,254) Z(I),PRESS(I),O3(I),H2O(I) 
 254    FORMAT(1PE9.3,3(E10.2))
 255   CONTINUE
       close(84)
      endif
C-KK  Need to find coldtrap by locating water mixing ratio minimum. 

	Jcold = 0

	DO J = 1, NZ
 	   IF (JCOLD .EQ. 0) THEN
           IF (T(J) .LT. T(J+1)) JCOLD = J
	   END IF
	END DO

	OPEN(unit=19,file= DIRIO//'/mixing_ratios.out')
c     Transfer surface mixing ratios
        WRITE(19,*) FAR
	WRITE(19,*) USOL(LCH4,1) 
	WRITE(19,*) FCO2
	WRITE(19,*) O2(1)
	WRITE(19,*) Jcold
        WRITE(19,*) O3COL
       close(19)

C   The following WRITEs out the o3 mixing ratio values so that the 
c   profiles can be compared
C
c      DO 333 I=1,NZ
c        WRITE(15,332) O3(I),Z(I)
c 332  FORMAT(E10.2,1X,1PE9.3)
c 333  CONTINUE 
C   The following WRITEs out the temperature and altitude so that 
C   the profiles can be compared.
C
c      DO 367 I=1,NZ
c	WRITE(14,366) T(I),Z(I)
c 366  FORMAT(1PE9.3,1X,1PE9.3)
c 367  CONTINUE     
c
c   MaKe file for plots
c-as Not needed for now
c      WRITE(10,699)
c 699  format(2x,'Alt',6x,'O',10x,'O3',9x,'CH4',8x,'CO',9x,'N2O',8x,
c     1  'CH3Cl',6x,'n(OH)',6x,'n(O1D)',5x,'PO2',8x,'PO3D')  
c      do 700 i=1,nz
c      zkm = z(i)/1.e5
c      po2tot = po2(i) + po2d(i)
c 700  WRITE(10,701) zkm,usol(2,i),usol(7,i),usol(10,i),usol(11,i),
c     1  usol(14,i),usol(23,i),sl(4,i),sl(31,i),po2tot,po3d(i)
c 701  format(f5.1,1x,1p10e11.4)
      GO TO 21
  20  write(90,300)I
 300  FORMAT(//1X,'NMAX EXCEEDED FOR SPECIES ',I3)
      GO TO 21
  25  write(90,301)IERR
 301  FORMAT(//1X,'ERROR IN REACTION ',I3)
C
  21  CONTINUE
c      close(82)
      STOP
      END 

c---------------------------------------------------------------


      SUBROUTINE OUTPUTP(N,NSTEPS,TIME,FLOW)
       INCLUDE 'INCLUDECHEM/parNZ.inc'
       INCLUDE 'INCLUDECHEM/parNQ_NQT.inc'
       INCLUDE 'INCLUDECHEM/parNR.inc'
       INCLUDE 'INCLUDECHEM/parNF.inc'
       INCLUDE 'INCLUDECHEM/parNSP_NSP1_NSP2.inc'
       INCLUDE 'INCLUDECHEM/parNMAX.inc'
      INCLUDE 'INCLUDECHEM/comDIRP.inc'
      INCLUDE 'INCLUDECHEM/comABLOK.inc'
      INCLUDE 'INCLUDECHEM/comBBLOK.inc'
      INCLUDE 'INCLUDECHEM/comCBLOK.inc'
      INCLUDE 'INCLUDECHEM/comDBLOK.inc'
      INCLUDE 'INCLUDECHEM/comFBLOK1.inc'
      INCLUDE 'INCLUDECHEM/comNBLOK.inc'
      INCLUDE 'INCLUDECHEM/comSULBLK.inc'
      INCLUDE 'INCLUDECHEM/comZBLOK.inc'
      INCLUDE 'INCLUDECHEM/comAERBLK.inc'
      INCLUDE 'INCLUDECHEM/comSATBLK.inc'
      INCLUDE 'INCLUDECHEM/comRRATS.inc'

      DIMENSION FUP(NQT),FLOW(NQT),CON(NQT),FLUXCH(NQT,NZ)
     2  ,ZF(NZ)
C
      ISKIP = 4
      JSKIP = ISKIP
      IF(N.EQ.NSTEPS) ISKIP = 2
      TIMEY = TIME/3600./24./365.
      write(90,100)TIME,TIMEY
 100  FORMAT(/1X,'TIME =',E11.4,5X,'TIMEY =',1pe13.4,1X,'YEARS')
      write(90,101)NPHOT
 101  FORMAT(/1X,'NPHOT =',I3)
C
      write(90,105)
 105  FORMAT(/1X,'MIXING RATIOS OF LONG-LIVED SPECIES'/)
      IROW = 12
      LR = NQ/IROW + 1
      RL = FLOAT(NQ)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 8 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
 110  FORMAT(/5X,'Z',8X,13(A8,1X))
      DO 20 I=1,3
  20  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
      DO 21 I=4,NZ,ISKIP
  21  write(90,120) Z(I),(USOL(K,I),K=K1,K2)
 120  FORMAT(1X,1P13E9.2)
c      OPEN(unit=86,file= DIRIO//'/h2omixing.out')
      DO I=1,NZ
c      IF ((LH2O.GE.K1).AND.(LH2O.LE.K2)) WRITE(86,121)Z(I),USOL(LH2O,I)
      END DO
c      close(86)  
 121  FORMAT(1PE10.3, 1PE10.3)
      IF (N.EQ.0) GO TO 8
      write(90,140)
 140  FORMAT(/1X,'TP, TL')
      write(90,145)(TP(K),K=K1,K2)
      write(90,145)(TL(K),K=K1,K2)
 145  FORMAT(10X,1P12E9.2)
   8  CONTINUE
C-AP
      write(90,106)
 106  FORMAT(//1X,'MIXING RATIO OF AEROSOL'/)
      write(90,185)
 185  FORMAT(5X,'Z',6X,'SO4AER')
      DO 18 J=1,3
  18  write(90,182) Z(J),SO4AER(J)
      DO 19 J=4,NZ,ISKIP
  19  write(90,182) Z(J),SO4AER(J)
 182  FORMAT(1X,1P2E9.2)
C-AP
      IF (N.EQ.0) RETURN
C
      write(90,183) TP(LSO4AER)
      write(90,184) TL(LSO4AER)
 183  FORMAT(/2X,'TP',6X,1P1E9.2)
 184  FORMAT(2X,'TL',6X,1P1E9.2)
C
      write(90,150) O3COL
 150  FORMAT(//1X,'OZONE COLUMN DEPTH = ',1PE11.4)
C-AP
      write(90,152) H2SCOL,SO2COL
 152  FORMAT(/1X,'SULFUR COLUMN DEPTHS:  H2S =',1PE10.3,2X,'SO2 =',
     2  E10.3,2X)
C-AP
	DO i = 1, NZ
	 IF (USOL(3,i) .LT. USOL(3,i+1)) THEN
		JCOLD = i
		GOTO 352
	 END IF
	END DO
 352  CONTINUE
      write(90,151) JCOLD, USOL(3,JCOLD)
 151  FORMAT(/1X,I3,' FH2O AT COLD TRAP =',1PE10.3)
      IF(N.LT.NSTEPS) RETURN
C
C ***** PRINT ON LAST ITERATION ONLY *****
      DO 1 I=1,NZ
   1  ZF(I) = Z(I) + 0.5*DZ
C
      DO 3 K=1,NQ
      DO 2 I=1,NZ
   2  SL(K,I) = USOL(K,I)*DEN(I)
      DO 4 I=1,NZ1
   4  FLUXCH(K,I) = - DK(I)*(USOL(K,I+1) - USOL(K,I))/DZ
   3  CONTINUE
C
      K = LSO4AER
      J = 1
      DO 5 I=1,NZ1
      FLUXCH(K,I) = -DK(I)*(SO4AER(I+1) - SO4AER(I))/DZ
     2  - 0.5*(WFALL(I,J)*DEN(I)*SO4AER(I) + WFALL(I+1,J)*DEN(I+1)
     3         *SO4AER(I+1))
   5  CONTINUE
C
      DO 15 K=1,NQT
      FLOW(K) = FLUXCH(K,1) - (YP(K,1) - YL(K,1)*SL(K,1))*DZ
      FUP(K) = FLUXCH(K,NZ1) + (YP(K,NZ) - YL(K,NZ)*SL(K,NZ))*DZ
      CON(K) = TP(K) - TL(K) + FLOW(K) - FUP(K)
 502  format('K=',I2,' FLOW(K)=',1PE10.3,' FUP(K)=',E10.3)
  15  CONTINUE

C-AP
      FLOW(3) = FLUXCH(3,11)
      CON(3) = TP(3) - TL(3) + FLOW(3) - FUP(3)
      DO 6 I=1,10
   6  FLUXCH(3,I) = 0.
C
      write(90,125)
 125  FORMAT(/1X,'NUMBER DENSITIES OF LONG-LIVED SPECIES'/)
      DO 9 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 22 I=1,NZ,ISKIP
       write(90,120) Z(I),(SL(K,I),K=K1,K2)
c      if ((LO3.GE.K1).AND.(LO3.LE.K2)) WRITE(83,355) Z(I),SL(LO3,I)
  22  CONTINUE
   9  CONTINUE
 355  FORMAT(1PE10.3,1PE12.3)
c     close(83)

      ISKIP = 4
      write(90,155)
 155  FORMAT(/1X,'FLUXES OF LONG-LIVED SPECIES'/)
      ZFL = 0.
      ZFT = ZF(NZ)
      DO 10 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQT
C-AP      IF (L.EQ.LR) K2 = NQ
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,120) ZFL,(FLOW(K),K=K1,K2)
      DO 23 I=1,NZ,ISKIP
  23  write(90,120) ZF(I),(FLUXCH(K,I),K=K1,K2)
      write(90,120) ZFT,(FUP(K),K=K1,K2)
  10  CONTINUE
C
      write(90,175)
 175  FORMAT(/1X,'RAINOUT RATE, PHIDEP, TLOSS AND LOWER B.C.'/)
      write(90,176)
 176  FORMAT(1X,'FOLLOWED BY TP, TL, FUP, FLOW, CON'/)
      DO 13 L=1,LR
      K1 = 1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NQT
      write(90,110) (ISPEC(K),K=K1,K2)
      write(90,145) (SR(K),K=K1,K2)
      write(90,145)(PHIDEP(K),K=K1,K2)
      write(90,145)(TLOSS(K),K=K1,K2)
      write(90,146) (LBOUND(K),K=K1,K2)
 146  FORMAT(14X,12(I1,8X))
      write(90,145)
      write(90,145) (TP(K),K=K1,K2)
      write(90,145) (TL(K),K=K1,K2)
      write(90,145) (FUP(K),K=K1,K2)
      write(90,145) (FLOW(K),K=K1,K2)
      write(90,145) (CON(K),K=K1,K2)     
  13  CONTINUE

C-AP 175  FORMAT(/1X,'RAINOUT OF LONG-LIVED SPECIES'/)
C-AP      DO 13 L=1,LR
C-AP      K1 = 1 + (L-1)*IROW
C-AP      K2 = K1 + IROW - 1
C-AP      IF (L.EQ.LR) K2 = NQ
C-AP      PRINT 110, (ISPEC(K),K=K1,K2)
C-AP  13  PRINT 145,(SR(K),K=K1,K2)
C-AP
C
C   COMPUTE CONSERVATION OF SULFUR
      SULDEP = - (FLOW(LHS) + FLOW(LSO) + FLOW(LH2SO4)
     2  + FLOW(LHSO) + FLOW(LSO4AER))
      IF (LBOUND(LSO2).EQ.0) SULDEP = SULDEP - FLOW(LSO2)
      IF (LBOUND(LH2S).EQ.0) SULDEP = SULDEP - FLOW(LH2S)
      SULRAN = SR(LH2S) + SR(LHS) + SR(LSO) + SR(LSO2) +
     2  SR(LH2SO4) + SR(LHSO) + SR(LSO4AER) 
      SULLOS = SULDEP + SULRAN
      SULPRO = 0.
      IF (LBOUND(LSO2).NE.0) SULPRO = SULPRO + FLOW(LSO2)
      IF (LBOUND(LH2S).NE.0) SULPRO = SULPRO + FLOW(LH2S)
      SO4LOS = TLOSS(LH2SO4) + TLOSS(LSO4AER)
      write(90,177) SULLOS,SULPRO,SO4LOS
 177  FORMAT(/1X,'CONSERVATION OF SULFUR:',/5X,'SULLOS =',1PE10.3,
     2  2X,'SULPRO =',E10.3,2X,'SO4LOS =',E10.3)
C-AP
      write(90,179)
  179 FORMAT(/1X,'INTEGRATED REACTION RATES'/)
      write(90,180) RAT
  180 FORMAT(1X,1P10E10.3)
C
      write(90,160)
 160  FORMAT(/1X,'ATMOSPHERIC PARAMETERS AND PH EQ SPECIES')
      NPE = NSP - NQ
      NQ1 = NQ + 1
      LR = NPE/IROW + 1
      RL = FLOAT(NPE)/IROW + 1
      DIF = RL - LR
      IF (DIF.LT.0.001) LR = LR - 1
C
      DO 12 L=1,LR
      K1 = NQ1 + (L-1)*IROW
      K2 = K1 + IROW - 1
      IF (L.EQ.LR) K2 = NSP
      write(90,110) (ISPEC(K),K=K1,K2)
      DO 24 I=1,NZ,ISKIP
  24  write(90,120) Z(I),(SL(K,I),K=K1,K2)
  12  CONTINUE
C
      write(90,190)
 190  FORMAT(/1X,'ATMOSPHERIC PARAMETERS')
      write(90,195)
 195  FORMAT(/4X,'Z',9X,'T',9X,'EDD',7X,'DEN')
      write(90,200)(Z(I),T(I),EDD(I),DEN(I),I=1,NZ,ISKIP)
 200  FORMAT(1X,1P4E10.3)
C
      write(90,230)
 230  FORMAT(/1X,'SULFATE AEROSOL PARAMETERS')
      write(90,235)
 235  FORMAT(/4X,'Z',8X,'AERSOL',5X,'RPAR',6X,'WFALL',5X,'FSULF',4X,
     2  'TAUSED',4X,'TAUEDD',4X,'TAUC',6X,'H2SO4S',4X,'H2SO4',5X,
     3  'CONSO4',4X,'CONVER')
      write(90,240) (Z(I),AERSOL(I,1),RPAR(I,1),
     & WFALL(I,1),FSULF(I),
     2  TAUSED(I,1),TAUEDD(I),TAUC(I,1),H2SO4S(I),USOL(LH2SO4,I),
     3  CONSO4(I),CONVER(I,1),I=1,NZ,ISKIP)
 240  FORMAT(1X,1P12E10.3)
C
      RETURN
      END
