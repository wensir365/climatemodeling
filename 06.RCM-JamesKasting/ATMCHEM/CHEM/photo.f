
      SUBROUTINE PHOTO(ZY,AGL,LTIMES,ISEASON,IZYO2,IO2,INO,N)
      
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'
      INCLUDE '../INCLUDECHEM/parNR.inc'
      INCLUDE '../INCLUDECHEM/parNSP_NSP1_NSP2.inc'
      INCLUDE '../INCLUDECHEM/parNMAX.inc'
 
      DIMENSION SO2(NZ),SO3(NZ),SIGL(NZ),SIGNOL(NZ),CL(NZ),
     2  CNO(NZ),D0(2),WAV(108),XL(10)
      dimension S1(10),SS(108)

      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comCBLOK.inc'
      INCLUDE '../INCLUDECHEM/comQBLOK.inc'
      INCLUDE '../INCLUDECHEM/comRBLOK.inc'
      INCLUDE '../INCLUDECHEM/comSBLOK.inc'
      INCLUDE '../INCLUDECHEM/comTBLOK.inc'
      INCLUDE '../INCLUDECHEM/comO3BLOK.inc'
      INCLUDE '../INCLUDECHEM/comFLUXPHOTO.inc'
c      common/fluxphoto/FLUX(108),SFX(10),EFLUX(108),     
c     &  GFLUX(108),PhEFLUX(108),PHGFLUX(108),ESFX(10),GSFX(10),
c     &  PhESFX(10),PhGSFX(10)        
     
       data XL/1725.,1675.,1625.,1575.,1525.,1475.,1425.,1375.,
     & 1325.,1216./   
C
      PM = 1.67E-24
      BK = 1.38E-16
      RGAS = 8.3143E7
      WT = O2(1)*32. + CO2(1)*44. + (1.-O2(1)-CO2(1))*28. + 0.4
      RMG = RGAS/(WT*G)
      PI = 3.14159
      ZYR = ZY*PI/180.
      U0 = COS(ZYR)
      AM = 1./U0
      NZ1 = NZ - 1
      RO3 = 0.9
      D0(1) = 1.2E-6
      D0(2) = 2.6E-6
      ALNO = 5.1E7
      QKNO = 1.5E-9
      DIJ = 1.65E9
      NPHOT = NPHOT + 1
      HC = 6.625E-27 * 3.00E10

      IF (N.EQ.1) write(90,119)
 119  FORMAT(/1X,'ENERGY FLUXES IN W/m^2/nm and photons/m^2/s/nm
     &  (NOT DIURNALLY AVERAGED)')
      IF (N.eq.1) write(90,118)
 118  FORMAT(//2X,'L',6X,'WAV',10X,'EFLUX',9X,'GFLUX',
     & 10X,'S(1)',8x,'PhEFLUX',8x,'PhGFLUX')
      
C ***** SET PHOTOLYSIS RATES TO ZERO *****
      DO 1 I=1,NZ
      PO2(I) = 0.
      PO3(I) = 0.
      PO3D(I) = 0.
      PH2O(I) = 0.
      PH2O2(I) = 0.
      PCO2(I) = 0.
      PCO2D(I) = 0.
      PHCO(I) = 0.
      PH2(I) = 0.
      PO2D(I) = 0.
      PHO2(I) = 0.
      PCH4(I) = 0.
      PMOOH(I) = 0.
      PN2O(I) = 0.
      PHNO3(I) = 0.
      PNO(I) = 0.
      PNO2(I) = 0.
      PHNO4(I) = 0.
      PCCL3F(I) = 0.
      PCCL2F2(I) = 0.
      PCCL4(I) = 0.
      PCH3CL(I) = 0.
      PMCCL3(I) = 0.
      PCL2(I) = 0.
      PCLO2(I) = 0.
      PHCL(I) = 0.
      PHOCL(I) = 0.
      PNOCL(I) = 0.
      PCLONO(I) = 0.
      PCLONO2(I) = 0.
      pno3(i) = 0.
      pn2o5(i) = 0.
      pcl2o2(i) = 0.
C-AP
      PSO2(I) = 0.
      PSO21(I) = 0.
      PSO23(I) = 0.
      PSO(I) = 0.
      PH2S(I) = 0.
      PH2SO4(I) = 0.
C-AP                                                           
   1  CONTINUE

C
C ***** CALCULATE COLUMN DEPTHS ABOVE EACH COLLOCATION POINT FOR O2
C       (TOO) AND O3 (TO3).  ADD OTHER SPECIES HERE IF YOU WANT
C       ADDITIONAL MAJOR ABSORBERS. *****
      HA = RMG*T(NZ)
      HAD = HA*DEN(NZ)
      TTOT(NZ) = DEN(NZ)*HA
      TO2(NZ) = O2(NZ)*HAD
      TO3(NZ) = O3(NZ)*HAD
      TCO2(NZ) = CO2(NZ)*HAD
      TH2O(NZ) = H2O(NZ)*HAD
C-AP
      TSO2(NZ) = FSO2(NZ) * HAD
      TH2S(NZ) = H2S(NZ) * HAD
      TCH4(NZ) = CH4(NZ) * HAD
C-AP                                                 
      DO 2 M=1,NZ1
      I = NZ - M
      DZ = Z(I+1) - Z(I)
      HA = RMG*0.5*(T(I) + T(I+1))
      EFAC = (1. - EXP(-DZ/HA))*DEN(I)*HA
      TTOT(I) = TTOT(I+1) + EFAC
      TO2(I) = TO2(I+1) + EFAC*SQRT(O2(I)*O2(I+1))
      TO3(I) = TO3(I+1) + EFAC*SQRT(O3(I)*O3(I+1))
      TCO2(I) = TCO2(I+1) + EFAC*SQRT(CO2(I)*CO2(I+1))
      TH2O(I) = TH2O(I+1) + EFAC*SQRT(H2O(I)*H2O(I+1))
C-AP
      TSO2(I) = TSO2(I+1) + EFAC*SQRT(FSO2(I)*FSO2(I+1))
      TH2S(I) = TH2S(I+1) + EFAC*SQRT(H2S(I)*H2S(I+1))
      TCH4(I) = TCH4(I+1) + EFAC*SQRT(CH4(I)*CH4(I+1))
c-AP                         
   2  CONTINUE
      IF(LTIMES.GT.0) GO TO 3          
C
C  SCALE FLUX(L) VALUES FOR SOLAR FLUX AT EARTH (FSCALE=1) OR
C  MARS (FSCALE=0.43).
C
      DO 800 L=1,108
  800  FLUX(L)=FLUX(L)*FSCALE
        
C   CALCULATE N2O5 CROSS SECTIONS AT LONG WAVELENGTHS
      DO 147 L=44,61
      WAVN = 0.5*(WAVL(L) + WAVU(L))/10.
      TEMP = 250.
 147  SN2O5(L) = 1.E-20 * EXP(2.735 + ((4728.5 - 17.127*WAVN)/TEMP))
C
C   CALCULATE NO2 QUANTUM YIELDS
      DO 47 L=47,60
      WAVN = 0.5*(WAVL(L) + WAVU(L))/10.
  47  RNO2(L) = 1. - 8.E-4*(WAVN - 275.)
C
C ***** CALCULATE ALLEN AND FREDERICK SCHUMANN-RUNGE BAND COEFS *****
   3  CONTINUE
      IF (LTIMES.GT.0 .AND. ISEASON.LT.3) GO TO 13
C
C   REPEAT THIS SECTION ONLY IF PRESSURE AND TEMPERATURE VARY WITH TIME
      DO 4 I=1,NZ
   4  PLOG(I) = ALOG10(DEN(I)*BK*T(I)/1.E3)
      CALL O3PHOT
C
      DO 5 L=1,17
      DO 6 I=1,NZ
   6  SIGL(I) = 0.
C
      IF (L.GE.15) GO TO 7
      DO 8 I=1,NZ
   8  BIGX(I) = PLOG(I)
      GO TO 9
   7  CONTINUE
      DO 10 I=1,NZ
  10  BIGX(I) = T(I)
C
   9  KMAX = KA(L)
      DO 11 K=1,KMAX
      DO 11 I=1,NZ
  11  SIGL(I) = SIGL(I) + A(L,K)*BIGX(I)**(K-1)
      DO 12 I=1,NZ
  12  SIG0(I,L) = 10.**SIGL(I)
   5  CONTINUE
C
C   COEFFICIENTS FOR NITROUS OXIDE (NO)
      DO 30 L=1,2
      DO 31 I=1,NZ
  31  SIGNOL(I) = 0.

      DO 32 K=1,9
      DO 32 I=1,NZ
  32  SIGNOL(I) = SIGNOL(I) + ANO(K,L)*PLOG(I)**(K-1)
      DO 33 I=1,NZ
  33  SIGNO0(I,L) = 10.**SIGNOL(I)
  30  CONTINUE
  13  CONTINUE
      IF(LTIMES.GT.0 .AND. IZYO2.LT.1) GO TO 14
C
C   REPEAT THIS SECTION ONLY IF SOLAR ZENITH ANGLE OR O2 VARIES
C   WITH TIME
      DO 39 I=1,NZ
  39  TO2L(I) = ALOG10(TO2(I))
      DO 15 L=1,17
      DO 16 I=1,NZ
  16  CL(I) = 0.
C
      KMAX = KB(L)
      DO 17 K=1,KMAX
      DO 17 I=1,NZ
  17  CL(I) = CL(I) + B(L,K)*TO2L(I)**(K-1)
C
      DO 18 I=1,NZ
      C = 10.**CL(I)
      SD = SIG0(I,L) * AM**(-C)
  18  SRO2(I,L) = AMIN1(SD,2.E-19)
  15  CONTINUE
C
C   COEFFICIENTS FOR NO
      DO 35 L=1,2
      DO 36 I=1,NZ
  36  CNO(I) = 0.
C
      DO 37 K=1,5
      DO 37 I=1,NZ
  37  CNO(I) = CNO(I) + BNO(K,L)*TO2L(I)**(K-1)
      DO 38 I=1,NZ
      SD = SIGNO0(I,L) * AM**CNO(I)
  38  SIGNO(I,L) = AMIN1(SD,1.E-15)
  35  CONTINUE
  14  CONTINUE
c      print*,'In photo.f before SHORTWAVE LOOP'
C
C ***** START SHORTWAVE LOOP (1754 - 2532 A) *****
      DO 19 L=1,35
      IF (L.LE.17) GO TO 27
      DO 25 I=1,NZ
      SO3(I) = SO31(L)
  25  SO2(I) = SO2HZ(L)
      GO TO 28
C
  27  CONTINUE
      DO 26 I=1,NZ
      SO3(I) = SO31(L)
  26  SO2(I) = SO2HZ(L) + SRO2(I,L)
  28  CONTINUE
      SSO2T = SSO2(L) + SSO21(L)
      KN = 1
      ALP = 1.
      IF (IO2.EQ.1 .AND. L.LE.17) KN = NK(L)
C
      M=0
C   LOOP OVER K'S AT LOW O2 LEVELS
      DO 19 K=1,KN 
      IF (IO2.EQ.0) GO TO 20
      IF (L.GT.17) GO TO 66
      ALP = ALPHA(L,K)
      DO 65 I=1,NZ
  65  SO2(I) = SO2HZ(L) + BETA(L,K)
  66  CONTINUE
C
C   SKIP MULTIPLE SCATTERING ROUTINE IF THE PURE ABSORPTION OPTICAL
C   DEPTH AT 20 KM IS MUCH GREATER THAN THE SCATTERING OPTICAL DEPTH
      WAV(L) = 0.5 * (WAVU(L) + WAVL(L))
      SIGR = SIGRAY(WAV(L)) * (1. + 1.5*CO2(1))
      TAURAY = SIGR * TTOT(20)
      TAUABS = SO2(20)*TO2(20) + SO3(20)*TO3(20) + SCO2(L)*TCO2(20)
     2  + SH2O(L)*TH2O(20) + SSO2T*TSO2(20) + SH2S(L)*TH2S(20)
      RATIO = TAUABS/TAURAY
C-AP      IF (RATIO.GT.5.) GO TO 20
      CALL TWOSTR(WAV(L),SIGR,U0,SO3,SO2,SCO2(L),SH2O(L),SSO2T,
     2 SH2S(L))
      M=1
      GO TO 21
C
  20  CONTINUE
      TAU = SO2(NZ)*TO2(NZ) + SO3(NZ)*TO3(NZ) + SH2O(L)*TH2O(NZ)
     2  + SCO2(L)*TCO2(NZ)
      S(NZ) = EXP(-AM*TAU)
C
      DO 22 J=1,NZ1
      I = NZ - J
      SIGO2 = 0.5 * (SO2(I) + SO2(I+1))
      DTAU = SIGO2*(TO2(I) - TO2(I+1)) + SO3(I)*(TO3(I) - TO3(I+1))
     2  + SH2O(L)*(TH2O(I) - TH2O(I+1)) + SCO2(L)*(TCO2(I)
     3  - TCO2(I+1))
  22  S(I) = S(I+1) * EXP(-AM*DTAU)

  21  FLX = FLUX(L)*AGL*ALP
      ss(L) =s(1)

      DO 23 I=1,NZ
      PO2(I) = PO2(I) + FLX*SO2(I)*S(I)
      PO3D(I) = PO3D(I) + FLX*SO3(I)*RO3*S(I)
      PO3(I) = PO3(I) + FLX*SO3(I)*(1.-RO3)*S(I)
      PH2O(I) = PH2O(I) + FLX*SH2O(L)*S(I)
      PH2O2(I) = PH2O2(I) + FLX*SH2O2(L)*S(I)
      PCO2(I) = PCO2(I) + FLX*SCO2(L)*S(I)
      PHO2(I) = PHO2(I) + FLX*SHO2(L)*S(I)
      PHNO3(I) = PHNO3(I) + FLX*SHNO3(L)*S(I)
      PN2O(I) = PN2O(I) + FLX*SN2O(L)*S(I)
      PMOOH(I) = PMOOH(I) + FLX*SMOOH(L)*S(I)
      PNO2(I) = PNO2(I) + FLX*SNO2(L)*RNO2(L)*S(I)
      PHNO4(I) = PHNO4(I) + FLX*SHNO4(L)*S(I)
      PCCL3F(I) = PCCL3F(I) + FLX*SCCL3F(L)*S(I)
      PCCL2F2(I) = PCCL2F2(I) + FLX*SCCL2F2(L)*S(I)
      PCCL4(I) = PCCL4(I) + FLX*SCCL4(L)*S(I)
      PCH3CL(I) = PCH3CL(I) + FLX*SCH3CL(L)*S(I)
      PMCCL3(I) = PMCCL3(I) + FLX*SMCCL3(L)*S(I)
      PCL2(I) = PCL2(I) + FLX*SCL2(L)*S(I)
      PCLO2(I) = PCLO2(I) + FLX*SCLO2(L)*S(I)
      PHCL(I) = PHCL(I) + FLX*SHCL(L)*S(I)
      PHOCL(I) = PHOCL(I) + FLX*SHOCL(L)*S(I)
      PNOCL(I) = PNOCL(I) + FLX*SCLNO(L)*S(I)
      PCLONO(I) = PCLONO(I) + FLX*SCLONO(L)*S(I)
      PCLONO2(I) = PCLONO2(I) + FLX*SCLONO2(L)*S(I)
      pn2o5(i) = pn2o5(i) + flx*sn2o5(l)*s(i)
      pcl2o2(i) = pcl2o2(i) + flx*scl2o2(l)*s(i)
      PSO2(I) = PSO2(I) + FLX*SSO2(L)*S(I)
      PSO21(I) = PSO21(I) + FLX*SSO21(L)*S(I)
c     PSO(I) = PSO(I) + FLX*SSO(L)*S(I)
      PH2S(I) = PH2S(I) + FLX*SH2S(L)*S(I)
      PH2SO4(I) = PH2SO4(I) + FLX*SHCL(L)*S(I)
C-AP
  23  CONTINUE
C-AP    
c PRinting the source function for each wavelength
       WAV(L) = 0.5 * (WAVU(L) + WAVL(L))
c     write(*,'(A,I4,2x,F8.1)')'wavelength ',L, WAV(L)
c      do i=1,NZ
c       write(*,'(I4,E10.3)')i,S(i)
c      enddo 

      IF (N.EQ.0) GO TO 72
      IF(L.LE.17 .AND. K.LT.KN) GO TO 72
       TAUR = SIGR*TTOT(1)
       DELWAV = WAVU(L) - WAVL(L)

c  Here the flux is converted from photons/cm^2/s to W/m^2/nm
c  cm^ converted to m^2 (1e4 factor)
c  HC/lambda in cgs converted to Joules (1e7 factor)
c  WAV converted from Angstrongs to cm (1e8 factor)
c  DELWAV converted from Angstrongs to nanometers (1e-1 factor)
c      EFLUX(L) = (HC*1.e8/WAV(L))*(1.e4/1.e7)*FLUX(L)/
c     & ((1.e-1*DELWAV))
c      GFLUX(L) = EFLUX(L)*S(1)
c Converting FLUX from photons/cm^2/s to photons/m^2/s/nm
c      PhEFLUX(L) = 1.e4*FLUX(L)/((1.e-1*DELWAV))
c      PhGFLUX(L) = PhEFLUX(L)*S(1)
c      write(90,120) L,WAV(L),TAUR,EFLUX(L),GFLUX(L),
c     & S(1),PhEFLUX(L),PhGFLUX(L)

  72  CONTINUE
C
C   NO PREDISSOCIATION IN THE D00 (1910 A) AND D10 (1830 A) BANDS
      NOL = LLNO(L)
      IF (NOL.EQ.0) GO TO 19
      IF (INO.EQ.1) GO TO 48
C
C   FREDERICK AND ALLEN METHOD (EFFECTIVE CROSS SECTIONS)
      DO 40 I=1,NZ
  40  PNO(I) = PNO(I) + FLX*SIGNO(I,NOL)*S(I)
      GO TO 19
C
C   OLD (CIESLIK AND NICOLET) METHOD WITH INTENSITIES UPDATED TO
C   FREDERICK AND HUDSON (1979)
  48  CONTINUE
      DO 49 I=1,NZ
      RN2 = DIJ/(ALNO + DIJ + QKNO*DEN(I))
  49  PNO(I) = PNO(I) + 0.5*D0(NOL)*S(I)*RN2*AGL*ALP 
    
  19  CONTINUE
C ***** END SHORTWAVE LOOP *****

      DO 60 I=1,NZ
  60  SO2(I) = 0.
c      print*,'In photo.f AFTER SHORTWAVE LOOP'
C ***** START LONGWAVE LOOP (2532 - 8550 A) *****
      DO 57 L=36,108
      IF (L.LT.39 .OR. L.GT.59) GO TO 61
C
      DO 29 I=1,NZ
      TI = AMAX1(T(I),203.)
      TI = AMIN1(T(I),273.)
      FR = (TI - 203.)/70.
  29  SO3(I) = FR*SO32(L) + (1.-FR)*SO31(L)
      GO TO 64
C
  61  CONTINUE
      DO 63 I=1,NZ
  63  SO3(I) = SO31(L)
C-AP
  64  SSO2T = 0.
      SH2SL = 0.
      IF (L.GT.68) GO TO 165
      SSO2T = SSO21(L) + SSO23(L)
      SH2SL = SH2S(L)
 165  CONTINUE
C
C-AP
C-AP  64  CONTINUE
C
C   SKIP MULTIPLE SCATTERING ROUTINE IF THE PURE ABSORPTION OPTICAL
C   DEPTH AT 20 KM IS MUCH GREATER THAN THE SCATTERING OPTICAL DEPTH
      WAV(L) = 0.5 * (WAVU(L) + WAVL(L))
      SIGR = SIGRAY(WAV(L)) * (1. + 1.5*CO2(1))
      TAURAY = SIGR*TTOT(20)
      TAUABS = SO3(20)*TO3(20) + SSO2T*TSO2(20) + SH2SL*TH2S(20)
      RATIO = TAUABS/TAURAY
C-AP Always do TWOSTR
C-AP      IF (RATIO.GT.5.) GO TO 51
      CALL TWOSTR(WAV(L),SIGR,U0,SO3,SO2,0.,0.,SSO2T,SH2SL)
      GO TO 52
C
  51  CONTINUE
      TAU = SO3(NZ)*TO3(NZ)
      S(NZ) = EXP(-AM*TAU)
C
      DO 53 J=1,NZ1
      I = NZ - J
      SIGO3 = 0.5 * (SO3(I) + SO3(I+1))
      DTAU = SIGO3 * (TO3(I) - TO3(I+1))
  53  S(I) = S(I+1) * EXP(-AM*DTAU)
C
  52  FLX = FLUX(L)*AGL
      SS(L) = s(1)
c   Calculate values for Whittet plots
      wav(l) = 0.5*(wavl(l) + wavu(l))

      DO 54 I=1,NZ
      BL = BT(I)*(0.1*WAV(L) - ALM0(I))
      RO3 = AT(I)*ATAN(BL) + CT(I)
      RO3 = AMAX1(RO3,0.)
      RO3 = AMIN1(RO3,0.9)
      PO3D(I) = PO3D(I) + FLX*SO3(I)*RO3*S(I)
      PO3(I) = PO3(I) + FLX*SO3(I)*(1.-RO3)*S(I)
 54   CONTINUE
C
      IF(L.NE.68) GO TO 170
C   Set NO3 photorate at this wavelength. Reduce the rate by 'fudge'
C   to simulate nighttime chemistry
      fudge = 0.25
      DO 169 I=1,NZ
 169  PNO3(I) = 0.18 * AGL * S(I)/S(NZ) * fudge
 170  CONTINUE
C
      IF (L.GT.68) GO TO 57
      DO 55 I=1,NZ
      PH2O2(I) = PH2O2(I) + FLX*SH2O2(L)*S(I)
      PHCO(I) = PHCO(I) + FLX*SH2CO(L)*RHCO(L)*S(I)
      PH2(I) = PH2(I) + FLX*SH2CO(L)*RH2(L)*S(I)
      PHNO3(I) = PHNO3(I) + FLX*SHNO3(L)*S(I)
      PMOOH(I) = PMOOH(I) + FLX*SMOOH(L)*S(I)
      PNO2(I) = PNO2(I) + FLX*SNO2(L)*RNO2(L)*S(I)
      PHNO4(I) = PHNO4(I) + FLX*SHNO4(L)*S(I)
      PCCL3F(I) = PCCL3F(I) + FLX*SCCL3F(L)*S(I)
      PCCL4(I) = PCCL4(I) + FLX*SCCL4(L)*S(I)
      PCL2(I) = PCL2(I) + FLX*SCL2(L)*S(I)
      PCLO2(I) = PCLO2(I) + FLX*SCLO2(L)*S(I)
      PHOCL(I) = PHOCL(I) + FLX*SHOCL(L)*S(I)
      PNOCL(I) = PNOCL(I) + FLX*SCLNO(L)*S(I)
      PCLONO(I) = PCLONO(I) + FLX*SCLONO(L)*S(I)
      PCLONO2(I) = PCLONO2(I) + FLX*SCLONO2(L)*S(I)
      pn2o5(i) = pn2o5(i) + flx*sn2o5(l)*s(i)
      pcl2o2(i) = pcl2o2(i) + flx*scl2o2(l)*s(i)
C-AP
      PSO21(I) = PSO21(I) + FLX*SSO21(L)*S(I)
      PSO23(I) = PSO23(I) + FLX*SSO23(L)*S(I)
      PH2S(I) = PH2S(I) + FLX*SH2S(L)*S(I)
C-AP
  55  CONTINUE

c PRinting the source function for each wavelength
c      write(*,'(A,I4,2x,F8.1)')'wavelength ',L,WAV(L)
c      do i=1,NZ
c       write(*,'(I4,E10.3)')i,S(i)
c      enddo

       IF (N.EQ.0) GO TO 57
       TAUR = SIGR*TTOT(1)
       DELWAV = WAVU(L) - WAVL(L)

c  Here the flux is converted from photons/cm^2/s to W/m^2/nm
c  cm^ converted to m^2 (1e4 factor)
c  HC/lambda in cgs converted to Joules (1e7 factor)
c  WAV converted from Angstrongs to cm (1e8 factor)
c  DELWAV converted from Angstrongs to nanometers (1e-1 factor)
c      EFLUX(L) = (HC*1.e8/WAV(L))*(1.e4/1.e7)*FLX/((1.e-1*DELWAV)*AGL)
c      GFLUX(L) = EFLUX(L)*S(1)
c Converting FLUX from photons/cm^2/s to photons/m^2/s/nm
c      PhEFLUX(L) = 1.e4*FLX/((1.e-1*DELWAV)*AGL)
c      PhGFLUX(L) = PhEFLUX(L)*S(1)

c      IF (L.EQ.1) write(90,118)
c      write(90,120) L,WAV(L),TAUR,EFLUX(L),
c     & GFLUX(L),S(1),PhEFLUX(L),PhGFLUX(L)
c 120  FORMAT(1X,I2,2X,1P7E14.5)
C-AP
  57  CONTINUE
c      print *,'In photo.f before far UV LOOP'
C ***** FAR UV (1216 - 1750 A) *****
      IF (LTIMES.GT.0) GO TO 62
C-AP
C-AP
C  CALCULATE CH4 OPTICAL DEPTH AT LYMAN ALPHA FOR UPPER LAYER
C-AP      TAUCH4 = TCH4(NZ)*SCH4(10)
C-AP      PRINT 121, TAUCH4
c 121  FORMAT(/'CH4 OPTICAL DEPTH OF UPPER LAYER (Lya):*',1PE10.3)

C  SCALE SFX(L) VALUES FOR SOLAR FLUX AT EARTH (FSCALE=1) OR
C  MARS (FSCALE=0.43)
      DO 900 L=1,10
  900 SFX(L)=SFX(L)*FSCALE
C-AP
  62  CONTINUE
C
c-as Methane has been added as major absorber (03/24/2004)
      DO 50 L=1,10
      DO 50 I=1,NZ
      QQ = SFX(L)*50.*EXP(-(SIGMA(1,L)*TO2(I) + SIGMA(2,L)*TCO2(I)
     &  +SIGMA(3,L)*TH2O(I)+SSO2A(L)*TSO2(I)+SCH4(L)*TCH4(I))*AM)
      if(i.eq.1) s1(L)=QQ/(SFX(L)*50.)
      PO2D(I) = PO2D(I) + QQ*SIGMA(1,L)*AGL
      PCO2D(I) = PCO2D(I) + QQ*SIGMA(2,L)*AGL
      PH2O(I) = PH2O(I) + QQ*SIGMA(3,L)*AGL
C-AP
      PSO2(I) = PSO2(I) + QQ*SSO2A(L)*AGL
c     PSO(I) = PSO(I) + QQ*SSOA(L)*AGL
      PCH4(I) = PCH4(I) + QQ*SCH4(L)*AGL    
  50  CONTINUE
C
C ***** FILL UP RATE MATRIX *****
      DO 56 I=1,NZ
      AR(23,I) = PO2D(I)
      AR(24,I) = PO2(I)
      AR(25,I) = PH2O(I)
      AR(26,I) = PO3D(I)
      AR(27,I) = PO3(I)
      AR(28,I) = PH2O2(I)
      AR(29,I) = PCO2(I)
      AR(38,I) = PH2(I)
      AR(39,I) = PHCO(I)
      AR(42,I) = PCO2D(I)
      AR(50,I) = PHO2(I)
      AR(51,I) = PCH4(I)
      AR(52,I) = PMOOH(I)
      AR(53,I) = PN2O(I)
      AR(54,I) = 1.7E-3
      AR(55,I) = PHNO3(I)
      AR(56,I) = PNO(I)
      AR(57,I) = PNO2(I)
      AR(95,I) = PHNO4(I)
      AR(101,I) = PCH3CL(I)
      AR(132,I) = PCL2(I)
      AR(133,I) = PCLO2(I)
      AR(134,I) = PHCL(I)
      AR(135,I) = PHOCL(I)
      AR(136,I) = PNOCL(I)
      AR(137,I) = PCLONO(I)
      AR(138,I) = PCLONO2(I)
      AR(142,I) = PCL2O2(I)
      AR(151,I) = PNO3(I)
      AR(153,I) = PN2O5(I)
C-AP**********************************
C-AP Adding Sulfur photochemistry
      AR(156,I) = 0.
      AR(157,I) = 0.7*PSO2(I)
      AR(158,I) = PH2S(I)
      AR(186,I) = PSO21(I)
      AR(187,I) = PSO23(I)
      AR(188,I) = PH2SO4(I)
      AR(189,I) = 0.
      AR(211,I) = PHO2(I)
C-KK	The following reaction has been renumbered from 218 to 182
      AR(182,I) = 0.3*PSO2(I)
C-AP**********************************
  56  CONTINUE
      if(N.eq.1) then
c Notice that for the 10 intervasl in the far UV SFX is in
c  photons/s/cm^2/A    
      do L1=10,1,-1      
      ESFX(L1) = (HC*1.e8/XL(L1))*(1.e4/1.e7)*SFX(L1)/1.e-1
      if(l1.eq.10) ESFX(L1)=(HC*1.e8/XL(L1))*(1.e4/1.e7)*
     & SFX(L1)*50.
      GSFX(L1) = ESFX(L1)*S1(L1)
c Converting FLUX from photons/cm^2/s to photons/m^2/s/nm
      PhESFX(L1) = 1.e4 *SFX(L1)/1.e-1
      if(l1.eq.10)PhESFX(L1)=1.e4 *SFX(L1)
      PhGSFX(L1) = PhESFX(L1)*S1(L1)
      write(90,120) L1,XL(L1),ESFX(L1),GSFX(L1),
     & S1(L1),PhESFX(L1),PhGSFX(L1)
      enddo
c      print*,'after the print of far UV'
      do L=1,108
      DELWAV = WAVU(L) - WAVL(L)
      WAV(L) = 0.5 * (WAVU(L) + WAVL(L))
c  Here the flux is converted from photons/cm^2/s to W/m^2/nm
c  cm^ converted to m^2 (1e4 factor)
c  HC/lambda in cgs converted to Joules (1e7 factor)
c  WAV converted from Angstrongs to cm (1e8 factor)
c  DELWAV converted from Angstrongs to nanometers (1e-1 factor)
c      print*,'WAV(L),DELWAV ',L,WAV(L), DELWAV
      EFLUX(L) = (HC*1.e8/WAV(L))*(1.e4/1.e7)*FLUX(L)/
     & (1.e-1*DELWAV)
      GFLUX(L) = EFLUX(L)*SS(L)
c Converting FLUX from photons/cm^2/s to photons/m^2/s/nm
      PhEFLUX(L) = 1.e4*FLUX(L)/(1.e-1*DELWAV)
      PhGFLUX(L) = PhEFLUX(L)*SS(L)
      write(90,120) L,WAV(L),EFLUX(L),GFLUX(L),
     & SS(L),PhEFLUX(L),PhGFLUX(L)
      enddo
 120  FORMAT(1X,I3,2x,F7.1,2X,1P5E14.5)
      endif
c     print*,'In photo.f AFTER print LOOP'
      LTIMES = LTIMES + 1
      RETURN
      END








