
      SUBROUTINE SEDMNT(FSULF,N)

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/parNQ_NQT.inc'

      DIMENSION FSULF(NZ),TAUTRN(NZ),RHOP(NZ)
      DIMENSION TAURAN(NZ,3),ALAM(NZ),TAUCPK(NZ,3),ETA(NZ)
      DIMENSION CUNING(NZ,3),amass(NZ,3),THERMSP(NZ,3),
     2  TAURELAXC(NZ,3),TAURELAX(NZ,3),AFPL(NZ,3),delta(NZ,3),
     3  BETAF(NZ,3),TAUCPKF(NZ,3)

      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comGBLOK.inc'
      INCLUDE '../INCLUDECHEM/comNBLOK.inc'
      INCLUDE '../INCLUDECHEM/comAERBLK.inc'

C
C   THIS SUBROUTINE CALCULATES FALL VELOCITIES AND ESTIMATES PARTICLE
C   SIZE BASED ON THEIR COAGULATION LIFETIMES
C
C   CONSTANTS from Pruppacher and Klett page 450
          A = 1.257
          B = 0.4
          C = 1.1
C-AP
      BK = 1.38*1.E-16
      PI = 3.14159
      NZ1 = NZ - 1
C
      DO J=1,NZ
      ALAM(J) = 1.63E14/DEN(J)
      ETA(J)=ABS((1.718 + 0.0049*(T(J)-273.) - 1.2
     & *(1.E-5)*(T(J)-273.)*(T(J)-273.))*1.E-4)
      ENDDO

C-AP Because we have only one aerosol we do not do DO-loop
      K = 1
C-AP      DO 10 K=1,3
C   (1 = SULFATE, 2 = S8, 3 = HYDROCARBON)
C
C-AP ESTIMATION OF THE AEROSOL FREE PATH LENGTH
      DO J = 1,NZ
       ALPH = A + B*EXP(-C*RPAR(J,K)/ALAM(J))
       CUNING(J,K) = 1 + ALPH*ALAM(J)/RPAR(J,K)

C-AP Here we assume that the density of aerosol is 1 g/cm3
C-AP Notation is similar Fusch 1964
       amass(J,K) = (4./3.)*PI*RPAR(J,K)**3
       THERMSP(J,K) = SQRT((8*BK*T(J))/(pi*amass(J,K)))
       TAURELAXC(J,K)=2*RPAR(J,K)*RPAR(J,K)/(9*ETA(J))
       TAURELAX(J,K) = TAURELAXC(J,K)*CUNING(J,K)
       AFPL(J,K) = THERMSP(J,K)*TAURELAX(J,K)
       delta(J,K) = (((2*RPAR(J,K)+AFPL(J,K))**3 - (4*
     & RPAR(J,K)*RPAR(J,K)+AFPL(J,K)*AFPL(J,K))**1.5)/
     & (6*RPAR(J,K)*AFPL(J,K)) - 2*RPAR(J,K))*SQRT(2.)
      ENDDO

C-AP Calculation of the correction to the coagulation kernel
      DO J = 1,NZ
       BETAF(J,K) = 1/(RPAR(J,K)/(RPAR(J,K)+delta(J,K)/2)
     & +PI*AFPL(J,K)/(2*SQRT(2.)*RPAR(J,K)))
      ENDDO
C-AP      PRINT *, 'TESTSEDIM3'
C-AP ******************************************************
C ESTIMATE COAGULATION AND SEDIMENTATION LIFETIMES (TOON AND FARLOW, 1981)
      DO 1 J=1,NZ
C-AP      PRINT *, 'AERSOLINSED=', AERSOL(J,K), 'LEVEL=',J
C-AP      PRINT *, 'RPARINSED=', RPAR(J,K), 'LEVEL=',J
      TAUC(J,K) = 1.E6/(AERSOL(J,K)*SQRT(RPAR(J,K)))
      TAUCPK(J,K) = 3*ETA(J)/(4*AERSOL(J,K)*BK*
     & T(J)*CUNING(J,K))
      TAUCPKF(J,K) = TAUCPK(J,K)/BETAF(J,K)
      TAUC(J,K) = TAUCPKF(J,K)
      TAURAN(J,K) = 1./(RAINGC(LH2SO4,J) + 1.E-20)
   1  TAUSED(J,K) = HSCALE(J)/WFALL(J,K)

C-AP
C   FIND MINIMUM OF DIFFUSION AND SEDIMENTATION LIFETIMES, THEN SCALE
C   PARTICLE SIZES
      DO 2 J=1,NZ
      TAUTRN(J) = AMIN1(TAUSED(J,K),TAUEDD(J))
      TAUTRN(J) = AMIN1(TAUTRN(J),TAURAN(J,K))
      RPAR(J,K) = RPAR(J,K) * (TAUTRN(J)/TAUC(J,K))**0.25
      IF (K.EQ.3) THEN
        RPAR(J,K) = AMAX1(RPAR(J,K),1.3E-7)
      ELSE
        RPAR(J,K) = AMAX1(RPAR(J,K),1.E-5)
      ENDIF
   2  CONTINUE
C
C Write TAUTRN on last iteration
C
C-AP      IF ( K.EQ.3) THEN
C-AP      WRITE(18,900)
C-AP 900  FORMAT(/'# HC(aerosol)',/'#  Z',3X,'TAUC',6X,'TAUCPK',
C-AP     2   4X,'TAUSED',4X,'TAUTRN',4X,'TAUEDD')
C-AP 900  FORMAT(/'# HC(aerosol)',/'#  Z',3X,'TAUC',6X,'TAURAN',
C-AP     2   4X,'TAUSED',4X,'TAUTRN',4X,'TAUEDD')
C-AP      DO 22 J=1,NZ
C-AP      WRITE(18,901) Z(J),TAUC(J,3),TAUCPK(J,3),TAUSED(J,3),
C-AP     2 TAUTRN(J), TAUEDD(J)
C-AP      WRITE(18,901) Z(J),TAUC(J,3),TAURAN(J,3),TAUSED(J,3),
C-AP     2 TAUTRN(J),
C-AP     3  TAUEDD(J)
C-AP 901  FORMAT(1P6E10.3)
C-AP  22  CONTINUE
C-AP      ENDIF
C
C   DON'T ALLOW PARTICLES TO DECREASE IN SIZE AT LOW ALTITUDES
      DO 3 I=1,NZ1
      J = NZ - I
   3  RPAR(J,K) = AMAX1(RPAR(J,K),RPAR(J+1,K))
C
C   COMPUTE PARTICLE-TO-GAS CONVERSION FACTORS AND DENSITIES
      IF (K.NE.1) GO TO 5
      DO 4 J=1,NZ
      R = RPAR(J,K)
      CONVER(J,K) = 4.6E7*FSULF(J) * (R/1.E-5)**3
   4  RHOP(J) = 1. + 0.8*FSULF(J)
      GO TO 7
C
   5  CONTINUE
      IF (K.NE.2) GO TO 9
      DO 6 J=1,NZ
      R = RPAR(J,K)
      CONVER(J,K) = 2.03E7 * (R/1.E-5)**3
   6  RHOP(J) = 2.07
      GO TO 7
C
   9  CONTINUE
      DO 11 J=1,NZ
      R = RPAR(J,K)
      CONVER(J,K) = 7.06E7 * (R/1.E-5)**3
  11  RHOP(J) = 1.4
   7  CONTINUE
C
C     NOW COMPUTE FALL VELOCITIES
      DO 8 J=1,NZ
       R = RPAR(J,K)
C-AP      ETA = 1.77E-4 * SQRT(T(J)/288.)
C-AP  From Prupacher & Klett
      F1 = 2./9. * RHOP(J)*R*R*G/ETA(J)
      ALPH = A + B*EXP(-C*R/ALAM(J))
   8  WFALL(J,K) = F1*(1. + ALAM(J)*ALPH/R)
C-AP****************************************************************
C-AP  10  CONTINUE
C
      RETURN
      END
