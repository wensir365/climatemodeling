
      SUBROUTINE TWOSTR(WAV,SIGR,U0,SO3,SO2,SCO2,SH2O,SSO2,SH2S)
C
C   This is my version of the Toon et al. 2-stream code.  (Ref.: JGR
C   94, 16287, 1989).  It vectorizes over height, rather than wavelength,
C   and is designed to work with PRIMS3 and its companion photochemical 
C   models.
C
C   For now, at least, it is hardwired as the quadrature approximation.
C     NP is the number of different types of particles
C
      INCLUDE '../INCLUDECHEM/parNZ.inc'
 
      PARAMETER(NZZ1=NZ+1, NZ2=2*NZ, NP=2)
      DIMENSION SO2(NZ),SO3(NZ)
      DIMENSION TAU(NZ), TAUC1(NZ),G1(NZ), GAM1(NZ), GAM2(NZ), GAM3(NZ),
     1  GAM4(NZ), ALAM(NZ), CGAM(NZ), E1(NZ), E2(NZ), E3(NZ), E4(NZ),
     2  CP0(NZ), CPB(NZ), CM0(NZ), CMB(NZ), Y1(NZ), Y2(NZ), W0(NZ), 
     3  TAUSG(NZ), TAUSP(NZ), DIRECT(NZZ1), AMEAN(NZZ1), FMT(NZ)
      DIMENSION A(NZ2), B(NZ2), D(NZ2), E(NZ2), Y(NZ2)
      DIMENSION W0P(NP),QEXT(NP)
C
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comAERBLK.inc'
      INCLUDE '../INCLUDECHEM/comCBLOK.inc'
      INCLUDE '../INCLUDECHEM/comTBLOK.inc'

C     U1 = 0.5  (Eddington value)
      SQ3 = SQRT(3.)
      PI = 3.14159
      U1 = 1./SQ3
      U0M = 1./U0
      U0M2 = U0M*U0M
      U1M = 1./U1
      GP = 0.8
      NZM1 = NZ - 1
      MZ2 = NZ2
C
C   Particle 1 is sulfate, 2 is S8
      W0P(1) = 1.
      W0P(2) = 0.5
      IF (WAV .GT. 3500.) W0P(2) = 1.
      DO 20 J=1,NP
  20  QEXT(J) = 2.
C   Above value of Qext is valid for large particles
C
C   Calculate the optical depths of the different layers.  TAUA is absorption,
C   TAUSG is scattering by gases, TAUSP is scattering by particles, TAU is
C   total extinction.  Note that the grid for this subroutine is numbered 
c   from top to bottom, whereas the main program is numbered from bottom to
c   top.
C   First do gases 
      DO 21 I=1,NZ
      N = NZZ1 - I
      TAUSP(N) = 0.
      TAUSG(N) = SIGR*DEN(I)*DZ
      TAUA = ( SO3(I)*O3(I) + SO2(I)*O2(I) + SCO2*CO2(I) + SH2O*H2O(I)
     2   ) * DEN(I)*DZ
  21  TAU(N) = TAUA + TAUSG(N)
C
C   Now do particles.  Must combine their scattering optical depths into a
C   single array in order to calculate W0 and G.  Don't need any arrays for
C   pure absorption
      J=1
      DO 1 I=1,NZ
      N = NZZ1 - I
      TAUP = QEXT(J)*PI*RPAR(I,J)*RPAR(I,J)*AERSOL(I,J)*DZ
      TAUSP(N) = TAUSP(N) + W0P(J)*TAUP
   1  TAU(N) = TAU(N) + TAUP

C   Calculate W0 and G by averaging over Rayleigh and Mie scatterers.  Avoid
C   letting W0 equal exactly 1.
      DO 22 N=1,NZ
      W0(N) = (TAUSG(N) + TAUSP(N))/TAU(N)
      W0(N) = AMIN1(W0(N), 0.999)
  22  G1(N) = TAUSP(N)*GP/(TAUSP(N) + TAUSG(N))
C
C-RP  A. Pavlov's modification of code, acc. to Joseph et al. 1976
      DO 23 N=1,NZ
      FMT(N) = G1(N)*G1(N)
      TAU(N) = TAU(N)*(1. - W0(N)*FMT(N))
      W0(N) = W0(N)*(1. - FMT(N))/(1. - W0(N)*FMT(N))
  23  G1(N) = G1(N)/(1. + G1(N))
C-RP  ***end*** 
C   Calculate the gamma's, lambda's, and e's
      DO 2 N=1,NZ
C     GAM1(N) = (7. - W0(N)*(4.+3.*G1(N)))/4.
C     GAM2(N) = - (1. - W0(N)*(4.-3.*G1(N)))/4.
C     GAM3(N) = (2. - 3.*G1(N)*U0)/4.
C   (Eddington values above; quadrature values below)
C
      GAM1(N) = SQ3*(2. - W0(N)*(1.+G1(N)))/2.
      GAM2(N) = SQ3*W0(N)*(1.-G1(N))/2.
      GAM3(N) = (1. - SQ3*G1(N)*U0)/2.
      GAM4(N) = 1. - GAM3(N)
C
      ALAM(N) = SQRT(GAM1(N)*GAM1(N) - GAM2(N)*GAM2(N))
      CGAM(N) = (GAM1(N) - ALAM(N))/GAM2(N)
      EMLT = EXP(-ALAM(N)*TAU(N))
C
      E1(N) = 1. + CGAM(N)*EMLT
      E2(N) = 1. - CGAM(N)*EMLT
      E3(N) = CGAM(N) + EMLT
   2  E4(N) = CGAM(N) - EMLT
C
C   Calculate A, B, and D, i.e. the coefficients of the tridiagonal
C      matrix
C   Top of atmosphere
      A(1) = 0.
      B(1) = E1(1)
      D(1) = -E2(1)
C
C   Odd coefficients
      DO 3 N=1,NZM1
      L = 2*N + 1
      A(L) = E2(N)*E3(N) - E4(N)*E1(N)
      B(L) = E1(N)*E1(N+1) - E3(N)*E3(N+1)
   3  D(L) = E3(N)*E4(N+1) - E1(N)*E2(N+1)
C
C   Even coefficients
      DO 4 N=1,NZM1
      L = 2*N
      A(L) = E2(N+1)*E1(N) - E3(N)*E4(N+1)
      B(L) = E2(N)*E2(N+1) - E4(N)*E4(N+1)
   4  D(L) = E1(N+1)*E4(N+1) - E2(N+1)*E3(N+1)
C
C   Bottom of atmosphere
      A(NZ2) = E1(NZ) - ALB*E3(NZ)
      B(NZ2) = E2(NZ) - ALB*E4(NZ)
      D(NZ2) = 0.
C
C   Now, set up the RHS of the equation:
C   TAUC1(N) is the optical depth above layer N
      TAUC1(1) = 0.
      DO 5 N=2,NZ
   5  TAUC1(N) = TAUC1(N-1) + TAU(N-1)
C
C   DIRECT(N) is the direct solar flux at the top of layer N.  Values
C   are normalized to unity.  DIRECT(NZZ1) is the direct flux at the ground.
      DIRECT(1) = 1.
      DO 6 N=1,NZ
      ET0 = EXP(-TAUC1(N)/U0)
      ETB = ET0 * EXP(-TAU(N)/U0)
      DIRECT(N+1) = ETB
      DENOM = ALAM(N)*ALAM(N) - U0M2
      FACP = W0(N) * ((GAM1(N)-U0M)*GAM3(N) + GAM4(N)*GAM2(N))
      FACM = W0(N) * ((GAM1(N)+U0M)*GAM4(N) + GAM2(N)*GAM3(N))
C
      CP0(N) = ET0*FACP/DENOM
      CPB(N) = ETB*FACP/DENOM
      CM0(N) = ET0*FACM/DENOM
   6  CMB(N) = ETB*FACM/DENOM
      SSFC = ALB*U0*DIRECT(NZZ1)
C
C   Odd coefficients
      E(1) = - CM0(1)
      DO 7 N=1,NZM1
      L = 2*N + 1
   7  E(L) = (CP0(N+1)-CPB(N))*E3(N) + (CMB(N)-CM0(N+1))*E1(N)
C
C   Even coefficients
      DO 8 N=1,NZM1
      L = 2*N
   8  E(L) = (CP0(N+1)-CPB(N))*E2(N+1) - (CM0(N+1)-CMB(N))*E4(N+1)
      E(NZ2) = SSFC - CPB(NZ) + ALB*CMB(NZ)
C
C   Call the tridiagonal solver (from LINPACK).  E is the RHS of the matrix
C   equation on input and is the solution vector Y on output
      CALL SGTSL(MZ2,A,B,D,E,NFLAG)
      IF (NFLAG .NE. 0) PRINT 100, NFLAG
 100  FORMAT(/1X,'Tridiagonal solver failed in TWOSTR, NFLAG =',I4)
C
      DO 9 N=1,NZ
      L = 2*N
      L1 = L-1
      Y1(N) = E(L1)
   9  Y2(N) = E(L)
C
C   Calculate the mean intensities, AMEAN(N), at the boundaries between
C   the layers.  AMEAN(N) is the intensity at the top of layer N.
      AMEAN(1) = U1M * (Y1(1)*E3(1) - Y2(1)*E4(1) + CP0(1)) + 1.
      DO 10 N=1,NZ
      I = NZZ1 - N
  10  AMEAN(N+1) = U1M * (Y1(N)*(E1(N)+E3(N)) + Y2(N)*(E2(N)+E4(N)) 
     1  + CPB(N) + CMB(N)) + DIRECT(N+1)
c-mm
c-mm  Reset any AMEAN values that may go negative.  Check error file
c-mm  to be sure this only happens near the ground where AMEAN ~0. 
c-mm  (allowing us to simply take the absolute value)
c-mm
      DO 12 N=1,NZZ1
      IF (AMEAN(N).LT.0.0) THEN
c        WRITE(82,103) WAV,N,AMEAN(N)
         AMEAN(N) = ABS(AMEAN(N))
      ENDIF
 12   CONTINUE
 103  FORMAT('WAVE =',F6.1,' AMEAN(',I3,')=',1PE11.3)
C
C   Convert back to main program grid.  S(I) is the mean intensity at the
C   midpoint of layer I.
      DO 11 I=1,NZ
      N = NZZ1 - I
 11   S(I) = SQRT(AMEAN(N)*AMEAN(N+1))
 
      RETURN
      END
