 
      SUBROUTINE PSATRAT(H2O)

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      DIMENSION HL(NZ),H2O(NZ),A(NZ)
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
      INCLUDE '../INCLUDECHEM/comSATBLK.inc'
C

      T0 = 273.15
C-AP      P0 = 6.103
      P0 = 6.103E-3
      AMV = 18.
      VAPL = 597.3
      SUBL = 677.
      R = 1.9872
      A0 = 0.553
      BK = 1.38E-16
C-AP      PS = 1.E-3 * DEN(1)*BK*T(1)
      PS = 1.E-6 * DEN(1)*BK*T(1)
C
C-AP Ask Jim why JTROP and not NZ
      DO 1 J=1,NZ
      HL(J) = SUBL
      A(J) = 0.
      IF (T(J).LT.T0) GO TO 1
      HL(J) = VAPL + A0*T0
      A(J) = A0
   1  CONTINUE

C	print*,"First JTROP loop complete. Line 2563:"
C	print*,"P0 = ", P0, " T0 = ", T0, " AMV = ", AMV, " R = ", R
C	do j = 1, NZ
C	  print*, T(J)
C	end do

C
      DO 2 J=1,NZ
      P1 = P0 * (T0/T(J))**(AMV*A(J)/R)
      P2 = EXP(AMV*HL(J)/R * (1./T0 - 1./T(J)))
      PV = P1 * P2
C-AP      P = 1.E-3 * DEN(J)*BK*T(J)
      P(J) = 1.E-6 * DEN(J)*BK*T(J)
      H2OSAT(J) = PV/P(J)
   2  CONTINUE


C-AP
C-AP CALCULATE TROPOSPHERIC H2O CONCENTRATIONS
      DO 3 J=1,JTROP
      REL = 0.77 * (P(J)/PS - 0.02)/0.98
      RELH(J) = REL
   3  H2O(J) = REL * H2OSAT(J) 
C
      RETURN
      END


