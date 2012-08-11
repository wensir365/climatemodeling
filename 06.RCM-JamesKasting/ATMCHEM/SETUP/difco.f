      SUBROUTINE DIFCO

      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
C-AP
      NZ1 = NZ - 1
C                                                                              
C ***** DK(I) = K*N AT GRID STEP I+1/2 *****
      DO 1 I=1,NZ1
      EDDAV = SQRT(EDD(I)*EDD(I+1))
      DENAV = SQRT(DEN(I)*DEN(I+1))
   1  DK(I) = EDDAV*DENAV
C
      RETURN
      END
