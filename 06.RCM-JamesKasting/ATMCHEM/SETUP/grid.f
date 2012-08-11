
      SUBROUTINE GRID
      INCLUDE '../INCLUDECHEM/parNZ.inc'
      INCLUDE '../INCLUDECHEM/comABLOK.inc'
C
C ***** SET UP THE VERTICAL GRID ZS *****
      DZ = 1.E5
      DO 1 I=1,NZ
       Z(I) = (I - 0.5)*DZ
   1  CONTINUE
C
      RETURN
      END
