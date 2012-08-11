C------------------------------------------------------------------------------
C-KK	Code added Sept 27 2001 for integration with the photochemical code
C-KK 	atm_chem.f   
C------------------------------------------------------------------------------
  
      SUBROUTINE OUTPUT_PHOTO(T, FI, water, ALT,O3save)
      INCLUDE '../CLIMA/INCLUDECLIM/parND.inc'
      INCLUDE '../ATMCHEM/INCLUDECHEM/parNZ.inc'
     
      DIMENSION T(ND), FI(4,ND), alt_new(NZ), T_new(NZ)
      DIMENSION water(NZ), ALT(ND)
      INCLUDE '../CLIMA/INCLUDECLIM/comCBLOK.inc'

	alt_new(NZ) = 0.5
	do i = (NZ-1), 1, -1
C	  interpolate temperatures at correct altitudes
	  alt_new(i) = alt_new(i+1) + 1.
	end do

	ISTART = ND
	
	DO J = NZ, 1, -1
	 DO i = ISTART,1,-1
	   IS = i
	   IF (ALT(IS) .GT. alt_new(J)) GOTO 350
	 END DO
  350  CONTINUE
	IS1 = IS+1
	ISTART = IS
	DZI = ALT(IS) - ALT(IS1)
	DZJ = ALT(IS) - alt_new(J)
	FR = DZJ/DZI
C-KK			begin temperature interpolation
	T1log=ALOG(T(IS))
	T2log=ALOG(T(IS1))
	T_temp = FR*T2log + (1-FR)*T1log
	T_new(J) = EXP(T_temp)
C-KK			begin water interpolation
	FH1log = ALOG(FI(1,IS))
	FH2log = ALOG(FI(1,IS1))
	water_temp = FR*FH2log + (1-FR)*FH1log
	water(J) = EXP(water_temp)
	END DO


C-KK	This statement reverses the indexing (has to be done here and not
C-KK	prior because earlier reversal would confuse the issue during the
C-KK	interpolation procedure. The interpolated altitudes, temperatures,
C-KK	water profile and ozone profile are recorded to file fromClima2Photo.dat

	JCOLD = 0

	DO J = NZ, 2, -1
	   WRITE(54,*) alt_new(J), T_new(J), water(J)
 	   IF (JCOLD .EQ. 0) THEN
		IF (T_new(J) .LT. T_new(J-1)) JCOLD = NZ - J
c               IF (water(J) .LT. water(J-1)) JCOLD = NZ - J
	   END IF
	END DO

	J=1
	WRITE(54,*) alt_new(J), T_new(J), water(J)
        CLOSE(54)

	FCH4 = FI(3,ND)
	FCO2 = FI(2,ND)	
        WRITE(55,102) FAR, FCH4, FCO2, FO2, JCOLD, O3save
  102	FORMAT(1PE10.3,10x,'!Argon'/,1PE10.3,10x,'!Methane'/
     &  ,1PE10.3,10x,'!Carbon Dioxide'/,1PE10.3,10x,'!Oxygen'/,
     &  I3,17x,'!Tropopause layer',/,1PE12.4,10x,'!O3 column depth')
        close(55)
	END
