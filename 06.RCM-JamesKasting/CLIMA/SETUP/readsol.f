      SUBROUTINE READSOL

c-as This subroutine reads the incoming flux from the Sun and the
c-as parameters for the subroutine SOLAR and CONVEC
      INCLUDE '../INCLUDECLIM/parND.inc'
      INCLUDE '../INCLUDECLIM/parNF.inc'
      INCLUDE '../INCLUDECLIM/parNS_NS4.inc'
      INCLUDE '../INCLUDECLIM/parNT.inc'
      INCLUDE '../INCLUDECLIM/parMT.inc'
      INCLUDE '../INCLUDECLIM/parNSOL_NSOLUV.inc'
       
      DIMENSION ALPHAZ(4,2),BETAZ(4,2),NPROB(2),
     &  NG(2),SIGG(4,2,NSOL)

      INCLUDE '../INCLUDECLIM/comABLOK1.inc'
      INCLUDE '../INCLUDECLIM/comSOLARBLK.inc'
      INCLUDE '../INCLUDECLIM/comCH4BLOCK.inc'
      INCLUDE '../INCLUDECLIM/comFBLOK.inc'
      INCLUDE '../INCLUDECLIM/comGBLOK1.inc'
      INCLUDE '../INCLUDECLIM/comPRESS.inc'
      INCLUDE '../INCLUDECLIM/comCBLOK.inc'
      INCLUDE '../INCLUDECLIM/comSIGUVS.inc'

      DATA ALPHACH4NEW/0.08566225,0.1803808,0.23395695,0.23395695,
     & 0.1803808,0.08566225/

C READ THE STEAM TABLE
      do i=1,4
       READ(33,*)
      enddo
      READ(33,300) (TTAB(I),PVAP(I),DPVAP(I),SVC(I),DSV(I),DSC(I),
     2  RHOV(I),DRHOV(I),I=1,NT)
 300  FORMAT(1P8E12.5)
      DSV(NT) = 0.

C   READ THE HIGH TEMPERATURE (UNSATURATED) STEAM TABLE
      do i=1,2
       READ(33,*)
      enddo
      DO 58 I=1,7
      M1 = 10*(I-1) + 1
      M2 = M1 + 9
      READ(33,302) (PCP(M),M=M1,M2)
 302  FORMAT(5X,10F11.0)
C
      DO 59 N=1,75
      READ(33,303) TCP(N),(DPDTL(M,N),M=M1,M2)
      READ(33,304) (DRDTL(M,N),M=M1,M2)
      READ(33,304) (BETAM(M,N),M=M1,M2)
  59  READ(33,304)
  58  CONTINUE 
 303  FORMAT(1X,F6.0,1P10E11.4)
 304  FORMAT(7X,1P10E11.4)
      close(33)

c    READ THE SATURATED CO2 TABLE
      do ii=1,5
       READ(34,*)
      enddo
      READ(34,306) (TCTAB(I),PCVAP(I),BETASC(I),DPCVAP(I),DRCVAP(I),
     2  SVSC(I),DSVC(I),DSCC(I),I=1,MT)
 306  FORMAT(1X,F7.2,7E10.3)
C
c  READ THE UNSATURATED CO2 TABLE
      READ(34,177)
      DO 76 L=1,6
      IS = 1 + 6*(L-1)
      IF = IS + 5
      READ(34,307) (PCC(I),I=IS,IF)
 307  FORMAT(///10X,6(F4.1,7X))
      DO 68 K=1,25
      READ(34,308) TCC(K),(BETAMC(K,I),I=IS,IF)
 308  FORMAT(/F5.1,2X,6E11.4)
      READ(34,309) (CPC(K,I),I=IS,IF)
 309  FORMAT(7X,6E11.4)
      READ(34,309) (DVDTC(K,I),I=IS,IF)
  68  READ(34,309) (DVDPC(K,I),I=IS,IF)
  76  CONTINUE
      close(34)

***INPUT TO SOLAR-TWO-STREAM*** (DMW)
C
C   TAUAER = OPTICAl DEPTH OF AEROSOL (IF NO AEROSOL, TAUAER = 0)
C   SIGERT = EXTINCTION COEFFICIENT FOR AEROSOL (IF NO AEROSOL, SIGERT = 0)
C   FMA = "F" (SECOND MOMENT/5) FOR AEROSOL (IF NO AEROSOL, FMA = 0)
C
      DO 128 I=1,NSOL
         TAUAER(I) = 0.
         SIGERT(I) = 0.
         FMA(I) = 0.
 128  CONTINUE
C
C  READING of new data (Added by Pat Kasting oct-2005)
C  I corresponds to interval 22-38
C  J corresponds to Temps 112,188,295 [Kelvin]
C  K corresponds to pressures 0.0001,0.001,0.01,0.1,1.0 [Bars]
	
	DO I=1,17
	 READ(40,*)
	 READ(40,*)
	  DO J=1,3
	   READ(40,*)
	   READ(40,*)
          DO K=1,5
	     READ(40,902) (BETACH4NEW(I,J,K,L),L = 1,6)
	    END DO
	  END DO
	END DO
	close(40)
902   FORMAT(15X,1PE11.5,2X,E11.5,2X,E11.5,2X,
     1 E11.5,2X,E11.5,2X,E11.5)
C
C  READING of Kathy Rages data (near IR CH4 exponential sums)
      do i=1,4
       READ(35,*)
      enddo
      DO 173 I = 1,17
       READ(35,182) GAMMAEXP188(I)
       READ(35,176) (ALPHACH4T188(K,I),K = 1,4)
       READ(35,179) (BETACH4T188(K,I),K = 1,4)
 173  CONTINUE
      do i=1,4
       READ(35,*)
      enddo 
      DO 174 I = 1,17
       READ(35,182) GAMMAEXP295(I)
       READ(35,176) (ALPHACH4T295(K,I),K = 1,4)
       READ(35,179) (BETACH4T295(K,I),K = 1,4)
 174  CONTINUE
     
c Reading Karkoshka data for CH4
      do i=1,5
        READ(35,*)
      enddo
      DO 172 I = 1,21
       READ(35,181) (ALPHACH4Kark(K,I),K = 1,4)
       READ(35,181) (BETACH4Kark(K,I),K = 1,4)
       READ(35,177)
 172  CONTINUE
     
 176  FORMAT(4X,F6.4,4x,F6.4,4x,F6.4,4x,F6.4)
 177  FORMAT(/)
 178  FORMAT(//)
 179  FORMAT(F9.5,2x,F9.4,1x,F9.3,1x,F12.2)
 181  FORMAT(2x,F8.5,2x,F8.5,2x,F8.5,2x,F8.5)
 182  FORMAT(4x,F8.5)

C   READ EXPONENTIAL SUM DATAFILES

**** H2O parameters for exponential sums
      DO 56 I=1,30
       READ(35,501)
  56   READ(35,500) ((BETH2O(K,L,I),K=1,4),L=1,4)
      DO 156 I=31,NF
       READ(35,501)
 156   READ(35,500) ((BETH2O(K,L,I),K=1,4),L=1,5)
   
C ***** TEMPORARY FILL FOR H2O EXP SUMS AT 10 BARS *****
      DO 57 K=1,4
      DO 57 I=1,30
  57  BETH2O(K,5,I) = BETH2O(K,4,I)

**** CO2 parameters for exponential sums
      READ(35,*)
      READ(35,*)
      DO 54 I=9,23
      READ(35,501)
  54  READ(35,500) ((BETCO2(K,L,I),K=1,4),L=1,5)
      DO 55 I=27,45
      READ(35,501)
  55  READ(35,500) ((BETCO2(K,L,I),K=1,4),L=1,5)
      READ(35,501)
      READ(35,500) ((BETCO2(K,L,48),K=1,4),L=1,5)

 500  FORMAT(12X,E14.8,5X,E14.8,5X,E14.8,5X,E14.8)
 501  FORMAT(////)
      close(35) 

C   CONVERT TO CM2/GM
      DO 34 K=1,4
      DO 34 L=1,5
      DO 34 I=1,NF
      BETH2O(K,L,I) = BETH2O(K,L,I) * 6.023E23/18.
  34  BETCO2(K,L,I) = BETCO2(K,L,I) * 6.023E23/44.
C
C   FILL UP SOLAR ABSORPTION MATRICES
      DO 35 K=1,4
      DO 35 L=1,5
      DO 35 I=14,NSOL
      J = 69 - I
      BETIR1(K,L,I) = BETH2O(K,L,J)
  35  BETIR2(K,L,I) = BETCO2(K,L,J)
C
C   Read Solar Data - formerly Ackerman's routine SOLAR (7-95 DMW)
C   First find FCO2.
c   SOLINT = SOLAR FLUX (ERG/CM^2/S)
C   NG = TAGS FOR THE GASES THAT ABSORB IN A PARTICULAR WAVELENGTH

      READ(32,*)
      do j=1,NSOLUV  
       read(32,*)XLL,XLU,SOLUV(j)
       XLAMBUV(j) = (XLL+XLU)*1.e-4/2.
      enddo
      read(32,*)
      DO 435 I=1,NSOL
      READ(32,*) IJUNK
      READ(32,*) XLAMBDA(I),(NPROB(L),L=1,2),SOLINT(I),(NG(L),L=1,2)
         DO 436 K=1,4
            READ(32,571) (ALPHAZ(K,L),L=1,2),(BETAZ(K,L),L=1,2)
 571        FORMAT(1X,2(F10.8,2X),2E14.7)
 436     CONTINUE
         DO 437 L=1,2
            NPR(L,I) = NPROB(L)
            NGAS(L,I) = NG(L)
            DO 438 K=1,4    
               SIGG(K,L,I) = BETAZ(K,L)
               WGHT(K,L,I) = ALPHAZ(K,L)
 438        CONTINUE
 437     CONTINUE  
 435  CONTINUE
      close(32)
      IF (FCO2 .LT. 0.1) THEN
         NPR(2,15) = 1
         NGAS(2,15) = 4
         WGHT(1,2,15) = 1.
         WGHT(2,2,15) = 0.
         SIGG(1,2,15) = 0.02
      END IF
      DO 440 K=1,4
         DO 441 L=1,5
            DO 442 I=1,13
               BETIR1(K,L,I)=SIGG(K,1,I)
 442        CONTINUE
            DO 443 I=1,20
               BETIR2(K,L,I)=SIGG(K,2,I)
 443        CONTINUE
 441     CONTINUE
 440  CONTINUE

*** Read the near UV absorption coefficients for O2, O3 and CO2
*** (addition by Antigona Segura 08/2005)
      do jj=1,3
        read(39,*)
      enddo
      do i=1,16
        read(39,*)x,x,x,siguvO3(i),siguvO2(i),siguvCO2(i)
      enddo
      close(39)
C
      RETURN
      END

