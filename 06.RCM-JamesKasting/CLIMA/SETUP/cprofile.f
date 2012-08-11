
      SUBROUTINE CPROFILE(TSTRAT,P,T,DZ,FH2O,FCO2V,BETA,Idry,JCOLD)
      
      INCLUDE '../INCLUDECLIM/parND.inc'
      DIMENSION P(ND),T(ND),DZ(ND),FH2O(ND),BETA(ND),
     & FLAGCONVEC(ND),FCO2V(ND)
      INCLUDE '../INCLUDECLIM/comCBLOK.inc'
      INCLUDE '../INCLUDECLIM/comEBLOK.inc'
      INCLUDE '../INCLUDECLIM/comCO2BLOK.inc'

C
C   THIS SUBROUTINE CALCULATES TROPOSPHERIC TEMPERATURES AND
C   HUMIDITIES BY CALLING CONVEC.
C
      BETA1 = 1.
      BETA2 = 1.
      ITROP = 1

c Settinbg all the temperatures at the stratospheric temperature
      do i=1, ND
       FLAGCONVEC(i)=0.
       T(i)=TSTRAT
       FCO2V(i)=FCO2
      enddo

C Calculating tropospheric  temperatures and water
C   SOLVE FROM THE GROUND UP
      PG = P(ND)
      CALL SATRAT(TG,PSAT)
      T(ND) = TG
      FH2O(ND) = RELHUM(PG) * PSAT/PG
      IF(PSAT.GT.PG) FH2O(ND) = POCEAN/PG
    
      DO 2 J1=ND,JCOLD+1,-1
       T1 = T(J1)
       F1 = FH2O(J1)
       P1 = P(J1)
       P2 = P(J1-1)
       DZP = DZ(J1)
       FC1 = FCO2
       CALL CONVEC(T1,T2,P1,P2,F1,FH2,FC1,FC2,DZP,ITROP,cflag,Idry)
       FLAGCONVEC(J1)=cflag
       BETA(J1) = BETA1
       BETA(J1-1) = BETA2
       JCOLD = J1
       IF(T2.LT.T(J1-1)) GOTO 6
       T(J1-1) = T2
       FCO2V(J1-1)= FC2
   2   FH2O(J1-1) = FH2

c Calculating stratospheric water contents

   6    IF(IMW.NE.2) GO TO 4

      DO 3 J=1,(JCOLD-1)
       TJ = T(J)
       PJ = P(J)
       CALL SATRAT(TJ,PSAT)
       FH2 = RELHUM(PJ) * PSAT/PJ
   3   FH2O(J) = FH2 
      RETURN
C
   4  CALL SATRAT(T(JCOLD),PSAT)
      FSAT= PSAT/P(JCOLD)
      DO 5 J=(JCOLD-1),1,-1
      CALL SATRAT(T(J),PSAT)
      FSATUR = (PSAT/P(J))
      FH2O(J) = FSATUR*RELHUM(P(J))
   5  FH2O(J) = AMIN1(FSAT,FH2O(J+1))

      RETURN
      END
C ----------------------------------------------------------------





