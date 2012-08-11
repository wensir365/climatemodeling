
      SUBROUTINE IR(T,P,CGAS,FI)

c  This subroutine calculates the infrared flux

      INCLUDE '../INCLUDECLIM/parND.inc'
      INCLUDE '../INCLUDECLIM/parNF.inc'
      INCLUDE '../INCLUDECLIM/parNGS.inc'
      
      INCLUDE '../INCLUDECLIM/comIRDATA.inc'
      INCLUDE '../INCLUDECLIM/comIRBLK.inc'
      INCLUDE '../INCLUDECLIM/comCBLOK.inc'
      INCLUDE '../INCLUDECLIM/comALTBLOK.inc'
      INCLUDE '../INCLUDECLIM/comCPART.inc'
      INCLUDE '../INCLUDECLIM/comCONS.inc'
      INCLUDE '../INCLUDECLIM/comWAVE.inc'

      DIMENSION P(ND),BPLANCK(ND),CPR(NF),TPR(NF),TPRIND(ND), 
     2 TAUGCO2(ND),TAUGH2O(ND),TAUGCH4(ND),TAUGIR(ND),FUP(ND),FDN(ND),
     3 FUPA(ND),FDNA(ND),T(ND),CGAS(ND,NGS),TAUCONTIN(ND)
      REAL KAPPALAYER(55,8,3,ND),kappa(55,8,3),KAPPALAYEROZ(8,ND),
     2 WEIGHTOZC(8),KAPPAOZC(8),TAUGOZ(ND),TAUCONTINT(NF),CNUT,FI(4,ND),
     3 TAUTOTAL(NF),TRANSLAYER(ND),TWGHTT(8),OMG0IR(ND-1),ASYIR(ND-1),
     4 TAULAMIR(ND-1),TAUAEXTIR(ND-1),TAUASIR(ND-1),TAUSIR(ND-1)
      
      DATA HP,SIGMA/6.63E-27, 5.67E-5/

C   PRESSURE-INDUCED CO2 ABSORPTION (FROM JIM POLLACK)
      DATA CPR/4.3E-5, 3.8E-5, 1.2E-5, 2.8E-6, 7.6E-7, 4.5E-7, 2.3E-7,
     2  5.4E-7, 1.6E-6, 10*0., 4.E-7, 4.E-6, 1.4E-5, 1.0E-5,
     3  1.2E-6, 2.0E-7, 5.0E-8, 3.0E-8, 28*0./
      DATA TPR/3.4, 2.2, 1.9, 52*1.7/

       INTEGER COUNTERIR
C
       COUNTERIR = 0
       SRFALBIR = 0.
       NLAYERS = ND - 1
       BCON = 2.*HP/C/C
       HK = HP/BK
C       np = 1.E+1

C Read the IR exponential sums
c      CALL IREXPSUMS(WEIGHT,xkappa)

      DO 7 IL = 1,NLAYERS
       CALL INTERP(T(IL),P(IL),xkappa,kappa)
       DO 8 I=1, 55
       DO 9 J=1, 8
       DO 10 K=1, 3
       KAPPALAYER(I,J,K,IL) = kappa(I,J,K)
  10  CONTINUE
  9   CONTINUE
  8   CONTINUE
  7   CONTINUE
      DO IL = 1,NLAYERS
       CALL INTERPOZONE(P(IL),WEIGHTOZC,KAPPAOZC)
        DO J = 1, 8
        KAPPALAYEROZ(J,IL) = KAPPAOZC(J)
        ENDDO
      ENDDO
c       PRINT*, T,P
       DO 99 I=1, ND
        FUPIR(I) = 0.
        FDNIR(I) = 0.
 99   CONTINUE
       DO I=1,55
       TAUTOTAL(I) = 0.0
       ENDDO
****** Loop over frequency      
       DO 1 I=1, 55
       DO 20 J=1,ND 
       VAC = AV1(I)
       BPLANCK(J) = PLANCK(VAC,T(J),HP,C,HK)
  20   CONTINUE
       DO 15 J=1,ND
         FUPA(J) = 0.0
         FDNA(J) = 0.0
         FUP(J) = 0.0
         FDN(J) = 0.0
         TRANSLAYER(J) = 0.0 
  15   CONTINUE
****** AEROSOLS
       DO IL=1,NLAYERS
        r = RAER(IL)
        np = PARTICLES(IL)
        TAUAEXTIR(IL) = QEXTIR(I,IL)*PI*r*r*DALT(IL)*np*1.E+5
        TAUASIR(IL) = OMG0AIR(I,IL)*TAUAEXTIR(IL)
       ENDDO

       IF (I.EQ.15) THEN
       TAUEXTIRTOTAL = 0.
       TAUASIRTOTAL = 0.
       DO IL=1,NLAYERS
       TAUEXTIRTOTAL = TAUEXTIRTOTAL + TAUAEXTIR(IL)
       TAUASIRTOTAL = TAUASIRTOTAL + TAUASIR(IL)
       ENDDO
c       PRINT *, '******************************'
c       PRINT *, 'TAUEXTIRTOTAL'
c       PRINT *, TAUEXTIRTOTAL
c       PRINT *, 'TAUSIRTOTAL'
c       PRINT *, TAUASIRTOTAL
       TAUAABSIRTOTAL = TAUEXTIRTOTAL - TAUASIRTOTAL
c       PRINT *, 'TAUAABSIRTOTAL'
c       PRINT *, TAUAABSIRTOTAL
c       PRINT *, '*******************************'
       ENDIF 
       
****  8 - 12 UM CONTINUUM
       IF ((I.GT.14).and.(I.LT.21)) THEN
C       PRINT *, 'TAUCONTIN'
       DO IL = 1,NLAYERS
       DP = ABS(P(IL+1)-P(IL))*1.E6
       TV = 1.25E-22 + 1.67E-19 * EXP(-2.62E-13*VAC)
       TF = EXP(1800. * (1./T(IL) - 1./296.))
       CNUT = TV*TF
       TAUCONTIN(IL) = CNUT*FI(1,IL)*FI(1,IL)*P(IL)*DP/(DM*SM*G)
C       PRINT 19, TAUCONTIN(IL) 
       ENDDO
       ELSE
       DO IL = 1,NLAYERS
       TAUCONTIN(IL) = 0.
       ENDDO
       ENDIF
C******
C
       DO 2 K1=1, 8 
       DO 3 K2=1, 8 
       DO 4 K3=1, 6 
        TWGHT = WEIGHT(K1,1)*WEIGHT(K2,2)*WEIGHT(K3,3)
        DO 11 IL = 1,NLAYERS
         PPE = (1. + 0.3*FI(2,IL))*P(IL)     !CO2
         TPE = (300./T(IL))**TPR(I)
         TAUGCO2(IL) = KAPPALAYER(I,K1,1,IL)*CGAS(IL,3)/(44.*SM)
         TAUGH2O(IL) = KAPPALAYER(I,K2,2,IL)*CGAS(IL,2)/(18.*SM)
         TAUGCH4(IL) = KAPPALAYER(I,K3,3,IL)*CGAS(IL,5)
     & *(4.46*DM)/(16.*SM)
****** Pressure induced absorption by CO2 
         CGAS1 = CGAS(IL,3)/1.963E-3  
         TPRIND(IL) = CPR(I)*PPE*TPE*CGAS1
*******
         TAUGIR(IL) = TAUGCH4(IL)+TAUGCO2(IL)+TAUGH2O(IL)+TPRIND(IL)
     & +TAUCONTIN(IL)
  11    CONTINUE
******* Ozone absorption 
        IF (I.EQ.18) THEN
      DO K4=1,8
        TWGHTT(K4) = TWGHT*WEIGHTOZC(K4)
      DO IL=1,NLAYERS 
      TAUGOZ(IL) =KAPPALAYEROZ(K4,IL)*CGAS(IL,4)*(4.46E-5*DM)/(48.*SM)
      TAUGIR(IL) = TAUGCH4(IL)+TAUGCO2(IL)+TAUGH2O(IL)+TPRIND(IL)
     & +TAUCONTIN(IL) + TAUGOZ(IL)
      ENDDO
      DO IL=1, NLAYERS
          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL)
          TAUSIR(IL) = TAUASIR(IL) 
          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL) 
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)
      ENDDO 
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK)
C
               J = 1
  13          IF (J .LE. ND) THEN
                  FUPA(J)=FUPA(J)+TWGHTT(K4)*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHTT(K4)*FDN(J)
                  J=J+1
                  GOTO 13
               END IF
        ENDDO
        GOTO 4
        ENDIF
***************
      DO IL=1, NLAYERS
          TAULAMIR(IL) = TAUAEXTIR(IL) + TAUGIR(IL)
          TAUSIR(IL) = TAUASIR(IL) 
          ASYIR(IL) = ASYAIR(I,IL) 
          OMG0IR(IL) = TAUSIR(IL)/TAULAMIR(IL) 
          OMG0IR(IL) = AMIN1(OMG0IR(IL),0.99999)
          OMG0IR(IL) = AMAX1(OMG0IR(IL),1.E-5)
      ENDDO 
      CALL DELTA2STRIR(SRFALBIR,ASYIR,TAULAMIR,OMG0IR,
     & FUP,FDN,BPLANCK)
C
C
               J = 1
  16          IF (J .LE. ND) THEN
                  FUPA(J)=FUPA(J)+TWGHT*FUP(J)
                  FDNA(J)=FDNA(J)+TWGHT*FDN(J)
                  J=J+1
                  GOTO 16
               END IF
C
  4     CONTINUE             
  3     CONTINUE            
  2     CONTINUE           
         DO 14 J=1,ND
             FDNIR(J)=FDNIR(J)+W(I)*FDNA(J)
             FUPIR(J)=FUPIR(J)+W(I)*FUPA(J)
             BPLANCK(J) = 0.
  14    CONTINUE
C        PRINT*, 'TRANSLAYER'
C        DO J=1,NLAYERS
C        PRINT 19, TRANSLAYER(J)
C         TAUTOTAL(I) = TAUTOTAL(I) - alog(TRANSLAYER(J))
C        ENDDO 
   1     CONTINUE                   !**END LOOP over frequency**
  19   FORMAT(/1X,1PE10.4)
 100  FORMAT(1X,1P10E12.5) 
C      PRINT*,'TAUTOTAL'
C      PRINT 100, TAUTOTAL
C
      RETURN
      END
