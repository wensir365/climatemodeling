
        SUBROUTINE GRIDAER 

        INCLUDE '../INCLUDECLIM/parND.inc'
        
        REAL  Z(18), NPART(18),RPART(18)

        INCLUDE '../INCLUDECLIM/comALTBLOK.inc'
        INCLUDE '../INCLUDECLIM/comCPART.inc'

        DATA Z/0.5, 4.5, 8.5, 12.5, 16.5, 20.5, 24.5, 28.5, 32.5,
     2  36.5, 40.5, 44.5, 48.5, 52.5, 56.5, 60.5, 64.5, 68.5/
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.003 FCH4=0.001 
C        DATA NPART/3.4E-1, 0.9, 1.9, 6.7, 25.7, 27.3, 25.3, 26.2,
C     2  35.3, 58.6, 102., 178.8, 313.6, 549.6, 962.3, 1683., 2879.,
C     3  4514./
C  Radius of particles in cm
C        DATA RPART/3.84E-5, 3.84E-5, 3.84E-5, 3.84E-5, 3.84E-5,
C     2  3.71E-5, 3.27E-5, 2.73E-5, 2.11E-5, 1.55E-5, 1.12E-5, 
C     3  8.14E-6, 5.9E-6, 4.3E-6, 3.1E-6, 2.25E-6, 1.56E-6, 9.1E-7/
C
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.001 FCH4=0.001 
        DATA NPART/7.7E-1, 2.05, 4.3, 14.5, 48.6, 47.6, 42.2, 40.3,
     2  47.1, 72.1, 124., 216.6, 379.6, 664.8, 1163., 2026., 3448.,
     3  5638./
C  Radius of particles in cm
        DATA RPART/5.3E-5, 5.3E-5, 5.3E-5, 5.3E-5, 5.3E-5,
     2  5.E-5, 4.5E-5, 3.8E-5, 3.1E-5, 2.3E-5, 1.7E-5, 
     3  1.2E-5, 8.8E-6, 6.4E-6, 4.6E-6, 3.3E-6, 2.3E-6, 1.5E-6/
C
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.002 FCH4=0.001 
C        DATA NPART/5.5E-1, 1.45, 3., 10.4, 37.2, 37.6, 33.9, 33.4,
C     2  41.4, 65.9, 114., 200.1, 350.7, 614.4, 1075., 1880., 3209.,
C     3  5167./
C  Radius of particles in cm
C        DATA RPART/4.6E-5, 4.6E-5, 4.6E-5, 4.6E-5, 4.6E-5,
C     2  4.4E-5, 3.9E-5, 3.3E-5, 2.6E-5, 1.95E-5, 1.4E-5, 
C     3  1.E-5, 7.4E-6, 5.4E-6, 3.9E-6, 2.8E-6, 1.96E-6, 1.2E-6/
C USE THIS PARTICLE DISTRIBUTION FOR FCO2=0.0025 FCH4=0.001 
C        DATA NPART/4.4E-1, 1.18, 2.5, 8.5, 31.5, 32.5, 29.6, 29.9,
C     2  38.4, 62.4, 108.5, 190., 333.1, 583.7, 1022., 1787., 3053.,
C     3  4851./
C  Radius of particles in cm
C        DATA RPART/4.3E-5, 4.3E-5, 4.3E-5, 4.3E-5, 4.3E-5,
C     2  4.1E-5, 3.6E-5, 3.E-5, 2.4E-5, 1.76E-5, 1.3E-5, 
C     3  9.E-6, 6.7E-6, 4.9E-6, 3.5E-6, 2.6E-6, 1.8E-6, 1.1E-6/
        DO I =1,18
         RPART(I) = RPART(I)
        ENDDO 
        DO IL = 1,ND
        I = 1
        IF (ALT(IL).LE.Z(I)) THEN
        PARTICLES(IL) = 0.
        RAER(IL) = RPART(I)
        GOTO 4
        ENDIF
    1   IF (ALT(IL).GT.Z(I)) THEN
        I = I+1
        GOTO 1
        ENDIF
        PARTICLES(IL) = NPART(I-1)+((NPART(I)-NPART(I-1))/
     &  (Z(I)-Z(I-1)))*(ALT(IL)-Z(I-1))
        RAER(IL) = RPART(I-1)+((RPART(I)-RPART(I-1))/
     &  (Z(I)-Z(I-1)))*(ALT(IL)-Z(I-1))
    4   ENDDO
        DO IL = 1,ND
        PARTICLES(IL) = PARTICLES(IL)/100
        ENDDO
        PRINT *, 'NPARTICLES'
        PRINT *, PARTICLES
        PRINT *, 'Radius'
        PRINT *, RAER
        RETURN
        END
