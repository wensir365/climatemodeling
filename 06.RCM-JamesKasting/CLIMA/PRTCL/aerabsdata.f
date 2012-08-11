
        SUBROUTINE AERABSDATA

        DIMENSION wavirst(55), wavsolst(38)
        CHARACTER :: DIRDATA*14
        
        INCLUDE '../INCLUDECLIM/comHYDROCARB.inc'
       
C  radius of particles is given in cm
        DATA radstand/0.001E-5, 0.002E-5, 0.003E-5, 0.004E-5,
     &  0.005E-5, 0.006E-5, 0.007E-5, 0.008E-5, 0.009E-5, 
     &  0.001E-4, 0.002E-4, 0.003E-4, 0.004E-4,
     &  0.005E-4, 0.006E-4, 0.007E-4, 0.008E-4, 0.009E-4, 
     &  0.01E-4, 0.02E-4, 0.03E-4, 0.04E-4, 0.05E-4,
     &  0.06E-4, 0.07E-4,
     &  0.08E-4, 0.09E-4, 0.1E-4, 0.2E-4, 0.3E-4, 0.4E-4, 0.5E-4,
     &  0.6E-4, 0.7E-4, 0.8E-4, 0.9E-4,
     &  1.E-4, 2.E-4, 3.E-4, 4.E-4, 5.E-4, 6.E-4, 7.E-4, 8.E-4,
     &  9.E-4, 1.E-3/
         DIRDATA ='../CLIMA/DATA/'
         OPEN(UNIT=37,FILE= DIRDATA//'/irtotal.DAT')
         OPEN(UNIT=38,FILE= DIRDATA//'/soltotal.DAT')
C READING IR particle absorption data
         DO J=1,46
          READ(37,100)
         DO i=1,55
         READ(37,*) wavirst(i),Qextirst(J,i),w0irst(J,i),
     &   girst(J,i)
         ENDDO
         ENDDO
C READING SOLAR particle absorption data
         DO J=1,46
          READ(38,100)
         DO i=1,38
         READ(38,*) wavsolst(i),Qextsolst(J,i),w0solst(J,i),
     &   gsolst(J,i)
         ENDDO
         ENDDO
  100    format(/)
         CLOSE (37)
         CLOSE (38)
         RETURN
         END
