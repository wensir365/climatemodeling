      SUBROUTINE CHOOSE_STAR(FLUXC,FLUXUVC)

c This subroutine reads the spectra from stars other than the SUN

      INCLUDE '../INCLUDECLIM/parNSOL_NSOLUV.inc'
      DIMENSION FLUXC(NSOL),FLUXUVC(NSOLUV)
      CHARACTER :: STARR*3,dirDATA*7
      INCLUDE '../INCLUDECLIM/comSTR.inc'
      
      dirDATA = '../DATA'       
      
      if(STARR.eq.'Sun') goto 1
      IF(STARR=="dMV")OPEN(UNIT=46,FILE= dirDATA//'/ADLeo_surf.pdat')
      IF(STARR=="F2V")OPEN(UNIT=46,FILE= dirDATA//'/F2V_surf.pdat')
      IF(STARR=="K2V")OPEN(UNIT=46,FILE= dirDATA//'/K2V_surf.pdat')

       READ(46,*)
       DO j=1,NSOLUV
         READ(46,*)xl,xf,x,x
         FLUXUVC(j) = xf
       ENDDO
       READ(46,*)
       DO j=1,NSOL
         READ(46,*)xl,xf,x,x
         FLUXC(j) = xf
       ENDDO
       close(46)
C

 1    continue
      RETURN
      END

