 
      SUBROUTINE READSTAR(FLUXFAC)

      CHARACTER :: STARR*3, dirDATA*4

      INCLUDE '../INCLUDECHEM/comSTR.inc'
      INCLUDE '../INCLUDECHEM/comFLUXPHOTO.inc'
      dirDATA = 'DATA'

      select case(STARR)
       case('Sun')
        OPEN(unit=76, file= dirDATA//'/faruv_sun.pdat')
      case('dMV')
        OPEN(unit=76, file= dirDATA//'/faruv_ADLeo.pdat')
        OPEN(unit=75,file=dirDATA//'/ADLeo_photo.pdat')
      case('F2V')
        OPEN(unit=76, file= dirDATA//'/faruv_F2V.pdat')
        OPEN(unit=75,file=dirDATA//'/F2V_photo.pdat')
      case('K2V')
        OPEN(unit=76, file= dirDATA//'/faruv_K2V.pdat')
        OPEN(unit=75,file=dirDATA//'/K2V_photo.pdat')
      end select

c Reading far UV data
        READ(76,100)
        do L=1,10
          read(76,*)X,SFX(L)
          SFX(L)=SFX(L)*FLUXFAC
        end do
      close(76) 
 100  FORMAT(/)        
  
****** Read UV flux data from stars other than the Sun
* The Sun UV flux is read in the subroutine READPHOTO
       if(STARR.eq.'Sun') RETURN
        read(75,*)
        do l=1,108
          read(75,*)ll,FLUX(l),x,x
          FLUX(ll) = FLUX(ll)*FLUXFAC
        enddo
      close(75)
      
      RETURN
      END

