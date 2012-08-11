C------------------------------------------------------------------
C   Code added 6/15/01, Kara Krelove, for integration of Mlawer's RRTM
C------------------------------------------------------------------

        SUBROUTINE TRANSLATEM(G,FI,T,PF,ND1,DM,BKM)

C    This subroutine is designed to correctly translate between SurfTem
C    and RRTM, including flipping of layer indices, and READing 
C    molecular abundances into the appropriate grid. 

      INCLUDE '../INCLUDECLIM/parND.inc'
      INCLUDE '../INCLUDECLIM/parNLAYERS.inc'
 
      PARAMETER (NS=7)

	REAL muH, muatm, Ncol, newalt

      DIMENSION T(ND),PF(ND),FI(4,ND), ALTF(ND), DALTF(ND-1), TF(ND)
      INCLUDE '../INCLUDECLIM/comCBLOK.inc'
      INCLUDE '../INCLUDECLIM/comALTBLOK.inc'

c      COMMON/CBLOK/FO2,FN2,FCO2,FAR,FCH4
c      COMMON/ALTBLOK/DALT(ND-1),ALT(ND)

      COMMON/ MLAWERI/  layers, numspec, newalt(ND), TempT(0:NLAYERS), 
     & 			Pres(0:NLAYERS), gasses(7, 0:NLAYERS), COLDEP(ND)


C        Assign mixing ratios to proper placement in array

        do i = 1, 7
	  select case (i)          
	   case (1) 		! water
             do n = 0, NLAYERS
        	gasses(1, n) = FI(1, ND-n)
             end do 
           case (2) 		! CO2
             do n = 0, NLAYERS
        	gasses(2, n) = FI(2, ND-n)
             end do
           case (3)		! ozone
             do n = 0, NLAYERS
        	gasses(3, n) = FI(4, ND-n)
             end do
           case (4) 
              !OR READ from a different file somewhere; or a data block. N2O.
              do n = 0, NLAYERS
        	gasses(4, n) = 1E-20
              end do
           case (5)
              !OR ditto from above. CO. 
              do n = 0, NLAYERS
        	gasses(5, n) = 1E-20
              end do
           case (6)		! methane
             do n = 0, NLAYERS
        	gasses(6, n) = FI(3, ND-n)
             end do
           case (7)
             do n = 0, NLAYERS
        	gasses(7, n) = FO2
             end do
        end select
        end do

C       Translate T grid:
C		T at Flux grid points
	do j = 1, ND-1
	  TF(j) = (T(j)+T(j+1))/2.
	end do
	TF(ND) = T(ND)

        do i = 0, NLAYERS
	  TempT(i) = TF(ND-i)
	end do


C        Translate P grid:
        do k = 0, NLAYERS
           Pres(k) = PF(ND-k)*1013.
	end do

C        Calculate ALT at flux points (currently at T points)

      ALTF(ND) = 0.
      DO JS=1,ND1
      J = ND - JS
      TAvg = 0.5*(TF(J) + TF(J+1))
      water = 0.5 * (FI(1,J) + FI(1,J+1))
      AM = 18.*water + DM
      BMG = BKM/AM
      ALTF(J) = ALTF(J+1) + BMG*TAvg*ALOG(PF(J+1)/PF(J))*1.E-5
      DALTF(J) = ALTF(J) - ALTF(J+1)
      ENDDO

C       Calculate total atmospheric column depth
C	and altitude at flux points (midway between T points).	
	do ij = 0, NLAYERS
           newalt(ij+1) = ALTF(ND-ij)
	end do

	muH = 1.67E-24
c molecular weigth in gr/molecule
	muatm = DM*muH
c   convert from gr/molecule to kg/molecule
	muatm = muatm*1E-3
        
        do j = 1, ND-1
c Nitrogen concentration
          FBROAD =1.-( FI(1,j)+FI(2,j)+FI(3,j)+FI(4,j)+FO2)
	  Ncol = FBROAD*((Pres(j-1)-Pres(j))*1.E2)/((G/100.)*muatm)
c   convert from m^-2 to cm^-2
	  COLDEP(j) = Ncol*1.E-4
c        print*,j,coldep(j)
        end do
          FBROAD = 1. - (FI(1,ND)+FI(2,ND)+FI(3,ND)+FI(4,Nd)+FO2)
          COLDEP(ND)= FBROAD*(Pres(52)*1.E2)*1.e-4/((G/100.)*muatm)

        layers = NLAYERS
	numspec = 7
	
	return
 
        END
