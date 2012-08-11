module data

	implicit none

	real, parameter :: PI=3.1415926
	real, parameter :: Ki=1.0, Kj=1.0, Kk=1.0
	integer, parameter :: Ni=100, Nj=100, Nk=50

	real, parameter :: length_i=100, length_j=100, length_k=50		! domain size in physical space
	real, parameter :: di=length_i/real(Ni), 				&
			   dj=length_j/real(Nj), 				&
			   dk=length_k/real(Nk)					! grid interval in physical space
	real, parameter :: dt=0.05						! should < CFL=0.1 (sec)
	integer, parameter :: totalsteps=4000

	logical, parameter :: switch_adve = .true.
	logical, parameter :: switch_diff = .false.

	!--- the order of dimension (i,j,k) is for generating netCDF file ---
	real, dimension(1:Ni,1:Nj,1:Nk) :: posi, posj, posk				! physical space
	real, dimension(1:Ni,1:Nj,1:Nk) :: u, v, w, n0, n1				! mathmatical space
	integer, dimension(1:Ni,1:Nj,1:Nk) :: pright,pleft,pnorth,psouth,pup,pdown	! pointers

contains
        function gaussian(x,sigma,mean)
                implicit none
                real, intent(in) :: x,sigma,mean
                real :: gaussian
                real :: power
                gaussian = 1.0/sqrt(2.0*pi*sigma*sigma)
                power = ((x-mean)**2)/-2.0/(sigma**2)
                gaussian = gaussian*exp(power)
        end function

	subroutine init
		implicit none
		integer :: i,j,k
		real :: r

		!--- setting u/v/w ---
		u = 10.0
		v = 0.0
		w = 0.0

		!--- pointers ---
		! i : cyclic		all available
		do k = 1, Nk
			do j = 1, Nj
				do i = 1, Ni
					!--- setting coordinates ---
					posi(i,j,k) = di*real(i)-di/2.0
					posj(i,j,k) = dj*real(j)-dj/2.0
					posk(i,j,k) = dk*real(k)-dk/2.0

					!--- setting pointers ---
					pright(i,j,k)	= i+1		! pright(Ni,:,:) is incorrect
					pleft(i,j,k)	= i-1		! pleft(1,:,:) is incorrect
					pnorth(i,j,k)	= j+1		! pnorth(:,Nj,:) is incorrect
					psouth(i,j,k)	= j-1		! psouth(:,1,:) is incorrect
					pup(i,j,k)	= k+1		! pup(:,:,Nk) is incorrect
					pdown(i,j,k)	= k-1		! pdown(:,:,1) is incorrect

					!--- setting IC ---
					r = sqrt((posi(i,j,k)-50.0)**2.0+(posj(i,j,k)-50.0)**2.0+(posk(i,j,k)-25.0)**2.0)
					n0(i,j,k) = 100*gaussian(r,10.0,0.0)
				end do
			end do
		end do
		pright(Ni,:,:)	= 1			! correct pright(Ni,:,:) ---> cyclic BC
		pleft(1,:,:)		= Ni		! correct pleft(1,:,:) -----> cyclic BC

		print *, "------------------------------"
		print *, " 3D Advection-Diffusion Model "
		print *, "    Climate Modeling 2010     "
		print *, "                              "
		print *, "    steps = ", totalsteps
		print *, "    dt    = ", dt
		print *, "    di    = ", di
		print *, "    dj    = ", dj
		print *, "    dk    = ", dk
		print *, "    Ni    = ", Ni
		print *, "    Nj    = ", Nj
		print *, "    Nk    = ", Nk
		print *, "                              "
		print *, "    advection term =", switch_adve
		print *, "    diffusion term =", switch_diff
		print *, "                              "
		print *, " For more info, please email  "
		print *, " Xinyu Wen <xwen@pku.edu.cn>  "
		print *, "------------------------------"
	end subroutine

end module data
