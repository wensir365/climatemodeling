program main

	implicit none
	integer, parameter :: N=18
	real, parameter :: S0=1368.0,A=204.0,B=2.17,k=3.81,Tc=-10.0,c=1.0,m=1.0
	real, dimension(1:N) :: Sf, lat
	real, parameter :: dt = 0.01

	real, dimension(1:N) :: t0,t1
	integer :: i

	Sf = (/ 0.5,0.531,0.624,0.77,0.892,1.021,1.12,1.189,1.219, &
	        1.219,1.189,1.12,1.021,0.892,0.77,0.624,0.531,0.5 /)
	lat= (/ -85,-75,-65,-55,-45,-35,-25,-15,-5,5,15,25,35,45,55,65,75,85 /)

	t0 = 0.0

	do i = 1, 10000
		call forward(t0,t1)
		t0 = t1
		print "(1i10.0,18f10.2)", i,t0
	end do

	do i = 1, N
		print *, lat(i),t0(i)
	end do

contains

	subroutine forward(t0,t1)
		implicit none
		real, dimension(1:N), intent(in) :: t0
		real, dimension(1:N), intent(out) :: t1
		real :: term1, term2, term3, tmean, albedo, tend
		integer :: i,j

		do i = 1, N
			!--- shortwave ---
			if (t0(i)>Tc) then
				albedo = 0.3
			else
				albedo = 0.62
			end if
			term1 = S0/4.0*(1-albedo)*Sf(i)

			!--- longwave ---
			term2 = A+B*t0(i)

			!--- transport ---
			tmean = 0.0
			do j = 1,N
				tmean = tmean + t0(j)
			end do
			tmean = tmean / N
			term3 = k*(t0(i)-tmean)

			!--- forward ---
			tend = term1-term2-term3
			t1(i) = t0(i) + dt*(tend/c/m)
		end do
	end subroutine

end program main
