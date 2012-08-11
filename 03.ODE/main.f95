program main
	use integrate
	implicit none

	real :: x0,x1, y0,y1, z0,z1
	integer :: i

	x0=1.0; y0=1.0; z0=1.0

	do i = 1, 10000
		call rk2(x0,y0,z0,x1,y1,z1)
		print *, i,x1,y1,z1
	
		x0=x1; y0=y1; z0=z1
	end do

end program main
