program lorenz

	use equation
	use integrate
	implicit none

	integer :: i		! for loop
	real :: x0, x1, y0, y1, z0, z1

	!--- initialize ---
	x0 = 30.0; y0 = 30.0; z0 = 30.0

	!--- main loop ---
	do i = 1, 9999
		call eu_forward(eqx(x0,y0,z0),x0,x1)	! integrate eq X
		call eu_forward(eqy(x0,y0,z0),y0,y1)	! integrate eq Y
		call eu_forward(eqz(x0,y0,z0),z0,z1)	! integrate eq Z

		print *, x1,y1,z1			! print latest xyz

		x0=x1; y0=y1; z0=z1			! replace old
	end do

end program lorenz
