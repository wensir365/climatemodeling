module integrate
	use equation
	implicit none

	real :: dt=0.01

contains

	subroutine rk2(x0,y0,z0, x1,y1,z1)
		implicit none
		real, intent(in) :: x0,y0,z0
		real, intent(out) :: x1,y1,z1
		real :: q1x,q1y,q1z, q2x,q2y,q2z

		!--- step 1
		q1x = dt*eqx(x0,y0,z0)
		q1y = dt*eqy(x0,y0,z0)
		q1z = dt*eqz(x0,y0,z0)

		!--- step 2
		q2x = dt*eqx(x0+q1x,y0+q1y,z0+q1z)-q1x
		q2y = dt*eqy(x0+q1x,y0+q1y,z0+q1z)-q1y
		q2z = dt*eqz(x0+q1x,y0+q1y,z0+q1z)-q1z

		!--- step 3
		x1 = x0+q1x+0.5*q2x
		y1 = y0+q1y+0.5*q2y
		z1 = z0+q1z+0.5*q2z
	end subroutine

	subroutine rk4(x0,y0,z0, x1,y1,z1)
		implicit none
		real, intent(in) :: x0,y0,z0
		real, intent(out) :: x1,y1,z1
		real :: q1x,q1y,q1z, q2x,q2y,q2z, q3x,q3y,q3z, q4x,q4y,q4z

		!--- step 1
		q1x = dt*eqx(x0,y0,z0)
		q1y = dt*eqy(x0,y0,z0)
		q1z = dt*eqz(x0,y0,z0)

		!--- step 2
		q2x = dt*eqx(x0+q1x*0.5,y0+q1y*0.5,z0+q1z*0.5)
		q2y = dt*eqy(x0+q1x*0.5,y0+q1y*0.5,z0+q1z*0.5)
		q2z = dt*eqz(x0+q1x*0.5,y0+q1y*0.5,z0+q1z*0.5)

		!--- step 3
		q3x = dt*eqx(x0+q2x*0.5,y0+q2y*0.5,z0+q2z*0.5)
		q3y = dt*eqy(x0+q2x*0.5,y0+q2y*0.5,z0+q2z*0.5)
		q3z = dt*eqz(x0+q2x*0.5,y0+q2y*0.5,z0+q2z*0.5)

		!--- step 4
		q4x = dt*eqx(x0+q3x,y0+q3y,z0+q3z)
		q4y = dt*eqy(x0+q3x,y0+q3y,z0+q3z)
		q4z = dt*eqz(x0+q3x,y0+q3y,z0+q3z)

		!--- step 5
		x1 = x0+(q1x+2.0*q2x+2.0*q3x+q4x)/6.0
		y1 = y0+(q1y+2.0*q2y+2.0*q3y+q4y)/6.0
		z1 = z0+(q1z+2.0*q2z+2.0*q3z+q4z)/6.0
	end subroutine

	subroutine abm(x0,y0,z0, x1,y1,z1, x2,y2,z2)
		implicit none
		real, intent(in) :: x0,y0,z0, x1,y1,z1
		real, intent(out) :: x2,y2,z2
		real :: xest,yest,zest

		!--- step 1
		xest = x1+0.5*dt*(3.0*eqx(x1,y1,z1)-eqx(x0,y0,z0))
		yest = y1+0.5*dt*(3.0*eqy(x1,y1,z1)-eqy(x0,y0,z0))
		zest = z1+0.5*dt*(3.0*eqz(x1,y1,z1)-eqz(x0,y0,z0))

		!--- step 2
		x2 = x1+dt*(5.0*eqx(xest,yest,zest)+8.0*eqx(x1,y1,z1)-eqx(x0,y0,z0))/12.0
		y2 = y1+dt*(5.0*eqy(xest,yest,zest)+8.0*eqy(x1,y1,z1)-eqy(x0,y0,z0))/12.0
		z2 = z1+dt*(5.0*eqz(xest,yest,zest)+8.0*eqz(x1,y1,z1)-eqz(x0,y0,z0))/12.0
	end subroutine

end module integrate
