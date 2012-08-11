module integrate

	use data
	use equation

	implicit none
	real, parameter :: dt=0.3

contains

	subroutine rk4(x0,x1)
		implicit none
		real, dimension(0:N-1,0:N-1), intent(in) :: x0
		real, dimension(0:N-1,0:N-1), intent(out) :: x1
		real, dimension(0:N-1,0:N-1) :: tend,q1,q2,q3,q4

		call advectiontend(x0,tend)
		q1 = dt*tend

		call advectiontend(x0+0.5*q1,tend)
		q2 = dt*tend

		call advectiontend(x0+0.5*q2,tend)
		q3 = dt*tend

		call advectiontend(x0+q3,tend)
		q4 = dt*tend

		x1 = x0+(q1+2.0*q2+2.0*q3+q4)/6.0
	end subroutine

end module integrate
