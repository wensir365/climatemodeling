module equation

	use data
	implicit none

contains

	subroutine advectiontend(x,tend)
		implicit none
		real, dimension(0:N-1,0:N-1), intent(in) :: x
		real, dimension(0:N-1,0:N-1), intent(out) :: tend

		tend = 0.0
		tend(1:N-2,1:N-2) = -1.0*u(1:N-2,1:N-2)*(x(1:N-2,2:N-1)-x(1:N-2,0:N-3))/(2.0*dx) &
				    -1.0*v(1:N-2,1:N-2)*(x(2:N-1,1:N-2)-x(0:N-3,1:N-2))/(2.0*dy)
	end subroutine

end module equation
