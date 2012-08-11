module integrate

	implicit none
	real :: dt=0.01

contains

	subroutine eu_forward(tend,x0,x1)
		implicit none
		real, intent(in) :: tend, x0
		real, intent(out) :: x1
		x1 = tend*dt+x0
	end subroutine

end module integrate
