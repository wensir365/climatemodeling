module equation

	implicit none

contains

	function eqx(x,y,z)
		implicit none
		real, parameter :: sigma = 10.0
		real, intent(in) :: x,y,z
		real :: eqx

		eqx = sigma*(y-x)
	end function

	function eqy(x,y,z)
		implicit none
		real, parameter :: rho = 28.0
		real, intent(in) :: x,y,z
		real :: eqy

		eqy = rho*x-y-x*z
	end function eqy

	function eqz(x,y,z)
		implicit none
		real, parameter :: beta = 8.0/3.0
		real, intent(in) :: x,y,z
		real :: eqz

		eqz = x*y-beta*z
	end function

end module equation
