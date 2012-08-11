module equation
	implicit none
	real, parameter :: a=10.0, b=28.0, c=8.0/3.0

contains

	function eqx(x,y,z)
		implicit none
		real, intent(in) :: x,y,z
		real :: eqx
		eqx = -1.0*a*(x-y)
	end function

	function eqy(x,y,z)
		implicit none
		real, intent(in) :: x,y,z
		real :: eqy
		eqy = b*x-y-x*z
	end function

	function eqz(x,y,z)
		implicit none
		real, intent(in) :: x,y,z
		real :: eqz
		eqz = x*y-c*z
	end function

end module

