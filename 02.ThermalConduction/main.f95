program main

	implicit none
	real :: tf,t1,t2,t3, t1new,t2new,t3new
	real, parameter :: k1=1, k2=1, k3=1, c=1000, m=1, dt=1
	integer :: i

	tf=100; t1=20; t2=15; t3=10

	do i = 1, 10000
		tf	= 100 + 20*sin(real(i)/300)

		t1new = forward(t1,tend1(tf,t1,t2,t3))
		t2new = forward(t2,tend2(tf,t1,t2,t3))
		t3new = forward(t3,tend3(tf,t1,t2,t3))

		print *, i, tf, t1, t2, t3
		t1=t1new; t2=t2new; t3=t3new
	end do

contains

	function tend1(tf,t1,t2,t3)
		implicit none
		real, intent(in) :: tf,t1,t2,t3
		real :: tend1
		tend1 = k1*(tf-t1)-k2*(t1-t2)
	end function

	function tend2(tf,t1,t2,t3)
		implicit none
		real, intent(in) :: tf,t1,t2,t3
		real :: tend2
		tend2 = k2*(t1-t2)-k3*(t2-t3)
	end function

	function tend3(tf,t1,t2,t3)
		implicit none
		real, intent(in) :: tf,t1,t2,t3
		real :: tend3
		tend3 = k3*(t2-t3)
	end function

	function forward(t,tend)
		implicit none
		real, intent(in) :: t, tend
		real :: forward
		forward = tend/c/m*dt+t
	end function

end program main
