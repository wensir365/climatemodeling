!---------------------------
! 1-D Advection Equation
!
!    p(N)        p(N)
!   ------ = -c ------
!    p(t)        p(x)
!
! demo script for CM2010
!
! By 
!
! Xinyu Wen, Peking Univ.
! xwen@pku.edu.cn
! Oct. 29, 2010
!---------------------------

program main

	implicit none
	integer, parameter :: N=100
	real, parameter :: Pi=3.1415926, sgm=10.0, A=100.0, initpos=50.0, c=10.0, dx=1.0, dt=0.02
	real, dimension(1:N) :: x0,x1
	integer, dimension(1:N) :: pleft,pright
	real :: idealx
	integer :: i

	!--- init ---
	do i = 1,N						! given N distribution
		x0(i) = A*gaussian(real(i)-0.5,sgm,initpos)
	end do

	do i = 1,N-1						! set the right pointer
		pright(i) = i+1
	end do
	pright(N) = 1

	do i = 2,N						! set the left pointer
		pleft(i) = i-1
	end do
	pleft(1) = N

	!--- integrate ---
	do i = 1, 1000
		call rk4(x0,x1)
		!call forward(x0,x1)
		x0 = x1
		idealx = mod(initpos+dt*real(i)*c,100.0)
		print "(1i10,1f10.4,100f10.4)", i, idealx, x0
	end do

contains

        function gaussian(x,sigma,mean)
                implicit none
                real, intent(in) :: x,sigma,mean
                real :: gaussian
                real :: power
                gaussian = 1.0/sqrt(2.0*Pi*sigma*sigma)
                power = ((x-mean)**2)/-2.0/(sigma**2)
                gaussian = gaussian*exp(power)
        end function

	subroutine advectiontend(x0,tend)
		implicit none
		real, dimension(1:N), intent(in) :: x0
		real, dimension(1:N), intent(out) :: tend
		integer :: i
		do i = 1,N
			tend(i) = -1.0*c*(x0(pright(i))-x0(pleft(i)))/(2.0*dx)
		end do
	end subroutine

	subroutine rk4(x0,x1)
		implicit none
		real, dimension(1:N), intent(in) :: x0
		real, dimension(1:N), intent(out) :: x1
		real, dimension(1:N) :: tend,q1,q2,q3,q4
                !--- step1
                call advectiontend(x0,tend)
                q1 = dt*tend
                !--- step2
                call advectiontend(x0+0.5*q1,tend)
                q2 = dt*tend
                !--- step3
                call advectiontend(x0+0.5*q2,tend)
                q3 = dt*tend
                !--- step4
                call advectiontend(x0+q3,tend)
                q4 = dt*tend
                !--- step5
                x1 = x0+(q1+2.0*q2+2.0*q3+q4)/6.0
	end subroutine

	subroutine forward(x0,x1)
		implicit none
		real, dimension(1:N), intent(in) :: x0
		real, dimension(1:N), intent(out) :: x1
		real, dimension(1:N) :: tend
		call advectiontend(x0,tend)
		x1 = x0+dt*tend
	end subroutine

end program main
