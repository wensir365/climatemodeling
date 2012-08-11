!----------------------------
! 1D Burger's Equation
! Demo Script for CM2010
!
!  p(u)       p(u)
! ------ = -u ----
!  p(t)       p(x)
!
! by Xinyu Wen, Peking Univ.
! xwen@pku.edu.cn
! Nov 17, 2010
!----------------------------

program main

	implicit none
	integer, parameter :: N=100
	real, parameter :: Pi=3.1415926, k=1.0, L=101.0, dx=1.0, dt=0.02
	real, dimension(1:N) :: x0,x1
	integer, dimension(1:N) :: pleft,pright
	real :: idealx
	integer :: i

	!--- init ---
	do i = 1,N
		!x0(i) = 1.0*sin(2*Pi*k*real(i)/L) +	&
		!        0.5*sin(2*Pi*2.0*real(i)/L) +	&
		!        0.2*sin(2*Pi*3.0*real(i)/L) +	&
		!        0.8*sin(2*Pi*4.0*real(i)/L) +	&
		!        0.2*sin(2*Pi*5.0*real(i)/L)
		x0(i) = cos(2*Pi*k*real(i)/L)
	end do

	do i = 1,N-1
		pright(i) = i+1
	end do
	pright(N) = 1

	do i = 2,N
		pleft(i) = i-1
	end do
	pleft(1) = N

	!--- integrate ---
	do i = 1, 10000
		call rk2(x0,x1)
		!call rk4(x0,x1)
		!call forward(x0,x1)
		x0 = x1
		!idealx = mod(initpos+dt*real(i)*c,100.0)
		!print "(1i10,1f10.4,100f10.4)", i, idealx, x0
		print "(1i10,100f10.4)", i, x0
	end do


contains

	subroutine advectiontend(x0,tend)
		implicit none
		real, dimension(1:N), intent(in) :: x0
		real, dimension(1:N), intent(out) :: tend
		integer :: i
		do i = 1,N
			tend(i) = -1.0*x0(i)*(x0(pright(i))-x0(pleft(i)))/(2.0*dx)
		end do
	end subroutine

        subroutine rk2(x0,x1)
                implicit none
                real, dimension(1:N), intent(in) :: x0
                real, dimension(1:N), intent(out) :: x1
                real, dimension(1:N) :: tend,q1,q2

                call advectiontend(x0,tend)
                q1 = dt*tend

                call advectiontend(x0+q1,tend)
                q2 = dt*tend-q1

                x1 = x0+q1+0.5*q2
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

end program main
