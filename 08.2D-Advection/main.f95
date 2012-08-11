!---------------------------------
! 2-D Advection Equation
!
!    p(N)        p(N)      P(N)
!   ------ = -u ------ -v ------
!    p(t)        p(x)      p(y)
!
! demo script for CM2010
!
! By
!
! Xinyu Wen, Peking Univ.
! xwen@pku.edu.cn
! Nov. 3, 2010
!---------------------------------

program main

	use data
	use integrate

	implicit none
	integer :: i
	character(len=4) :: id		! 4-digit expression of i
	character(len=100) :: filename

	call init_uvx
	idealx = xhillcenter
	idealy = yhillcenter
	filename = "result.0000.txt"
	call print_current_uvx(filename)

	do i = 1, 2000
		write(id,"(i4.4)") i	! make id, which is the 4-digit format of i
		filename = "result."//id//".txt"
		call rk4(x0,x1)
		x0 = x1

		idealx = xcenter-25.0*cos(omega*dt*real(i))
		idealy = ycenter-25.0*sin(omega*dt*real(i))
		call print_current_uvx(filename)
	end do

end program main
