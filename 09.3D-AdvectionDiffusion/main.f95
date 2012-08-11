!-------------------------------------------
! 3D Advection-Diffusion Model
! Demo Script for CM2010
!
!  p(N)      p(uN)   P(vN)   p(wN)
! ------ = - ----- - ----- - -----
!  p(t)      p(x)    p(y)    p(z)
!
!               p(N)2      p(N)2      p(N)2
!          + Kx ----- + Ky ----- + Kz -----
!               P2(x)      p2(y)      p2(z)
!
! by Xinyu Wen, Peking Univ.
! xwen@pku.edu.cn
! Nov 4, 2010
!-------------------------------------------

program main

	use io
	use data
	use integrate
	implicit none

	integer :: t
	logical :: constantwind=.true.

	call init
	call output_netcdf(0)

	do t = 1, totalsteps
		if (constantwind) then
			! do nothing
		else
			! call upgradewind(u,v,w,time)
		end if

		call rk4(u,v,w,n0,n1)
		n0 = n1

		call output_netcdf(t)
	end do

end program main
