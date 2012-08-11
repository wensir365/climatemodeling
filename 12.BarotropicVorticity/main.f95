!----------------------------------------
!       Barotropic Vorticity Equation
!
!       Xinyu Wen
!       E-mail: xwen@pku.edu.cn
!
!       this program was developed
!       for the class of CM2010
!       Nov. 2010
!----------------------------------------

program main

        use data
        use integrate
	use io
        implicit none

        integer :: i

        call init
	call input
	call output_netcdf(0,0)
	call output_netcdf(0,1)

        do i = 1, totalstep
		call rk2(phi)
                !call rk4(phi)
                print *, i, "      max/min of h: ", maxval(phi), minval(phi)
		call output_netcdf(i,1)
        end do

	call output_netcdf(0,2)

        print *, "time step (sec) = ", dt
        print *, "number of total steps = ", totalstep
	print *, "dx=", dx
	print *, "dy=", dy

end program main
