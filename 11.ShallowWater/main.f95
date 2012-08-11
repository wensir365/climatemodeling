!----------------------------------------
!       Shallow Water Equations
!
!       Xinyu Wen
!	E-mail: xwen@pku.edu.cn
!
!	this program was developed
!	for the class of CM2010
!	Nov. 2010
!----------------------------------------

program main

        use data
        use integrate
	use io
        implicit none

        integer :: i, totalstep=2000

        call init
	call output_netcdf(0,0)

        do i = 1, totalstep
                call forcing(u,v,i)
                call rk4(u,v,h)
                print *, i,"... ok!"
                print *, "      max/min of u: ", maxval(u), minval(u)
                print *, "      max/min of v: ", maxval(v), minval(v)
                print *, "      max/min of h: ", maxval(h), minval(h)
		call output_netcdf(i,1)
        end do
	call output_netcdf(0,2)

        print *, "time step (sec) = ", dt
        print *, "number of total steps = ", totalstep

end program main
