module data
        implicit none

        real, parameter :: Pi=3.1415926
	real :: g, Re, Omega

	character(len=100) :: description, relaxation
        logical :: betaplane
	real :: criteria

	real :: lat0
        real, parameter :: deltalat=20.0
        integer, parameter :: Nx=360, Ny=2*deltalat

        real :: Lx, Ly
        real, dimension(Nx) :: lon, x
        real, dimension(Ny) :: lat, y, f
        real :: dx, dy, dt
	integer :: stepoutinterval, totalstep

        real, dimension(Nx,Ny)   :: phi
        integer, dimension(Nx) :: pr, pl

contains

        function rad(degree)
                implicit none
                real, intent(in) :: degree
                real :: rad
                rad = degree/180.0*Pi
        end function

        subroutine init
                implicit none
                integer :: i,j

		!--- read parameters from namelist
		namelist /physics/ g, Re, Omega, lat0
		namelist /bve/ description, relaxation, betaplane, criteria, dt, stepoutinterval, totalstep
		read(*,nml=physics)
		read(*,nml=bve)

		print *, "--- Barotropic Vorticity Equation ---"
		print *, "g = ", g
		print *, "Re = ", Re
		print *, "Omega = ", Omega
		print *, "lat0 = ", lat0
		print *, "deltalat = ", deltalat
		print *, "description = ", trim(description)
		print *, "relaxation = ", trim(relaxation)
		print *, "betaplane = ", betaplane
		print *, "criteria = ", criteria
		print *, "dt = ", dt
		print *, "stepoutinterval =", stepoutinterval
		print *, "totalstep =", totalstep
		print *, "--- B E G I N  H E R E !!! ---"

                Lx = 2.0*Pi*Re*cos(rad(lat0))
                Ly = 2.0*Re*rad(deltalat)

                !--- init coordinates and pointers
                do i = 1, Nx
                        lon(i) = real(i)-0.5
                        x(i) = rad(lon(i))*Re*cos(rad(lat0))

                        pr(i) = i+1
                        pl(i) = i-1
                end do
                pr(Nx) = 1
                pl(1) = Nx

                do j = 1, Ny
                        lat(j) = lat0-deltalat+real(j)-0.5
                        y(j) = Re*rad(lat(j)-lat0)
			if (betaplane) then
				f(j) = 2.0*Omega*sin(rad(lat(j)))
			else
				f(j) = 2.0*Omega*sin(rad(lat0))
			end if
                end do

                dx = Re*cos(rad(lat0))*rad(1.0)
                dy = Re*rad(1.0)
		dx = (dx+dy)/2.0
		dy = (dx+dy)/2.0

		if (dt<=0) then
                	dt = min(dx,dy)/310.0/2.0
		else
			print *, "dt was given by YOU!   dt=", dt
		end if

        end subroutine

end module data
