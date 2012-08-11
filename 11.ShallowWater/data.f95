module data
        implicit none

        real, parameter :: Pi=3.1415926
	real :: g, Re, Omega
	real :: Hscale

	character(len=100) :: description
        logical :: term_advection, term_coriolis, term_diffusion, beta_effect
	real :: Kx, Ky

	real :: lat0
        real, parameter :: deltalat=30.0
        integer, parameter :: Nx=360, Ny=2*deltalat

        real :: Lx, Ly
        real, dimension(Nx) :: lon, x
        real, dimension(Ny) :: lat, y
        real, dimension(Ny) :: f

        real :: dx, dy, dt

        real, dimension(Nx,Ny)   :: h
        real, dimension(Nx,Ny)   :: u
        real, dimension(Nx,Ny+1) :: v

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
		namelist /physics/ g, Re, Omega, Hscale, lat0
		namelist /swe/ description, term_advection, term_coriolis, term_diffusion, beta_effect, Kx, Ky
		read(*,nml=physics)
		read(*,nml=swe)

		print *, "--- Shallow Water Equations ---"
		print *, "g = ", g
		print *, "Re = ", Re
		print *, "Omega = ", Omega
		print *, "lat0 = ", lat0
		print *, "deltalat = ", deltalat
		print *, "Kx = ", Kx
		print *, "Ky = ", Ky
		print *, "description = ", trim(description)
		print *, "term_advection = ", term_advection
		print *, "term_coriolis = ", term_coriolis
		print *, "term_diffusion = ", term_diffusion
		print *, "beta_effect = ", beta_effect
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
			if (beta_effect) then
                        	f(j) = 2.0*Omega*sin(rad(lat0))+2.0*Omega*cos(rad(lat0))*rad(lat(j)-lat0)*20.0
			else
				f(j) = 2.0*Omega*sin(rad(lat0))
			end if
                end do

                dx = x(2)-x(1)
                dy = y(2)-y(1)
                dt = min(dx,dy)/310.0/2.0

                !--- init u,v,h
                u = 0.0; v=0.0; h=Hscale
        end subroutine

        subroutine forcing(u,v,tt)
                implicit none
                real, dimension(Nx,Ny), intent(inout)   :: u
                real, dimension(Nx,Ny+1), intent(inout) :: v
                integer, intent(in) :: tt
                integer :: xx,yy
                integer, parameter :: periodstep=200
                real, parameter :: amp=1.0

                xx = Nx/2
                yy = Ny/2

		if (tt<=periodstep) then
                	!v(xx,yy) = v(xx,yy)+amp*sin(2*Pi*real(tt)/real(periodstep))
                	u(xx,yy) = u(xx,yy)+amp*sin(2*Pi*real(tt)/real(periodstep))
		end if
        end subroutine

end module data
