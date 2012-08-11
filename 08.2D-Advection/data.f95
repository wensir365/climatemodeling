module data

	implicit none

	integer, parameter :: N=101
	real, parameter :: xlength=100.0, ylength=100.0
	real, parameter :: PI=3.1415926
	real, parameter :: omega=0.1		! rad/s

	real, parameter :: dx=xlength/real(N-1), dy=ylength/real(N-1)
	real, parameter :: xcenter=50.0, ycenter=50.0
	real, dimension(0:N-1,0:N-1) :: u,v, x0,x1

	real :: idealx, idealy

	real, parameter :: xhillcenter=25.0,yhillcenter=50.0, Reffect=5.0, Hscale=50.0

contains

	function gaussian(x,sigma,mean)
		implicit none
		real, intent(in) :: x, sigma, mean
		real :: gaussian
		real :: power
		power = ((x-mean)**2.0)/-2.0/(sigma**2.0)
		gaussian = 1.0/sqrt(2.0*PI*sigma*sigma)
		gaussian = gaussian*exp(power)
	end function

	subroutine init_uvx
		implicit none
		integer :: i,j
		real :: r,windspeed,theta
		real :: xx, yy

		do j = 0,N-1
			do i = 0,N-1
				!--- Setting U/V IC ---
				xx = dx*real(i)
				yy = dy*real(j)
				r = sqrt((xx-xcenter)**2.0+(yy-ycenter)**2.0)

				if (r>0.01) then
					windspeed = omega*r
					if (xx>xcenter) then	! right region
						theta = atan((yy-ycenter)/(xx-xcenter))
						u(j,i) = windspeed*sin(theta)*-1.0
						v(j,i) = windspeed*cos(theta)
					end if
					if (xx<xcenter) then	! left region
						theta = atan((yy-ycenter)/(xcenter-xx))
						u(j,i) = windspeed*sin(theta)*-1.0
						v(j,i) = windspeed*cos(theta)*-1.0
					end if
					if (abs(xx-xcenter)<0.01) then	! central line
						v(j,i) = 0.0
						if (ycenter>yy) then
							u(j,i) = windspeed
						else
							u(j,i) = -1.0*windspeed
						end if
					end if
				else
					u(j,i) = 0
					v(j,i) = 0
				end if

				!--- Setting N IC ---
				r =sqrt((xx-xhillcenter)**2.0+(yy-yhillcenter)**2.0)
				x0(j,i) = Hscale*gaussian(r,Reffect,0.0)
			end do
		end do
	end subroutine

	subroutine print_current_uvx(fn)
		implicit none
		character(len=100), intent(in) :: fn
	
		print *, "writting ... "//trim(fn)
		open(unit=99,file=trim(fn),action="write")
			write(99,*) idealx, idealy
			write(99,"(101f10.4)") u
			write(99,"(101f10.4)") v
			write(99,"(101f10.4)") x0
		close(unit=99)
	end subroutine

end module data
