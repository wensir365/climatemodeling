module integrate

	use data
	use equation

	implicit none

contains

	subroutine forward(u,v,w,n0,n1)
		implicit none
		real, dimension(1:Ni,1:Nj,1:Nk), intent(in) :: u,v,w,n0
		real, dimension(1:Ni,1:Nj,1:Nk), intent(out) :: n1
		real, dimension(1:Ni,1:Nj,1:Nk) :: tend
		call advectiondiffusion(u,v,w,n0,tend)
		n1 = n0+dt*tend
	end subroutine

	subroutine rk2(u,v,w,n0,n1)
		implicit none
		real, dimension(1:Ni,1:Nj,1:Nk), intent(in) :: u,v,w,n0
		real, dimension(1:Ni,1:Nj,1:Nk), intent(out) :: n1
		real, dimension(1:Ni,1:Nj,1:Nk) :: tend,q1,q2

		call advectiondiffusion(u,v,w,n0,tend)
		q1 = dt*tend

		call advectiondiffusion(u,v,w,n0+q1,tend)
		q2 = dt*tend-q1

		n1 = n0+q1+0.5*q2
        end subroutine

	subroutine rk4(u,v,w,n0,n1)
		implicit none
		real, dimension(1:Ni,1:Nj,1:Nk), intent(in) :: u,v,w,n0
		real, dimension(1:Ni,1:Nj,1:Nk), intent(out) :: n1
		real, dimension(1:Ni,1:Nj,1:Nk) :: tend,q1,q2,q3,q4

		call advectiondiffusion(u,v,w,n0,tend)
		q1 = dt*tend

		call advectiondiffusion(u,v,w,n0+0.5*q1,tend)
		q2 = dt*tend

		call advectiondiffusion(u,v,w,n0+0.5*q2,tend)
		q3 = dt*tend

		call advectiondiffusion(u,v,w,n0+q3,tend)
		q4 = dt*tend

		n1 = n0+(q1+2.0*q2+2.0*q3+q4)/6.0
	end subroutine

end module integrate
