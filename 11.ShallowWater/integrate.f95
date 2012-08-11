module integrate
        use data
        implicit none

contains

        subroutine swe_u(u,v,h,tend)
                implicit none
                real, dimension(Nx,Ny), intent(in)   :: h
                real, dimension(Nx,Ny), intent(in)   :: u
                real, dimension(Nx,Ny+1), intent(in) :: v
                real, dimension(Nx,Ny), intent(out)  :: tend

                real :: term_a, term_pgf, term_c, term_d
                integer :: i,j
                real :: localv

                tend = 0.0
                do i = 1, Nx
                        do j = 1, Ny
                                localv = (v(pl(i),j+1)+v(pl(i),j)+v(i,j+1)+v(i,j))/4.0
                                term_c = f(j)*localv
                                term_pgf = g*(h(i,j)-h(pl(i),j))/dx
                                if (j==1) then
                                        term_a = u(i,j)*(u(pr(i),j)-u(pl(i),j))/(2.0*dx)
                                        term_a = term_a+localv*(u(i,j+1)-u(i,j))/dy	! reflective bc
                                end if
                                if (j==Ny) then
                                        term_a = u(i,j)*(u(pr(i),j)-u(pl(i),j))/(2.0*dx)
                                        term_a = term_a+localv*(u(i,j)-u(i,j-1))/dy	! reflective bc
                                end if
                                if (1<j .and. j<Ny) then
                                        term_a = u(i,j)*(u(pr(i),j)-u(pl(i),j))/(2.0*dx)
                                        term_a = term_a+localv*(u(i,j+1)-u(i,j-1))/(2.0*dy)
                                end if
                                term_d = Kx*(u(pl(i),j)-2.0*u(i,j)+u(pr(i),j))/(dx*dx)

                                tend(i,j) = -1.0*term_pgf
                                if (term_advection) then
                                        tend(i,j) = tend(i,j)-term_a
                                end if
                                if (term_coriolis) then
                                        tend(i,j) = tend(i,j)+term_c
                                end if
                                if (term_diffusion) then
					if ((j==Ny) .or. (j==1)) then
						tend(i,j) = tend(i,j)+term_d*10000
					else
                                        	tend(i,j) = tend(i,j)+term_d
					end if
                                end if
                        end do
                end do
        end subroutine

        subroutine swe_v(u,v,h,tend)
                implicit none
                real, dimension(Nx,Ny), intent(in) :: h
                real, dimension(Nx,Ny), intent(in) :: u
                real, dimension(Nx,Ny+1), intent(in) :: v
                real, dimension(Nx,Ny+1), intent(out) :: tend

                real :: term_a, term_pgf, term_c, term_d
                integer :: i,j
                real :: localu

                tend = 0.0
                do i = 1, Nx
                        do j = 2, Ny
                                localu = (u(i,j-1)+u(i,j)+u(pr(i),j-1)+u(pr(i),j))/4.0
                                term_c = f(j)*localu
                                term_pgf = g*(h(i,j)-h(i,j-1))/dy
                                term_a = localu*(v(pr(i),j)-v(pl(i),j))/(2.0*dx)+       &
                                         v(i,j)*(v(i,j+1)-v(i,j-1))/(2.0*dy)
                                term_d = Ky*(v(i,j-1)-2.0*v(i,j)+v(i,j+1))/(dy*dy)

                                tend(i,j) = -1.0*term_pgf
                                if (term_advection) then
                                        tend(i,j) = tend(i,j)-term_a
                                end if
                                if (term_coriolis) then
                                        tend(i,j) = tend(i,j)-term_c
                                end if
                                if (term_diffusion) then
					if ((j==2) .or. (j==Ny)) then
						tend(i,j) = tend(i,j)+term_d*10000
					else
                                        	tend(i,j) = tend(i,j)+term_d
					end if
                                end if
                        end do
                end do
        end subroutine

        subroutine swe_h(u,v,h,tend)
                implicit none
                real, dimension(Nx,Ny), intent(in) :: h
                real, dimension(Nx,Ny), intent(in) :: u
                real, dimension(Nx,Ny+1), intent(in) :: v
                real, dimension(Nx,Ny), intent(out) :: tend

                real :: term_a, term_div, term_d
                integer :: i,j
                real :: localu,localv

                tend = 0.0
                do i = 1, Nx
                        do j = 1, Ny
                                localu = (u(i,j)+u(pr(i),j))/2.0
                                localv = (v(i,j+1)+v(i,j))/2.0
                                term_div = h(i,j)*((u(pr(i),j)-u(i,j))/dx+(v(i,j+1)-v(i,j))/dy)
                                if (j==1) then
                                        term_a = localu*(h(pr(i),j)-h(pl(i),j))/(2.0*dx)+        &
                                                 localv*(h(i,j+1)-h(i,j))/dy

					term_d = Kx*(h(pr(i),j)-2.0*h(i,j)+h(pl(i),j))/(dx*dx)+  &
					         Ky*(h(i,j+1)-0)/(dy*dy)
                                end if
                                if (j==Ny) then
                                        term_a = localu*(h(pr(i),j)-h(pl(i),j))/(2.0*dx)+        &
                                                 localv*(h(i,j)-h(i,j-1))/dy

					term_d = Kx*(h(pr(i),j)-2.0*h(i,j)+h(pl(i),j))/(dx*dx)+  &
					         Ky*(0-h(i,j-1))/(dy*dy)
                                end if
                                if (1<j .and. j<Ny) then
                                        term_a = localu*(h(pr(i),j)-h(pl(i),j))/(2.0*dx)+        &
                                                 localv*(h(i,j+1)-h(i,j-1))/(2.0*dy)

					term_d = Kx*(h(pr(i),j)-2.0*h(i,j)+h(pl(i),j))/(dx*dx)+  &
					         Ky*(h(i,j+1)-2.0*h(i,j)+h(i,j-1))/(dy*dy)
                                end if

                                tend(i,j) = -1.0*term_div
                                if (term_advection) then
                                        tend(i,j) = tend(i,j)-term_a
                                end if
                                !if (term_diffusion) then
                                !        tend(i,j) = tend(i,j)+term_d
                                !end if
                        end do
                end do
        end subroutine

        subroutine rk4(u,v,h)
                implicit none
                real, dimension(Nx,Ny),   intent(inout) :: h
                real, dimension(Nx,Ny),   intent(inout) :: u
                real, dimension(Nx,Ny+1), intent(inout) :: v

                real, dimension(Nx,Ny)   :: tendh,q1h,q2h,q3h,q4h
                real, dimension(Nx,Ny)   :: tendu,q1u,q2u,q3u,q4u
                real, dimension(Nx,Ny+1) :: tendv,q1v,q2v,q3v,q4v

                !--- step 1
                call swe_u(u,v,h,tendu)
                call swe_v(u,v,h,tendv)
                call swe_h(u,v,h,tendh)
                q1u = dt*tendu
                q1v = dt*tendv
                q1h = dt*tendh

                !--- step 2
                call swe_u(u+0.5*q1u,v+0.5*q1v,h+0.5*q1h,tendu)
                call swe_v(u+0.5*q1u,v+0.5*q1v,h+0.5*q1h,tendv)
                call swe_h(u+0.5*q1u,v+0.5*q1v,h+0.5*q1h,tendh)
                q2u = dt*tendu
                q2v = dt*tendv
                q2h = dt*tendh

                !--- step 3
                call swe_u(u+0.5*q2u,v+0.5*q2v,h+0.5*q2h,tendu)
                call swe_v(u+0.5*q2u,v+0.5*q2v,h+0.5*q2h,tendv)
                call swe_h(u+0.5*q2u,v+0.5*q2v,h+0.5*q2h,tendh)
                q3u = dt*tendu
                q3v = dt*tendv
                q3h = dt*tendh

                !--- step 4
                call swe_u(u+q3u,v+q3v,h+q3h,tendu)
                call swe_v(u+q3u,v+q3v,h+q3h,tendv)
                call swe_h(u+q3u,v+q3v,h+q3h,tendh)
                q4u = dt*tendu
                q4v = dt*tendv
                q4h = dt*tendh

                !--- update
                u = u+(q1u+2.0*q2u+2.0*q3u+q4u)/6.0
                v = v+(q1v+2.0*q2v+2.0*q3v+q4v)/6.0
                h = h+(q1h+2.0*q2h+2.0*q3h+q4h)/6.0
        end subroutine

end module integrate
