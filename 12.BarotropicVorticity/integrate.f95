module integrate
        use data
        implicit none

contains

	function xcendiffpoint(m,i,j)
		implicit none
		real, dimension(Nx,Ny), intent(in) :: m
		integer, intent(in) :: i,j
		real :: xcendiffpoint
		xcendiffpoint = (m(pr(i),j)-m(pl(i),j))/2.0
	end function

	subroutine xcendiff(m,mx)
		implicit none
		real, dimension(Nx,Ny), intent(in)  :: m
		real, dimension(Nx,Ny), intent(out) :: mx
		integer :: i, j
		!$omp parallel do
		do i = 1, Nx
			do j = 1, Ny
				mx(i,j) = xcendiffpoint(m,i,j)
			end do
		end do
		!$omp end parallel do
	end subroutine

	function ycendiffpoint(m,i,j)
		implicit none
		real, dimension(Nx,Ny), intent(in) :: m
		integer, intent(in) :: i,j
		real :: ycendiffpoint
		if (j==1) then
			ycendiffpoint = (m(i,j+1)-m(i,j))
		end if
		if ((1<j).and.(j<Ny)) then
			ycendiffpoint = (m(i,j+1)-m(i,j-1))/2.0
		end if
		if (j==Ny) then
			ycendiffpoint = (m(i,j)-m(i,j-1))
		end if
	end function

	subroutine ycendiff(m,my)
		implicit none
		real, dimension(Nx,Ny), intent(in)  :: m
		real, dimension(Nx,Ny), intent(out) :: my
		integer :: i, j
		!$omp parallel do
		do i = 1, Nx
			do j = 1, Ny
				my(i,j) = ycendiffpoint(m,i,j)
			end do
		end do
		!$omp end parallel do
	end subroutine

	function LaplacePoint(m,i,j)
		implicit none
		real, dimension(Nx,Ny), intent(in) :: m
		integer, intent(in) :: i, j
		real :: LaplacePoint
		if ((j==1).or.(j==Ny)) then
			LaplacePoint = (m(pr(i),j)-2.0*m(i,j)+m(pl(i),j))/(dx*dx)
		else
			LaplacePoint = (m(pr(i),j)-2.0*m(i,j)+m(pl(i),j))/(dx*dx) + &
			               (m(i,j+1)-2.0*m(i,j)+m(i,j-1))/(dy*dy)
		end if
	end function

	function LaplacePoint2(m,i,j)
                implicit none
                real, dimension(Nx,Ny), intent(in) :: m
                integer, intent(in) :: i, j
                real :: LaplacePoint2
                if ((j==1).or.(j==Ny)) then
                        LaplacePoint2 = (m(pr(i),j)-2.0*m(i,j)+m(pl(i),j)) + 0.0
                else
                        LaplacePoint2 = (m(pr(i),j)-2.0*m(i,j)+m(pl(i),j)) + &
                                        (m(i,j+1)-2.0*m(i,j)+m(i,j-1))
                end if
        end function

	subroutine LaplaceOperator(m,ml)
		implicit none
		real, dimension(Nx,Ny), intent(in)  :: m
		real, dimension(Nx,Ny), intent(out) :: ml
		integer :: i, j
		!$omp parallel do
		do i = 1, Nx
			do j = 1, Ny
				ml(i,j) = LaplacePoint(m,i,j)
			end do
		end do
		!$omp end parallel do
	end subroutine

	subroutine reverseHelmholtz(m,mh)	! LaplaceOperator(mh) should = m !!!
		implicit none
		real, dimension(Nx,Ny), intent(in)  :: m
		real, dimension(Nx,Ny), intent(out) :: mh
		real, dimension(Nx,Ny) :: nowm, R
		integer :: i, j, N

		!--- criteria adjustable???
		!if (N>1000) then
			!criteria = criteria*10
		!end if

		N=1

		!--- first guess
		mh = 0.0
		call LaplaceOperator(mh,nowm)
		!print *, "nowm (M/m) =", maxval(nowm), minval(nowm)
		!print *, "m (M/m) =", maxval(m), minval(m)
		R = nowm-m
		!print *, "old R (M/m) =", maxval(R), minval(R)
		!print *, "Repsilon (M/m) =", maxval(Repsilon), minval(Repsilon)

		do
			!--- simultaneous relaxation
			if (relaxation=="simultaneous") then
				!$omp parallel do
				do i = 1, Nx
				do j = 1, Ny
					mh(i,j) = mh(i,j)+R(i,j)/4.0
				end do
				end do
				!$omp end parallel do
				!$omp parallel do
				do i = 1, Nx
				do j = 1, Ny
					R(i,j) = LaplacePoint2(mh,i,j)-m(i,j)
				end do
				end do
				!$omp end parallel do
			end if

			!---sequential relaxation
			if (relaxation=="sequential") then
				do i = 1, Nx
				do j = 1, Ny
					mh(i,j) = mh(i,j)+R(i,j)/4.0
					R(i,j) = LaplacePoint2(mh,i,j)-m(i,j)
				end do
				end do
			end if

			!print *, "new R (M/m) =", maxval(R), minval(R)
			if (maxval(abs(R))<criteria) then
				exit
			end if
			N = N+1
			!print *, N, "mh (M/m) =", maxval(mh), minval(mh)
			!print *, N, "nowm (M/m) =", maxval(nowm), minval(nowm)
		end do
	end subroutine

	subroutine bve(z,ztend)
		implicit none
		real, dimension(Nx,Ny), intent(in) :: z
		real, dimension(Nx,Ny), intent(out) :: ztend
		real, dimension(Nx,Ny) :: zeta, a1, a2, b1, b2
		real, dimension(Nx,Ny) :: J1, J2, J3, J4, J5, J6, J7
		integer :: j

		!--- calculate relative and absolute vorticity
		call LaplaceOperator(z,zeta)
		do j = 1, Ny
			zeta(:,j) = zeta(:,j)+f(j)
		end do

		!--- calculate jacobian
		call ycendiff(z,a1)		! u
		call xcendiff(zeta,a2)
		call xcendiff(z,b1)		! v
		call ycendiff(zeta,b2)
		J1 = a1*a2-b1*b2

		call xcendiff(zeta,a1)
		call ycendiff(a1*z,a2)
		call ycendiff(zeta,b1)
		call xcendiff(b1*z,b2)
		J2 = a2-b2

		call ycendiff(z,a1)
		call xcendiff(a1*zeta,a2)
		call xcendiff(z,b1)
		call ycendiff(b1*zeta,b2)
		J3 = a2-b2

		!J4 = (J1+J2)/2.0
		!J5 = (J2+J3)/2.0
		!J6 = (J1+J3)/2.0
		J7 = (J1+J2+J3)/3.0

		!--- calculate Helmholtz reversely
		call reverseHelmHoltz(J7,ztend)

		!--- southmost and northmost tendency = 0
		ztend(:,1) = 0.0
		ztend(:,Ny) = 0.0

		!print *, "a1(u) (M/m)=", maxval(-a1), minval(-a1)
		!print *, "a2(vorx) (M/m)=", maxval(a2), minval(a2)
		!print *, "b1(v) (M/m)=", maxval(b1), minval(b1)
		!print *, "b2(vory) (M/m)=", maxval(b2), minval(b2)
		!print *, "z(M/m)=", maxval(z), minval(z)
		!print *, "zeta(M/m)=", maxval(zeta), minval(zeta)
		!print *, "tend(M/m)=", maxval(tend), minval(tend)
		!print *, "ztend(M/m)=", maxval(ztend), minval(ztend)

	end subroutine

        subroutine rk2(z)
                implicit none
                real, dimension(Nx,Ny), intent(inout) :: z
                real, dimension(Nx,Ny) :: tend,q1,q2

                call bve(z,tend)
                q1 = dt*tend

                call bve(z+q1,tend)
                q2 = dt*tend-q1

                z = z+q1+0.5*q2
        end subroutine

        subroutine rk4(z)
                implicit none
                real, dimension(Nx,Ny), intent(inout)  :: z 
                real, dimension(Nx,Ny)   :: tend,q1,q2,q3,q4

                !--- step 1
                call bve(z,tend)
                q1 = dt*tend

                !--- step 2
                call bve(z+0.5*q1,tend)
                q2 = dt*tend

                !--- step 3
                call bve(z+0.5*q2,tend)
                q3 = dt*tend

                !--- step 4
                call bve(z+q3,tend)
                q4 = dt*tend

                !--- update
                z = z+(q1+2.0*q2+2.0*q3+q4)/6.0
        end subroutine

end module integrate
