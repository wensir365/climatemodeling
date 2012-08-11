module equation

	use data
	implicit none

contains

	subroutine advectiondiffusion(u,v,w,n,tend)
		implicit none
		real, dimension(1:Ni,1:Nj,1:Nk), intent(in) :: u,v,w,n
		real, dimension(1:Ni,1:Nj,1:Nk), intent(out) :: tend
		real, dimension(1:Ni,1:Nj,1:Nk) :: uN, vN, wN
		integer :: i,j,k
		real :: Adveci,Advecj,Adveck, Diffi,Diffj,Diffk

		uN = u*n
		vN = v*n
		wN = w*n
		tend = 0.0000

		do k = 2, Nk-1
			do j = 2, Nj-1
				do i = 1, Ni
					!=== ADVECTION ===
					!--- flux form ---
					Adveci = uN(pright(i,j,k),j,k)-uN(pleft(i,j,k),j,k)
					Adveci = Adveci/(2.0*di)

					!--- regular form ---
					!Adveci = n(pright(i,j,k),j,k)-n(pleft(i,j,k),j,k)
					!Adveci = Adveci/(2.0*di)
					!Adveci = Adveci*u(i,j,k)

					!--- flux form ---
					Advecj = vN(i,pnorth(i,j,k),k)-vN(i,psouth(i,j,k),k)
					Advecj = Advecj/(2.0*dj)
					!--- regular form ---
					!Advecj = n(i,pnorth(i,j,k),k)-n(i,psouth(i,j,k),k)
					!Advecj = Advecj/(2.0*dj)
					!Advecj = Advecj*v(i,j,k)

					!--- flux form ---
					Adveck = wN(i,j,pup(i,j,k))-wN(i,j,pdown(i,j,k))
					Adveck = Adveck/(2.0*dk)
					!--- regular form ---
					!Adveck = w(i,j,pup(i,j,k))-w(i,j,pdown(i,j,k))
					!Adveck = Adveck/(2.0*dk)
					!Adveck = Adveck*w(i,j,k)

					!=== DIFFUSION ===
					Diffi = n(pright(i,j,k),j,k)-2*n(i,j,k)+n(pleft(i,j,k),j,k)
					Diffi = Diffi/(di**2.0)
					Diffi = Ki*Diffi

					Diffj = n(i,pnorth(i,j,k),k)-2*n(i,j,k)+n(i,psouth(i,j,k),k)
					Diffj = Diffj/(dj**2.0)
					Diffj = Kj*Diffj

					Diffk = n(i,j,pup(i,j,k))-2*n(i,j,k)+n(i,j,pdown(i,j,k))
					Diffk = Diffk/(dk**2.0)
					Diffk = Kk*Diffk

					if (switch_diff) then
						tend(i,j,k) = tend(i,j,k)+Diffi+Diffj+Diffk
					end if
					if (switch_adve) then
						tend(i,j,k) = tend(i,j,k)-Adveci-Advecj-Adveck
					end if
				end do
			end do
		end do

	end subroutine

end module equation
