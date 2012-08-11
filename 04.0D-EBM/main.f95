program main
	implicit none
	real, parameter :: g=0.4,dt=0.01
	real, parameter :: S0=1368.0,alpha=0.3,SBC=5.67e-8
	real, parameter :: k=0.025,c=1.0,m=1.0

	integer :: i
	real :: Ts=10.0,Ta=450.0, Tsnew,Tanew

	do i = 1,5000
		call forward(Ts,Ta,Tsnew,Tanew)
		print *, i, Ts, Ta
		Ts=Tsnew
		Ta=Tanew
	end do

contains
	subroutine forward(ts,ta,tsnew,tanew)
	  real, intent(in) :: ts,ta
	  real, intent(out) :: tsnew,tanew
	  real :: tendts,tendta

	  tendts = S0*(1-alpha)/4.0-SBC*Ts**4.0+g*SBC*Ta**4.0-k*(Ts-Ta)
	  tendta = SBC*Ts**4.0-SBC*Ta**4.0+k*(Ts-Ta)

	  tsnew = ts+dt*(tendts/m/c)
	  tanew = ta+dt*(tendta/m/c)
	end subroutine

end program main
