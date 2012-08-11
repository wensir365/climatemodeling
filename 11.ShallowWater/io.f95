module io

	use netcdf
	use data
	implicit none

contains

	subroutine check(status)
		integer, intent(in) :: status
		if(status/=nf90_noerr) then 
		print *, trim(nf90_strerror(status))
		stop "Stopped due to NETCDF I/O"
		end if
	end subroutine

	subroutine output_netcdf(t,mode)	! mode: 0=create; 1=append; 2=close
		implicit none
		integer, intent(in) :: t, mode
		character(len=100) :: id, fon, path
		integer, save :: ncid,timeid,jid,iid,vartimeid,varjid,variid,varuid,varvid,varhid
		character (len=8) :: date
		character (len=10) :: time
		character (len=5) :: zone
		integer, dimension(8) :: values

		real, dimension(Nx,Ny) :: outu, outv, outh
		integer :: i,j

		!--- regrid to mass-point mesh
		do i = 1, Nx
			do j = 1, Ny
				if (i<Nx) then
					outu(i,j) = (u(i,j)+u(i+1,j))/2.0
				else
					outu(i,j) = (u(i,j)+u(1,j))/2.0
				end if
				outv(i,j) = (v(i,j)+v(i,j+1))/2.0
				outh(i,j) = h(i,j)
			end do
		end do

		!if (mod(t,10)>0) then		! output every 10 timesteps
		!	return
		!end if


		if (mode==0) then

		!--- generate filename ---
		write(id,"(i8.8)") t
		path = "./"
		fon = trim(path)//trim(description)//".output.nc"

		!--- open ---
		call check(NF90_CREATE(trim(fon), NF90_CLOBBER, ncid))

		!--- define dims ---
		call check(NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, timeid))
		call check(NF90_DEF_DIM(ncid, "lat",    Ny,             jid))
		call check(NF90_DEF_DIM(ncid, "lon",    Nx,             iid))

		!--- define vars ---
		call check(NF90_DEF_VAR(ncid, "time", NF90_FLOAT, (/timeid/), vartimeid))
		call check(NF90_DEF_VAR(ncid, "lat",  NF90_FLOAT, (/jid/),    varjid))
		call check(NF90_DEF_VAR(ncid, "lon",  NF90_FLOAT, (/iid/),    variid))
		call check(NF90_DEF_VAR(ncid, "U",    NF90_FLOAT, (/iid, jid, timeid/), varuid))
		call check(NF90_DEF_VAR(ncid, "V",    NF90_FLOAT, (/iid, jid, timeid/), varvid))
		call check(NF90_DEF_VAR(ncid, "H",    NF90_FLOAT, (/iid, jid, timeid/), varhid))
		!call check(NF90_ENDDEF(ncid))

		!--- attributes ---
		!call check(NF90_REDEF(ncid))
		call check(NF90_PUT_ATT(ncid, timeid, "long_name", "time"))
		call check(NF90_PUT_ATT(ncid, timeid, "units",     "seconds since 0000-00-00 00:00:00"))
		call check(NF90_PUT_ATT(ncid, varjid, "long_name", "latitude"))
		call check(NF90_PUT_ATT(ncid, varjid, "units",     "degrees_north"))
		call check(NF90_PUT_ATT(ncid, variid, "long_name", "longitude"))
		call check(NF90_PUT_ATT(ncid, variid, "units",     "degrees_east"))
		call check(NF90_PUT_ATT(ncid, varuid, "long_name", "zonal wind"))
		call check(NF90_PUT_ATT(ncid, varuid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varvid, "long_name", "meridional wind"))
		call check(NF90_PUT_ATT(ncid, varvid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varhid, "long_name", "height of free surface"))
		call check(NF90_PUT_ATT(ncid, varhid, "units",     "m"))

		!--- global ---
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "time_step_in_seconds", dt))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "term_advection", logical2string(term_advection)))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "term_coriolis",  logical2string(term_coriolis)))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "term_diffusion", logical2string(term_diffusion)))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "beta_effect",    logical2string(beta_effect)))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "Model",   "Shallow Water Equations"))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "History", "Xinyu Wen, Peking University, 2010"))
		call date_and_time(date,time,zone,values)
		call check(NF90_PUT_ATT(ncid,NF90_GLOBAL,"Stamp","date(yyyymmdd):"//date//" time(hhmmss.sss):"//time//" timezone(hhmm):"//zone))
		call check(NF90_ENDDEF(ncid))

		call check(NF90_PUT_VAR(ncid, variid,    lon))
		call check(NF90_PUT_VAR(ncid, varjid,    lat))
		end if

		if (mode==1) then
		!--- write ---
		call check(NF90_PUT_VAR(ncid, vartimeid, real(t)*dt, start=(/t/)     ))
		call check(NF90_PUT_VAR(ncid, varuid,    outu,       start=(/1,1,t/) ))
		call check(NF90_PUT_VAR(ncid, varvid,    outv,       start=(/1,1,t/) ))
		call check(NF90_PUT_VAR(ncid, varhid,    outh,       start=(/1,1,t/) ))
		end if

		if (mode==2) then
		!--- close ---
		call check(NF90_CLOSE(ncid))
		print *, "writting to ... "//trim(fon)
		end if

	end subroutine

	function logical2string(l)
		implicit none
		logical, intent(in) :: l
		character(len=10) :: logical2string
		if (l) then
			logical2string = "True"
		else
			logical2string = "False"
		end if
	end function

end module io
