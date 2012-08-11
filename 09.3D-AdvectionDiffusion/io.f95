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

	subroutine output_netcdf(t)
		implicit none
		integer, intent(in) :: t
		character(len=100) :: id, description, fon, path
		integer :: ncid,timeid,kid,jid,iid,vartimeid,varkid,varjid,variid,varuid,varvid,varwid,varnid
		character (len=8) :: date
		character (len=10) :: time
		character (len=5) :: zone
		integer, dimension(8) :: values

		if (mod(t,10)>0) then		! output every 10 timesteps
			return
		end if

		!--- generate filename ---
		write(id,"(i8.8)") t
		description = "advection-diffusion.output."
		path = "result/"
		fon = trim(path)//trim(description)//trim(id)//".nc"

		!--- open ---
		call check(NF90_CREATE(trim(fon), NF90_CLOBBER, ncid))

		!--- define dims ---
		call check(NF90_DEF_DIM(ncid, "time", NF90_UNLIMITED, timeid))
		call check(NF90_DEF_DIM(ncid, "k",    Nk,             kid))
		call check(NF90_DEF_DIM(ncid, "j",    Nj,             jid))
		call check(NF90_DEF_DIM(ncid, "i",    Ni,             iid))

		!--- define vars ---
		call check(NF90_DEF_VAR(ncid, "time", NF90_FLOAT, (/timeid/), vartimeid))
		call check(NF90_DEF_VAR(ncid, "k",    NF90_FLOAT, (/kid/),    varkid))
		call check(NF90_DEF_VAR(ncid, "j",    NF90_FLOAT, (/jid/),    varjid))
		call check(NF90_DEF_VAR(ncid, "i",    NF90_FLOAT, (/iid/),    variid))
		call check(NF90_DEF_VAR(ncid, "U",    NF90_FLOAT, (/iid, jid, kid, timeid/), varuid))
		call check(NF90_DEF_VAR(ncid, "V",    NF90_FLOAT, (/iid, jid, kid, timeid/), varvid))
		call check(NF90_DEF_VAR(ncid, "W",    NF90_FLOAT, (/iid, jid, kid, timeid/), varwid))
		call check(NF90_DEF_VAR(ncid, "N",    NF90_FLOAT, (/iid, jid, kid, timeid/), varnid))
		call check(NF90_ENDDEF(ncid))

		!--- write ---
		call check(NF90_PUT_VAR(ncid, vartimeid, real(t)*dt))
		call check(NF90_PUT_VAR(ncid, variid,    posi(:,1,1)))
		call check(NF90_PUT_VAR(ncid, varjid,    posj(1,:,1)))
		call check(NF90_PUT_VAR(ncid, varkid,    posk(1,1,:)))
		call check(NF90_PUT_VAR(ncid, varuid,    u))
		call check(NF90_PUT_VAR(ncid, varvid,    v))
		call check(NF90_PUT_VAR(ncid, varwid,    w))
		call check(NF90_PUT_VAR(ncid, varnid,    n0))

		!--- attributes ---
		call check(NF90_REDEF(ncid))
		call check(NF90_PUT_ATT(ncid, timeid, "long_name", "time"))
		call check(NF90_PUT_ATT(ncid, timeid, "units",     "seconds since 0000-00-00 00:00:00"))
		call check(NF90_PUT_ATT(ncid, varkid, "long_name", "Z"))
		call check(NF90_PUT_ATT(ncid, varkid, "units",     "m"))
		call check(NF90_PUT_ATT(ncid, varjid, "long_name", "Y"))
		call check(NF90_PUT_ATT(ncid, varjid, "units",     "m"))
		call check(NF90_PUT_ATT(ncid, variid, "long_name", "X"))
		call check(NF90_PUT_ATT(ncid, variid, "units",     "m"))
		call check(NF90_PUT_ATT(ncid, varuid, "long_name", "zonal wind"))
		call check(NF90_PUT_ATT(ncid, varuid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varvid, "long_name", "meridional wind"))
		call check(NF90_PUT_ATT(ncid, varvid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varwid, "long_name", "vertical wind"))
		call check(NF90_PUT_ATT(ncid, varwid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varnid, "long_name", "tracer concentration"))
		call check(NF90_PUT_ATT(ncid, varnid, "units",     "kg/m3"))

		!--- global ---
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "Model",   "PKU Dispersion Model"))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "History", "Xinyu Wen, Peking University, 2010"))
		call date_and_time(date,time,zone,values)
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "Stamp",   "date(yyyymmdd):"//date//" time(hhmmss.sss):"//time//" timezone(hhmm):"//zone))
		call check(NF90_ENDDEF(ncid))

		!--- close ---
		call check(NF90_CLOSE(ncid))

		print *, "writting to ... "//trim(fon)
	end subroutine

end module io
