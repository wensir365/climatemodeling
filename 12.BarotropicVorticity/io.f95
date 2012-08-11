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

        subroutine input
                implicit none
                character(len=100) :: fin
                integer :: ncid,varid,j
                fin = "ic/ic.nc"
                call check( NF90_OPEN(trim(fin),NF90_NOWRITE,ncid) )
                call check( NF90_INQ_VARID(ncid,"z500",varid) )
                call check( NF90_GET_VAR(ncid,varid,phi,start=(/1,1,1/),count=(/Nx,Ny,1/)) )
                call check( NF90_CLOSE(ncid) )

		!--- convert geopotential height to streamfunction ---
		do j = 1, Ny
			phi(:,j) = phi(:,j)*g/(2.0*Omega*sin(lat(j)*Pi/180.0))
		end do
        end subroutine

	subroutine output_netcdf(t,mode)	! mode: 0=create; 1=append; 2=close
		use integrate
		implicit none
		integer, intent(in) :: t, mode
		character(len=100) :: id, fon, path
		integer, save :: ncid,timeid,jid,iid,vartimeid,varjid,variid,varuid,varvid,varhid
		character (len=8) :: date
		character (len=10) :: time
		character (len=5) :: zone
		integer, dimension(8) :: values
		integer :: i,j
		integer :: current
		real, dimension(Nx,Ny) :: z500, u500, v500

		if (mod(t,stepoutinterval)>0) then
			!print *, "     writting omitted"
			return
		end if

		current = t/stepoutinterval

		!--- convert back to geopotential height (m)
		do j = 1, Ny
			z500(:,j) = phi(:,j)*2.0*Omega*sin(lat(j)*Pi/180.0)/g
		end do

		!--- compute u and v
		call ycendiff(phi,u500)
		call xcendiff(phi,v500)
		u500 = -u500

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
		call check(NF90_DEF_VAR(ncid, "Z",    NF90_FLOAT, (/iid, jid, timeid/), varhid))
		call check(NF90_DEF_VAR(ncid, "U",    NF90_FLOAT, (/iid, jid, timeid/), varuid))
		call check(NF90_DEF_VAR(ncid, "V",    NF90_FLOAT, (/iid, jid, timeid/), varvid))
		!call check(NF90_ENDDEF(ncid))

		!--- attributes ---
		!call check(NF90_REDEF(ncid))
		call check(NF90_PUT_ATT(ncid, timeid, "long_name", "time"))
		call check(NF90_PUT_ATT(ncid, timeid, "units",     "hours since 0000-00-00 00:00:00"))
		call check(NF90_PUT_ATT(ncid, varjid, "long_name", "latitude"))
		call check(NF90_PUT_ATT(ncid, varjid, "units",     "degrees_north"))
		call check(NF90_PUT_ATT(ncid, variid, "long_name", "longitude"))
		call check(NF90_PUT_ATT(ncid, variid, "units",     "degrees_east"))
		call check(NF90_PUT_ATT(ncid, varhid, "long_name", "geopotential height"))
		call check(NF90_PUT_ATT(ncid, varhid, "units",     "m"))
		call check(NF90_PUT_ATT(ncid, varuid, "long_name", "zonal wind"))
		call check(NF90_PUT_ATT(ncid, varuid, "units",     "m/s"))
		call check(NF90_PUT_ATT(ncid, varvid, "long_name", "meridional wind"))
		call check(NF90_PUT_ATT(ncid, varvid, "units",     "m/s"))

		!--- global ---
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "time_step_in_seconds", dt))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "beta_plane",    logical2string(betaplane)))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "Model",   "Shallow Water Equations"))
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "History", "Xinyu Wen, Peking University, 2010"))
		call date_and_time(date,time,zone,values)
		call check(NF90_PUT_ATT(ncid, NF90_GLOBAL, "Stamp",   "date(yyyymmdd):"//date//" time(hhmmss.sss):"//time//" timezone(hhmm):"//zone))
		call check(NF90_ENDDEF(ncid))

		call check(NF90_PUT_VAR(ncid, variid,    lon))
		call check(NF90_PUT_VAR(ncid, varjid,    lat))
		end if



		if (mode==1) then
		!--- write ---
		call check(NF90_PUT_VAR(ncid, vartimeid, real(current), start=(/current+1/)     ))
		call check(NF90_PUT_VAR(ncid, varhid,    z500,          start=(/1,1,current+1/) ))
		call check(NF90_PUT_VAR(ncid, varuid,    u500,          start=(/1,1,current+1/) ))
		call check(NF90_PUT_VAR(ncid, varvid,    v500,          start=(/1,1,current+1/) ))
		print *, "                   writting ...", current
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
