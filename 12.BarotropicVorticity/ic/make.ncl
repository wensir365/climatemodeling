load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; You need download daily mean NCEP/NCAR Reanalysis first
; Example:
;   - hgt.2010.nc
;   - uwnd.2010.nc
;   - vwnd.2010.nc

begin

	lat0	= 55
	deltalat= 20
	days	= 5

	;--- generating grid ---
	lon		= fspan(0.5,359.5,360)
	lon@long_name	= "longitude"
	lon@units	= "degrees_east"

	lat		= fspan(lat0-deltalat+0.5,lat0+deltalat-0.5,2*deltalat)
	lat@long_name	= "latitude"
	lat@units	= "degrees_north"

	;--- input ---
	f		= addfile("hgt.2010.nc","r")
	time		= f->time
	N		= dimsizes(time)
	z		= short2flt(f->hgt(N-days*4:,{500},::-1,:))

	f2		= addfile("uwnd.2010.nc","r")
	u		= short2flt(f2->uwnd(N-days*4:,{500},::-1,:))
	f3		= addfile("vwnd.2010.nc","r")
	v		= short2flt(f3->vwnd(N-days*4:,{500},::-1,:))

	;--- regridding ---
	newz = linint2_Wrap(z&lon,z&lat,z,True,lon,lat,0)
	newu = linint2_Wrap(u&lon,u&lat,u,True,lon,lat,0)
	newv = linint2_Wrap(v&lon,v&lat,v,True,lon,lat,0)
	printVarSummary(newz)

	;--- output ---
	fon = "ic.nc"
	system("rm -f "+fon)
	fo = addfile(fon,"c")
	fo->z500 = newz
	fo->u500 = newu
	fo->v500 = newv

	;--- checking ---
	z0	= newz(0,:,:)
	phi	= 9.8*z0

	dims	= dimsizes(phi)
	do j = 0, dims(0)-1
		phi(j,:) = phi(j,:)/(2*7.3e-5*sin(lat(j)/180*3.1415926))
	end do
	;printMinMax(phi,True)

	uu	= center_finite_diff_n(phi,111194.8,False,0,0)*-1
	vv	= center_finite_diff_n(phi,63778.81,False,0,1)
	vor	= center_finite_diff_n(v,63778.81,False,0,1)-center_finite_diff_n(u,111194.8,False,0,0)

	printMinMax(uu,True)
	printMinMax(vv,True)
	printMinMax(vor,True)


	vorticity = uv2vr_cfd(u(0,{35:75},:),v(0,{35:75},:),f2->lat({35:75}),f2->lon,0)

	print("--- M/m of u/v---")
	printMinMax(newu(0,:,:),True)
	printMinMax(newv(0,:,:),True)
	print("--- true vorticity ---")
	printMinMax(vorticity,True)

	fon2	= "try.nc"
	system("rm -f "+fon2)
	fo2	= addfile(fon2,"c")
	fo2->u	= newu(0,:,:)
	fo2->uu = uu
end
