load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
	N	= 101
	fo	= getenv("FO")
	fn	= getenv("FN")
	f	= asciiread(fn,-1,"float")

	idealx	= f(0)
	idealy	= f(1)

	block1	= f(2:2+N*N-1)
	block2	= f(2+N*N:2+N*N*2-1)
	block3	= f(2+N*N*2:)

	utmp	= onedtond(block1,(/N,N/))
	vtmp	= onedtond(block2,(/N,N/))
	xtmp	= onedtond(block3,(/N,N/))

	utmp!0	= "i"
	utmp!1	= "j"
	u	= utmp(j|:,i|:)

	vtmp!0	= "i"
	vtmp!1	= "j"
	v	= vtmp(j|:,i|:)

	xtmp!0	= "i"
	xtmp!1	= "j"
	x	= xtmp(j|:,i|:)

	;--- plotting ---
	res						= True
	res@gsnDraw					= False
	res@gsnFrame					= False

	res@cnLinesOn					= True
	res@cnLineLabelsOn				= False
	res@cnInfoLabelOn				= False
	res@cnFillOn					= True
	res@gsnSpreadColors				= True
	res@gsnContourZeroLineThicknessF		= 2.0

	res@cnLevelSelectionMode			= "ManualLevels"
	res@cnMinLevelValF				= 0
	res@cnMaxLevelValF				= 4
	res@cnLevelSpacingF				= 0.2

	res@lbLabelBarOn                                = True
	;res@lbLabelAutoStride                           = True
	res@lbLabelStride                              = 2
	res@pmLabelBarWidthF                            = 0.5
	res@pmLabelBarHeightF                           = 0.05
	res@pmLabelBarOrthogonalPosF                    = 0.04
	res@lbLabelFontHeightF                          = 0.01
	res@lbTitleFontHeightF                          = 0.01

	res@tiMainString				= fo
	;res@tiXAxisString				= "X (0-100m)"
	;res@tiYAxisString				= "Y (0-100m)"

	wks	= gsn_open_wks("eps",fo)
	gsn_define_colormap(wks,"BlAqGrYeOrRe")
	plot	= gsn_csm_contour(wks,x,res)

	;--- adding analytical solution ---
	resm						= True
	resm@gsMarkerIndex				= 16
	resm@gsMarkerColor				= "white"
	resm@gsMarkerSizeF				= 0.005

	dum1	= gsn_add_polymarker(wks,plot,idealx,idealy,resm)
	print("adding ideal point at "+idealx+" / "+idealy)

	resl						= True
	resl@gsLineColor				= "white"
	dum2	= gsn_add_polyline(wks,plot,(/idealx,50/),(/idealy,50/),resl)

	draw(plot)
	frame(wks)
	
end
