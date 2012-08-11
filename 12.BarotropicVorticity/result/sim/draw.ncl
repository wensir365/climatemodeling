load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"

begin
	;--- open file ---
	fn	= "../../nh.2010-11-17-00.output.nc"
	f	= addfile(fn,"r")
	x	= f->Z

	time	= f->time
	time@units	= "hours since 2010-11-17 00:00:00"
	Ntime	= dimsizes(time)

	;--- s3: plotting ---
	res							= True
	res@gsnDraw						= False
	res@gsnFrame						= False

	res@mpPerimOn						= True
	res@gsnPolar						= "NH"
	res@mpMinLatF						= 30
	res@mpMaxLatF						= 90
	res@mpCenterLonF					= 90
	res@mpGridPolarLonSpacingF				= 90
	res@mpFillOn						= False

	res@cnFillOn						= True
	res@cnLinesOn						= False
	res@cnLineLabelsOn					= False
	res@gsnSpreadColors					= True
	res@gsnAddCyclic					= True

	res@mpPerimLineThicknessF				= 3
	res@mpOutlineBoundarySets				= "National"
	res@mpGeophysicalLineThicknessF				= 1
	res@mpGridLatSpacingF					= 90
	res@mpGridLonSpacingF					= 90
	res@mpGridLineThicknessF				= 0.2
	res@mpGridLineColor					= "black"
	res@mpGridLineDashPattern				= 0
	res@mpNationalLineThicknessF				= 0.5

	res@tmXBLabelFontHeightF				= 0.0065

	res@cnLevelSelectionMode				= "ManualLevels"
	res@cnMinLevelValF					= 5100
	res@cnMaxLevelValF					= 5800
	res@cnLevelSpacingF					= (res@cnMaxLevelValF - res@cnMinLevelValF)/40

	res@lbLabelBarOn					= True
	res@lbBoxLinesOn					= False
        res@lbLabelAutoStride					= True
	res@lbLabelStride					= 4
        res@pmLabelBarOrthogonalPosF				= 0.05
        res@pmLabelBarWidthF					= 0.3
        res@pmLabelBarHeightF					= 0.05
        res@lbLabelFontHeightF					= 0.01
        res@lbTitleFontHeightF					= 0.014
	res@lbTitleOn						= True
	res@lbTitleString					= "Z~B~500hPa~N~ (gpm)"
	res@lbTitlePosition					= "Top"

	res@gsnLeftString					= ""
	res@gsnRightString					= ""
	res@gsnCenterString					= "NH"

	;--- add Beijing marker
	bjres				= True
	bjres@gsMarkerIndex		= 12
	bjres@gsMarkerSizeF		= 7
	bjres@gsMarkerColor		= "red"
	bjres@gsMarkerThicknessF	= 2
	
        ;--- add z500
        res2                            = True
        res2@cnFillOn                   = False
        res2@cnLinesOn                  = True
        res2@cnLineLabelsOn             = False
        res2@cnInfoLabelOn              = False

	res2@cnLineDashPattern		= 1
        res2@cnLineThicknessF           = 2
        res2@cnLineColor                = "black"
        res2@cnLineLabelFontHeightF     = 0.008
        res2@cnLineLabelFontColor       = "black"

        ;res2@cnLevelSelectionMode       = "ManualLevels"
        ;res2@cnMinLevelValF             = 5500
        ;res2@cnMaxLevelValF             = 5530
        ;res2@cnLevelSpacingF            = 40

	res2@cnLevelSelectionMode	= "ExplicitLevels"
	res2@cnLevels			= (/ 5400 /)

        res2@gsnDraw                    = False
        res2@gsnFrame                   = False
        res2@gsnTickMarksOn             = False
        res2@gsnLeftString              = ""
        res2@gsnRightString             = ""
        res2@tiXAxisString              = ""
        res2@tiYAxisString              = ""

do i = 0, Ntime-1
	res@tiMainString	= ut_string(time(i),"%Y/%N/%d %HZ")
	id			= sprinti("%04i",i)
	fon	= "plot.nh500."+id
	wks	= gsn_open_wks("eps",fon)
	gsn_define_colormap(wks,"BlAqGrYeOrRe");"gui_default") ;"GreenYellow")

	plot = gsn_csm_contour_map(wks,x(i,:,:),res)
	dum_bj	= gsn_add_polymarker(wks,plot,116.3163,39.9924,bjres)
        dum_z500 = gsn_csm_contour(wks,x(i,:,:),res2)

	overlay(plot,dum_z500)

	draw(plot)
	frame(wks)

	system("convert -trim -density 120 "+fon+".eps "+fon+".png")
	system("rm -f "+fon+".eps")

	print(res@tiMainString+" ... done!   file:"+fon)
end do

end
