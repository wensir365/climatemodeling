load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
	fin	= getenv("FIN")
	fon	= getenv("FON")
	lev	= stringtofloat(getenv("LEV"))

	f	= addfile(fin,"r")
	n	= f->N(0,{lev},:,:)
	time	= f->time

        ;--- plotting ---
        res                                             = True
        res@gsnDraw                                     = False
        res@gsnFrame                                    = False

        res@cnLinesOn                                   = True
        res@cnLineLabelsOn                              = False
        res@cnInfoLabelOn                               = False
        res@cnFillOn                                    = True
        res@gsnSpreadColors                             = True
        res@gsnContourZeroLineThicknessF                = 2.0

        res@cnLevelSelectionMode                        = "ManualLevels"
        res@cnMinLevelValF                              = -2
        res@cnMaxLevelValF                              = 4
        res@cnLevelSpacingF                             = 0.2

        res@lbLabelBarOn                                = True
        res@lbLabelStride                               = 2
        res@pmLabelBarWidthF                            = 0.5
        res@pmLabelBarHeightF                           = 0.05
        res@pmLabelBarOrthogonalPosF                    = 0.04
        res@lbLabelFontHeightF                          = 0.01
        res@lbTitleFontHeightF                          = 0.01

        res@tiMainString                                = "Time="+time+"(s)"
	res@gsnCenterString				= lev+" (m)"
	res@gsnLeftString				= ""
	res@gsnRightString				= ""
        res@tiXAxisString                               = "X"
        res@tiYAxisString                               = "Y"

        wks     = gsn_open_wks("eps",fon)
        gsn_define_colormap(wks,"BlAqGrYeOrRe")
        plot    = gsn_csm_contour(wks,n,res)

        draw(plot)
        frame(wks)
end
