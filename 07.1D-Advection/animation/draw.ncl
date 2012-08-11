load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
	N		= 100
	TotalLines	= numAsciiRow("../result")
	f	= asciiread("../result",(/TotalLines,N+2/),"float")
	idx	= f(:,0)
	idealx	= f(:,1)
	x	= f(:,2:)
	x!0	= "time"
	x!1	= "x"
	;printVarSummary(x)

	xaxis	= fspan(0.5,99.5,N)

	;--- plotting ---
	res				= True
	res@gsnDraw                     = False
	res@gsnFrame                    = False
	res@vpWidthF                    = 0.7
	res@vpHeightF                   = 0.3

	res@xyLineThicknessF		= 4.0
	res@xyLineColor			= "blue"

	res@tiXAxisString		= "X"
	res@tiYAxisString		= "Tracer"

	res@trXMinF			= 0
	res@trXMaxF			= 100
	res@trYMinF			= -2
	res@trYMaxF			= 6

        resl                                            = True
        resl@gsLineColor                                = "red"

	do i = 0, TotalLines-1
		fon	= sprintf("%05.0f",idx(i))
		res@tiMainString	= fon ;"Time = "+fon
		wks	= gsn_open_wks("eps",fon)
		plot	= gsn_csm_xy(wks,xaxis,x(i,:),res)
		dum	= gsn_add_polyline(wks,plot,(/idealx(i),idealx(i)/),(/-100,100/),resl)
		draw(plot)
		frame(wks)
		print("plotting ... "+fon)
		system("convert -trim -density 150 "+fon+".eps "+fon+".jpg")
		;system("rm -f "+fon+".eps")
	end do
end
