load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
	TotalLines	= 900
	N		= 100
	f	= asciiread("../result",(/TotalLines,N+1/),"float")
	idx	= f(:,0)
	x	= f(:,1:)
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

	res@xyLineThicknessF		= 2.0
	res@xyLineColor			= "blue"

	res@tiXAxisString		= "X"
	res@tiYAxisString		= "U Wind (m/s)"

	res@trXMinF			= 0
	res@trXMaxF			= 100
	res@trYMinF			= -5
	res@trYMaxF			= 5

	do i = 0, TotalLines-1
		fon	= sprintf("%05.0f",idx(i))
		res@tiMainString	= "Time = "+fon
		wks	= gsn_open_wks("eps",fon)
		plot	= gsn_csm_xy(wks,xaxis,x(i,:),res)
		draw(plot)
		frame(wks)
		print("plotting ... "+fon)
		system("convert -trim -density 150 "+fon+".eps "+fon+".jpg")
		system("rm -f "+fon+".eps")
	end do
end
