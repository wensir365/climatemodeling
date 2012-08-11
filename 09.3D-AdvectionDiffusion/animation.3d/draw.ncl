load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
	fon	= getenv("FON")
	fin	= getenv("FIN")
	f	= addfile(fin,"r")

	x	= f->N(time|0,i|:,j|:,k|:)
	time	= f->time	

	dims	= dimsizes(x)

	xaxis	= fspan(0.,1.,dims(0))
	yaxis	= fspan(0.,1.,dims(1))
	zaxis	= fspan(0.,0.5,dims(2))

	;--- plotting ---
	rho	= 3.0
	theta	= 36.0
	phi	= 67.0
	xlb = "0 20 40 60 80 100"
	ylb = "0.0 0.2 0.4 0.6 0.8 1.0"
	zlb = "0.0 0.1 0.2 0.3 0.4 0.5"
	wks	= gsn_open_wks("ps",fon)

	tdsetp ("FOV",25.)

	tdez3d(wks,xaxis,yaxis,zaxis,x,0.2,rho,theta,phi,-6)
	tdgrds (wks,(/0,0,0.0/),(/1,1,0.5/),(/1,1,1/),12,0)
	tdgrds (wks,(/0,0,0.0/),(/1,1,0.5/),(/1,1,1/),12,1)
	tdlbls (wks,(/0,0,0.0/),(/1,1,0.5/),(/xlb,ylb,zlb/),(/"Dimension X","Dimension Y","Dimension Z"/),1)

	txres               = True
	txres@txFontHeightF = 0.03
	gsn_text_ndc(wks,"time = "+time+" (sec)",0.5,0.88,txres)

	frame(wks)
	system("convert -trim -density 120 "+fon+".ps "+fon+".jpg")
	system("rm -f "+fon+".ps")
end
