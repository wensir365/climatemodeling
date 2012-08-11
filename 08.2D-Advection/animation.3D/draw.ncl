load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
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

	u	= onedtond(block1,(/N,N/))
	v	= onedtond(block2,(/N,N/))
	x	= onedtond(block3,(/N,N/))

	x	= x/10.0
	printMinMax(x,True)

	xaxis	= fspan(0.,1.,N)
	yaxis	= fspan(0.,1.,N)


	;--- plotting ---
	rho	= 3.0
	theta	= 36.0
	phi	= 67.0
	xlb = "0.0 0.2 0.4 0.6 0.8 1.0"
	ylb = "0.0 0.2 0.4 0.6 0.8 1.0"
	zlb = "-0.2 0.0 0.2 0.4 0.6"
	wks	= gsn_open_wks("ps",fo)

	tdsetp ("FOV",25.)

	tdez2d(wks,xaxis,yaxis,x,rho,theta,phi,-6)
	tdgrds (wks,(/0,0,-0.2/),(/1,1,0.6/),(/1,1,1/),12,0)
	tdgrds (wks,(/0,0,-0.2/),(/1,1,0.6/),(/1,1,1/),12,1)
	tdlbls (wks,(/0,0,-0.2/),(/1,1,0.6/),(/xlb,ylb,zlb/),(/"Dimension X","Dimension Y","Tracer Concentration"/),1)

	txres               = True
	txres@txFontHeightF = 0.03
	gsn_text_ndc(wks,fo,0.5,0.88,txres)

	frame(wks)
	system("convert -trim -density 120 "+fo+".ps "+fo+".jpg")
	;system("rm -f "+fo+".ps")
	;system("convert "+fo+".ps "+fo+".eps")
end
