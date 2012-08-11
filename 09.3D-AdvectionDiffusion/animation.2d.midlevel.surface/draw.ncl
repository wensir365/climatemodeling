load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
        fon     = getenv("FON")
        fin     = getenv("FIN")
        f       = addfile(fin,"r")

        x       = f->N(time|0,i|:,j|:,{k|25})
        time    = f->time

	x	= x/10.0

        dims    = dimsizes(x)

        xaxis   = fspan(0.,1.,dims(0))
        yaxis   = fspan(0.,1.,dims(1))

	;--- plotting ---
	rho	= 3.0
	theta	= 36.0
	phi	= 67.0
	xlb = "0.0 0.2 0.4 0.6 0.8 1.0"
	ylb = "0.0 0.2 0.4 0.6 0.8 1.0"
	zlb = "-0.2 0.0 0.2 0.4 0.6"
	wks	= gsn_open_wks("ps",fon)

	tdsetp ("FOV",25.)

	tdez2d(wks,xaxis,yaxis,x,rho,theta,phi,-6)
	tdgrds (wks,(/0,0,-0.2/),(/1,1,0.6/),(/1,1,1/),12,0)
	tdgrds (wks,(/0,0,-0.2/),(/1,1,0.6/),(/1,1,1/),12,1)
	tdlbls (wks,(/0,0,-0.2/),(/1,1,0.6/),(/xlb,ylb,zlb/),(/"Dimension X","Dimension Y","Z=25m Concentration"/),1)

	txres               = True
	txres@txFontHeightF = 0.03
	gsn_text_ndc(wks,fon,0.5,0.88,txres)

	frame(wks)
	system("convert -trim -density 120 "+fon+".ps "+fon+".jpg")
	system("rm -f "+fon+".ps")
end
