# make file for coupled run

###############################
# Set variables
###############################
FC = pgf90
#PLATFORM = linux
FCFLAG = -c
FFLAG1 =

OUTPUT = runclima
CDOTOC = CLIMA.o
CDOTOR = RRTM.o

# Main files
CSUR = clima

# Subroutines for coupling
CPATH0 = COUPLE
COUTP = output_photo
CINPI = input_interp

# Subdirectories and subroutines for the CLIMATE model
CPATHR = RRTM
CRTM = rrtm
CREG = rtreg
CRTR = rtr
CATM = rrtatm
CSET = setcoef
CTAU = taumol
CRGC = rtregcld
CRTC = rtrcld
CEXT = extra
CRTX = rtrcldmr
CRGX = rtregcldmr
CKGS = k_g
CCLD = cldprop
#CUTL = util_$(PLATFORM)
CTRNT = translatem

# These are the climate code files
CPATH1 = SETUP
CPROF = cprofile
CREAD = readsol
CGRID = cgrid
CIRES = irexpsums
CINTP = interp
CINOZ = interpozone
COZON = ozone
CSTAR = choose_star
CGASC = gascon

CPATH2 = CONVEC
CCONV = convec
CSATH = satrat
CRELH = relhum
CSATC = satco2

CPATH3 = RADTRANS
CSOL  = solar
CRAYL = rayley
CDSTI = delta2strir
CDSTS = delta2str
CMTRX = tridag
CPLAN = planck

CPATH4 = PRTCL
CGRDA = gridaer
CADAT = aerabsdata
CINPA = interpar1


OBPATH = $(CDOTOC)/$(COUTP).o \
         $(CDOTOC)/$(CINPI).o \
         $(CDOTOC)/$(CSUR).o \
	 $(CDOTOR)/$(CRTM).o \
	 $(CDOTOR)/$(CREG).o \
	 $(CDOTOR)/$(CRTR).o \
	 $(CDOTOR)/$(CATM).o \
	 $(CDOTOR)/$(CSET).o \
	 $(CDOTOR)/$(CTAU).o \
	 $(CDOTOR)/$(CRGC).o \
	 $(CDOTOR)/$(CRTC).o \
	 $(CDOTOR)/$(CEXT).o \
	 $(CDOTOR)/$(CRTX).o \
	 $(CDOTOR)/$(CRGX).o \
	 $(CDOTOR)/$(CCLD).o \
	 $(CDOTOR)/$(CKGS).o \
         $(CDOTOR)/$(CTRNT).o \
	 $(CDOTOC)/$(CPROF).o \
	 $(CDOTOC)/$(CREAD).o \
	 $(CDOTOC)/$(CSTAR).o \
	 $(CDOTOC)/$(CGRID).o \
	 $(CDOTOC)/$(CIRES).o \
	 $(CDOTOC)/$(CINTP).o \
	 $(CDOTOC)/$(CINOZ).o \
	 $(CDOTOC)/$(COZON).o \
	 $(CDOTOC)/$(CCONV).o \
	 $(CDOTOC)/$(CSATH).o \
	 $(CDOTOC)/$(CRELH).o \
	 $(CDOTOC)/$(CSATC).o \
	 $(CDOTOC)/$(CGASC).o \
	 $(CDOTOC)/$(CSOL).o  \
	 $(CDOTOC)/$(CRAYL).o \
	 $(CDOTOC)/$(CDSTI).o \
	 $(CDOTOC)/$(CDSTS).o \
	 $(CDOTOC)/$(CMTRX).o \
	 $(CDOTOC)/$(CPLAN).o \
	 $(CDOTOC)/$(CGRDA).o \
	 $(CDOTOC)/$(CADAT).o \
	 $(CDOTOC)/$(CINPA).o 
###############################
# Load line
###############################

$(OUTPUT) : $(OBPATH) 
	$(FC) $(FFLAG1) -o $(OUTPUT) $(OBPATH)

###############################
# Object compile lines
###############################


### Compiling programs for coupling

$(CDOTOC)/$(COUTP).o : $(CPATH0)/$(COUTP).f
	$(FC) $(FCFLAG) $(CPATH0)/$(COUTP).f
	\mv $(COUTP).o $(CDOTOC)

$(CDOTOC)/$(CINPI).o : $(CPATH0)/$(CINPI).f
	$(FC) $(FCFLAG) $(CPATH0)/$(CINPI).f
	\mv $(CINPI).o $(CDOTOC)


## COMPILING THE CLIMATE CODE

$(CDOTOC)/$(CSUR).o : $(CSUR).f
	$(FC) $(FCFLAG) $(CSUR).f
	\mv $(CSUR).o $(CDOTOC)

$(CDOTOC)/$(CPROF).o : $(CPATH1)/$(CPROF).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CPROF).f
	\mv $(CPROF).o $(CDOTOC)

$(CDOTOC)/$(CREAD).o : $(CPATH1)/$(CREAD).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CREAD).f
	\mv $(CREAD).o $(CDOTOC)

$(CDOTOC)/$(CSTAR).o : $(CPATH1)/$(CSTAR).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CSTAR).f
	\mv $(CSTAR).o $(CDOTOC)

$(CDOTOC)/$(CGRID).o : $(CPATH1)/$(CGRID).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CGRID).f
	\mv $(CGRID).o $(CDOTOC)

$(CDOTOC)/$(CIRES).o : $(CPATH1)/$(CIRES).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CIRES).f
	\mv $(CIRES).o $(CDOTOC)

$(CDOTOC)/$(CINTP).o : $(CPATH1)/$(CINTP).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CINTP).f
	\mv $(CINTP).o $(CDOTOC)

$(CDOTOC)/$(CINOZ).o : $(CPATH1)/$(CINOZ).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CINOZ).f
	\mv $(CINOZ).o $(CDOTOC)

$(CDOTOC)/$(COZON).o : $(CPATH1)/$(COZON).f
	$(FC) $(FCFLAG) $(CPATH1)/$(COZON).f
	\mv $(COZON).o $(CDOTOC)

$(CDOTOC)/$(CGASC).o : $(CPATH1)/$(CGASC).f
	$(FC) $(FCFLAG) $(CPATH1)/$(CGASC).f
	\mv $(CGASC).o $(CDOTOC)

$(CDOTOC)/$(CCONV).o : $(CPATH2)/$(CCONV).f
	$(FC) $(FCFLAG) $(CPATH2)/$(CCONV).f
	\mv $(CCONV).o $(CDOTOC)

$(CDOTOC)/$(CSATH).o : $(CPATH2)/$(CSATH).f
	$(FC) $(FCFLAG) $(CPATH2)/$(CSATH).f
	\mv $(CSATH).o $(CDOTOC)

$(CDOTOC)/$(CRELH).o : $(CPATH2)/$(CRELH).f
	$(FC) $(FCFLAG) $(CPATH2)/$(CRELH).f
	\mv $(CRELH).o $(CDOTOC)

$(CDOTOC)/$(CSATC).o : $(CPATH2)/$(CSATC).f
	$(FC) $(FCFLAG) $(CPATH2)/$(CSATC).f
	\mv $(CSATC).o $(CDOTOC)

$(CDOTOC)/$(CSOL).o : $(CPATH3)/$(CSOL).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CSOL).f
	\mv $(CSOL).o $(CDOTOC)

$(CDOTOC)/$(CRAYL).o : $(CPATH3)/$(CRAYL).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CRAYL).f
	\mv $(CRAYL).o $(CDOTOC)

$(CDOTOC)/$(CDSTI).o : $(CPATH3)/$(CDSTI).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CDSTI).f
	\mv $(CDSTI).o $(CDOTOC)

$(CDOTOC)/$(CDSTS).o : $(CPATH3)/$(CDSTS).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CDSTS).f
	\mv $(CDSTS).o $(CDOTOC)

$(CDOTOC)/$(CMTRX).o :$(CPATH3)/$(CMTRX).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CMTRX).f
	\mv $(CMTRX).o $(CDOTOC)

$(CDOTOC)/$(CPLAN).o : $(CPATH3)/$(CPLAN).f
	$(FC) $(FCFLAG) $(CPATH3)/$(CPLAN).f
	\mv $(CPLAN).o $(CDOTOC)

$(CDOTOC)/$(CGRDA).o : $(CPATH4)/$(CGRDA).f
	$(FC) $(FCFLAG) $(CPATH4)/$(CGRDA).f
	\mv $(CGRDA).o $(CDOTOC)

$(CDOTOC)/$(CADAT).o : $(CPATH4)/$(CADAT).f
	$(FC) $(FCFLAG) $(CPATH4)/$(CADAT).f
	\mv $(CADAT).o $(CDOTOC)

$(CDOTOC)/$(CINPA).o :$(CPATH4)/$(CINPA).f
	$(FC) $(FCFLAG) $(CPATH4)/$(CINPA).f
	\mv $(CINPA).o $(CDOTOC)

# RRTM files
$(CDOTOR)/$(CRTM).o : $(CPATHR)/$(CRTM).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRTM).f
	\mv $(CRTM).o $(CDOTOR)

$(CDOTOR)/$(CRTR).o : $(CPATHR)/$(CRTR).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRTR).f
	\mv $(CRTR).o $(CDOTOR)

$(CDOTOR)/$(CREG).o : $(CPATHR)/$(CREG).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CREG).f
	\mv $(CREG).o $(CDOTOR)

$(CDOTOR)/$(CATM).o : $(CPATHR)/$(CATM).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CATM).f
	\mv $(CATM).o $(CDOTOR)

$(CDOTOR)/$(CSET).o : $(CPATHR)/$(CSET).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CSET).f
	\mv $(CSET).o $(CDOTOR)

$(CDOTOR)/$(CTAU).o : $(CPATHR)/$(CTAU).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CTAU).f
	\mv $(CTAU).o $(CDOTOR)

$(CDOTOR)/$(CRGC).o : $(CPATHR)/$(CRGC).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRGC).f
	\mv $(CRGC).o $(CDOTOR)

$(CDOTOR)/$(CRTC).o : $(CPATHR)/$(CRTC).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRTC).f
	\mv $(CRTC).o $(CDOTOR)

$(CDOTOR)/$(CEXT).o :$(CPATHR)/$(CEXT).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CEXT).f
	\mv $(CEXT).o $(CDOTOR)

$(CDOTOR)/$(CRTX).o : $(CPATHR)/$(CRTX).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRTX).f
	\mv $(CRTX).o $(CDOTOR)

$(CDOTOR)/$(CRGX).o : $(CPATHR)/$(CRGX).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CRGX).f
	\mv $(CRGX).o $(CDOTOR)

$(CDOTOR)/$(CCLD).o : $(CPATHR)/$(CCLD).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CCLD).f
	\mv $(CCLD).o $(CDOTOR)

#$(CDOTOR)/$(CUTL).o : $(CPATHR)/$(CUTL).f
#	$(FC) $(FCFLAG) $(CPATHR)/$(CUTL).f
#	\mv $(CUTL).o $(CDOTOR)

$(CDOTOR)/$(CKGS).o : $(CPATHR)/$(CKGS).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CKGS).f
	\mv $(CKGS).o $(CDOTOR)

$(CDOTOR)/$(CTRNT).o : $(CPATHR)/$(CTRNT).f
	$(FC) $(FCFLAG) $(CPATHR)/$(CTRNT).f
	\mv $(CTRNT).o $(CDOTOR)


