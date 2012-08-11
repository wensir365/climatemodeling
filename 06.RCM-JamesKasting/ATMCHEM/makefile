# make file for running the photochemical model

###############################
# Set variables
###############################
#FC = f90
#PLATFORM = sun
#FCFLAG = -c -xtypemap=real:64,integer:64 -xtarget=ultra -xarch=v9
#FFLAG1 = -xtypemap=real:64,integer:64 -xtarget=ultra -xarch=v9

FC = ifort
#PLATFORM = linux
FCFLAG = -c -r8 -i8  
FFLAG1 = -r8 -i8 

OUTPUT = runchem
CDOTO  = ATMCHEM.o

# Main files
CCHEM = atm_chem

# Subroutines for coupling
CPATH0 = COUPLE
CINPI = input_interp

# Subdirectories  and subroutines from the PHOTOCHEMICAL model
PPATH1 = SETUP
PCGRID = grid
PCRATE = rates
PCRPHT = readphoto
PCDENS = densty
PCDFCO = difco
PCSATR = psatrat
CGASC = gascon
PCRSTR = readstar

PPATH2 = CHEM
PCTBDY = tbdy
PCPHOT = photo
PCO3PH = o3phot
PCDOCH = dochem
PCCHPL = chempl
PCLTNG = ltning
PCTSTR = twostr

PPATH3 = PRTCL
PCSDMT = sedmnt
PCAERT = aertab
PCAERC = aercon

PPATH4 = MATRIX
PSGBSL = sgbsl
PSGBFA = sgbfa
PSGTSL = sgtsl
PSAXPY = saxpy
PSSCAL = sscal


OBPATH = $(CDOTO)/$(CCHEM).o \
	 $(CDOTO)/$(CINPI).o \
	 $(CDOTO)/$(PCGRID).o \
	 $(CDOTO)/$(PCRATE).o \
	 $(CDOTO)/$(PCRPHT).o \
	 $(CDOTO)/$(PCRSTR).o \
	 $(CDOTO)/$(PCTBDY).o \
	 $(CDOTO)/$(PCPHOT).o \
	 $(CDOTO)/$(PCO3PH).o \
	 $(CDOTO)/$(PCDENS).o \
	 $(CDOTO)/$(PCDFCO).o \
	 $(CDOTO)/$(PCSDMT).o \
	 $(CDOTO)/$(PCSATR).o \
	 $(CDOTO)/$(PCAERT).o \
	 $(CDOTO)/$(PCAERC).o \
	 $(CDOTO)/$(PCDOCH).o \
	 $(CDOTO)/$(PCCHPL).o \
	 $(CDOTO)/$(PCLTNG).o \
	 $(CDOTO)/$(PCTSTR).o \
         $(CDOTO)/$(PSGTSL).o \
         $(CDOTO)/$(PSGBSL).o \
         $(CDOTO)/$(PSGBFA).o \
         $(CDOTO)/$(PSAXPY).o \
         $(CDOTO)/$(PSSCAL).o 
###############################
# Load line
###############################

$(OUTPUT) : $(OBPATH) 
	$(FC) $(FFLAG1) -o $(OUTPUT) $(OBPATH)

###############################
# Object compile lines
###############################
$(CDOTO)/$(CCHEM).o : $(CCHEM).f
	$(FC) $(FCFLAG)  $(CCHEM).f
	\mv $(CCHEM).o $(CDOTO)

### Compiling programs for coupling

$(CDOTO)/$(CINPI).o : $(CPATH0)/$(CINPI).f
	$(FC) $(FCFLAG) $(CPATH0)/$(CINPI).f
	\mv $(CINPI).o $(CDOTO)

### COMPILING THE PHOTOCHEMICAL CODE

$(CDOTO)/$(PCGRID).o : $(PPATH1)/$(PCGRID).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCGRID).f
	\mv $(PCGRID).o $(CDOTO)

$(CDOTO)/$(PCRATE).o : $(PPATH1)/$(PCRATE).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCRATE).f
	\mv $(PCRATE).o $(CDOTO)

$(CDOTO)/$(PCRPHT).o : $(PPATH1)/$(PCRPHT).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCRPHT).f
	\mv $(PCRPHT).o $(CDOTO)

$(CDOTO)/$(PCRSTR).o : $(PPATH1)/$(PCRSTR).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCRSTR).f
	\mv $(PCRSTR).o $(CDOTO)

$(CDOTO)/$(PCDFCO).o : $(PPATH1)/$(PCDFCO).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCDFCO).f
	\mv $(PCDFCO).o $(CDOTO)

$(CDOTO)/$(PCDENS).o : $(PPATH1)/$(PCDENS).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCDENS).f
	\mv $(PCDENS).o $(CDOTO)

$(CDOTO)/$(PCSATR).o : $(PPATH1)/$(PCSATR).f
	$(FC) $(FCFLAG) $(PPATH1)/$(PCSATR).f
	\mv $(PCSATR).o $(CDOTO)

$(CDOTO)/$(PCTBDY).o : $(PPATH2)/$(PCTBDY).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCTBDY).f
	\mv $(PCTBDY).o $(CDOTO)

$(CDOTO)/$(PCPHOT).o :$(PPATH2)/$(PCPHOT).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCPHOT).f
	\mv $(PCPHOT).o $(CDOTO)

$(CDOTO)/$(PCO3PH).o : $(PPATH2)/$(PCO3PH).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCO3PH).f
	\mv $(PCO3PH).o $(CDOTO)

$(CDOTO)/$(PCDOCH).o : $(PPATH2)/$(PCDOCH).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCDOCH).f
	\mv $(PCDOCH).o $(CDOTO)

$(CDOTO)/$(PCCHPL).o :$(PPATH2)/$(PCCHPL).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCCHPL).f
	\mv $(PCCHPL).o $(CDOTO)

$(CDOTO)/$(PCLTNG).o : $(PPATH2)/$(PCLTNG).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCLTNG).f
	\mv $(PCLTNG).o $(CDOTO)

$(CDOTO)/$(PCTSTR).o : $(PPATH2)/$(PCTSTR).f
	$(FC) $(FCFLAG) $(PPATH2)/$(PCTSTR).f
	\mv $(PCTSTR).o $(CDOTO)

$(CDOTO)/$(PCSDMT).o : $(PPATH3)/$(PCSDMT).f
	$(FC) $(FCFLAG) $(PPATH3)/$(PCSDMT).f
	\mv $(PCSDMT).o $(CDOTO)

$(CDOTO)/$(PCAERT).o : $(PPATH3)/$(PCAERT).f
	$(FC) $(FCFLAG) $(PPATH3)/$(PCAERT).f
	\mv $(PCAERT).o $(CDOTO)

$(CDOTO)/$(PCAERC).o : $(PPATH3)/$(PCAERC).f
	$(FC) $(FCFLAG) $(PPATH3)/$(PCAERC).f
	\mv $(PCAERC).o $(CDOTO)

$(CDOTO)/$(PSGTSL).o : $(PPATH4)/$(PSGTSL).f
	$(FC) $(FCFLAG) $(PPATH4)/$(PSGTSL).f
	\mv $(PSGTSL).o $(CDOTO)

$(CDOTO)/$(PSGBSL).o : $(PPATH4)/$(PSGBSL).f
	$(FC) $(FCFLAG) $(PPATH4)/$(PSGBSL).f
	\mv $(PSGBSL).o $(CDOTO)

$(CDOTO)/$(PSGBFA).o : $(PPATH4)/$(PSGBFA).f
	$(FC) $(FCFLAG) $(PPATH4)/$(PSGBFA).f
	\mv $(PSGBFA).o $(CDOTO)

$(CDOTO)/$(PSAXPY).o : $(PPATH4)/$(PSAXPY).f
	$(FC) $(FCFLAG) $(PPATH4)/$(PSAXPY).f
	\mv $(PSAXPY).o $(CDOTO)

$(CDOTO)/$(PSSCAL).o : $(PPATH4)/$(PSSCAL).f
	$(FC) $(FCFLAG) $(PPATH4)/$(PSSCAL).f
	\mv $(PSSCAL).o $(CDOTO)


