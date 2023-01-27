#############################################################################
# Fortran compiler options and directives

SYSTYPE="Cascade"
#SYSTYPE="CCF"
#SYSTYPE="lyocral"
#SYSTYPE="lyocralP"
#SYSTYPE="mac"
#SYSTYPE="macP"
#SYSTYPE="ccage"
#SYSTYPE="horizon"

OPTIONS = -cpp # -fopenmp
OPTIONS += -DDISPLAY_PROGRESS_PERCENT


ifeq ($(SYSTYPE),"Cascade")
F90      = mpif90  -ffree-line-length-none -ffree-form -fallow-argument-mismatch
#FFLAGS  =  -g -fcheck=all -ftrapv -fbounds-check -fbacktrace
FFLAGS   = -O3
FFLAGS  += $(OPTIONS)
LDFLAGS  = -lmpi
endif

ifeq ($(SYSTYPE),"PSMN")
F90         = mpif90
#FFLAGS      = -g -check all -traceback
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif

ifeq ($(SYSTYPE),"lyocral")
F90         = ifort
#FFLAGS      = -g -traceback -fpp -check all -debug -warn all -ftrapuv 
##FFLAGS      = -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS      = -O3 -fpp -ftz -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
##FFLAGS      = -p -g
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif

ifeq ($(SYSTYPE),"lyocralP")
F90         = mpif90 -ffree-line-length-none -ffree-form
FFLAGS      = #-g -fcheck=all -ftrapv -fbacktrace -Wall
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = -lmpi
endif


ifeq ($(SYSTYPE),"CCF")
F90         = ifort
##FFLAGS      = -g -fpp -check all -debug -warn all -ftrapuv -traceback
FFLAGS      = -O3 -fpp -ftz  -xavx -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
FFLAGS += $(OPTIONS)
LDFLAGS     =  -lmpi 
endif

ifeq ($(SYSTYPE),"mac")
F90         = gfortran -ffree-line-length-none -ffree-form  ###-fopenmp
FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -Wall -fbounds-check
#FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/
endif


ifeq ($(SYSTYPE),"macP")
F90         = mpif90 #-openmpi-gcc8
FFLAGS      = -g -ffree-line-length-none -fbounds-check -fbacktrace -fcheck=all -Wall -ftrapv -ffree-form
#FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = #-lmpi
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/ #-lmpi
endif


ifeq ($(SYSTYPE),"ccage")
F90         = mpif90
#FFLAGS      = -g -fcheck=all -ffree-line-length-none -ffree-form -ftrapv -fbounds-check -fbacktrace
FFLAGS      = -O3 -ffree-line-length-none -ffree-form
FFLAGS += $(OPTIONS)
LDFLAGS     = 
endif

ifeq ($(SYSTYPE),"horizon")
F90         = mpif90
FFLAGS      = -fopenmp -g -traceback -O3 -fpp -ftz  -xavx -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all -heap-arrays10
FFLAGS += $(OPTIONS)
LDFLAGS     =  #-lmpi
endif


#############################################################################
# All objects
COMPOBJS = module_dust_model.o module_scatterer_model.o module_idealised_model.o module_gas_composition.o


GASOBJS = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) module_select.o module_mesh.o

OBJSPHO = $(GASOBJS) module_mesh.o module_mock.o module_photon.o 

OBJSMPI = $(OBJSPHO) module_parallel_mpi.o module_worker.o module_master.o

OBJLYAFG = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o

OBJFESC	= $(GASOBJS) module_mesh.o module_gray_ray.o


#############################################################################
#.SUFFIXES: .o .f90
#.f90.o:
#	$(F90) -c  $(FFLAGS) $<

%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDump PhotonsFromSourceModel PhotonsFromStars LyaPhotonsFromGas rascas

PhotonsFromSourceModel: $(OBJSPHO) module_ssp_lib.o PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) module_ssp_lib.o PhotonsFromSourceModel.o

PhotonsFromStars: $(OBJSPHO) module_ssp_lib.o PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) module_ssp_lib.o PhotonsFromStars.o

LyaPhotonsFromGas: $(OBJLYAFG)  LyaPhotonsFromGas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJLYAFG)  LyaPhotonsFromGas.o

CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

CreateDomDumpAMRModel: $(OBJSDOM) CreateDomDumpAMRModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDumpAMRModel.o

#ExtractSubvol: $(OBJSDOM) ExtractSubvol.o
#	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) ExtractSubvol.o

#FescFromStars: $(OBJFESC) FescFromStars.o
#	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJFESC) FescFromStars.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
