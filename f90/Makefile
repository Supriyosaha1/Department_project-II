#############################################################################
# Fortran compiler options and directives

# Compiler flavor: GNU or INTEL
COMPILER = GNU
# compiler name
F90 = gfortran
# Do we want a debug build? 1=Yes, 0=No
DEBUG = 0
# Use of OpenMP? 1=Yes, 0=No
OPENMP = 1
# Use of MPI? 1=Yes, 0=No
MPI=1

# MacOSX with strange path to internal libraries
MACOS11 = 0

#############################################################################
OPTIONS = -cpp
OPTIONS += -DDISPLAY_PROGRESS_PERCENT 
ifeq ($(OPENMP),1)
   OPTIONS += -fopenmp
endif

#ifeq ($(MPI),1)
#   LDFLAGS = -lmpi
#endif

# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -ffree-line-length-none -ffree-form
   ifeq ($(DEBUG),1)
      FFLAGS += -O0 -fbounds-check -Wuninitialized -Wall -g -ftrapv -fbacktrace -fimplicit-none -fcheck=all
      #FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan
   else
      FFLAGS += -O3
   endif
endif

# Intel compiler (ifort)
ifeq ($(COMPILER),INTEL)
   #FFLAGS += -fp-model source
   ifeq ($(DEBUG),1)
      FFLAGS += -warn all -O0 -g -traceback -check bounds -check all -debug -ftrapuv -fpp
      #FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      FFLAGS += -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
   endif
endif

FFLAGS += $(OPTIONS)
FFLAGST = $(FFLAGS) -DTRANSP

ifeq ($(MACOS11),1)
   LDFLAGS += -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/
endif

#############################################################################
# All objects

OBJSBASE = module_constants.o module_random.o module_utils.o module_voigt.o module_domain.o coolrates_module.o module_ramses.o module_uparallel.o

OBJSGAS = $(OBJSBASE) module_dust_model.o module_scatterer_model.o module_idealised_model.o module_gas_composition.o

OBJSDOM = $(OBJSGAS) module_select.o module_mesh.o

OBJSPHO = $(OBJSGAS) module_ssp_lib.o module_mesh.o module_mock.o module_photon.o 

OBJSMPI = $(OBJSPHO) module_parallel_mpi.o module_worker.o module_master.o


#############################################################################
#.SUFFIXES: .o .f90
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

%.ot:%.f90
	$(F90) $(FFLAGST) -c $^ -o $@
#############################################################################
ifeq ($(MPI),1)
all: clean GenerateAMRmodel CreateDomDump PhotonsFromStars LyaPhotonsFromGas HaPhotonsFromGas PhotonsFromSourceModel rascas rascas-serial
else
all: clean GenerateAMRmodel CreateDomDump PhotonsFromStars LyaPhotonsFromGas HaPhotonsFromGas PhotonsFromSourceModel rascas-serial
endif

PhotonsFromSourceModel: $(OBJSPHO) PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromSourceModel.o

PhotonsFromStars: $(OBJSPHO) PhotonsFromStars.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromStars.o

LyaPhotonsFromGas: $(OBJSBASE)  LyaPhotonsFromGas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSBASE) LyaPhotonsFromGas.o

HaPhotonsFromGas: $(OBJSBASE)  HaPhotonsFromGas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSBASE) HaPhotonsFromGas.o

CreateDomDump: $(OBJSDOM) CreateDomDump.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDump.o

GenerateAMRmodel: $(OBJSDOM) GenerateAMRmodel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) GenerateAMRmodel.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

rascas_transp: $(OBJSMPI) module_photon.ot rascas.ot
	$(F90) $(FFLAGST) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.ot

clean:
	rm -f *.o *.mod *.ot
#############################################################################
