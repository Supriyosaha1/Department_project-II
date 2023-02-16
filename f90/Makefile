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

# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -ffree-line-length-none -ffree-form
   LDMPIFLAGS =     #-lmpi
   ifeq ($(DEBUG),1)
      FFLAGS += -O0 -fbounds-check -Wuninitialized -Wall -g -ftrapv -fbacktrace -fimplicit-none -fcheck=all
      #FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan
   else
      FFLAGS += -O3
   endif
endif

# Intel compiler (ifort)
ifeq ($(COMPILER),INTEL)
   FFLAGS =
   LDMPIFLAGS = -lmpi
   #FFLAGS += -fp-model source
   ifeq ($(DEBUG),1)
      FFLAGS += -warn all -O0 -g -traceback -check bounds -check all -debug -ftrapuv -fpp
      #FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      FFLAGS += -O3 -fpp -ftz  -xavx -sox -ipo -fno-alias -fno-fnalias -no-prec-div -no-prec-sqrt -align all
   endif
endif

FFLAGS += $(OPTIONS)

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

#############################################################################
ifeq ($(MPI),1)
all: clean CreateDomDumpAMRModel PhotonsFromSourceModel rascas rascas-serial
else
all: clean CreateDomDumpAMRModel PhotonsFromSourceModel rascas-serial
endif

PhotonsFromSourceModel: $(OBJSPHO) PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) PhotonsFromSourceModel.o

CreateDomDumpAMRModel: $(OBJSDOM) CreateDomDumpAMRModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDumpAMRModel.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

clean:
	rm -f *.o *.mod
#############################################################################
