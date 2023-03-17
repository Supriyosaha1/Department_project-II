#############################################################################
# Fortran compiler options and directives

# Compiler flavor: GNU or INTEL
COMPILER = GNU
# compiler name
F90 = mpif90
# Do we want a debug build? 1=Yes, 0=No
DEBUG = 0
# Use of OpenMP? 1=Yes, 0=No
OPENMP = 1
# Use of MPI? 1=Yes, 0=No
MPI=1

RASCASPATH = ../../f90/

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
   FFLAGS =
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

OBJSBASE = $(RASCASPATH)module_constants.o $(RASCASPATH)module_random.o $(RASCASPATH)module_utils.o $(RASCASPATH)module_voigt.o \
	   $(RASCASPATH)module_domain.o $(RASCASPATH)coolrates_module.o $(RASCASPATH)module_ramses.o $(RASCASPATH)module_uparallel.o

OBJSGAS = $(OBJSBASE) $(RASCASPATH)module_dust_model.o $(RASCASPATH)module_scatterer_model.o module_idealised_model.o \
	  $(RASCASPATH)module_gas_composition.o

OBJSDOM = $(OBJSGAS) $(RASCASPATH)module_select.o $(RASCASPATH)module_mesh.o

OBJSPHO = $(OBJSGAS) $(RASCASPATH)module_ssp_lib.o $(RASCASPATH)module_mesh.o $(RASCASPATH)module_mock.o $(RASCASPATH)module_photon.o 

OBJSMPI = $(OBJSPHO) $(RASCASPATH)module_parallel_mpi.o $(RASCASPATH)module_worker.o $(RASCASPATH)module_master.o


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

PhotonsFromSourceModel: $(OBJSPHO) $(RASCASPATH)PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) $(RASCASPATH)PhotonsFromSourceModel.o

CreateDomDumpAMRModel: $(OBJSDOM) $(RASCASPATH)CreateDomDumpAMRModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) $(RASCASPATH)CreateDomDumpAMRModel.o

rascas-serial: $(OBJSPHO) $(RASCASPATH)rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) $(RASCASPATH)rascas-serial.o

rascas: $(OBJSMPI) $(RASCASPATH)rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) $(RASCASPATH)rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
