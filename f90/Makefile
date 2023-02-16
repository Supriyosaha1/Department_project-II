#############################################################################
# Fortran compiler options and directives

#SYSTYPE="mac"
SYSTYPE="macP"

OPTIONS = -cpp # -fopenmp
OPTIONS += -DDISPLAY_PROGRESS_PERCENT


ifeq ($(SYSTYPE),"mac")
F90         = gfortran -ffree-line-length-none -ffree-form  ###-fopenmp
#FFLAGS      = -g -fcheck=all -ftrapv -fbacktrace -Wall -fbounds-check
FFLAGS      = -O2
FFLAGS += $(OPTIONS)
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/
endif


ifeq ($(SYSTYPE),"macP")
F90         = mpif90-openmpi-gcc12 -ffree-line-length-none -ffree-form #mpif90 #-openmpi-gcc8
#FFLAGS      = -g -ffree-line-length-none -fbounds-check -fbacktrace -fcheck=all -Wall -ftrapv -ffree-form
FFLAGS      = -O3
FFLAGS += $(OPTIONS)
LDFLAGS     = #-lmpi
LDFLAGS     = #-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib/ #-lmpi
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
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDumpAMRModel PhotonsFromSourceModel rascas rascas-serial

PhotonsFromSourceModel: $(OBJSPHO) module_ssp_lib.o PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) module_ssp_lib.o PhotonsFromSourceModel.o

CreateDomDumpAMRModel: $(OBJSDOM) CreateDomDumpAMRModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) CreateDomDumpAMRModel.o

rascas-serial: $(OBJSPHO) rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) rascas-serial.o

rascas: $(OBJSMPI) rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
