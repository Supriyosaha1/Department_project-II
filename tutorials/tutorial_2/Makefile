#############################################################################
# Fortran compiler options and directives

#SYSTYPE="mac"
SYSTYPE="macP"

RASCASPATH = ../../f90/

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
COMPOBJS = $(RASCASPATH)module_dust_model.o $(RASCASPATH)module_scatterer_model.o module_idealised_model.o $(RASCASPATH)module_gas_composition.o

GASOBJS = $(RASCASPATH)module_constants.o $(RASCASPATH)module_random.o $(RASCASPATH)module_utils.o $(RASCASPATH)module_voigt.o $(RASCASPATH)module_domain.o $(RASCASPATH)coolrates_module.o $(RASCASPATH)module_ramses.o $(RASCASPATH)module_uparallel.o $(COMPOBJS)

OBJSDOM = $(GASOBJS) $(RASCASPATH)module_select.o $(RASCASPATH)module_mesh.o

OBJSPHO = $(GASOBJS) $(RASCASPATH)module_mesh.o $(RASCASPATH)module_mock.o $(RASCASPATH)module_photon.o 

OBJSMPI = $(OBJSPHO) $(RASCASPATH)module_parallel_mpi.o $(RASCASPATH)module_worker.o $(RASCASPATH)module_master.o

OBJLYAFG = $(RASCASPATH)module_constants.o $(RASCASPATH)module_random.o $(RASCASPATH)module_utils.o $(RASCASPATH)module_voigt.o $(RASCASPATH)module_domain.o $(RASCASPATH)coolrates_module.o $(RASCASPATH)module_ramses.o $(RASCASPATH)module_uparallel.o

OBJFESC	= $(GASOBJS) $(RASCASPATH)module_mesh.o $(RASCASPATH)module_gray_ray.o


#############################################################################
#.SUFFIXES: .o .f90
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@

#############################################################################
all: clean CreateDomDumpAMRModel PhotonsFromSourceModel rascas rascas-serial 

PhotonsFromSourceModel: $(OBJSPHO) $(RASCASPATH)module_ssp_lib.o $(RASCASPATH)PhotonsFromSourceModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) $(RASCASPATH)module_ssp_lib.o $(RASCASPATH)PhotonsFromSourceModel.o

CreateDomDumpAMRModel: $(OBJSDOM) $(RASCASPATH)CreateDomDumpAMRModel.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSDOM) $(RASCASPATH)CreateDomDumpAMRModel.o

rascas-serial: $(OBJSPHO) $(RASCASPATH)rascas-serial.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSPHO) $(RASCASPATH)rascas-serial.o

rascas: $(OBJSMPI) $(RASCASPATH)rascas.o
	$(F90) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJSMPI) $(RASCASPATH)rascas.o

clean:
	rm -f $(OBJSDOM) $(OBJSPHO) $(OBJSMPI) *.o *.mod
#############################################################################
