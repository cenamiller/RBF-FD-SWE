ifndef MAKE_DEFS

# swe model top directory
ARCH_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Directories
SOURCE_DIR          =	$(TOP_DIR)/swe_code
MAIN_DIR            =	$(SOURCE_DIR)/main
IO_DIR	            =	$(SOURCE_DIR)/io
RCM_DIR	            =	$(SOURCE_DIR)/rcm
LAYOUT_DIR	    =	$(SOURCE_DIR)/layout
MPI_DIR		    =	$(SOURCE_DIR)/mpi
OCL_DIR		    =	$(SOURCE_DIR)/ocl
COMMON_DIR	    =	$(MAIN_DIR)

# header includes
INCLUDE_COMMON	    =	-I$(MAIN_DIR)/include
INCLUDE_IO	    =	-I$(IO_DIR)/include
INCLUDE_RCM	    =	-I$(RCM_DIR)/include
INCLUDE_LAYOUT	    =	-I$(LAYOUT_DIR)/include
INCLUDE_MPI	    =	-I$(MPI_DIR)/include
INCLUDE_MAIN	    =	-I$(MAIN_DIR)/include
INCLUDE_OCL	    =	-I$(OCL_DIR)/include

INCLUDE_ALL	    =	$(INCLUDE_COMMON) $(INCLUDE_IO) $(INCLUDE_RCM) $(INCLUDE_LAYOUT) $(INCLUDE_MPI) $(INCLUDE_MAIN) $(INCLUDE_OCL)

# libs
SWE_LIBS	    =   $(MPI_DIR)/swe_mpi.a $(IO_DIR)/swe_io.a $(RCM_DIR)/swe_rcm.a $(LAYOUT_DIR)/swe_dl.a $(MAIN_DIR)/swe_main.a $(OCL_DIR)/swe_ocl.a

# linux utilities
MKDIR_P		    =	mkdir -p
LD_R		    =	ld -r

# define user model configuration options
ifndef CONFIG_DIR
CONFIG_DIR          =   intel
endif
include $(ARCH_DIR)arch/$(CONFIG_DIR)/config.swe
TOP_DIR             =   $(dir $(abspath $(dir $(abspath $(ARCH_DIR)))))

### Parallelization Defs
# MPI parallelization
ifeq ($(MPI),1)
	CONFIG          +=  -DUSE_MPI
	CC              =   $(MPICC)
endif

# OpenMP/OpenACC
ifeq ($(OPENMP),1)
	SMPAR_FLAGS     =   $(OMP_FLAGS)
endif

ifeq ($(OPENACC),1)
	SMPAR_FLAGS     =   $(ACC_FLAGS)
endif

# OpenCL
ifeq ($(OPENCL),1)
	CONFIG 		+=  -DUSE_OCL
ifeq ($(shell uname -s), Darwin)
	OCL_LIBS 	=   -framework OpenCL
	CONFIG		+=  -DOCL_FOLDER
endif
	SWE_LIBS	+=	$(OCL_LIBS)
	CONFIG		+=	$(OCL_FLAGS)
endif

ifeq ($(SPLIT_DEV),1)
	CONFIG		+=  -DSPLITDEV
endif

### I/O Options
# I/O method
ifeq ($(NCIO),1)
ifndef NETCDF
	NETCDF = /opt/local
endif
	CONFIG		+=	-DUSE_NCIO
	INCLUDE_NETCDF	+=	-I$(NETCDF)/include
	NETCDF_LIBS	+=	-L$(NETCDF)/lib -lnetcdf
	INCLUDE_ALL 	+=	$(INCLUDE_NETCDF)
	SWE_LIBS	+=	$(NETCDF_LIBS)
endif

ifeq ($(BINIO),1)
    CONFIG          +=  -DUSE_BINIO
endif

# History output
ifeq ($(HISTORY),1)
    CONFIG          +=  -DUSE_HIST
endif

### Layout Options
# DM/SVM Data Layout Options
ifeq ($(CFDL),1)
    CONFIG          +=  -DUSE_CFDL
endif

ifeq ($(SFDL),1)
    CONFIG          +=  -DUSE_SFDL
endif

# Padding/simd sizes
ifeq ($(PADDED_NUM),0)
    CONFIG          +=  -DPADDED_NUM=$(SIMD_LENGTH)
else
    CONFIG          +=  -DPADDED_NUM=$(PADDED_NUM)
endif
CONFIG              +=  -DSIMD_LENGTH=$(SIMD_LENGTH)

# differentiation evaluation loop vectorization method
ifeq ($(RHS_SIMD_METHOD),1)
    CONFIG          +=  -DRHS_OUTER_SIMD
endif
ifeq ($(RHS_SIMD_METHOD),2)
    CONFIG          +=  -DRHS_INNER_SIMD
endif

# check whether compiler is gcc or something else (small change needed in rk function)
ifeq ($(GCC),1)
    CONFIG          +=-DUSE_GCC
endif

# executable
ifndef EXEC_DIR
    EXEC_DIR        =   $(TOP_DIR)/run/default
endif

ifndef EXEC
    EXEC            =   swe_*.exe
endif

# standard compiler flags
CFLAGS	            +=	$(C99_FLAGS) $(ARCH_FLAGS) $(OPT_FLAGS) $(SMPAR_FLAGS) $(CONFIG)

# define MAKE_DEFS condition to avoid including twice
MAKE_DEFS           =   true

endif
