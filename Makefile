
#-------------------------------------------------------------------------
# Variables
#-------------------------------------------------------------------------

ifndef DEBUG
DEBUG=0
endif

# installation path
PREFIX=/usr/local

# source code and build main dirs
SRCDIR=src
OUTDIR=build

# adjust according to DEBUG
ifeq ($(DEBUG),1)
EXEC=wasa_dbg
    ifeq ($(OS),Windows_NT)
        OBJDIR=$(OUTDIR)\obj\debug
    else
        OBJDIR=$(OUTDIR)/obj/debug
    endif
else
    EXEC=wasa
    ifeq ($(OS),Windows_NT)
        OBJDIR=$(OUTDIR)\obj\release
    else
        OBJDIR=$(OUTDIR)/obj/release
    endif
endif

# script for updating revision number according to platform
ifeq ($(OS),Windows_NT)
    UPDATE_SCRIPT=$(SRCDIR)\update_revision_no.bat
else
    UPDATE_SCRIPT=$(SRCDIR)/update_revision_no.sh
endif

# source code will be searched for in this variable (necessary for build rule)
VPATH=$(SRCDIR)/General:$(SRCDIR)/Hillslope:$(SRCDIR)/Reservoir:$(SRCDIR)/River

# source files
FILES=Hillslope/erosion_h.f90 \
General/time_h.f90 \
General/utils_h.f90 \
General/common_h.f90 \
General/params_h.f90 \
Hillslope/hymo_h.f90 \
River/routing_h.f90 \
General/climo_h.f90 \
Hillslope/snow_h.f90 \
General/calcyear.f90 \
General/climo.f90 \
General/petcalc.f90 \
Reservoir/lake_h.f90 \
Reservoir/reservoir_h.f90 \
General/readgen.f90 \
Hillslope/check_climate.f90 \
Hillslope/erosion.f90 \
Hillslope/etp_max.f90 \
Hillslope/etp_soil.f90 \
Hillslope/etp_soil_hour.f90 \
Hillslope/model_state_io.f90 \
Hillslope/hymo_all.f90 \
Hillslope/readhymo.f90 \
Hillslope/sedi_yield.f90 \
Hillslope/sedi_yield_subbas.f90 \
Hillslope/soilwat.f90 \
Hillslope/snow_compute.f90 \
Reservoir/change_sec.f90 \
Reservoir/eq1_wu.f90 \
Reservoir/eq2_ashida.f90 \
Reservoir/eq3_tsinghua.f90 \
Reservoir/eq4_ackers.f90 \
Reservoir/hydraul_res.f90 \
Reservoir/lake.f90 \
Reservoir/lake_routing.f90 \
Reservoir/reservoir.f90 \
Reservoir/reservoir_routing.f90 \
Reservoir/sedbal.f90 \
Reservoir/sedbal_lake.f90 \
Reservoir/semres.f90 \
Reservoir/vert_dist.f90 \
River/bedload.f90 \
River/muskingum.f90 \
River/route_sediments.f90 \
River/routing.f90 \
River/routing_coefficients.f90 \
River/routing_new.f90 \
General/wasa.f90

# add main source code directory as prefix
SRC=$(addprefix $(SRCDIR)/,$(FILES))

# object files with their directory
OBJ_FILES=$(notdir $(SRC))
OBJ=$(OBJ_FILES:%.f90=$(OBJDIR)/%.o)


#-------------------------------------------------------------------------
# Compiler and Linker Flags
# (Un)comment according to your needs!
#-------------------------------------------------------------------------

#******gfortran*********
# IMPORTANT: GNU compiler since version 5.x changed behaviour: https://gcc.gnu.org/gcc-5/changes.html
# affects:
#   - use '-Wno-tabs' instead of '-Wtabs' for use of tab characters in source code
#
# compiler flag explanations
# COMMON
# -ffree-line-length-none: free line length in source code
# -c: compiling without linking -> produces object files (not needed to specify in CFLAGS!)
# -fimplicit-none: no implicit typing allowed; even if "implicit none" was not specified in a function, subroutine etc.
# DEBUG only (slow compiling and execution)
# -g: produce debugging symbols
# -Wtabs: tab characters in source code allowed; use -Wno-tabs instead since GNU v. 5.x
# -Wall: enable all common compiler warnings
# -Wextra: extra warnings
# -fcheck=all: enable all run-time tests (e.g. array subscripts etc.)
# -fbacktrace: output backtrace of run-time errors
# -Wno-maybe-uninitialized
# 	- -Wall turns on -Wmaybe-uninitialized: detects variables that are potentially used uninitialized
# 	- may cause spurious warnings, so turn it off with -Wno-maybe-uninitialized
# -Og: optimizations that do/should not interfere with debugging, gfortran 4.8+
# -ffpe-trap=list: enable certain floeating poinbt exception traps with list a comma separated list with: invalid, zero, overflow, underflow, inexact, and/or denormal
# RELEASE only (fast execution but warnings and errors omitted as far as possible)
# -O1 [-O2, -O3]: optimize for speed
# -Os: optimize for size of compiled executable

FC=gfortran
CDFLAGS=-g -fcheck=all -Wall -Wextra -ffree-line-length-none -fimplicit-none -Wno-maybe-uninitialized -Wno-tabs -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow
CRFLAGS=-ffree-line-length-none -fimplicit-none -Wno-maybe-uninitialized -Wno-tabs -O3 -s
LFLAGS=

ifeq ($(DEBUG),1)
CFLAGS=$(CDFLAGS)
else
CFLAGS=$(CRFLAGS)
endif


#******PGI*********
# FC=pgfortran
# CFLAGS=-Mfree
# LFLAGS=


#-------------------------------------------------------------------------
# Build-Rules
#-------------------------------------------------------------------------

$(OBJDIR)/%.o: %.f90
	$(FC) -J$(OBJDIR) -c $(CFLAGS) $< -o $@


#-------------------------------------------------------------------------
# Main-Targets
#-------------------------------------------------------------------------

help:
	@echo ""
	@echo "Building the WASA-SED model: make [OPTION] target"
	@echo ""
	@echo "Targets  : help      - Shows this text."
	@echo "           all       - Creates the model executable."
	@echo "           install   - Installs the model executable into PREFIX/bin (directory has to be accessible!)."
	@echo "           uninstall - Removes the model executable from PREFIX/bin (directory has to be accessible!)."
	@echo "           clean     - Deletes object files and local executables (retains build directories)."
	@echo "           distclean - Deletes the whole build directory with subdirectories."
	@echo ""
	@echo "Options  : DEBUG=1   - Creates a debug version."
	@echo "         : PREFIX    - Model executable installation path (default: /usr/local)."
	@echo ""
	@echo "Examples : make DEBUG=1 all"
	@echo "           make DEBUG=1 PREFIX=/my/path install"
	@echo "           make all"
	@echo "           make clean"
	@echo ""

all: prepare update_rev $(OUTDIR)/bin/$(EXEC)
	@echo ""
	@echo The generated WASA-SED executable can be found in $(OUTDIR)/bin.
	@echo FINISHED!

install:
	@echo Installing the model into $(PREFIX)/bin ...
	-@mkdir -p $(PREFIX)/bin
	-@cp $(OUTDIR)/bin/$(EXEC) $(PREFIX)/bin/$(EXEC)
	@echo FINISHED!

uninstall:
	@echo Removing $(EXEC) from $(PREFIX)/bin ...
	-@rm -f $(PREFIX)/bin/$(EXEC)
	@echo FINISHED!

clean:
	@echo "Cleaning the workspace (directories will be retained) ..."
	-@rm -f $(OBJDIR)/*.o
	-@rm -f $(OBJDIR)/*.mod
	-@rm -f $(OUTDIR)/bin/$(EXEC)
	@echo FINISHED!

distclean:
	@echo "Cleaning the workspace including all build (sub-)directories ..."
	-@rm -rf $(OUTDIR)
	@echo FINISHED!
	
#-------------------------------------------------------------------------
# Sub-Targets
#-------------------------------------------------------------------------

prepare:
	@echo "Output directories will be created (if they do not exist) ..."
	-@mkdir $(OUTDIR) 
	-@mkdir $(OBJDIR)
	ifeq ($(OS),Windows_NT)
        -@mkdir $(OUTDIR)\bin
    else
        -@mkdir $(OUTDIR)/bin
    endif
	
update_rev:
	@echo "Updating revision number ..."
	@echo "Compiling model source code ..."
	@$(UPDATE_SCRIPT) $(SRCDIR)

$(OUTDIR)/bin/$(EXEC): $(OBJ)
	@echo "Linking code ..."
	$(FC) $(LFLAGS) -o $(OUTDIR)/bin/$(EXEC) $(OBJ)
