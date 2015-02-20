# ATTENTION:
# at first run 'update_revision_no.sh' in your svn directory and copy General/svn_rev.var into your actual working tree

CC=gfortran
CFLAGS=-c

# DEBUG
FLAGS= -g -fcheck=all -Wall -Wextra -fimplicit-none -ffree-line-length-none -Wno-maybe-uninitialized -Wtabs -fbacktrace

# compiler flag explanations
# COMMON
# -ffree-line-length-none: free line length in source code
# -c: compiling without linking -> produces object files
# -fimplicit-none: no implicit typing allowed; even if "implicit none" was not specified in a function, subroutine etc.
# DEBUG only (slow compiling and execution)
# -g: produce debugging symbols
# -Wtabs: tab characters in source code allowed
# -Wall: enable all common compiler warnings
# -Wextra: extra warnings
# -fcheck=all: enable all run-time tests (e.g. array subscripts etc.)
# -fbacktrace: output backtrace of run-time errors
# -Wno-maybe-uninitialized
# 	- -Wall turns on -Wmaybe-uninitialized: detects variables that are potentially used uninitialized
# 	- may cause spurious warnings, so turn it off with -Wno-maybe-uninitialized
# -Og: optimizations that do/should not interfere with debugging, gfortran 4.8+
# RELEASE only (fast execution but warnings and errors omitted as far as possible)
# -O1 [-O2, -O3]: optimize for speed
# -Os: optimize for size of compiled executable

# CC=pgfortran
# CFLAGS=-c -Mfree


SOURCES=./Hillslope/erosion_h.f90 \
./General/common_h.f90 \
./General/params_h.f90 \
./Hillslope/hymo_h.f90 \
./River/routing_h.f90 \
./General/time_h.f90 \
./General/utils_h.f90 \
./General/calcyear.f90 \
./General/climo_h.f90 \
./General/climo.f90 \
./General/petcalc.f90 \
./Reservoir/lake_h.f90 \
./Reservoir/reservoir_h.f90 \
./General/readgen.f90 \
./Hillslope/check_climate.f90 \
./Hillslope/erosion.f90 \
./Hillslope/etp_max.f90 \
./Hillslope/etp_soil.f90 \
./Hillslope/etp_soil_hour.f90 \
./Hillslope/model_state_io.f90 \
./Hillslope/hymo_all.f90 \
./Hillslope/readhymo.f90 \
./Hillslope/sedi_yield.f90 \
./Hillslope/sedi_yield_subbas.f90 \
./Hillslope/soildistr.f90 \
./Hillslope/soilwat.f90 \
./Reservoir/change_sec.f90 \
./Reservoir/eq1_wu.f90 \
./Reservoir/eq2_ashida.f90 \
./Reservoir/eq3_tsinghua.f90 \
./Reservoir/eq4_ackers.f90 \
./Reservoir/hydraul_res.f90 \
./Reservoir/lake.f90 \
./Reservoir/lake_routing.f90 \
./Reservoir/reservoir.f90 \
./Reservoir/reservoir_routing.f90 \
./Reservoir/sedbal.f90 \
./Reservoir/sedbal_lake.f90 \
./Reservoir/semres.f90 \
./Reservoir/vert_dist.f90 \
./River/bedload.f90 \
./River/muskingum.f90 \
./River/route_sediments.f90 \
./River/routing.f90 \
./River/routing_coefficients.f90 \
./River/routing_new.f90 \
./General/wasa.f90

EXECUTABLE=wasa.lin

# top-level rule to compile the whole program.
all: $(SOURCES) $(EXECUTABLE)
# program is made of several source files.
$(EXECUTABLE): $(SOURCES)
	#./update_revision_no.sh
	$(CC) $(FLAGS) $(SOURCES) -o $(EXECUTABLE)

.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm *.mod *.o $(EXECUTABLE)