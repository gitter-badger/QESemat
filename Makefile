# the compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -Wall -g -fPIC
#FCFLAGS+= -freal-4-real-8
FCFLAGS+= -fdefault-real-8
FCFLAGS+= -ffpe-trap=invalid,zero,overflow
#FCFLAGS+= -O2
# tell the compiler where your .mod files are located
#FCFLAGS+= -I/usr/include
# tell the compiler where to put compiled modules
FCFLAGS+= -J./bin

# libraries needed for linking, unused in the examples
#LDFLAGS = -L./CERNlib.SLC5 -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = Test_Section QESemat
# Test_Flux Test_FluxTable
# QESemat
LIBS = libQES_flux.so libQES_event.so libQES_sect.so
DOCS = "doc"

# "make" builds all
all:  $(PROGRAMS)
libs: $(LIBS)
docs:
	$(MAKE) -C $(DOCS)
# DO NOT MOVE PREVIOUS 2 LINES LOWER THAN $^ DESCRIPTION! - otherwise will be compiled not what you expect

QESemat.o: +PhysMathConstants.o
EventRate.o: +PhysMathConstants.o
QESNuc_dQ2.o: +PhysMathConstants.o
QESkin_SM.o: +PhysMathConstants.o
d3sQES_dQ2dnudkF_SM.o: +PhysMathConstants.o
NucQESFF.o: +PhysMathConstants.o
MassNucleus.o: +PhysMathConstants.o
QESFree_dQ2.o: +PhysMathConstants.o
dsQESCC_dQ2.o: +PhysMathConstants.o
QESemat: +PhysMathConstants.o \
GeM.o MuL.o spline1.o \
DZEROX.o LambdaFunc.o DMINFC.o \
GeM_FV_SM.o MuL_Funs.o QESkin_SM.o QESkin.o NucQESFF.o rho_SM.o MassNucleus.o FactorPauli.o \
d3sQES_dQ2dnudkF_SM.o QESNuc_dQ2.o QESFree_dQ2.o dsQESCC_dQ2.o \
MA_QES_eff.o fui.o Flux.o EventRate.o

CrossSection.o: +PhysMathConstants.o
Test_Section: +PhysMathConstants.o \
GeM.o MuL.o spline1.o \
DZEROX.o LambdaFunc.o DMINFC.o \
GeM_FV_SM.o MuL_Funs.o QESkin_SM.o QESkin.o NucQESFF.o rho_SM.o MassNucleus.o FactorPauli.o \
d3sQES_dQ2dnudkF_SM.o QESNuc_dQ2.o QESFree_dQ2.o dsQESCC_dQ2.o \
MA_QES_eff.o fui.o Flux.o CrossSection.o

Test_Flux: Flux.o spline1.o
Test_FluxTable: Flux.o spline1.o

libQES_flux.so: Flux.o spline1.o
libQES_sect.so: QESNuc_dQ2.o QESFree_dQ2.o dsQESCC_dQ2.o
libQES_event.so: EventRate.o

BIN=bin/
SRC=src/
LIB=lib/
# IF YOU MOVE NEXT 2 LINES HIGHER THAN $^ DESCRIPTION nothing will be wrong. It seems...
#SOURCES=$(addprefix $(SRC),$(^:.o=.for))
SOURCES=$(addprefix $(SRC),$^)
OBJECTS=$(addprefix $(BIN),$^)
LIBRARIES=$(addprefix $(LIB),$^)
GARBAGE=$(foreach dir,$(BIN) ./,$(wildcard $(dir)*.o $(dir)*.mod  $(dir)*.MOD $(dir)fort.*))

# General rule for building shared library
%.so:
	ld -shared  $(OBJECTS) -o $(LIB)$@ 
	# $(FC) -shared $(FCFLAGS) -fPIC -o $@ $^ $(LDFLAGS)

# General rule for building prog from prog.o; $^ (GNU extension) is used in order to list additional object files on which the executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is used in order to list only the first prerequisite (the source file) and not the additional prerequisites such as module or include files
%.o: $(SRC)%.for
	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@
%.o: $(SRC)%.f90
	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@

# Utility targets
clean:
	$(RM) $(GARBAGE)
	$(RM) $(PROGRAMS)
	$(RM) $(LIBRARIES)

# DO NOT DELETE THIS LINE! - make depend depends on it
