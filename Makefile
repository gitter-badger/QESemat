# the compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -Wall -g
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
PROGRAMS = QESemat
#
#testflux test_fluxtable
# "make" builds all
all: $(PROGRAMS)
# DO NOT MOVE PREVIOUS 2 LINES LOWER THAN $^ DESCRIPTION! - otherwise will be compiled not what you expect

QESemat.o: +PhysMathConstants.o
EventRate.o: +PhysMathConstants.o
fui.o: +PhysMathConstants.o
QESNuc_dQ2.o: +PhysMathConstants.o
MuL_Funs.o: +PhysMathConstants.o
QESkin_SM.o: +PhysMathConstants.o
d3sQES_dQ2dnudkF_SM.o: +PhysMathConstants.o
rho_SM.o: +PhysMathConstants.o
GeM_FV_SM.o: +PhysMathConstants.o
NucQESFF.o: +PhysMathConstants.o
MassNucleus.o: +PhysMathConstants.o
QESFree_dQ2.o: +PhysMathConstants.o
QESkin.o: +PhysMathConstants.o
dsQESCC_dQ2.o: +PhysMathConstants.o
QESemat: +PhysMathConstants.o \
GeM.o MuL.o spline1.o \
DZEROX.o LambdaFunc.o DMINFC.o \
GeM_FV_SM.o MuL_Funs.o QESkin_SM.o QESkin.o NucQESFF.o rho_SM.o MassNucleus.o FactorPauli.o \
d3sQES_dQ2dnudkF_SM.o QESNuc_dQ2.o QESFree_dQ2.o dsQESCC_dQ2.o \
MA_QES_eff.o fui.o Flux.o EventRate.o

testflux: Flux.o spline1.o
test_fluxtable: Flux.o spline1.o

BIN=bin/
SRC=src/
# IF YOU MOVE NEXT 2 LINES HIGHER THAN $^ DESCRIPTION nothing will be wrong. It seems...
#SOURCES=$(addprefix $(SRC),$(^:.o=.for))
OBJECTS=$(addprefix $(BIN),$^)
GARBAGE=$(foreach dir,$(BIN) ./,$(wildcard $(dir)*.o $(dir)*.mod  $(dir)*.MOD $(dir)fort.*))

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

# DO NOT DELETE THIS LINE! - make depend depends on it
