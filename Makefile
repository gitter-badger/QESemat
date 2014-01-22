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
PROGRAMS = testflux test_fluxtable qesemat
# "make" builds all
all: $(PROGRAMS)
# DO NOT MOVE PREVIOUS 2 LINES LOWER THAN $^ DESCRIPTION! - otherwise will be compiled not what you expect

qesemat.o: +PhysMathConstants.o +InpOutUnits.o
fui.o: +PhysMathConstants.o
dsQESCC_dQ2_SM.o: +PhysMathConstants.o
FunMuL_SM.o: +PhysMathConstants.o
QESkin_SM.o: +PhysMathConstants.o
d3sQES_dQ2dnudkF_SM.o: +PhysMathConstants.o
rho_SM.o: +PhysMathConstants.o
FunGeM_SM.o: +PhysMathConstants.o
NucQESFF.o: +PhysMathConstants.o
MassNucleus.o: +PhysMathConstants.o
dsQESCC_dQ2_fN.o: +PhysMathConstants.o
QESkin.o: +PhysMathConstants.o
dsQESCC_dQ2.o: +PhysMathConstants.o
dsQESCC_dQ2_FP.o: +PhysMathConstants.o
setEds.o: +PhysMathConstants.o
qesemat: +PhysMathConstants.o +InpOutUnits.o \
GeM.o MuL.o spline1.o \
DZEROX.o LambdaFUNCTION.o DMINFC.o \
FunGeM_SM.o FunMuL_SM.o QESkin_SM.o QESkin.o NucQESFF.o rho_SM.o MassNucleus.o FactorPauli.o \
d3sQES_dQ2dnudkF_SM.o dsQESCC_dQ2_SM.o dsQESCC_dQ2_fN.o dsQESCC_dQ2_FP.o dsQESCC_dQ2.o \
setEds.o MA_QES_EFF.o fui.o InitFlux.o

testflux: InitFlux.o spline1.o
test_fluxtable: InitFlux.o spline1.o

BIN=bin/
SRC=src/
# IF YOU MOVE NEXT 2 LINES HIGHER THAN $^ DESCRIPTION nothing will be wrong. It seems...
SOURCES=$(addprefix $(SRC),$(^:.o=.for))
OBJECTS=$(addprefix $(BIN),$^)
GARBAGE=$(foreach dir,$(BIN) ./,$(wildcard $(dir)*.o $(dir)*.mod  $(dir)*.MOD $(dir)fort.*))

# General rule for building prog from prog.o; $^ (GNU extension) is used in order to list additional object files on which the executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is used in order to list only the first prerequisite (the source file) and not the additional prerequisites such as module or include files
%.o: $(SRC)%.for
	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@

# Utility targets
clean:
	$(RM) $(GARBAGE)

# DO NOT DELETE THIS LINE! - make depend depends on it
