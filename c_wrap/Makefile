# FC = gfortran
# FCFLAGS = -Wall -g
# FCFLAGS+= -ffpe-trap=invalid,zero,overflow

CC = g++
CFLAGS=-Wall -g 
#CFLAGS-std=c++11
CFLAGS+=-O3 
LDFLAGS= -L ../lib/ -lQES_flux -lgfortran
TARGET=Flux

all: $(TARGET)

%: %.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

%.o: $(SRC)%.cxx
	$(CC) -c $(CFLAGS) $< -o $@ $(filter *.o, $^)

# %.f.o: $(SRC)%.f90
# 	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@

clean:
	rm *.o $(TARGET)


