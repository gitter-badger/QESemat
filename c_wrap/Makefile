# FC = gfortran
# FCFLAGS = -Wall -g
# FCFLAGS+= -ffpe-trap=invalid,zero,overflow

CC = g++
# CFLAGS=-Wall -std=c++11 -g
CFLAGS+=-O3 
LDFLAGS= -L../ -lFlux -lgfortran
TARGET=Flux

all: $(TARGET)


%: %.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

%.o: $(SRC)%.cxx
	$(CC) -c $(CFLAGS) $< -o $@ $(filter *.o, $^)

# %.f.o: $(SRC)%.f90
# 	$(FC) $(FCFLAGS) -c $< -o $(BIN)$@

clean:
	rm *.o $(TARGET)


