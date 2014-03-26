# FC = gfortran
# FCFLAGS = -Wall -g
# FCFLAGS+= -ffpe-trap=invalid,zero,overflow

CC = g++
<<<<<<< HEAD
CFLAGS=-Wall -g #-std=c++11
=======
# CFLAGS=-Wall -std=c++11 -g
>>>>>>> 451dfb5190928a8688f201b6a119f66385bc8907
CFLAGS+=-O3 
LDFLAGS= -L ../lib/ -lQES_flux -lgfortran
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


