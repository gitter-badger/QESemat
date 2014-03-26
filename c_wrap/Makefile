# FC = gfortran
# FCFLAGS = -Wall -g
# FCFLAGS+= -ffpe-trap=invalid,zero,overflow

CC = g++
CFLAGS=-Wall -g #-std=c++11
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


