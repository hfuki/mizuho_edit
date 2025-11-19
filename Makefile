TARGET = wPFAW3D

## for intel oneAPI
FC = mpiifx
#FFLAGS = -O0 -check bounds -traceback -fpp
FFLAGS = -O2 -fpp -g -O0
IFLAGS = -I $(HOME)/local/include
LFLAGS = -L $(HOME)/local/lib
LIBS = -lpetsc
##[note] remeber to add $(HOME)/local/lib to LD_LIBRARY_PATH.

## for GNU compiler
#FC = mpif90
#FFLAGS = -O2 -cpp -fallow-argument-mismatch
#IFLAGS = -I $(HOME)/local/include
##FFLAGS = -O0 -cpp -ggdb -fbounds-check -fbacktrace -Wuninitialized -ffpe-trap=invalid,zero,overflow -fallow-argument-mismatch
#LFLAGS = -L $(HOME)/local/lib
#LIBS   = -lpetsc
##[note] remeber to add $(HOME)/local/lib to PATH.

##----------------------------------------------------------------------------
OBJS = wPFAW3D.o parallel.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(OBJS) $(LFLAGS) $(LIBS) -o $@

%.o: %.f90
	$(FC) -c $(FFLAGS) $(IFLAGS) $< -o $@

wPFAW3D.o: parallel.o

clean:
	rm -f *.o *.mod *~ $(TARGET)
