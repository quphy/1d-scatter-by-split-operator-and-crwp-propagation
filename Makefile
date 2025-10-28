FC = gfortran
FFLAGS = -O3
LIBS = -llapack -lblas

TARGET = scartt
SRC = scatter.f90 main.f90

$(TARGET): $(SRC)
		$(FC) $(FFLAGS) $(SRC) $(LIBS) -o $(TARGET)

clean:
		rm -f $(TARGET) *.o *.mod
