F90 = gfortran
OBJ = table_1d_module.o singlestar_module.o testtable.o
LIBS = 
FFLAGS = -frecord-marker=4 -O1 -ffree-line-length-none -x f95-cpp-input
VPATH = ../

testtable: $(OBJ)
	$(F90) $(OBJ) -o testtable $(LIBS)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
clean:	
	rm *.mod *.o
