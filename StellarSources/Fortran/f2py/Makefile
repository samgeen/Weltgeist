#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@cea.fr
#############################################################################
# Compilation time parameters
EXEC = singlestar
#############################################################################
DEFINES = 
FFLAGS = -x f95-cpp-input $(DEFINES) -frecord-marker=4 -fbacktrace -ffree-line-length-none
F90 = gfortran -O2 -fPIC $(FFLAGS)
FC = gcc -O2 -fPIC
#############################################################################
LIBS = -lgfortran
#############################################################################
# Sources directories are searched in this exact order
VPATH = ../
#############################################################################
# All objects
MODOBJ = table_1d_module.o singlestar_module.o
GLUEOBJ = singlestar_f2py.o
ALLLIB = $(GLUEOBJ)
#############################################################################
singlestar:	$(MODOBJ) singlestar_f2py.o
	f2py -m $(EXEC) ../singlestar_f2py.f90 -h $(EXEC).pyf --overwrite-signature
	f2py -lgfortran -c *.o $(EXEC).pyf
objects: $(MODOBJ) $(ALLLIB)
#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#############################################################################
%.o:%.c
	$(FC) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.mod
#############################################################################
