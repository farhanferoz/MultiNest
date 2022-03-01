NESTLIBDIR = ../
LIBS = -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
OBJFILES = params.o like.o nestwrap.o main.o

all: himmelblau 

%.o: %.f90
	$(FC) $(FFLAGS) -I$(NESTLIBDIR) -c $*.f90
 
himmelblau: $(OBJFILES)
	$(FC) -o ../himmelblau $(OBJFILES) \
	$(FFLAGS) $(LIBS)

clean:
	rm -f *.o *.mod ../himmelblau
