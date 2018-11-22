FC = $(F90C)

NESTLIBDIR = ./
CAMBDIR = ../camb
 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o symmetricbetadistribution.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -I$(CAMBDIR) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -I$(CAMBDIR) -c -o $@ $^ 

default: libnest3.a

all: libnest3.a 
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean

