#FC = ifort -lmpi
#CC = icc -lmpi
#CXX = icpc -lmpi
#FFLAGS += -O3 -DMPI
#CFLAGS += -O3 -DMPI

FC = gfortran
CC = gcc
CXX = g++
FFLAGS += -ffree-line-length-none -O3 #-DMPI
CFLAGS += -O3 #-DMPI

LAPACKLIB = -llapack

NESTLIBDIR = ./

export FC CC CXX FFLAGS CFLAGS LAPACKLIB

 
AR = ar r  
LINKLIB = ld -shared
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

default: libnest3.a

all: libnest3.a obj_detect eggbox eggboxC eggboxC++ gaussian gauss_shell \
rosenbrock himmelblau ackley
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
 
obj_detect:
	make -C example_obj_detect
 
gaussian:
	make -C example_gaussian
 
rosenbrock:
	make -C example_rosenbrock
 
ackley:
	make -C example_ackley
 
himmelblau:
	make -C example_himmelblau
 
eggbox:
	make -C example_eggbox
 
gauss_shell:
	make -C example_gauss_shell
	
eggboxC:
	make -C example_eggbox_C
	
eggboxC++:
	make -C example_eggbox_C++

clean: 
	-rm $(NESTLIBDIR)/libnest3.*  *.o *.mod
	
cleanall: clean_exec clean clean_obj_detect clean_gaussian clean_gauss_shell clean_eggbox \
clean_example_eggbox_C clean_example_eggbox_C++ clean_rosenbrock clean_himmelblau \
clean_ackley

clean_exec:
	-rm obj_detect gaussian rosenbrock ackley himmelblau gauss_shell eggbox eggboxC eggboxC++

clean_obj_detect:
	make -C example_obj_detect clean
	
clean_gaussian:
	make -C example_gaussian clean
	
clean_rosenbrock:
	make -C example_rosenbrock clean
	
clean_ackley:
	make -C example_ackley clean
	
clean_himmelblau:
	make -C example_himmelblau clean
	
clean_eggbox:
	make -C example_eggbox clean
	
clean_gauss_shell:
	make -C example_gauss_shell clean
	
clean_example_eggbox_C:
	make -C example_eggbox_C clean
	
clean_example_eggbox_C++:
	make -C example_eggbox_C++ clean
