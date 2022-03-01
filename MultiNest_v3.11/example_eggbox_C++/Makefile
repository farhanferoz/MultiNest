NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB) -lstdc++
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.cc
	$(CXX) $(CFLAGS) -c $*.cc
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -nofor_main -o ../eggboxC++ $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC++
