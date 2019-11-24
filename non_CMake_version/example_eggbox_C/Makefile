NESTLIBDIR = ../
LIBS =  -L$(NESTLIBDIR) -lnest3 $(LAPACKLIB)
 
OBJFILES = eggbox.o

all: eggbox 

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c
 
eggbox: $(OBJFILES)
	$(FC) $(FFLAGS) -nofor_main -o ../eggboxC $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod ../eggboxC
