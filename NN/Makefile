.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include ../Makefile.build

OBJECTS_NN = \
delaun.o \
del_sub.o \
nn.o \
nnplot.o \
qhullf_dummy.o \
stack.o \
stackpair.o \
volume.o

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $*.f

.c.o:
	$(CC) $(CFLAGS) $*.c

nn:	$(OBJECTS_NN)
	ar rcv libnn.a $(OBJECTS_NN)

clean:
	rm *.o

distclean:
	-rm -f *.o 
	-rm -f libnn.a
