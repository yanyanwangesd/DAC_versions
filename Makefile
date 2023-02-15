.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.build

OBJECTS = \
module_definitions.o \
add_and_remove_nodes_2.o \
DivideAndCapture.o \
floader.o \
initialize_geometry.o \
reinit_geometry.o \
calculate_delaunay.o \
calculate_centroid.o\
find_network.o \
find_order.o \
find_discharge.o \
find_surface.o \
find_polygon_surface.o \
find_hack.o \
find_strahler.o \
timer.o \
erode.o \
initialize_parameters.o \
uplift_and_advect.o \
captures_and_divides.o \
find_precipitation.o \
find_catchment.o \
output_geometry.o \
output_z_tau.o \
write_ascii.o \
update_erosion_rate_history.o \
finetri.o \
VTKfine.o \
VTK.o \
reboundary.o \
clean_dacmem.o \
set_params_as_function.o 
#csigfun.o \



# These are always built
.PHONY: all clean distclean resclean rundir

# Default. Builds everything if outdated
all: NN/libnn.a DivideAndCapture rundir

NN/libnn.a:
	$(MAKE) -C NN

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $*.f

.c.o:
	$(CC) $(CFLAGS) $*.

#.ONESHELL:
rundir: 
	if [ ! -d "RUN5" ]; then mkdir RUN5; fi 
	if [ ! -d "RESTART" ]; then mkdir RESTART; fi 
	if [ ! -d "ASCII" ]; then mkdir ASCII; fi 

DivideAndCapture:	$(OBJECTS)
	$(F90)  $(OPTIONS) $(OBJECTS) $(LIBS) -o DAC 


clean:
	-rm *.o *.mod 

# Cleans everything including libraries, executable, and results
distclean:
	-rm -f *.o *.mod *.vtk
	-rm -f DAC
	-rm -f massbalance.txt massbalance_ib.txt
	-cd NN; rm -f *.o *.a; cd -
	-cd RUN5; rm -f *.vtk ; cd -
	-cd RESTART; rm -f GEOFILE* ;cd -
	-cd ASCII; rm -f * ; cd -
	
	
# Cleans results only	
resclean: 
	-cd RUN5; rm -f *.vtk ; cd -
	-cd RESTART; rm -f GeoFile* ;cd -
	-cd ASCII; rm -f * ;cd -
