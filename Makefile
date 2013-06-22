#vortex/Makefile

#source files
SOURCES =dapif1.f dapif2.f dapif3.f dcfft.f matplot.f \
	 prini.f pvmulti.f vortex_motion.f
#object files for fmm testing  
OBJECTS = $(SOURCES:.f=.o)

#Current Directory
FMMDIREC =/home/visitor/nchan/fortran/vortex/fmmlib2d-1.2/src 
DIREC=/home/visitor/nchan/fortran/vortex/
GMRESDIREC=/home/visitor/nchan/fortran/vortex/GMRES
#libraries required
#libfmm2d.a is the compiled library of the latest fmm code.
#GMRES has the required files for the GMRES method 
LIBFMM =   libhfmm2d.a 
LIBGMRES = libgmres.a
LIB = $(LIBFMM) $(LIBGMRES)
#Set your fortran compiler here
FC  = gfortran

#Flags
FFLAGS = -fno-automatic -std=legacy
FLAGS1 = -c
FLAGS2 = -o

#type `make clean` to rm the object files
.PHONY: clean help 

lib:
	cd $(FMMDIREC); $(MAKE); \
	cp $(LIBFMM) $(DIREC)
	
	cd $(GMRESDIREC) ; \
	$(FC) $(FLAGS1) *.f ; \
	ar rvs $(LIBGMRES) *.o ; \
	rm -f *.o ; \
	cp $(LIBGMRES) $(DIREC)	
	
		

test: test.exe
	./test.exe

test.exe: $(OBJECTS)
	$(FC) $(FLAGS2)  $@ $(OBJECTS) $(LIB)

$(OBJECTS): $(SOURCES)
	$(FC) $(FLAGS1) $(SOURCES)

clean: 
	rm -f *.exe $(OBJECTS) text.exe

help:
	@echo make lib --- To get compiled libraries
	@echo make test --- To solve the vortex problem
	 
