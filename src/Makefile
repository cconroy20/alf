#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
F90 = gfortran
#---------------------
# Intel
#F90 = ifort 
#---------------------

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU
F90FLAGS = -O -fno-strict-overflow #-fopenmp
#---------------------
# Intel 
#F90FLAGS = -O3 -assume byterecl
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP
#---------------------

PROGS = alf write_a_model

COMMON =  nrtype.o sfvars.o nr.o nrutil.o ran_state.o sfutils.o getmodel.o \
	sfsetup.o linterp.o velbroad.o tsum.o add_response.o \
	contnormspec.o locate.o gasdev.o setup_params.o str2arr.o \
	ran1.o read_data.o lfit.o gaussj.o covsrt.o getm2l.o myran.o \
	powell.o linmin.o mnbrak.o brent.o func.o init_random_seed.o \
	getmass.o maskemlines.o getvelz.o emcee_advance.o update_runtot.o \
	linterp3.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

alf : alf.o $(COMMON) 
	$(F90) $(F90FLAGS) -o alf.exe alf.o $(COMMON)

write_a_model : write_a_model.o $(COMMON) 
	$(F90) -o write_a_model.exe write_a_model.o $(COMMON)

alf.o : sfvars.o nrtype.o nr.o sfutils.o ran_state.o nrutil.o

write_a_model.o : sfvars.o nrtype.o nr.o sfutils.o ran_state.o nrutil.o

ran_state.o: nrutil.o

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


