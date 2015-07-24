#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
#F90 = gfortran
#---------------------
# Intel
#F90 = ifort 
#---------------------
#MPI
F90 = mpifort 

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU
F90FLAGS = -O -Wall
#---------------------
# Intel 
#F90FLAGS = -O3 -assume byterecl
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP
#---------------------

PROGS = alf write_a_model

COMMON =  nrtype.o alf_vars.o nr.o nrutil.o ran_state.o alf_utils.o getmodel.o \
	setup.o linterp.o velbroad.o tsum.o add_response.o \
	contnormspec.o set_pinit_priors.o str2arr.o locate.o gasdev.o \
	ran1.o lfit.o gaussj.o covsrt.o powell.o linmin.o \
	mnbrak.o brent.o getm2l.o read_data.o func.o init_random_seed.o \
	getmass.o maskemlines.o getvelz.o emcee_advance_mpi.o update_runtot.o \
	linterp3.o free_workers.o function_parallel_map.o  myran.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

alf : alf.o $(COMMON) 
	$(F90) $(F90FLAGS) -o alf.exe alf.o $(COMMON)

write_a_model : write_a_model.o $(COMMON) 
	$(F90) -o write_a_model.exe write_a_model.o $(COMMON)

alf.o : alf_vars.o nrtype.o nr.o alf_utils.o ran_state.o nrutil.o

write_a_model.o : alf_vars.o nrtype.o nr.o alf_utils.o ran_state.o nrutil.o

ran_state.o: nrutil.o

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


