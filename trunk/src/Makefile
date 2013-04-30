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
F90FLAGS = -O3 -fno-strict-overflow
#---------------------
# Intel 
#F90FLAGS = -O3 -assume byterecl
# Intel (EM64T/AMD64 processor)
#F90FLAGS = -O3 -funroll-loops -xP
#---------------------

PROGS = specfit

COMMON =  nrtype.o sfvars.o nr.o nrutil.o ran_state.o sfutils.o getmodel.o \
	sfsetup.o linterp.o velbroad.o tsum.o set_wave_limits.o \
	contnormspec.o locate.o gasdev.o setup_params.o str2arr.o \
	ran1.o read_data.o lfit.o gaussj.o covsrt.o getm2l.o myran.o \
	powell.o linmin.o mnbrak.o brent.o func.o init_random_seed.o \
	imf.o getmass.o qromb.o polint.o trapzd.o maskemlines.o getvelz.o


all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

specfit : specfit.o $(COMMON) 
	$(F90) -o specfit.exe specfit.o $(COMMON)

specfit.o : sfvars.o nrtype.o nr.o sfutils.o ran_state.o nrutil.o

sps_utils.o: sfvars.o

ran_state.o: nrutil.o

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<

%.o : %.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


