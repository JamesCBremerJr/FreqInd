FCOMP     =  gfortran
FOPTS     =  -Ofast -march=native -w -fallow-argument-mismatch $(PREC)  
#FOPTS    +=  -freal-8-real-10 -fdefault-integer-8
#FOPTS    +=  -fdefault-real-8 -fdefault-integer-8

# specify BLAS and LAPACK library
LDOPT     =  -lblas -llapack


# Set the list of programs to compile

PROGRAMS = test_freqind experiment1 experiment3 experiment4 experiment5 experiment6

# Compile all of the test programs and the library
all	      	             : clean $(PROGRAMS) 


# List the dependencies for each module's test program

EXPERIMENT1_FILES              = experutils.o                                             \
                                 chebyshev.o                                              \
	                         expokit.o                                                \
                                 magnus.o                                                 \
                                 freqind.o

EXPERIMENT3_FILES              = experutils.o                                             \
                                 chebyshev.o                                              \
                                 freqind.o

EXPERIMENT4_FILES              = experutils.o                                             \
                                 chebyshev.o                                              \
                                 freqind.o

EXPERIMENT5_FILES              = experutils.o                                             \
                                 chebyshev.o                                              \
                                 freqind.o                                                \
                                 odetwo.o

EXPERIMENT6_FILES              = experutils.o                                             \
                                 chebyshev.o                                              \
                                 freqind.o                                                \
                                 odetwo.o

FREQIND_FILES                  = experutils.o                                             \
                                 chebyshev.o                                              \
                                 freqind.o


experiment1.o                  : $(EXPERIMENT1_FILES) experiment1.f90
experiment1                    : $(EXPERIMENT1_FILES) experiment1.o

experiment3.o                  : $(EXPERIMENT3_FILES) experiment3.f90
experiment3                    : $(EXPERIMENT3_FILES) experiment3.o

experiment4.o                  : $(EXPERIMENT4_FILES) experiment4.f90
experiment4                    : $(EXPERIMENT4_FILES) experiment4.o

experiment5.o                  : $(EXPERIMENT5_FILES) experiment5.f90
experiment5                    : $(EXPERIMENT5_FILES) experiment5.o

experiment6.o                  : $(EXPERIMENT6_FILES) experiment6.f90
experiment6                    : $(EXPERIMENT6_FILES) experiment6.o

test_freqind.o                 : $(FREQIND_FILES) test_freqind.f90
test_freqind                   : $(FREQIND_FILES) test_freqind.o


# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 


%.o		: %.cpp
	$(CPP) -c $(CPPOPTS)  $<

%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<


%.o		: %.c
	$(CC) -c $(COPTS)  $<


clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.py 
	rm -f *.pdf


