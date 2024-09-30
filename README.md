This package contains a frequency-independent solver for oscillatory second order
linear ordinary differential equations. It is an implementation of the algorithm 
described in the following preprint:

   Tara Stojimirovic and James Bremer, ``An accelerated frequency-indepedent solver 
   for oscillatory differential equations.''  arXiv:2409.18487

The code is written in Fortran and the implementation of our algorithm, which can be 
found in the file freqind.f90, has as its only dependency the file chebyshev.f90
which was written by the authors and is included in this package.  Test code
is provided in the aptly-named file ``test_freqind.f90.''  It can be compiled
and executed via the command:

   ``make test_freqind''

The code for all but one of the numerical experiments described in the preprint mentioned 
above is contained files experiment?.f90.  Most of these have no external dependences.  
However, the code in experiment1.f90 depends on LAPACK, BLAS and the EXPOKIT package for 
performing matrix exponentiation.  However, each experiment products python
scripts for which generate plots, and these scripts depend on numpy and mathplotlib.

The command ``make'' will cause all of the experiments to be complied and executed.  
The python scripts can then be executed to generate the plots in the preprint.
