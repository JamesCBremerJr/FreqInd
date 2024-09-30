!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EXPERIMENT 3:  Compute the accuracy with which the phase functions are computed
!                 by measuring the accuracy on the derivative of a trigonometric
!                 phase function for Legendre's differential equation on the interval
!                 [0,1-10^(-7)].
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment3_funs

use experutils
use freqind

type(freqind_vars_t)           :: vars
integer, parameter             :: nruns  = 100
integer                        :: ncheck 

contains


subroutine qfun_lege(n,ts,qs,qders,pars)
implicit double precision (a-h,o-z)
double precision           :: ts(n)
double precision           :: qs(n)
double precision           :: qders(n)
double precision, pointer  :: pars(:)

dnu   = pars(1)
qs    = (dnu)*(dnu+1)/(1-ts**2) + 1/(1-ts**2)**2
qders = (dnu)*(dnu+1)*2*ts/(1-ts**2)**2  + 4*ts/(1-ts**2)**3

end subroutine

end module

program experiment3


use experutils
use freqind
use experiment3_funs

implicit double precision (a-h,o-z)
double precision, allocatable          :: times(:,:), errs(:,:), dnus(:)

double precision, allocatable          :: ab(:,:)
double precision, allocatable          :: alpha(:,:)
double precision, allocatable          :: alphader(:,:)
double precision, allocatable          :: alphader2(:,:)
double precision, pointer              :: pars(:)

pi = acos(-1.0d0)

epsdisc = 1.0d-12
epssol  = 1.0d-12
k       = 16
call freqind_init(vars,k,epsdisc,epssol)


a      = 0.0000d0
b      = 1.0d0-1.0d-7

ifleft = 1
aval0  = 0
ncheck = 1000
jj1    = 7
jj2    = 21


nn = jj2-jj1+1

allocate(times(1,jj1:jj2),errs(1,jj1:jj2),dnus(jj1:jj2))
allocate(pars(1))

do jj=jj1,jj2

dnu     = 2**jj
n       = dnu
pars(1) = dnu


call elapsed(t1)
do irun=1,nruns
call freqind_phase(vars,ier,a,b,ifleft,aval0,qfun_lege,pars,nints,ab, &
  alpha,alphader,alphader2)
end do
call elapsed(t2)
time1 = (t2-t1)/nruns

errmax=0

do i=1,ncheck-1

call legeder(n,t,pval,pder)
call legeqder(n,t,qval,qder)

pval = pval * sqrt(1-t**2)*sqrt(pi/2)
qval = qval * sqrt(1-t**2)*sqrt(2/pi)
apval0 = 1/(pval**2+qval**2)


call freqind_eval(vars,nints,ab,alphader,t,apval)
errrel = abs(apval-apval0)/abs(apval0)
errmax = max(errrel,errmax)
end do


dnus(jj)    = dnu
errs(1,jj)  = errmax
times(1,jj) = time1

print *,dnu,time1,errmax

end do

close(iw)


call plot_functions2("experiment3-times.py","experiment3-times.pdf",1,"$n$","Time (seconds)",                      &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-7,10.0d0,                                                                         &
   nn,dnus,times,"","", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

call plot_functions2("experiment3-errs.py","experiment3-errs.pdf",1,"$n$","Max Relative Error", &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-16,1.0d0,                                                 &
   nn,dnus,errs,"","", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")


end program
