!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EXPERIMENT 1:   Comparison with the modified Magnus method 
!
!  We compare our method with the modified Magnus method of 
!
!    A. Iserles, “On the global error of discretization methods for highly-oscillatory 
!      ordinary differential equations,” BIT 42 (2002), 561-599
!
!  We do so by evaluating Legendre functions on the interval [0,0.9] and measuring the 
!  accuracy of the methods and the time taken for values of the degree n between 2^6
!  and 2^20.
!
!  We use an equispaced discretization grid for the modified Magnus method and allow
!  the stepsize to scale with frequency.  We execute the method twice; in the first
!  run, the stepsize scales as dnu^(3/4) and in the second we let it scale as  dnu^(1/2).  
!  The stepsize for the lowest frequency considered was set so that approximately 10^(-10) 
!  accuracy was achieved.  Lowering the stepsize too much leads to numerical instability,
!  hence our choice to aim for 10 digit accuracy.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment1_funs

use experutils
use freqind
use magnus

type(freqind_vars_t)           :: vars
integer, parameter             :: nruns  = 100
integer, parameter             :: ncheck = 1000
double precision               :: b 

contains


subroutine magnus_legendre(n,nsteps,time,errmax)
implicit double precision (a-h,o-z)
!
!  Calculate the solution
!
!    P_n(x) +  2/pi Q_n(x)
!
!  over the interval [0,0.9] by applying a modified Magnus method to Legendre's
!  differential equation.  The modified Magnus method uses equispaced steps
!  with the number of steps specified as an input parameter.
!
!  The subroutine returns the time taken as well as the maximum observed 
!  relative error.
!
!  Input parameters:
!    n - the order of the Legendre functions
!    nsteps - the number of steps to take
!
!  Output parameters:
!    time - the time in milliseconds
!    errmax - the maximum observed RELATIVE error
!

double precision, allocatable    :: ts(:)
double complex, allocatable      :: ys(:,:)
double complex                   :: y0(2), ima
double complex, pointer          :: pars(:)
double complex                   :: val0, val

data pi / 3.14159265358979323846264338327950288d0 /


ima = (0.0d0,1.0d0)

a   = 0.000d0

allocate(pars(1))
allocate(ys(2,nsteps+1),ts(nsteps+1))
pars(1) = n

call legeder(n,a,pval,pder)
call legeqder(n,a,qval,qder)
qval = qval*2/pi
qder = qder*2/pi

y0(1) = pval + ima*qval
y0(2) = pder + ima*qder

call elapsed(t1)
do irun=1,nruns
call  magnus_modified(a,b,nsteps,y0,ts,ys,gfun_lege,pars)
end do
call elapsed(t2)
time = (t2-t1)/nruns

errmax = 0

if (nsteps .lt. ncheck) then

  do ii=1,nsteps+1
  i    = ii
  t    = ts(i)
  val  = ys(1,i)/sqrt(1-t**2)
   
  call legeder(n,t,pval,pder)
  call legeqder(n,t,qval,qder)
  qval = qval*2/pi
  qder = qder*2/pi
   
  val0   = (pval+ima*qval)
  derr   = abs(val-val0)/abs(val0)
  errmax = max(derr,errmax)   
  end  do

else

  do ii=1,ncheck
  i    = 1 + (nsteps-1)*(ii-1.0d0)/(ncheck-1.0d0)
  t    = ts(i)
  
  val  = ys(1,i)/sqrt(1-t**2)
   
  call legeder(n,t,pval,pder)
  call legeqder(n,t,qval,qder)
  qval = qval*2/pi
  qder = qder*2/pi
   
  val0   = (pval+ima*qval)
  derr   = abs(val-val0)/abs(val0)
  errmax = max(derr,errmax)   
  end  do

endif

end subroutine

subroutine gfun_lege(t,val,der,pars)
implicit double precision (a-h,o-z)
double complex, pointer         :: pars(:)
dnu    = pars(1)
val    = 1/(1-t**2)**2 + dnu*(dnu+1)/(1-t**2)
der    = (0.4d1*t)/(1-t**2)**3+(0.2d1*dnu*(1+dnu)*t)/(1-t**2)**2
end subroutine


subroutine phase_legendre(n,time,errmax,condmax)
implicit double precision (a-h,o-z)
!
!  Calculate the solution
!
!    P_n(x) +  2/pi Q_n(x)
!
!  over the interval [0,0.9] by applying a phase function method to Legendre's
!  differential equation.
!
!  The subroutine returns the time taken as well as the maximum observed 
!  relative error.
!
!  Input parameters:
!    n - the order of the Legendre functions
!    nsteps - the number of steps to take
!
!  Output parameters:
!    time - the time in milliseconds
!    errmax - the maximum observed RELATIVE error

double precision, allocatable   :: ab(:,:)
double precision, allocatable   :: alpha(:,:)
double precision, allocatable   :: alphader(:,:)
double precision, allocatable   :: alphader2(:,:)
double precision, pointer       :: pars(:)
double complex                  :: ima, val, val0

data pi / 3.14159265358979323846264338327950288d0 /


k      = 16
ifleft = 1
aval0  = -pi/2*(n+0)
ima   = (0.0d0,1.0d0)
dnu   = n
eps0  = epsilon(0.0d0)

a     = 0.000d0


allocate(pars(1))
pars(1) = dnu

call elapsed(t1)
do irun=1,nruns
! call  freqind_positive(vars,ier,a,b,ifleft,aval0,qfun_lege,pars,nints,ab, &
!   alpha,alphader,alphader2)

call  freqind_phase(vars,ier,a,b,ifleft,aval0,qfun_lege,pars,nints,ab, &
  alpha,alphader,alphader2)

end do
call elapsed(t2)
time = (t2-t1)/nruns

if (ier .ne. 0) then
print *,"after freqind_positive, ier = ",ier
stop
endif

errmax  = 0
condmax = 0

do i=1,ncheck
t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
call freqind_eval(vars,nints,ab,alpha,alphader,alphader2,t,aval,apval,appval)
val    = exp(ima*aval)/sqrt(apval)*1/sqrt(1-t**2)*sqrt(2/pi)
cond   = abs(apval-appval/(2*apval))
cond   = abs(t)*cond*eps0

call legeder(n,t,pval,pder)
call legeqder(n,t,qval,qder)

qval   = qval*2/pi
qder   = qder*2/pi
val0   = pval+ima*qval

derr    = abs(val-val0)/abs(val0)
errmax  = max(derr,errmax)
condmax = max(condmax,cond)

end do


end subroutine


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

program experiment1


use experutils
use freqind
use magnus
use experiment1_funs

implicit double precision (a-h,o-z)

double precision, allocatable :: times(:,:), errs(:,:), dnus(:)

epsdisc = 1.0d-12
epssol  = 1.0d-12
k       = 16
call freqind_init(vars,k,epsdisc,epssol)

b       = 0.9d0

j1      = 6
j2      = 20

nsteps0 = 4150          
nn      = j2-j1+1

allocate(times(4,j1:j2),errs(4,j1:j2),dnus(j1:j2))

do j=j1,j2
dnu     = 2**j
n       = dnu
nsteps1 = nsteps0*(dnu/2**j1)**(0.5d0)
nsteps2 = nsteps0*(dnu/2**j1)**(0.75d0)

call magnus_legendre(n,nsteps1,time1,errmax1)
call magnus_legendre(n,nsteps2,time2,errmax2)
call phase_legendre(n,time3,errmax3,condmax)

dnus(j)    = dnu

times(1,j) = time1
times(2,j) = time2
times(3,j) = time3

errs(1,j)  = errmax1
errs(2,j)  = errmax2
errs(3,j)  = errmax3
errs(4,j)  = condmax

write (*,"(I2.2,2X,I7.7,4(D10.2),4(D15.6))") j,n,errmax1,errmax2,errmax3,time1,time2,time3,time2/time3

end do

call plot_functions2("experiment1-times.py","experiment1-times.pdf",3,"$n$","Time (seconds)",                      &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-7,10.0d0,                                                                         &
   nn,dnus,times(1:3,:),"best","Magnus $h=O(n^{-0.50})$*Magnus $h=O(n^{-0.75})$*Phase*Phase 10 threads*", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

call plot_functions2("experiment1-errs.py","experiment1-errs.pdf",4,"$n$","Max Relative Error", &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-16,1.0d0,                                                 &
   nn,dnus,errs,"best","Magnus $h=O(-n^{0.50})$*Magnus $h=O(-n^{0.75})$*Phase*Cond #*",     &
  "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

end program
