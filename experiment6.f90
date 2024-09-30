!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EXPERIMENT 6:  Solve a boundary value problem.
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment6_funs

use experutils
use freqind
use chebyshev
use odetwo

type(odetwo_vars_t)            :: odetwo_vars
type(freqind_vars_t)           :: vars
integer, parameter             :: nruns  = 100
integer, parameter             :: ncheck = 1000

double precision               :: omega0
double precision               :: dd
double precision               :: dd2

contains


subroutine phase(omega,domega,time,timeode,errmax)
implicit double precision (a-h,o-z)
!
!
!  Input parameters:
!    nthreads - the number of threads to use
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
double precision, allocatable   :: vals(:), vals0(:), ts(:)

double precision              :: amatr(2,2),cs(2),rhs(2)
double precision, allocatable :: ab1(:,:)
double complex, allocatable   :: y1s(:,:)
double complex, allocatable   :: y1ders(:,:)


double complex                :: yc, ypc, val1, val2, der1, der2, val0, val
double complex                :: bmatr(2,2), detz, c1,c2
type(c_ptr)                   :: userptr


data pi / 3.14159265358979323846264338327950288d0 /
eps0     = epsilon(0.0d0)


dnu      = n
a        =  0.0d0
b        =  1.0d0
aval0    = 0
ifleft   = 1


omega0   = omega
dd       = domega


allocate(vals(ncheck),vals0(ncheck),ts(ncheck))


call elapsed(t1)
do irun=1,nruns

call freqind_phase(vars,ier,a,b,ifleft,aval0,qfun,pars,&
  nints,ab,alpha,alphader,alphader2)


if (ier .ne. 0 ) then
print *,"!!!!!!!!!!!!!! ier = ",ier
stop
endif


c   = 0.0d0
val = 1.0d0
der = 0.0d0

call freqind_eval(vars,nints,ab,alpha,alphader,c,aval,apval)
amatr(1,1) = sin(aval)/sqrt(apval)
amatr(1,2) = cos(aval)/sqrt(apval)

call freqind_eval(vars,nints,ab,alpha,alphader,alphader2,c,aval,apval,appval)
amatr(2,1) = cos(aval)*sqrt(apval) - 0.5d0 * sin(aval)*appval/apval**1.5d0
amatr(2,2) =-sin(aval)*sqrt(apval) - 0.5d0 * cos(aval)*appval/apval**1.5d0

det = amatr(1,1)*amatr(2,2) - amatr(1,2)*amatr(2,1)
d1  = ( amatr(2,2) * val - amatr(1,2) * der ) / det
d2  = (-amatr(2,1) * val + amatr(1,1) * der ) / det

end do

eps = 1.0d-14
call elapsed(t2)
time = (t2-t1)/nruns

omega0 = omega

call elapsed(t1)
yc         = val
ypc        = der
call odetwo_nonlinear1(odetwo_vars,ier1,eps,a,b,c,yc,ypc,nints1,ab1,y1s,y1ders,odefun,userptr)

errmax = 0

call elapsed(t2)
timeode = t2-t1

do i=1,ncheck

t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
call odetwo_eval(odetwo_vars,nints1,ab1,y1s,t,val0)

call freqind_eval(vars,nints,ab,alpha,alphader,alphader2,t,aval,apval,appval)
call freqind_eval(vars,nints,ab,alpha,alphader,t,aval,apval)

val = sin(aval)/sqrt(apval)*d1 + cos(aval)/sqrt(apval)*d2

ts(i)    = t
vals0(i) = val0


errabs = abs(val0-val)
errmax = max(errabs,errmax)

end do



end subroutine



subroutine odefun(t,y,yp,f,dfdy,dfdyp,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
type(c_ptr)                 :: userptr
double precision            :: ts(1),qs(1),qders(1)
double precision, pointer   :: pars(:)

ts(1) = t
call qfun(1,ts,qs,qders,pars)

f     = -qs(1)*y
dfdy  = -qs(1)
dfdyp = 0

end subroutine


subroutine qfun(n,ts,qs,qders,pars)
implicit double precision (a-h,o-z)
double precision           :: ts(n)
double precision           :: qs(n)
double precision           :: qders(n)
double precision, pointer  :: pars(:)

omega   = omega0

qs    = omega**2 * (1000+ts**2   + sin(omega**dd*ts)/omega**dd2                      )
qders = omega**2 * (2*ts         + cos(omega**dd*ts)/omega**dd2*omega**dd            )

end subroutine



end module

program experiment6

use experutils
use freqind
use experiment6_funs

implicit double precision (a-h,o-z)

double precision, allocatable :: times(:,:), errs(:,:), omegas(:)


call odetwo_init(odetwo_vars)

epsdisc = 1.0d-14
epssol  = 1.0d-14
k       = 16
call freqind_init(vars,k,epsdisc,epssol)

jj1 = 6
jj2 = 20
nn  = jj2-jj1+1

allocate(errs(3,jj1:jj2),times(3,jj1:jj2), omegas(jj1:jj2))

do jj=jj1,jj2
omega   = 2**jj

dd2     = 3.0d0

call phase(omega,0.50d0,time1,time0,errmax1)
call phase(omega,1.00d0,time2,time0,errmax2)
call phase(omega,1.50d0,time3,time0,errmax3)

omegas(jj)  = omega
times(1,jj) = time1
times(2,jj) = time2
times(3,jj) = time3

errs(1,jj)  = errmax1
errs(2,jj)  = errmax2
errs(3,jj)  = errmax3

print "(I5.5,6(D13.5,1X))",jj,time1,time2,time3,errmax1,errmax2,errmax3

end do


call plot_functions2("experiment6-times.py","experiment6-times.pdf",3,"$\\omega$","Time (seconds)",                      &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-8,10.0d0,                                                                         &
   nn,omegas,times,"best","$\\beta=0.50$*$\\beta=1.00$*$\\beta=1.50$*", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

call plot_functions2("experiment6-errs.py","experiment6-errs.pdf",3,"$\\omega$","Max Absolute Error", &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-16,1.0d0,                                                 &
   nn,omegas,errs,"best","$\\beta=0.50$*$\\beta=1.00$*$\\beta=1.50$*", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")


end program
