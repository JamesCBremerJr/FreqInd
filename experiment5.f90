!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EXPERIMENT 5:  Test the phase function method by solving a boundary value problem
!                 for an oscillatory ODE.
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment5_funs

use experutils
use freqind
use chebyshev
use odetwo

type(odetwo_vars_t)            :: odetwo_vars
type(freqind_vars_t)           :: vars
integer, parameter             :: nruns  = 100
integer, parameter             :: ncheck = 1000

double precision               :: omega0
contains


subroutine phase(omega,time,timeode,errmax)
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

double precision                :: amatr(2,2)

double precision, allocatable   :: ab1(:,:)
double complex, allocatable     :: y1s(:,:)
double complex, allocatable     :: y1ders(:,:)

double precision, allocatable   :: ab2(:,:)
double complex, allocatable     :: y2s(:,:)
double complex, allocatable     :: y2ders(:,:)

double complex                  :: yc, ypc, val1, val2, der1, der2
double complex                  :: bmatr(2,2), detz, c1,c2, rhs(2), cs(2)
type(c_ptr)                     :: userptr


data pi / 3.14159265358979323846264338327950288d0 /
eps0     = epsilon(0.0d0)


dnu      = n
a        =  0.0d0
b        =  1.0d0
aval0    = 0
ifleft   = 1


z1    = 0.0d0
valz1 = 1.0d0

z2    = 1.0d0
valz2 = 1.0d0

omega0   = omega


allocate(vals(ncheck),vals0(ncheck),ts(ncheck))


call elapsed(t1)
do irun=1,nruns

call freqind_phase(vars,ier,a,b,ifleft,aval0,qfun,pars,&
  nints,ab,alpha,alphader,alphader2)

if (ier .ne. 0 ) then
print *,"!!!!!!!!!!!!!! ier = ",ier
stop
endif


call freqind_eval(vars,nints,ab,alpha,alphader,z1,aval,apval)
amatr(1,1) = sin(aval)/sqrt(apval)
amatr(1,2) = cos(aval)/sqrt(apval)

call freqind_eval(vars,nints,ab,alpha,alphader,z2,aval,apval)
amatr(2,1) = sin(aval)/sqrt(apval)
amatr(2,2) = cos(aval)/sqrt(apval)

det = amatr(1,1)*amatr(2,2) - amatr(1,2)*amatr(2,1)
d1  = ( amatr(2,2) * valz1 - amatr(1,2) * valz2 ) / det
d2  = (-amatr(2,1) * valz1 + amatr(1,1) * valz2 ) / det

end do


eps= 1.0d-13
call elapsed(t2)
time = (t2-t1)/nruns

omega0 = omega

call elapsed(t1)

c   = 0.0d0
yc  = 1
ypc = 0
call odetwo_nonlinear1(odetwo_vars,ier1,eps,a,b,c,yc,ypc,nints1,ab1,y1s,y1ders,odefun,userptr)

c   = 1.0d0
yc  = 1
ypc = 0
call odetwo_nonlinear1(odetwo_vars,ier2,eps,a,b,c,yc,ypc,nints2,ab2,y2s,y2ders,odefun,userptr)

if (ier1 +ier2 .gt. 0) then
print *,"@@@@@@@@@@@@@@@@@@@@@@@@@@@ ","ier1 = ",ier1," ier2 = ",ier2
stop
endif


!
!  find the constants
! 
call odetwo_eval(odetwo_vars,nints1,ab1,y1s,z1,bmatr(1,1))
call odetwo_eval(odetwo_vars,nints2,ab2,y2s,z1,bmatr(1,2))

call odetwo_eval(odetwo_vars,nints1,ab1,y1s,z2,bmatr(2,1))
call odetwo_eval(odetwo_vars,nints2,ab2,y2s,z2,bmatr(2,2))

detz = bmatr(1,1)*bmatr(2,2) - bmatr(1,2)*bmatr(2,1)
c1   = ( bmatr(2,2) * valz1 - bmatr(1,2) * valz2 ) / detz
c2   = (-bmatr(2,1) * valz1 + bmatr(1,1) * valz2 ) / detz


errmax =0 

call elapsed(t2)
timeode = t2-t1


do i=1,ncheck

t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)

call odetwo_eval(odetwo_vars,nints1,ab1,y1s,t,val1)
call odetwo_eval(odetwo_vars,nints2,ab2,y2s,t,val2)
val0 = val1*c1 + val2*c2

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

qs =omega**2*((0.2d1*exp(-ts))/(0.1d0+ts**2)+(1+omega*ts**2+  &
0.3d1*omega**2*ts**2)/(1+omega**2-0.1d1*omega*(1+ts**2)))

qders  =  omega**2*((-0.4d1*exp(-ts)*ts)/(0.1d0+ts**2)**2-  &
(0.2d1*exp(-ts))/(0.1d0+ts**2)+(0.2d1*omega*ts*(1+  &
omega*ts**2+3*omega**2*ts**2))/(1+omega**2-omega*(1+  &
ts**2))**2+(0.2d1*omega*ts+0.6d1*omega**2*ts)/(1+omega**2-  &
0.1d1*omega*(1+ts**2)))

end subroutine



end module

program experiment5

use experutils
use freqind
use experiment5_funs

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

allocate(errs(1,jj1:jj2),times(1,jj1:jj2), omegas(jj1:jj2))

do jj=jj1,jj2
omega   = 2**jj

call phase(omega,time1,time2,errmax1)

omegas(jj)  = omega
times(1,jj) = time1
errs(1,jj)  = errmax1

print "(I5.5,3(D13.5,1X))",jj,time1,time2,errmax1

end do


call plot_functions2("experiment5-times.py","experiment5-times.pdf",1,"$\\omega$","Time (seconds)",                      &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-6,1000.0d0,                                                                         &
   nn,omegas,times,"","", & 
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

call plot_functions2("experiment5-errs.py","experiment5-errs.pdf",1,"$\\omega$","Max Absolute Error", &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-16,1.0d0,                                                 &
   nn,omegas,errs,"","",&
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")


end program
