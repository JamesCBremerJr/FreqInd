!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  EXPERIMENT 4:   Test the phase function code by evaluating the Wronskian normalized
!                  Gegenbauer polynomials on the interval [0,0.999].
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module experiment4_funs

use experutils
use freqind

type(freqind_vars_t)           :: vars
integer, parameter             :: nruns  = 100
integer, parameter             :: ncheck = 1000

contains

subroutine phase_geg(n,da,time,errmax)
implicit double precision (a-h,o-z)
!
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
double precision, allocatable   :: vals(:), vals0(:), ts(:)

data pi / 3.14159265358979323846264338327950288d0 /

dnu      = n
a        = 0
b        = 0.999d0
aval0    = -pi/2*(n+1)
ifleft   = 1




allocate(vals(ncheck),vals0(ncheck),ts(ncheck))
allocate(pars(2))

pars(1) = dnu
pars(2) = da

call elapsed(t1)
do irun=1,nruns
call freqind_phase(vars,ier,a,b,ifleft,aval0,qfun_geg,pars,&
  nints,ab,alpha,alphader,alphader2)
end do
call elapsed(t2)
time = (t2-t1)/nruns



  do i=1,ncheck
  t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
  call freqind_eval(vars,nints,ab,alpha,alphader,t,aval,apval)

  val     = -sin(aval)/sqrt(apval)*(1-t**2)**(-(0.5d0+da)/2)
  ts(i)   = t
  vals(i) = val
  end do



!  !$OMP PARALLEL DEFAULT(SHARED)
!  !$OMP DO PRIVATE(i,t,der)
  do i=1,ncheck
  t = ts(i)
  call gegpolw(n,da,t,vals0(i),der)
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

 errmax = maxval(abs(vals-vals0))


end subroutine


subroutine qfun_geg(n,ts,qs,qders,pars)
implicit double precision (a-h,o-z)
double precision           :: ts(n)
double precision           :: qs(n)
double precision           :: qders(n)
double precision, pointer  :: pars(:)
dnu   = pars(1)
da    = pars(2)-0.5d0
qs    = (dnu+da)*(dnu+da+1)/(1-ts**2) + (1-da**2)/(1-ts**2)**2
qders = (dnu+da)*(dnu+da+1)*2*ts/(1-ts**2)**2  + 4*(1-da**2)*ts/(1-ts**2)**3

end subroutine

end module

program experiment4

use experutils
use freqind
use experiment4_funs

implicit double precision (a-h,o-z)

double precision, allocatable :: times(:,:), errs(:,:), dnus(:)


epsdisc = 1.0d-12
epssol  = 1.0d-12
k       = 16
call freqind_init(vars,k,epsdisc,epssol)


jj1 = 6
jj2 = 20
nn  = jj2-jj1+1

allocate(errs(3,jj1:jj2),times(3,jj1:jj2), dnus(jj1:jj2))

do jj=jj1,jj2
n   = 2**jj
dnu = n
da1  =-0.499d0
da2  = 0.25d0
da3  = 1.00d0

call phase_geg(n,da1,time1,errmax1)
call phase_geg(n,da2,time2,errmax2)
call phase_geg(n,da3,time3,errmax3)

times(1,jj) = time1
times(2,jj) = time2
times(3,jj) = time3

errs(1,jj) = errmax1
errs(2,jj) = errmax2
errs(3,jj) = errmax3

dnus(jj)   = dnu

print "(I5.5,7(D13.5,1X))",jj,time1,time2,time3,errmax1,errmax2,errmax3

end do


call plot_functions2("experiment4-times.py","experiment4-times.pdf",3,"$n$","Time (seconds)",                      &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-7,10.0d0,                                                                         &
   nn,dnus,times,"best","$\\beta=-0.499$*$\\beta=0.25$*$\\beta=1.00$**", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")

call plot_functions2("experiment4-errs.py","experiment4-errs.pdf",3,"$n$","Max Absolute Error", &
  2,1,2.0d0**j1,2.0d0**j2,1.0d-16,1.0d0,                                                 &
   nn,dnus,errs,"best","$\\beta=-0.499$*$\\beta=0.25$*$\\beta=1.00$**", &
   "solid*dashdot*dashed*:*","#116611*#992255*#2222aa*#111111*")


end program
