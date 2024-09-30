module test_freqind_funs

use chebyshev
use experutils

contains


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

! qs    = -da**2+(da+dnu)*(1+da+dnu)/cosh(ts)**2
! qders = -2*(da+dnu)*(1+da+dnu)*tanh(ts)/cosh(ts)**2
end subroutine


subroutine qfun_alf(n,ts,qs,qders,pars)
implicit double precision (a-h,o-z)
double precision           :: ts(n)
double precision           :: qs(n)
double precision           :: qders(n)
double precision, pointer  :: pars(:)
dnu   = pars(1)
dmu   = pars(2)

! qs    = (dnu+da)*(dnu+da+1)/(1-ts**2) + (1-da**2)/(1-ts**2)**2
! qders = (dnu+da)*(dnu+da+1)*2*ts/(1-ts**2)**2  + 4*(1-da**2)*ts/(1-ts**2)**3

! qs    = -da**2+(da+dnu)*(1+da+dnu)/cosh(ts)**2
! qders = -2*(da+dnu)*(1+da+dnu)*tanh(ts)/cosh(ts)**2

qs    = -dmu**2 + (dnu+1)*dnu/cosh(ts)**2
qders = -2*dnu*(dnu+1)/cosh(ts)**2*tanh(ts)

end subroutine


end module

program test_freqind

use experutils
use freqind
use test_freqind_funs

implicit double precision (a-h,o-z)
type(freqind_vars_t)            :: vars

double precision, allocatable   :: ab(:,:)
double precision, allocatable   :: alpha(:,:)
double precision, allocatable   :: alphader(:,:)
double precision, allocatable   :: alphader2(:,:)
double precision, allocatable   :: alphader0(:,:)

double precision, allocatable   :: vals(:), vals0(:), ts(:)
double precision, pointer       :: pars(:)

nruns    = 1000
ncheck   = 1000
pi       = acos(-1.0d0)
eps0     = epsilon(0.0d0)

k        = 16
epsdisc  = 1.0d-12
epssol  = 1.0d-12
call freqind_init(vars,k,epsdisc,epssol)

allocate(vals(ncheck),vals0(ncheck),ts(ncheck))
allocate(pars(2))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Test freqind_fast
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n        = 1000
dnu      = n
da       = -0.25d0

aval0    = -pi/2*(n+1)
ifleft   = 1

pars(1) = dnu
pars(2) = da

nints = 10
allocate(ab(2,nints))

do i=1,nints
! ab(1,i) = a + (b-a)*(i-1.0d0)/(nints-0.0d0)
! ab(2,i) = a + (b-a)*(i-0.0d0)/(nints-0.0d0)
int = nints-i+1
ab(1,int) = 1.0d0-1.8d0**(-nints+i)
ab(2,int) = 1.0d0-1.8d0**(-nints+i-1)
end do
a=ab(1,1)
b=ab(2,nints)

write (*,"(A)") "before freqind_high, ab = "
write (*,"(6(D15.7))") ab
allocate(alpha(k,nints), alphader(k,nints), alphader2(k,nints))

call elapsed(t1)
do irun=1,nruns
call freqind_high(vars,ier,a,b,ifleft,aval0,qfun_geg,pars,&
  nints,ab,alpha,alphader,alphader2)
end do
call elapsed(t2)
dtime = (t2-t1)/nruns*1000


if (ier .ne. 0) then
  print *,"after freqind_high, ier = ",ier
else
  write (*,"(A)")     "freqind_high time (ms) = "
  write (*,"(D15.7)") dtime

  call elapsed(t1)
  do irun=1,nruns
  do i=1,ncheck
  t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
  call freqind_eval(vars,nints,ab,alpha,alphader,t,aval,apval)
  val     = -sin(aval)/sqrt(apval)*(1-t**2)**(-(0.5d0+da)/2)
  ts(i)   = t
  vals(i) = val
  end do
  end do
  call elapsed(t2)
  dtime = (t2-t1)/nruns *1.0d0/(ncheck+0.0d0)
  write (*,"(A)")     "average eval time (secs) = "
  write (*,"(D15.7)") dtime

  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(i,t,der)
  do i=1,ncheck
  t = ts(i)
  call gegpolw(n,da,t,vals0(i),der)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
   
  errmax = maxval(abs(vals-vals0))

  write (*,"(A)") "Gegenbauer errmax = "
  write (*,"(D15.7)") errmax

endif

write(*,"(A)") ""
write(*,"(A)") ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Test the freqind_positive routine on a somewhat larger interval
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


n        = 1000
dnu      = n
da       =-0.492132232d0
a        = 0
b        = 0.9999d0
aval0    = -pi/2*(n+1)
ifleft   = 1

pars(1) = dnu
pars(2) = da

call elapsed(t1)
do irun=1,nruns
call freqind_positive(vars,ier,a,b,ifleft,aval0,qfun_geg,pars,&
  nints,ab,alpha,alphader,alphader2)
end do
call elapsed(t2)
dtime = (t2-t1)/nruns*1000

if (ier .ne. 0) then
   print *,"after freqind_positive ier = ",ier
else
  write (*,"(A)")     "freqind_positive time (ms) = "
  write (*,"(D15.7)") dtime

  write (*,"(A)") "after freqind_positive, ab  = "
  write (*,"(6(D15.7))") ab 


  call elapsed(t1)
  do irun=1,nruns
  do i=1,ncheck
  t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
  call freqind_eval(vars,nints,ab,alpha,alphader,t,aval,apval)

  val     = -sin(aval)/sqrt(apval)*(1-t**2)**(-(0.5d0+da)/2)
  ts(i)   = t
  vals(i) = val
  end do
  end do
  call elapsed(t2)
  dtime = (t2-t1)/nruns *1.0d0/(ncheck+0.0d0)
  write (*,"(A)")     "average eval time (secs) = "
  write (*,"(D15.7)") dtime


  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(i,t,der)
  do i=1,ncheck
  t = ts(i)
  call gegpolw(n,da,t,vals0(i),der)
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  errmax = maxval(abs(vals-vals0))

  write (*,"(A)") "Gegenbauer errmax = "
  write (*,"(D15.7)") errmax
endif

write(*,"(A)") ""
write(*,"(A)") ""



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Use the most robust phase routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


n        = 100
dnu      = n
da       =-0.492132232d0
a        = 0
b        = 0.99999d0
aval0    = -pi/2*(n+1)
ifleft   = 1

pars(1) = dnu
pars(2) = da

call elapsed(t1)
do irun=1,nruns
call freqind_phase(vars,ier,a,b,ifleft,aval0,qfun_geg,pars,&
 nints,ab,alpha,alphader,alphader2)
end do
call elapsed(t2)
dtime = (t2-t1)/nruns*1000

if (ier .ne. 0) then
   print *,"after freqind_phase ier = ",ier
else
  write (*,"(A)")     "freqind_phase time (ms) = "
  write (*,"(D15.7)") dtime

  write (*,"(A)") "after freqind_positive, ab  = "
  write (*,"(6(D15.7))") ab 


  call elapsed(t1)
  do irun=1,nruns
  do i=1,ncheck
  t = a + (b-a) * (i-1.0d0)/(ncheck-1.0d0)
  call freqind_eval(vars,nints,ab,alpha,alphader,t,aval,apval)
  val     = -sin(aval)/sqrt(apval)*(1-t**2)**(-(0.5d0+da)/2)
  ts(i)   = t
  vals(i) = val
  end do
  end do
  call elapsed(t2)
  dtime = (t2-t1)/nruns *1.0d0/(ncheck+0.0d0)
  write (*,"(A)")     "average eval time (secs) = "
  write (*,"(D15.7)") dtime


  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(i,t,der)
  do i=1,ncheck
  t = ts(i)
  call gegpolw(n,da,t,vals0(i),der)
  end do
  !$OMP END DO
  !$OMP END PARALLEL
   
  errmax = maxval(abs(vals-vals0))
  write (*,"(A)") "Gegenbauer errmax = "
  write (*,"(D15.7)") errmax
endif




end program
