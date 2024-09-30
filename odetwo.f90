!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains a highly-robust but fairly slow solver for second order ordinary
!  differential equations.  It operates via an adaptive integration method and solutions
!  are represented as piecewise Chebyshev expansions. 
!
!  The following routines should be regarded as public:
! 
!    odetwo_init - initialize the solver by populating a structure which is passed to
!      the other subroutines
!
!    odetwo_nonlinear - adaptively solve a problem of the form
!
!        {  y'(t) = f(t,y(t))         a <= t <= b
!        {  y(c)  = yc
!        {  y'(c) = ypc
! 
!      where  a <= c <= b and f:R x C^2 --> C^2 is a smooth function supplied by 
!      an external subroutine
!
!    odetwo_eval - evaluate a solution produced by odetwo_nonlinear
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module odetwo

use chebyshev
use iso_c_binding
use ieee_arithmetic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The structure which stores any data needed by the procedures in this modules and
!  which is populated by the initialization routine.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type       odetwo_vars_t

integer                                    :: kcheb
integer                                    :: ntest
integer                                    :: maxstack
integer                                    :: maxints
integer                                    :: maxiters
integer                                    :: maxsteps

double precision, allocatable              :: xscheb(:)
double precision, allocatable              :: whtscheb(:)
double precision, allocatable              :: acoefs(:,:)
double precision, allocatable              :: aintl(:,:)
double precision, allocatable              :: aintr(:,:)
double precision, allocatable              :: aintl2(:,:)
double precision, allocatable              :: aintr2(:,:)
double precision, allocatable              :: dsigns(:)

end type   odetwo_vars_t


interface


subroutine odetwo_fun1(t,y,yp,f,dfdy,dfdyp,userptr)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
type(c_ptr)                 :: userptr
end subroutine

subroutine odetwo_fun2(t,y,yp,f,dfdy,dfdyp,pars)
use iso_c_binding
implicit double precision (a-h,o-z)
double complex              :: y, yp, f, dfdy, df, dfdyp
double complex, pointer     :: pars(:)
end subroutine

end interface

interface          odetwo_nonlinear
module procedure   odetwo_nonlinear1
end interface      odetwo_nonlinear

contains


subroutine odetwo_init(vars,kcheb0,ntest0,maxstack0,maxints0,maxsteps0,maxiters0)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)                         :: vars
integer, optional                           :: kcheb0, ntest0, maxstack0, maxsteps0
integer, optional                           :: maxints0, maxiters0
!
!  Initialize the structure containing any data needed by the other procedures 
!  in this module.
!
!  If vars is the only parameter passed to this routine, then reasonable defaults
!  are chosen.
!
!  Input parameters:
!    kcheb0 - the number of terms in the local Chebyshev expansions used to represent
!      solutions
!    ntest0 - the number of trailing terms in the local Chebyshev expansions to 
!      check when testing the fitness of a purported solution
!    maxstack0 - the length of the stack used by the adaptive procedure
!    maxints0 - the maximum number of intervals for the discretization scheme used
!      to represent solutions
!    maxsteps0 - the maximum number of iterations for the trapezoid rule
!    maxiters0 - the maximum number of iterations for Newton's method
!
!  Output parameters:
!    vars - the structure containing all of the data needed by the other procedures in
!      this module
!

if (.not. present(kcheb0) ) then
kcheb    = 16
ntest    = 2
maxstack = 30000000
maxints  = 30000000
maxsteps = 8
maxiters = 8
else
kcheb    = kcheb0
ntest    = ntest0
maxstack = maxstack0
maxints  = maxints0
maxsteps = maxsteps0
maxiters = maxiters0
endif

vars%kcheb    = kcheb
vars%ntest    = ntest
vars%maxstack = maxstack
vars%maxints  = maxints
vars%maxsteps = maxsteps
vars%maxiters = maxiters

call chebyshev_quad(kcheb,vars%xscheb,vars%whtscheb)
call chebyshev_coefsmatrix(kcheb,vars%acoefs)
call chebyshev_intlmatrix(kcheb,vars%aintl)
call chebyshev_intrmatrix(kcheb,vars%aintr)

allocate(vars%aintl2(kcheb,kcheb))
allocate(vars%aintr2(kcheb,kcheb))

vars%aintl2 = matmul(vars%aintl,vars%aintl)
vars%aintr2 = matmul(vars%aintr,vars%aintr)


! An array which helps in the calculation of barcyentric interpolation
! weights
allocate(vars%dsigns(kcheb))

vars%dsigns(1) = 1
do i=2,kcheb
vars%dsigns(i) = -vars%dsigns(i-1)
end do
vars%dsigns(1)     = vars%dsigns(1)/2
vars%dsigns(kcheb) = vars%dsigns(kcheb)/2

end subroutine


subroutine odetwo_fit(vars,vals,dcoefs)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)         :: vars
double complex              :: vals(:), coefs(vars%kcheb)

kcheb  = vars%kcheb
ntest  = vars%ntest
coefs  = matmul(vars%acoefs,vals)
dd1    = norm2(abs(coefs))
dd2    = norm2(abs(coefs(kcheb-ntest+1:kcheb)))
if (dd1 .eq. 0) dd1=1
dcoefs = dd2/dd1

end subroutine



subroutine odetwo_nonlinear1(vars,ier,eps,a,b,c,yc,ypc,nints,ab,ys,yders,odefun,userptr)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)                                 :: vars
double complex, allocatable, intent(out)            :: ys(:,:), yders(:,:)
double precision, allocatable, intent(out)          :: ab(:,:)
double complex                                      :: yc, ypc
procedure(odetwo_fun1)                              :: odefun
type(c_ptr)                                         :: userptr
!  
!  Adaptively solve the problem
!
!      {  y'(t) = f(t,y(t)),         a <= t <= b,
!      {  y(c)  = yc
!      {  y'(c) = ypc
! 
!  where a <= c <= b and f:R x C^2 --> C^2  is a smooth function suppplied
!  via an external subroutine.  The solution is represented via a piecewise
!  Chebyshev expansion given over a collection of subintervals.
!
!  Input parameters:
!    vars - the structure populated by odetwo_init
!    eps - the desired precision for the calculations
!    (a,b) - the interval over which the problem is given
!    c - the point c at which the values of y and y' are given
!    yc - the desired value of y(c)
!    ypc - the desired value of y'(c)
!
!  Output parameters:
!    nints - the number of discretization intervals
!    ab - a (2,nints) array whose jth column gives the endpoints of the jth
!      discretization interval
!    ys - the values of the obtained solution at the discretization nodes
!    ys - the values of the derivative of the obtained solution at the discretization 
!      nodes
!

double precision, allocatable                       :: stack(:,:), ab1(:,:), ab2(:,:), ts0(:)
double complex, allocatable                         :: ps0(:), qs0(:), fs0(:)
double complex, allocatable                         :: hs(:), hders(:), hder2s(:)
double complex, allocatable                         :: ys1(:,:), yders1(:,:), yder2s1(:,:)
double complex, allocatable                         :: ys2(:,:), yders2(:,:), yder2s2(:,:)

double complex                                      :: f, dfdy, dfdyp, ima

ier     = 0
epstrap = eps
epsnewt = eps
ima     = (0.0d0,1.0d0)

!
!  Read the parameters from the vars structure
! 

k             = vars%kcheb
ntest         = vars%ntest
maxstack      = vars%maxstack
maxints       = vars%maxints
maxsteps      = vars%maxsteps
maxiters      = vars%maxiters

!
!  Allocate temporary arrays
!


allocate(stack(2,maxstack))
 
allocate(ts0(k),ps0(k),qs0(k),fs0(k))
allocate(hs(k),hders(k),hder2s(k))

allocate(ab1(2,maxints), ys1(k,maxints), yders1(k,maxints), yder2s1(k,maxints) )
allocate(ab2(2,maxints), ys2(k,maxints), yders2(k,maxints), yder2s2(k,maxints) )

!
!  Solve going backward from the point c
!


nints1 = 0
ind    = 0

if (a < c) then

nstack     = 1
stack(1,1) = a
stack(2,1) = c

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints1 .eq. 0) then
ys1(k,nints1+1)    = yc
yders1(k,nints1+1) = ypc
else
ys1(k,nints1+1)    = ys1(1,nints1)
yders1(k,nints1+1) = yders1(1,nints1)
endif


! Build an initial solution using the trapezoid rule
call odetwo_traptvp1(jer,epstrap,maxsteps,k,ts0,odefun,ys1(:,nints1+1),yders1(:,nints1+1),yder2s1(:,nints1+1),userptr)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders1(:,nints1+1)  = matmul(vars%aintr*(b0-a0)/2,yder2s1(:,nints1+1)) + yders1(k,nints1+1)
ys1(:,nints1+1)     = matmul(vars%aintr*(b0-a0)/2,yders1(:,nints1+1))  + ys1(k,nints1+1)

!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys1(i,nints1+1),yders1(i,nints1+1),f,dfdy,dfdyp,userptr)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s1(i,nints1+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!

hs(k)    = 0
hders(k) = 0
call odetwo_lineartvp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)


if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys1(:,nints1+1)))) exit

!
!  Update the solution
!

ys1(:,nints1+1)     = ys1(:,nints1+1) + hs
yders1(:,nints1+1)  = yders1(:,nints1+1) + hders
yder2s1(:,nints1+1) = yder2s1(:,nints1+1) + hder2s

end do


 call odetwo_fit(vars,ys1(:,nints1+1),dy)
!call odetwo_fit(vars,yders1(:,nints1+1),dyp)
if (dy .lt. eps ) ifaccept = 1
endif

! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints1+1 .gt. maxints) then
ier    = 4
return
endif

nints1        = nints1+1
ab1(1,nints1) = a0
ab1(2,nints1) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

endif


end do

endif

!
!  Solve going forward from the point c
!

nints2 = 0
ind    = 0

if (b > c) then

nstack     = 1
stack(1,1) = c
stack(2,1) = b

do while (nstack > 0 )

a0         = stack(1,nstack)
b0         = stack(2,nstack)
nstack     = nstack-1

ts0        = min(b0,max(a0,(b0-a0)/2 * vars%xscheb + (b0+a0)/2))
ifaccept   = 0
dcoefs     = -1
dy         = -1
dy0        = -1
dyp        = -1

! Evaluate the coefficients 

!
!  Use the implicit trapezoidal method to construct an approximation
!  to use as an initial guess for Newton's method.
!

if (nints2 .eq. 0) then
ys2(1,nints2+1)    = yc
yders2(1,nints2+1) = ypc
else
ys2(1,nints2+1)    = ys2(k,nints2)
yders2(1,nints2+1) = yders2(k,nints2)
endif



! Build an initial solution using the trapezoid rule
call odetwo_trapivp1(jer,epstrap,maxsteps,k,ts0,odefun,ys2(:,nints2+1),yders2(:,nints2+1),yder2s2(:,nints2+1),userptr)

if (jer .eq. 0) then

!  The values of y, y' and y'' calculated using the implicit trapezoidal method satisfy
!  the differential equation at the specified nodes, but they are not consistent with eachother.  
!  That is, y' is not the derivative of y in any reasonable sense.
!
!  We integrate the second derivative twice in order to produce consistent
!  approximations.
!
!  Skipping this step is fatal.
!

yders2(:,nints2+1)  = matmul(vars%aintl*(b0-a0)/2,yder2s2(:,nints2+1)) + yders2(1,nints2+1)
ys2(:,nints2+1)     = matmul(vars%aintl*(b0-a0)/2,yders2(:,nints2+1))  + ys2(1,nints2+1)


!
!  Perform Newton iterations.
!

do inewt = 1, maxiters

!
!  Form the coefficients for the linearized problem and compute
!  the error in the current solution.
!
do i=1,k
call odefun(ts0(i),ys2(i,nints2+1),yders2(i,nints2+1),f,dfdy,dfdyp,userptr)
ps0(i)    = -dfdyp
qs0(i)    = -dfdy
fs0(i)    = f-yder2s2(i,nints2+1)
end do

!
!  Call the code_linear_ivp routine to solve the linearized system.
!


hs(1)    = 0
hders(1) = 0
call odetwo_linearivp1(vars,a0,b0,k,ps0,qs0,fs0,hs,hders,hder2s)

if ( norm2(abs(hs)) .lt. eps*norm2(abs(ys2(:,nints2+1)))) exit

!
!  Update the solution
!

ys2(:,nints2+1)     = ys2(:,nints2+1) + hs
yders2(:,nints2+1)  = yders2(:,nints2+1) + hders
yder2s2(:,nints2+1) = yder2s2(:,nints2+1) + hder2s

end do

call odetwo_fit(vars,ys2(:,nints2+1),dy)
!call odetwo_fit(vars,yders2(:,nints2+1),dyp)

if (dy .lt. eps) ifaccept = 1

endif


! ind = ind+1
! write (*,"(2(I6),3(D20.10,2X),1X,I1)")    ind,jer,a0,b0,dyp,ifaccept
! write (13,"(2(I6),3(D20.10,2X),1X,I1)")   ind,jer,a0,b0,dyp,ifaccept

if (ifaccept .eq. 1) then

if (nints2+1 .gt. maxints) then
ier    = 4
return
endif

nints2        = nints2+1
ab2(1,nints2) = a0
ab2(2,nints2) = b0

else

if (nstack+2 .gt. maxstack) then
ier = 8
return
endif
 
nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2
endif


end do

endif


!
!  Copy out the solution.
!

nints = nints1+nints2

allocate(ab(2,1:nints))
ab(:,1:nints1)       = ab1(:,nints1:1:-1)
ab(:,nints1+1:nints) = ab2(:,1:nints2)


allocate(ys(k,nints), yders(k,nints) )

ys(:,1:nints1)          = ys1(:,nints1:1:-1)
ys(:,nints1+1:nints )   = ys2(:,1:nints2)

yders(:,1:nints1)       = yders1(:,nints1:1:-1)
yders(:,nints1+1:nints) = yders2(:,1:nints2)


end subroutine


subroutine odetwo_eval(vars,nints,ab,vals,x,val)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)              :: vars
double precision                 :: ab(:,:)
double complex                   :: vals(:,:), val
!
!  Evaluate the solution of 
!
!  Input parameters:
!    vars - the structure initialized by odetwo_init
!    nints - the number of intervals in the piecewise Chebyshev scheme representing
!      the solution
!    ab - an (2,nints) array whose jth column gives the endpoints of the
!     jth discretization interval
!    vals - the values of the solution at the discretization nodes,
!      as returned by odetwo_nonlinear
!    x- the point at which to evaluate the solution
!
!  Output parameters:
!    val - the value of the solution at the point x
!
!
double precision                 :: whts(vars%kcheb)

eps0  = epsilon(0.0d0)
kcheb = vars%kcheb

!
!  Find the interval containing x ... unless the number of intervals is 
!  extremely large, it is fastest to use a "brute force" linear search.
!

do int=1,nints-1
if (x .lt. ab(2,int)) exit
end do
a = ab(1,int)
b = ab(2,int)

!
!  Map the point x from [a,b] to [-1,1]
!
xx = (2*x - (b+a) ) /(b-a)

!
!  Perform barcyentric interpolation
!
do i=1,kcheb
  dd      = xx - vars%xscheb(i)
  if ( abs(dd) .le. eps0) then
  val = vals(i,int)
  return
  endif
  whts(i) = vars%dsigns(i)/dd
end do

val    = dot_product(whts,vals(:,int)) / sum(whts)

end subroutine



subroutine odetwo_linearivp1(vars,a,b,k,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)           :: vars
integer, intent(in)           :: k
double complex, intent(in)    :: ps(k),qs(k),fs(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)
!
!  Solve an initial value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                              (3)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!    (a,b) - the interval on which the ODE (3) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (3)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (3)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(1) - the value of y(a)
!    yders(1) - the value of y'(a)
!
!  Output parameters:
!
!    ys - the values of the solution y of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (3) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!

double complex, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)
double complex              :: alpha,beta

!
!  Allocate memory for the procedure and setup some parameters.
!

allocate(amatr(k,k),xs(k),sigma(k),rhs(k))


xs       = max(a,min(b,(b-a)/2 *vars%xscheb + (b+a)/2))
alpha    = ys(1)
beta     = yders(1)

!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-a) + \int_a^t (t-s) sigma(s) ds,
!
!  insert this representation into (1), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintl(i,:)*(b-a)/2
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintl2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-a))
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call odetwo_solve_c(k,amatr,sigma)


!
!  Calculate y(t) and y'(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(vars%aintl,sigma)
ys     = ((b-a)/2)**2*matmul(vars%aintl2,sigma)

do i=1,k
ys(i)     = ys(i) + alpha + beta*(xs(i)-a)
yders(i)  = yders(i) + beta
end do

end subroutine



subroutine odetwo_lineartvp1(vars,a,b,k,ps,qs,fs,ys,yders,yder2s)
implicit double precision (a-h,o-z)
type(odetwo_vars_t)           :: vars
integer, intent(in)           :: k
double complex, intent(in)    :: ps(k),qs(k),fs(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)

!
!  Solve a terminal value for the ordinary differential equation
!
!     y''(t) + p(t) y'(t) + q(t) y(t) = f(t)                                              (4)
!
!  on the interval [a,b] using a standard Chebyshev spectral method.
!
!  Input parameters:
!
!    (a,b) - the interval on which the ODE (4) is given
!    k - the number of Chebyshev points on the interval [a,b]
!    ps - an array specifying the values of the function p(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!    qs - an array specifying the values of the function q(t) appearing in (4)
!      at the k Chebyshev node on [a,b]
!    fs - an array speciying the values of the function f(t) appearing in (4)
!      at the k Chebyshev nodes on [a,b]
!
!    ys(k) - the value of y(b)
!    yders(k) - the value of y'(b)
!
!  Output parameters:
!
!    ys - the values of the solution y of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yders - the values of the solution y' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!    yder2s - the values of the solution y'' of (4) at the nodes of the k-point
!      Chebyshev grid on [a,b]
!

double complex, allocatable :: amatr(:,:),xs(:),sigma(:),rhs(:)
double complex              :: alpha,beta
integer, allocatable        :: ipiv(:)

!
!  Allocate memory for the procedure and setup some parameters.
!
allocate(amatr(k,k),xs(k),sigma(k))

xs       = max(a,min(b,(b-a)/2 * vars%xscheb + (b+a)/2))
alpha    = ys(k)
beta     = yders(k)

!
!  We represent the solution in the form
!
!      y(t) = alpha + beta (t-b) + \int_b^t (t-s) sigma(s) ds,
!
!  insert this representation into (2), and solve the resulting system of
!  linear equations in order to obtain the values of sigma.  
!    

amatr  = 0
do i=1,k
amatr(i,i) = 1.0d0
sigma(i)   = 0.0d0
end do

!
! Handle the p(t) * y'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps(i) * vars%aintr(i,:)*(b-a)/2
sigma(i)   = sigma(i) - ps(i)*beta
end do

!
!  Handle the q(t) y(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs(i) * vars%aintr2(i,:)*((b-a)/2)**2
sigma(i)   = sigma(i) - qs(i) * (alpha + beta*(xs(i)-b))
end do


!
!  Form the right-hand side.
!
do i=1,k
sigma(i) = sigma(i) + fs(i)
end do

!
!  Use a QR decomposition to invert the linear system
!

call odetwo_solve_c(k,amatr,sigma)


!
!  Calculate y(t) and y'(t) from sigma.
!

yder2s = sigma
yders  = (b-a)/2*matmul(vars%aintr,sigma)
ys     = ((b-a)/2)**2*matmul(vars%aintr2,sigma)

do i=1,k
ys(i)     = ys(i) + alpha + beta*(xs(i)-b)
yders(i)  = yders(i) + beta
end do

end subroutine



subroutine odetwo_trapivp1(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,userptr)
implicit double precision (a-h,o-z)
integer, intent(out)          :: ier
integer, intent(in)           :: k
double precision, intent(in)  :: ts(k)
double complex, intent(out)   :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun1)       :: odefun
type(c_ptr)                   :: userptr

!
!  Use the implicit trapezodial method to crudely approximate the solution of an
!  initial value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    eps - precision for the calculation
!    maxsteps - the maximum number of Newton steps
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (4) is to be solved
!    odefun - the user-supplied external routine conforming to the interface odetwo_fun1
!
!    ys(1) - the first entry of this array is the initial value for the solution y
!    yders(1) - the first entry of this array is the initial value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!      ier = 0    indicates successful execution
!      ier = 16   means that the Newton procedure did not converge
!      ier = 64   means that NaN or Inf was encountered
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier      = 0
! eps0     = epsilon(0.0d0)
! eps      = eps0*10

!
!  Evaluate the second derivative at the left-hand endpoint of the interval.
!

call odefun(ts(1),ys(1),yders(1),yder2s(1),dfdy,dfdyp,userptr)


do i=2,k
t0 = ts(i-1)
t  = ts(i)
h  = t-t0

y0   = ys(i-1)
yp0  = yders(i-1)
ypp0 = yder2s(i-1)

!
!  Set the initial guess.
!

!yp1 = yp0 + h *ypp0
yp1 = yp0 
y1  = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,userptr)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Conduct Newton iterations in an attempt to improve the siutation.
!

do iter=1,maxsteps

!
!  Record the current approximation.
!
! dd0       = dd
! ys(i)     = y1
! yders(i)  = yp1
! yder2s(i) = ypp1

!
!  Take a Newton step.
!

val   = dd
der   = 1.0d0-h/2*(dfdy*h/2+dfdyp)
delta = val/der

yp1   = yp1-delta
y1    = y0 + h/2*(yp0+yp1)
call odefun(t,y1,yp1,ypp1,dfdy,dfdyp,userptr)
dd    = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp1)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do


if (ieee_is_nan(real(y1)) .OR. ieee_is_nan(imag(y1)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y1)) .OR. .not. ieee_is_finite(imag(y1))) then
ier = 64
return
endif

ys(i)     = y1
yders(i)  = yp1
yder2s(i) = ypp1


if (iter .gt. maxsteps) then
ier = 16
return
endif

end do




end subroutine


subroutine odetwo_traptvp1(ier,eps,maxsteps,k,ts,odefun,ys,yders,yder2s,userptr)
implicit double precision (a-h,o-z)
integer, intent(out)           :: ier
integer, intent(in)            :: k
double precision, intent(in)   :: ts(k)
double complex, intent(out)    :: ys(k),yders(k),yder2s(k)
procedure (odetwo_fun1)        :: odefun
type(c_ptr)                    :: userptr
!
!  Use the implicit trapezoidal method to crudely approximate the solution of a
!  terminal value problem for the nonlinear ordinary differential equation
!
!    y''(t) = f(t,y(t),y'(t))                                                             (5)
! 
!  on the interval [a,b].  The user specifies the nodes on which the solution
!  of is to be computed.
!
!  Input parameters:
!
!    k - the number of nodes at which the solution of (5) is to be computed
!    ts - an array of length k which supplies the a sorted list of nodes at which 
!     (5) is to be solved
!    odefun - a user-supplied external procedure of type "codefunction" which
!     supplied the values of the function f as well as the derivatives of
!     f w.r.t. y' and y  (see the definition of codefunction above)
!
!    ys(k) - the last entry of this array is the terminal value for the solution y
!    yders(k) - the first entry of this array is the terminal value for the solution y'
!    userptr - a "void *" pointer which is passed as an argument to the user-specified
!     external function odefun
!
!  Output parameters:
!
!    ier - an error return code;
!       ier = 0    indicates successful execution
!       ier = 1024 means NaN or Infinity was encountered; this usually means
!                  that a spectacular overflow took place
!
!    ys - the values of the obtained approximation of the solution of (5) at
!     the specified nodes
!    yders - the values of the obtained approximation of the derivative of the solution of 
!     (5) at the specified nodes
!    yder2s - the values of the obtained approximation of the second derivative of the solution 
!     of  (5) at the specified nodes
!  

double complex :: y0,yp0,ypp0,ypp1,yp1,y1,dd,dd0,dfdy,dfdyp,val,der,delta

ier        = 0

call odefun(ts(k),ys(k),yders(k),yder2s(k),dfdy,dfdyp,userptr)

do i=k-1,1,-1
t  = ts(i)
h  = ts(i+1)-ts(i)

y1   = ys(i+1)
yp1  = yders(i+1)
ypp1 = yder2s(i+1)

!
!  Set the initial guess.
!

yp0 = yp1
y0  = y1 - h/2*(yp0+yp1)

call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,userptr)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)

!
!  Try to improve it via Newton iterations
!

do iter=1,maxsteps

!
!  Record the current guess
!

! dd0       = dd
! ys(i)     = y0
! yders(i)  = yp0
! yder2s(i) = ypp0

!
!  Take a Newton step
!
val   = dd
der   = -1.0d0-h/2*(-dfdy*h/2+dfdyp)
delta = val/der

yp0   = yp0 -delta
y0    = y1 - h/2*(yp0+yp1)
call odefun(t,y0,yp0,ypp0,dfdy,dfdyp,userptr)
dd  = yp1 - yp0 - h/2*(ypp0+ypp1)


ddd = abs(yp0)
if (ddd .eq. 0) ddd=1
if (abs(delta) .lt. eps * ddd)  exit

end do



if (iter .gt. maxsteps) then
ier = 16
return
endif

if (ieee_is_nan(real(y0)) .OR. ieee_is_nan(imag(y0)) ) then
ier = 64
return
endif

if (.not. ieee_is_finite(real(y0)) .OR. .not. ieee_is_finite(imag(y0))) then
ier = 64
return
endif

ys(i)     = y0
yders(i)  = yp0
yder2s(i) = ypp0

end do

end subroutine




subroutine odetwo_solve_c(n,a,rhs)
implicit complex *16 (a-h,o-z)
double complex         :: a(n,n),u(2,2),rhs(n)
! 
!  This subroutine uses a version of QR-decomposition to solve
!  the user-supplied system of linear algebraic equations
!  Ax = y.
!
!  THE MATRIX A IS DESTROYED BY THIS ROUTINE AND THE VECTOR RHS IS
!  OVERWRITTEN.
! 
!  Input parameters:
!    n - the dimensionality of the system being solved
!    a - the (n,n) complex-valued coefficient matrix WHICH WILL
!      BE DESTROYED BY THIS ROUUTINE
!    rhs - the right-hand side of the system to be solved, which
!      is overwritten by this subroutine
! 
!  Output parameters:
!    rhs - the solution of the linear system
!

double precision       ::  dmax,dmin,size22,dd,eps0
double complex         ::  aa(2)
integer, allocatable   ::  ipiv(:)

!
! Transpose the input matrix a
! 

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)*conjg(a(j,i))
size22=size22+a(i,j)*conjg(a(i,j))
end do
end do

! 
! Eliminate the upper right triangle
! 
do i=1,n-1
! 
do j=n,i+1,-1
! 
aa(1)=a(i,j-1)
aa(2)=a(i,j)
u21=-aa(2)
u22=aa(1)
dd=u22*conjg(u22)+u21*conjg(u21)
if(dd .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
dd=sqrt(dd)
u(2,2)=u22/dd
u(2,1)=u21/dd
u(1,1)=-conjg(u(2,2))
u(1,2)=conjg(u(2,1))
endif


do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j)
a(ii,j-1) = d1
a(ii,j)   = d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix to rhs
!

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

end subroutine


end module
