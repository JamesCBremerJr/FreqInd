!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for rapidly constructing slowly-varying phase functions
!  which represent the solutions second order differential equations of the form
!
!    y''(t) + q(t) y(t) = 0, a < t < b,                                                    (1)
!
!  where q is real-valued and slowly-varying on [a,b].  There are several different
!  routines for constructing phase functions, the most general of which is fully
!  adaptive and applies in the case in which q(t) is large and positive on at least
!  some portion of [a,b] and has no turning points of EVEN order.  Other routines 
!  require q(t) to be positive on all of [a,b].
!
!  A phase function for (1) is, of course, a function such that 
!
!    sin( alpha(t) )                 cos(alpha(t))
!    ---------------      and       -----------------                                      (2)
!    sqrt(alpha'(t))                sqrt(alpha'(t))
!
!  comprise a basis in the space of solutions of (1).  Their raison d'être is that
!  highly oscillatory and other rapidly varying solutions of equations of the form 
!  (1) which are costly to represent via polynomial expansions and other standard
!  techniques can be represented via slowly-varying phase functions.
!
!  This code could be easily modified to handle the ostensibly more general case of 
!  equations of the form
!
!    y''(t) + p(t) y'(t) + q(t) y(t) = 0, a < t < b,                                       (3)
!
!  but, of course, all equations of this form can be put into the form (1) with
!  a simple transformation. 
!
!  The following subroutines should be regarded as publicly callable:
!
!    freqind_init - set the various parameters and construct any other data needed 
!      by the other routines in this module.  This routine must be called before 
!      any other subroutine in this module.
!
!    freqind_phase - construct a phase function for (1) in the case in which q(t) is 
!      large and positive somewhere on the interval [a,b] and has only turning points of 
!      ODD order.  In particular, q(t) can be negative on some portion of [a,b].
!
!      This is the most general and robust of the phase function routines in this
!      module.  It works in cases in which q(t) is negative on parts of [a,b], as 
!      long as it is large and positive somewhere and all of its turning points are
!      of ODD order.  Moreover, It adaptively discretizes the solution very carefully,
!      even checking the solution's accuracy in the high-frequency regime, where doing 
!      so verges on paranoia.
!
!      No attempt to parallelize this code has been made.
!
!    freqind_positive - construct a phase function for an equation of the form (1)
!      in the case in which q(t) is positive on all of [a,b].
!
!      This code constructs an adaptive piecewise Chebyshev discretization of the 
!      Liouville-Green phase and uses it to represent the phase function.  Difficulties 
!      will be  encountered if the coefficient q(t) is too close to 0, in which case
!      the procedure will fail with an error code.
!
!    freqind_high - construct a phase function for an equation of the form (1) in
!      the case in which q(t) is of high-frequency throughout the interval [a,b].
!
!      This is the fastest but least robust of the phase function construction 
!      routines provided in this module.  It does no adaptive discretization
!      whatsoever, but rather requires the user to specify a piecewise
!      Chebyshev discretization scheme for representing the phase function.
!      It also assumes that the coefficient q(t) is large and positive through
!      the interval [a,b].  On the other hand, it is very fast and because
!      of these assumptions the calculations can easily be conducted in parallel.
!
!    freqind_eval - evaluate a phase function constructed and up to two of its
!      derivatives at an arbitrary point in the solution domain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module freqind

use chebyshev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This structure stores algorithm parameters as well as all necessary spectral
!  integration and differentiation matrices.  It needs to be initialized by 
!  freqind_init and passed to all of the other routines in this module.
!
!  Parameters:
!
!    k -  the numbr of terms in the piecewise Chebyshev expansions used to represent 
!      phase functions 
!    ntest - number of trailing Chebyshev coefficients to examine for ``goodness of fit''
!    maxints - the maximum number of allowable intervals in a piecewise Chebyshev 
!      expansion used to represent a phase function
!    maxstack - the maximum number of entries in the stacks used in adaptive
!      discretization
!    maxnewt - the maximum number of Newton iterations performed during the solution
!      of the Riccati equation
!    maxjac - the maximum number of Jacobi iteration performed when solving the linearized
!      equations which arise when the Newton method is used to solve the Riccati equation
!    thresh - the threshold for determining if an interval is in the "high-frequency"
!      regime
!    epsdisc - the desired precision for adaptive discretization of the piecewise
!      Chebyshev expansions used to represent phase functions
!    epssol - the desired precision for the solver which is applied to the
!      Riccati equation; this should be somewhat smaller than eps
!
!
!  Variables related to spectral integration / differentiation:
!
!    xs - the nodes of the k-point Chebyshev extremal grid on the interval [-1,1]
!    whts - the weights of the k-point Cheb
!    acoefs - the (k,k) matrix taking values at the grid points to coefficients
!      in the Chebyshev expansion
!    adiff - the (k,k) spectral differentiation matrix
!    aintl? - the ?th "left" spectral integration matrix 
!    aintr? - the ?th "right" spectral integration matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type     freqind_vars_t
integer                           :: k
integer                           :: ntest
integer                           :: maxints
integer                           :: maxstack
integer                           :: maxnewt
integer                           :: maxjac
double precision                  :: thresh
double precision                  :: threshlow
double precision                  :: epsdisc
double precision                  :: epssol
double precision                  :: dknorm      

double precision, allocatable     :: xs(:)
double precision, allocatable     :: whts(:)
double precision, allocatable     :: adiff(:,:)
double precision, allocatable     :: aintl(:,:)
double precision, allocatable     :: aintl2(:,:)
double precision, allocatable     :: aintl3(:,:)
double precision, allocatable     :: aintr(:,:)
double precision, allocatable     :: aintr2(:,:)
double precision, allocatable     :: aintr3(:,:)
double precision, allocatable     :: acoefs(:,:)
double precision, allocatable     :: acoefs2(:,:)
double precision, allocatable     :: ainterp(:,:)
double precision, allocatable     :: dsigns(:)

end type freqind_vars_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The interface for the user-supplied external subroutines which specify the values of 
!  the coefficient q(t) in (1).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
interface

subroutine freqind_fun(n,ts,qs,qders,pars)
double precision               :: ts(n)
double precision               :: qs(n)
double precision               :: qders(n)
double precision, pointer      :: pars(:)
end subroutine

end interface

interface        freqind_eval
module procedure freqind_eval1
module procedure freqind_eval2
module procedure freqind_eval3
end interface    freqind_eval

contains


subroutine freqind_init(vars,k,epsdisc,epssol)
implicit double precision (a-h,o-z)
type(freqind_vars_t), intent(out)    :: vars
!
!  Sets parameters for the algorithms implemented in this module and construct
!  any needed spectral integration and differential matrices.
!
!  See the comments pertaining the the freqind_vars_t structure above for a description 
!  of the function of each parameter.  
!
!  Input parameters:
!    k - the number of points used in the Chebyshev expansions
!
!  Output parameters:
!    vars - a structure containing all of the data needed by the routines
!      in this module
!

double precision, allocatable  :: xsout(:), whtsout(:)

double precision, allocatable :: amatr(:,:)
double precision              :: wr(k),wi(k), work(10000), sings(k)

eps0            = epsilon(0.0d0)
ntest           = 2
maxints         = 10000
maxstack        = 10000
maxnewt         = 10
maxjac          = 10
thresh          = 10.0d0
threshlow       = 1.0d0

! epsdisc         = max(1.0d-20,epsilon(0.0d0)*100)
! epssol          = max(1.0d-20,epsilon(0.0d0)*100)

! if (eps0 .lt. 1.0d-20) then
! epsdisc         = 1.0d-19
! epssol          = 1.0d-20
! else
! epsdisc         = 1.0d-13
! epssol          = 1.0d-13
! endif

vars%k          = k
vars%ntest      = ntest
vars%maxints    = maxints
vars%maxstack   = maxstack
vars%maxnewt    = maxnewt
vars%maxjac     = maxjac
vars%thresh     = thresh
vars%threshlow  = threshlow
vars%epsdisc    = epsdisc
vars%epssol     = epssol

call chebyshev_quad(k,vars%xs,vars%whts)
call chebyshev_diffmatrix(k,vars%adiff)
call chebyshev_intlmatrix(k,vars%aintl)
call chebyshev_intrmatrix(k,vars%aintr)
call chebyshev_coefsmatrix(k,vars%acoefs)




allocate(vars%acoefs2(ntest+1,k))
vars%acoefs2(1,:)         = vars%acoefs(1,:)
vars%acoefs2(2:ntest+1,:) = vars%acoefs(k-ntest+1:k,:)

allocate(xsout(2*k))

a0 = -1.0d0
b0 =  0.0d0
c0 =  1.0d0

xsout(1:k)     = (b0-a0)/2*vars%xs + (b0+a0)/2
xsout(k+1:2*k) = (c0-b0)/2*vars%xs + (c0+b0)/2
call chebyshev_interpmatrix(k,2*k,xsout,vars%ainterp)

allocate(vars%aintl2(k,k))
allocate(vars%aintl3(k,k))
vars%aintl2 = matmul(vars%aintl,vars%aintl)
vars%aintl3 = matmul(vars%aintl,vars%aintl2)

allocate(vars%aintr2(k,k))
allocate(vars%aintr3(k,k))
vars%aintr2 = matmul(vars%aintr,vars%aintr)
vars%aintr3 = matmul(vars%aintr,vars%aintr2)

! An array which helps in the calculation of barcyentric interpolation
! weights
allocate(vars%dsigns(k))

vars%dsigns(1) = 1
do i=2,k
vars%dsigns(i) = -vars%dsigns(i-1)
end do
vars%dsigns(1) = vars%dsigns(1)/2
vars%dsigns(k) = vars%dsigns(k)/2



end subroutine


subroutine freqind_high(vars,ier,a,b,ifleft,aval0,qfun,pars,nints,ab, &
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                         :: vars
procedure(freqind_fun)                       :: qfun
double precision, pointer                    :: pars(:)
double precision                             :: ab(:,:)
double precision                             :: alpha(:,:)
double precision                             :: alphader(:,:)
double precision                             :: alphader2(:,:)
!
!  Construct a phase function for an equation of the form (1) in the case in which
!  q(t) is large and positive on the entire solution interval [a,b].  That is to say,
!  the solutions of equation (1) are highly-oscillatory throughout [a,b].
!
!  This routine performs no adaptive discretization whatsoever and instead requires
!  that the user specify a discretization scheme for representing the phase
!  function.
!
!  The procedure will fail with an error code if the Newton iterations used to
!  construct the phase function do not converge, which happens when q(t) is 
!  smaller than needed.  
!
!  Input parameters:
!    vars - the structure returned by freqind_init
!    ifcheck - an integer flag indicating whether (ifcheck=1) or not (ifcheck=0) the accuracy 
!      of the solution should be 
!    [a,b] - the domain over which the ODE to solve is given
!    ifleft - an integer indicating whether the value of the phase function
!      alpha is specified at a (ifleft=1) or b (ifleft=0)
!    aval0 - the initial or terminal value of the phase function
!    nints - the number of intervals in the piecewise Chebyshev discretization
!      scheme used to represent the phase function
!    ab - an (2,nints) array whose jth column gives the left and right endpoints
!      of the jth discretization interval --- note that we must have
!      a_j < b_j = a_{j+1} < ..
!
!    qfun - an external subroutine conforming to the freqind_fun interface
!      which specifies the coefficient q(t) in (1)
!    pars - a user-supplied array of doubles which is passed to qfun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0     indicates successful execution
!      ier = 256   indicates that the Newton iterations used to solve the Riccati
!                  Riccati equation in the high-frequency regime did not converge
!      ier = 512   indicates the obtained phase function was not well-represented by the 
!                  discretization scheme on one or more intervals --- the subroutine
!                  only checks for this condition when ifcheck = 1
!
!    alpha -  this user-allocated array will contain the values of the phase function 
!      at the discretization nodes
!    alphader -  this user-allocated array will contain the  values of the derivative of 
!      the phase function at the discretization nodes
!    alphader2 -  this user-allocated array will contain the  values of the second 
!      derivative of the phase function at the discretization nodes


double precision              :: qs(vars%k), qders(vars%k), ts(vars%k), aint(vars%k,vars%k)
double complex                :: rs(vars%k), rders(vars%k), coefs1(vars%ntest+1)

ier        = 0
ima        = (0.0d0,1.0d0)
k          = vars%k
epssol     = vars%epssol
maxnewt    = vars%maxnewt
maxjac     = vars%maxjac
ntest      = vars%ntest

if(ifleft .eq. 1) then
aint=vars%aintl
else
aint=vars%aintr
endif

!
!  Solve the Riccati equation
!

ier1 = 0
ier2 = 0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(a0,b0,int,jer,rs,rders,ts,qs,qders,coefs1,dd1,dd2,dcoefs)
!$OMP DO
do int=1,nints
  a0   = ab(1,int)
  b0   = ab(2,int)

  ts = vars%xs * (b0-a0)/2 + (b0+a0)/2
  call qfun(k,ts,qs,qders,pars)
  ! qders = 2/(b0-a0)/2 * matmul(vars%adiff,qs)  
  
  call freqind_highfreq(vars,jer,epssol,maxnewt,maxjac,a0,b0,k,qs,qders, &
    rs,rders)

  alphader(:,int)    = imag(rs) 
  alphader2(:,int)   = imag(rders)
  alpha(:,int)       = (b0-a0)/2 * matmul(aint,alphader(:,int))

  ! UNCOMMENT TO TEST THE PHASE FUNCTION FOR "GOODNESS OF FIT"
!   coefs1   = matmul(vars%acoefs2,rs)
!   dd1     = abs(coefs1(1))
!   dd2     = maxval(abs(coefs1(2:ntest)))
!   dcoefs  = dd2/dd1
!   if (dcoefs .gt. epsdisc) then
! !$OMP ATOMIC 
!      ier2 = ier2+1
!   endif

  if (jer .ne. 0) then
!$OMP ATOMIC 
     ier1 = ier1+1
  endif
  
end do
!$OMP END DO
!$OMP END PARALLEL

if (ier1 .gt. 0 ) then
ier = 256
return
endif

if (ier2 .gt. 0 ) then
ier = 512
return
endif

!
!  Add the user-specified constant to the phase function
!

if (ifleft .eq. 1) then
  aval = aval0
  do int=1,nints
  alpha(:,int)       = aval + alpha(:,int)
  aval               = alpha(k,int)
  end do
else
  aval = aval0
  do int=nints,1,-1
  alpha(:,int)       = aval + alpha(:,int)
  aval               = alpha(1,int)
  end do
endif

end subroutine

! R^{-1} 2/(b-a) D + I

subroutine freqind_positive(vars,ier,a,b,ifleft,aval0,qfun,pars,nints,ab, &
  alpha,alphader,alphader2)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                         :: vars
procedure(freqind_fun)                       :: qfun
double precision, pointer                    :: pars(:)
double precision, allocatable, intent(out)   :: ab(:,:)
double precision, allocatable, intent(out)   :: alpha(:,:)
double precision, allocatable, intent(out)   :: alphader(:,:)
double precision, allocatable, intent(out)   :: alphader2(:,:)
!
!  Construct a phase function for an equation of the form (1) in the case in which
!  q(t) is positive on the interval [a,b].
!
!  This version of the code adaptively discretizes the Liouville-Green phase to
!  construct a piecewise Chebyshev discretization scheme for representing the phase 
!  function and then constructs the phase function on the interval in parallel.
!
!  The Liouville-Green phase becomes singular as q(t) goes to 0, and this code
!  will encounter difficulty in discretizing it if q(t) gets too close to 0.
!  In this event, it will fail with an error code.
!
!  Input parameters:
!    vars - the structure returned by freqind_init
!    [a,b] - the domain of the differential equation
!    ifleft - an integer indicating whether the value of the phase function
!      alpha is specified at a (ifleft=1) or b (ifleft=0)
!    aval0 - the initial or terminal value of the phase function
!    qfun - an exteranl subroutine conforming to the freqind_fun interface
!      which specifies the coefficient q(t) in (1) as well as its derivative
!    pars - a user-supplied array of doubles which is passed to qfun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0     indicates successful execution
!      ier = 4     means that the maximum number of intervals was exceeded
!                  while forming an adaptive discretization 
!      ier = 8     means that the internal stack used in adaptive discretization
!                  overflowed
!      ier = 128   means that there were no intervals in the high-frequency
!                  regime
!      ier = 256   indicates that the Newton iterations used to solve the 
!                  Riccati equation in the high-frequency regime did not converge
!    nints - the number of intervals in the adaptively determined discretization
!      scheme used to represent the phase function
!    ab - a (2,nints) array whose jth column gives the endpoints of the jth
!      discretization interval
!    alpha -  the values of the phase function at the discretization nodes
!    alphader -  the values of the derivative of the phase function at the
!      discretization nodes
!    alphader2 - the values of the second derivative of the phase function at
!      the discretization nodes
!

double precision, allocatable :: stack(:,:), ab0(:,:), freqs(:)
double precision, allocatable :: ts(:), qs(:,:), qders(:,:), aint(:,:), coefs2(:)
double precision, allocatable :: ws(:,:,:),wders(:,:,:),wder2s(:,:,:),w(:),wder(:),wder2(:)
double precision, allocatable :: w0(:)

integer, allocatable          :: isolve(:)


double complex, allocatable   :: rs(:), rders(:), coefs(:)
double complex                :: ima, rb, rpb
double precision              :: amatr(3,3), rhs(3)
integer                       :: OMP_GET_THREAD_NUM


ier        = 0
ima        = (0.0d0,1.0d0)
k          = vars%k
ntest      = vars%ntest
maxints    = vars%maxints
maxstack   = vars%maxstack
maxnewt    = vars%maxnewt
maxjac     = vars%maxjac
thresh     = vars%thresh
threshlow  = vars%threshlow
epsdisc    = vars%epsdisc
epssol     = vars%epssol
dknorm     = vars%dknorm

allocate(stack(2,maxstack),ts(k), qs(k,maxints), qders(k,maxints))
allocate(rs(k), rders(k), ab0(2,maxints), freqs(maxints)  )
allocate(ws(k,3,maxints),wders(k,3,maxints),wder2s(k,3,maxints),w(k),wder(k))
allocate(aint(k,k),coefs(ntest+1),coefs2(ntest+1),wder2(k),w0(k))

if(ifleft .eq. 1) then
aint=vars%aintl
else
aint=vars%aintr
endif

!
!  Form a discretization scheme for the phase function.
!

nints       = 0
stack(1,1)  = a
stack(2,1)  = b
nstack      = 1
 
do while (nstack .gt. 0) 
  a0       = stack(1,nstack)
  b0       = stack(2,nstack)
  nstack   = nstack-1
  ifaccept = 1
  dr       = 1d300
  dq       = 1d300

  ts = vars%xs * (b0-a0)/2 + (b0+a0)/2
  call qfun(k,ts,qs(:,nints+1),qders(:,nints+1),pars)
  freq  =  (b0-a0) * sqrt(minval(qs(:,nints+1))) 


    !
    !  Check the Liouville-Green solution for "goodness of fit" when the frequency
    !  is high.
    !
    if (freq .gt. thresh) then
      rs      = ima*sqrt(qs(:,nints+1)) - 0.25*qders(:,nints+1)/qs(:,nints+1)
      coefs   = matmul(vars%acoefs2,rs)
      dd1     = abs(coefs(1))
      dd2     = maxval(abs(coefs(2:ntest)))
      dr      = dd2/dd1
      if (dr .gt. epsdisc) ifaccept = 0
    else
      !
      !  In the low-frequency regime, require that the frequency of the interval
      !  is bounded above by threshlow and check the q(t) is well-discretized
      !
      if (freq .gt. threshlow) then
        ifaccept = 0
      else
        coefs2  = matmul(vars%acoefs2,qs(:,nints+1))
        dd1     = abs(coefs2(1))
        dd2     = maxval(abs(coefs2(2:ntest)))
        dq      = dd2/dd1
        if ( dq .gt. epsdisc) ifaccept = 0
      endif

    endif


  if (ifaccept .eq. 0) then
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
   
  else
    if (nints+1 .ge. maxints) then
    ier = 4
    return
    endif
   
    nints           = nints+1
    ab0(1,nints)    = a0
    ab0(2,nints)    = b0
    freqs(nints)    = freq   
  endif
 
end do


!
!  Solve over each interval separately 
!

allocate(ab(2,nints),alpha(k,nints),alphader(k,nints),alphader2(k,nints))
allocate(isolve(nints))

isolve = 0
ab     = ab0(:,1:nints)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(int,a0,b0,freq,jer,ts,rs,rders,coefs,dd1,dd2,dcoefs)
!$OMP DO
do int=1,nints
  a0   = ab(1,int)
  b0   = ab(2,int)
  freq = freqs(int)
  jer  = 0
  
  if(freq .gt. thresh) then

    call freqind_highfreq(vars,jer,epssol,maxnewt,maxjac,a0,b0,k,qs(:,int),qders(:,int), &
      rs,rders)

    alphader(:,int)    = imag(rs) 
    alphader2(:,int)   = imag(rders)
    alpha(:,int)       = (b0-a0)/2 * matmul(aint,alphader(:,int))
    isolve(int)        = 1

     ! coefs   = matmul(vars%acoefs2,alphader(:,int))
     ! dd1     = abs(coefs(1))
     ! dd2     = maxval(abs(coefs(2:ntest)))
     ! dcoefs  = dd2/dd1
     ! print *,a0,b0,dcoefs

  else

    ts = (b0-a0)/2*vars%xs + (b0+a0)/2
    call freqind_appell3(vars,a0,b0,k,ts,qs(:,int),qders(:,int),ws(:,:,int),&
      wders(:,:,int),wder2s(:,:,int))

     ! coefs   = matmul(vars%acoefs2,1/ws(:,2,int))
     ! dd1     = abs(coefs(1))
     ! dd2     = maxval(abs(coefs(2:ntest)))
     ! dcoefs  = dd2/dd1
     ! print *,a0,b0,dcoefs

  endif

  
  if (jer .ne. 0) then
!$OMP ATOMIC 
     ier = ier+1
  endif

end do
!$OMP END DO
!$OMP END PARALLEL


if (ier .gt. 0 ) then
ier = 256
return
endif

if (sum(isolve) == 0) then
ier = 128
return
endif

!
!  Sweep from left to right setting the phase functions for any intervals on which
!  the
!

do int=2,nints

  if (isolve(int) .eq. 0 .AND. isolve(int-1) .eq. 1) then
    a0   = ab(1,int)
    b0   = ab(2,int)

    qval   = qs(k,int-1)
    apval  = alphader(k,int-1)
    appval = alphader2(k,int-1)
    call freqind_atow(qval,apval,appval,wval,wpval,wppval)
    
    ! do ii=1,3
    ! amatr(1,ii) = ws(1,ii,int)
    ! amatr(2,ii) = wders(1,ii,int)
    ! amatr(3,ii) = wder2s(1,ii,int)
    ! end do

    ! rhs(1) = wval
    ! rhs(2) = wpval
    ! rhs(3) = wppval

    ! call freqind_qrsolve(3,amatr,rhs)
    ! c1 = rhs(1)
    ! c2 = rhs(2)
    ! c3 = rhs(3)

    ! w    = c1*ws(:,1,int)    + c2*ws(:,2,int)    + c3*ws(:,3,int)
    ! wder = c1*wders(:,1,int) + c2*wders(:,2,int) + c3*wders(:,3,int)
    
    w    = wval*ws(:,1,int)    + wpval*ws(:,2,int)    + wppval*ws(:,3,int)
    wder = wval*wders(:,1,int) + wpval*wders(:,2,int) + wppval*wders(:,3,int)

    alphader(:,int)    = 1/w
    alphader2(:,int)   = -wder/w**2
    alpha(:,int)       = (b0-a0)/2 * matmul(aint,alphader(:,int))
    isolve(int)        = 1
  endif

end do


do int=nints-1,1,-1

  if (isolve(int) .eq. 0 .AND. isolve(int+1) .eq. 1) then
    a0   = ab(1,int)
    b0   = ab(2,int)

    qval   = qs(1,int+1)
    apval  = alphader(1,int+1)
    appval = alphader2(1,int+1)
    call freqind_atow(qval,apval,appval,wval,wpval,wppval)

    ! ts = (b0-a0)/2*vars%xs+(b0+a0)/2
    ! call freqind_appell_tvp2(vars,a0,b0,k,ts,qs(:,int),qders(:,int),w,wder,wder2,wval,wpval,wppval)

    do ii=1,3
    amatr(1,ii) = ws(k,ii,int)
    amatr(2,ii) = wders(k,ii,int)
    amatr(3,ii) = wder2s(k,ii,int)
    end do

    rhs(1) = wval
    rhs(2) = wpval
    rhs(3) = wppval

    call freqind_qrsolve(3,amatr,rhs)
    c1 = rhs(1)
    c2 = rhs(2)
    c3 = rhs(3)

    w    = c1*ws(:,1,int)    + c2*ws(:,2,int)    + c3*ws(:,3,int)
    wder = c1*wders(:,1,int) + c2*wders(:,2,int) + c3*wders(:,3,int)

    alphader(:,int)    = 1/w
    alphader2(:,int)   = -wder/w**2
    alpha(:,int)       = (b0-a0)/2 * matmul(aint,alphader(:,int))

    isolve(int)        = 1
  endif

end do


!
!  Add the appropriate constants to the phase function
!

if (ifleft .eq. 1) then
  aval = aval0
  do int=1,nints
  alpha(:,int)       = aval + alpha(:,int)
  aval               = alpha(k,int)
  end do
else
  aval = aval0
  do int=nints,1,-1
  alpha(:,int)       = aval + alpha(:,int)
  aval               = alpha(1,int)
  end do
endif

end subroutine



subroutine freqind_phase(vars,ier,a,b,ifleft,aval0,qfun,pars,nints,ab,alpha,alphader, &
  alphader2)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                       :: vars
procedure(freqind_fun)                     :: qfun
double precision, pointer                  :: pars(:)
double precision, allocatable, intent(out) :: ab(:,:)
double precision, allocatable, intent(out) :: alpha(:,:)
double precision, allocatable, intent(out) :: alphader(:,:)
double precision, allocatable, intent(out) :: alphader2(:,:)
!
!  Construct a phase function for an equation of the form (1) in the case in which q(t)
!  is large on at least some portion of [a,b] and any turning points are of ODD order.
!  
!  This is the most robust routine for constructing phase functions.  It is fully
!  adaptive, even testing the accuracy of the phase function in the high-frequency regime,
!  which is almost certainly unncessary.  Moreover, it can handle the case in which
!  q(t) is negative on some part of [a,b], as long q(t) is large and positive somewhere
!  on [a,b] and all turning points are of ODD order.
!
!  Input parameters:
!    vars - the structure returned by freqind_init
!    [a,b] - the solution domain of the differential equation
!    ifleft - an integer flag indicating whether the value of the phase function
!      alpha is specified at the point a (ifleft = 1) or b (ifleft = 0)
!    aval0 - the desired value of alpha at the initial or terminal point of the
!      interval
!    qfun - an exteranl subroutine conforming to the freqind_fun interface
!      which specifies the coefficient q(t) in (1)
!    pars - a user-supplied array of doubles which is passed to qfun
!
!  Output parameters:
!    ier - an error return code;
!      ier = 0     indicates successful execution
!      ier = 4     means that the maximum number of intervals was exceeded
!      ier = 8     means that the internal stack used in adaptive discretization of
!                  q(t) overflowed
!      ier = 256   indicates that the Newton iterations used to solve the 
!                  Riccati equation in the high-frequency regime did not converge
!    nints - the number of intervals in the adaptively determined discretization
!      scheme used to represent the phase function
!    ab - a (2,nints) array whose jth column gives the endpoints of the jth
!      discretization interval
!    alpha -  the values of the phase function at the discretization nodes
!    alphader -  the values of the derivative of the phase function at the
!      discretization nodes
!    alphader2 - the values of the second derivative of the phase function at
!      the discretization nodes
!


integer, allocatable          :: istack(:)
double precision, allocatable :: stack(:,:)
double precision, allocatable :: ts0(:), coefs0(:), qs0(:), qders0(:)
double complex, allocatable   :: coefs1(:), rs0(:)

double precision, allocatable :: ab0(:,:), freqs0(:)
double complex, allocatable   :: rs(:,:), rders(:,:)

double complex, allocatable   :: rs1(:,:), rders1(:,:)
double precision, allocatable :: ab1(:,:), freqs1(:)

double complex, allocatable   :: rs2(:,:), rders2(:,:)
double precision, allocatable :: ab2(:,:), freqs2(:)

double complex                :: ima, rb, rpb, ra, rpa


ier        = 0
eps0       = epsilon(0.0d0)
ima        = (0.0d0,1.0d0)
k          = vars%k
ntest      = vars%ntest
maxints    = vars%maxints
epssol     = vars%epssol
epsdisc    = vars%epsdisc
thresh     = vars%thresh
maxnewt    = vars%maxnewt
maxjac     = vars%maxjac

allocate(stack(2,maxints), istack(maxints))
allocate(ab0(2,maxints),freqs0(maxints))
allocate(rs1(k,maxints),rders1(k,maxints),ab1(2,maxints),freqs1(maxints))
allocate(rs2(k,maxints),rders2(k,maxints),ab2(2,maxints),freqs2(maxints))
allocate(ts0(k),coefs0(ntest+1),qs0(k),qders0(k),coefs1(ntest+1),rs0(k))

!
!  Adaptively discretize q(t) and (possibly also) q'(t) in order to constuct an initial
!  collection of intervals
!

nints0      = 0
stack(1,1)  = a
stack(2,1)  = b
nstack      = 1

do while (nstack .gt. 0) 
a0       = stack(1,nstack)
b0       = stack(2,nstack)
nstack   = nstack-1
ifaccept = 1

ts0 = vars%xs * (b0-a0)/2 + (b0+a0)/2
call qfun(k,ts0,qs0,qders0,pars)

coefs0  = matmul(vars%acoefs2,qs0)
dd1     = maxval(abs(coefs0(2:ntest+1)))
dd2     = abs(coefs0(1))
dcoefs  = dd1/dd2
if (dcoefs .gt. epsdisc) ifaccept = 0


if (ifaccept .eq. 0) then
if (nstack+2 .gt. maxints) then
ier = 8
return
endif

nstack          = nstack+1
stack(1,nstack) = (a0+b0)/2
stack(2,nstack) = b0

nstack          = nstack+1
stack(1,nstack) = a0
stack(2,nstack) = (a0+b0)/2

else

if (nints0+1 .ge. maxints) then
ier = 4
return
endif


nints0             = nints0+1
ab0(1,nints0)      = a0
ab0(2,nints0)      = b0

endif

end do

!
!  Sweep left to right, discretizing q adaptively, solving in all the high-frequency intervals 
!  and over any intervals to the right of a high-frequency interval.
!


nints1  = 0
iright  = -1
ileft   = maxints+1

 
do int0=1,nints0
 a0   = ab0(1,int0)
 b0   = ab0(2,int0)

  nstack             = 1
  stack(1,1)         = a0
  stack(2,1)         = b0
  
  do while( nstack .gt. 0) 
    a0       = stack(1,nstack)
    b0       = stack(2,nstack)
    nstack   = nstack-1

    ts0 = (b0-a0)/2*vars%xs + (b0+a0)/2 
    call qfun(k,ts0,qs0,qders0,pars)
    call freqind_freq(a0,b0,k,qs0,freq)
    ! qders0 = 2/(b0-a0)*matmul(vars%adiff,qs0)

    ifprocessed = 0

    if (freq .gt. thresh) then
       
       call freqind_highfreq(vars,jer,epssol,maxnewt,maxjac,a0,b0,k,qs0,qders0, &
         rs1(:,nints1+1),rders1(:,nints1+1))

       if (jer .eq. 0) ifprocessed=1

    endif
    
    if ( ifprocessed == 0) then


      if ( iright .gt. 0) then

         call freqind_appell_ivp(vars,a0,b0,k,ts0,qs0,qders0,rs1(:,nints1+1), &
          rders1(:,nints1+1),rs1(k,nints1),rders1(k,nints1))

        ifprocessed=1
      endif
    endif

    if (ifprocessed ==0) then
        nints1         = nints1+1
        ab1(1,nints1)  = a0
        ab1(2,nints1)  = b0
        freqs1(nints1) = 0
        cycle
    endif

    ifaccept = 0
    coefs1   = matmul(vars%acoefs2,rs1(:,nints1+1))
    dd1      = maxval(abs(coefs1(2:ntest+1)))    
    dd2      = abs(coefs1(1))
    if (dd2 .eq. 0) dd2=1
    dcoefs   = dd1/dd2

    if (dcoefs .lt. epsdisc) ifaccept = 1


    if (ifaccept .eq. 0) then

      if (nstack+2 .gt. maxints) then
      ier = 8
      return
      endif
      
      nstack          = nstack+1
      stack(1,nstack) = (a0+b0)/2
      stack(2,nstack) = b0
   
      nstack          = nstack+1
      stack(1,nstack) = a0
      stack(2,nstack) = (a0+b0)/2

    else


      if (nints1+1 .gt. maxints) then
      ier = 4
      return
      endif
   
      nints1         = nints1+1
      ab1(1,nints1)  = a0
      ab1(2,nints1)  = b0
      iright = max(iright,nints1)
      ileft  = min(ileft,nints1)

      ! print *,"[",a0,b0,"] accepted",iright,ileft

    endif


  end do

end do


!
!  If there are no high-frequency intervals, abort.
!
if (ileft .gt. maxints) then
ier = 128
return
endif


!
!  Sweep from right to left, starting from the lefttmost unprocessed interval, computing
!  the phase function as we go.
!


nints2 = 0
rb     = rs1(1,ileft)
rpb    = rders1(1,ileft)

do int1=ileft-1,1,-1
  nstack             = 1
  stack(1,1)         = ab1(1,int1)
  stack(2,1)         = ab1(2,int1)

  do while( nstack .gt. 0) 
    a0       = stack(1,nstack)
    b0       = stack(2,nstack)  
    nstack   = nstack-1
    
    ts0 = (b0-a0)/2*vars%xs + (b0+a0)/2 
    call qfun(k,ts0,qs0,qders0,pars)
    ! qders0 = 2/(b0-a0)*matmul(vars%adiff,qs0)
    call freqind_freq(a0,b0,k,qs0,freq)

         ! print *,"----[ ",a0,b0,freq," tvp ] --------------"

    call freqind_appell_tvp(vars,a0,b0,k,ts0,qs0,qders0,rs2(:,nints2+1), &
      rders2(:,nints2+1),rb,rpb)

    ifaccept = 1
    coefs1   = matmul(vars%acoefs2,rs2(:,nints2+1))
    dd1      = maxval(abs(coefs1(2:ntest+1)))
    dd2      = abs(coefs1(1))
    dcoefs   = dd1/dd2
    if (dcoefs .gt. epsdisc) ifaccept = 0
    ! print *,dcoefs,ifaccept

    if (ifaccept .eq. 0) then

    if (nstack+2 .gt. maxints) then
    ier = 8
    return
    endif

    nstack          = nstack+1
    stack(1,nstack) = a0
    stack(2,nstack) = (a0+b0)/2

    nstack          = nstack+1
    stack(1,nstack) = (a0+b0)/2
    stack(2,nstack) = b0
    else

    if (nints1+nints2+1 .gt. maxints) then
    ier = 4
    return
    endif

    nints2         = nints2+1
    ab2(1,nints2)  = a0
    ab2(2,nints2)  = b0

    rb  = rs2(1,nints2)
    rpb = rders2(1,nints2)

    endif

  end do


end do


!
!  Output the phase function
!

nn1   = nints1-ileft+1
nints = nn1+nints2

allocate(ab(2,nints), alpha(k,nints), alphader(k,nints), alphader2(k,nints) )

ab(:,1:nints2)                = ab2(:,nints2:1:-1)
ab(:,nints2+1:nints)          = ab1(:,ileft:nints1)

alphader(:,1:nints2)          = imag(rs2(:,nints2:1:-1))
alphader(:,nints2+1:nints)    = imag(rs1(:,ileft:nints1))

alphader2(:,1:nints2)         = imag(rders2(:,nints2:1:-1))
alphader2(:,nints2+1:nints)   = imag(rders1(:,ileft:nints1))

!
!  Add the appropriate constants to the phase function
!

if (ifleft .eq. 1) then
  aval = aval0
  do int=1,nints
  a0                 = ab(1,int)
  b0                 = ab(2,int)
  alpha(:,int)       = aval + (b0-a0)/2 * matmul(vars%aintl,alphader(:,int))
  aval               = alpha(k,int)
  end do
else
  aval = aval0
  do int=nints,1,-1
  a0                 = ab(1,int)
  b0                 = ab(2,int)
  alpha(:,int)       = aval + (b0-a0)/2 * matmul(vars%aintr,alphader(:,int))
  aval               = alpha(1,int)
  end do
endif

end subroutine


subroutine freqind_freq(a,b,k,qs,freq)
implicit double precision (a-h,o-z)
double precision          :: qs(k)

dmin = minval(qs)

if (dmin .gt. 0) then
  freq   =  (b-a) * sqrt(dmin)
else
   freq  = -1
endif

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Phase function evaluation routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine freqind_eval1(vars,nints,ab,alpha,x,aval)
implicit double precision (a-h,o-z)
type(freqind_vars_t)             :: vars
double precision                 :: ab(:,:)
double precision                 :: alpha(:,:)
!
!  Evaluate a phase function alpha constructed by one of the routines in this module
!  at a point x on [a,b].
!
!  Input parameters:
!    vars - the structure initialized by freqind_init
!    nints - the number of intervals in the piecewise Chebyshev scheme representing
!      the phase function
!    ab - an (2,nints) array whose jth column gives the endpoints of the
!     jth discretization interval
!    alpha - the values of the phase function at the discretization nodes
!    x - the point at which to evaluate the functions
!
!  Output parameters:
!    aval - the value of alpha(x)
!
!
!
double precision                 :: whts(vars%k)

eps0 = epsilon(0.0d0)
k    = vars%k

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
do i=1,k
  dd      = xx - vars%xs(i)
  if ( abs(dd) .le. eps0) then
  aval = alpha(i,int)
  return
  endif
  whts(i) = vars%dsigns(i)/dd
end do


aval    = dot_product(whts,alpha(:,int)) / sum(whts)

!aval    = sum(whts*alpha(:,int)) / sum(whts)

end subroutine


subroutine freqind_eval2(vars,nints,ab,alpha,alphader,x,aval,apval)
implicit double precision (a-h,o-z)
type(freqind_vars_t)             :: vars
double precision                 :: ab(:,:)
double precision                 :: alpha(:,:)
double precision                 :: alphader(:,:)

!
!  Evaluate a phase function alpha constructed by one of the routines in this module
!  and its first derivative at a point x on [a,b].
!
!  Input parameters:
!    vars - the structure initialized by freqind_init
!    nints - the number of intervals in the piecewise Chebyshev scheme representing
!      the phase function
!    ab - an (2,nints) array whose jth column gives the endpoints of the
!     jth discretization interval
!    alpha - the values of alpha(x) at the discretization nodes
!    alphader - the values of alpha'(x) at the discretization nodes
!    x - the point at which to evaluate the functions
!
!  Output parameters:
!    aval - the value of alpha(x)
!    apval - the value of alpha'(x)
!
!
double precision                 :: whts(vars%k)

eps0 = epsilon(0.0d0)
k    = vars%k

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
do i=1,k
  dd      = xx - vars%xs(i)
  if ( abs(dd) .le. eps0) then
  aval  = alpha(i,int)
  apval = alphader(i,int)
  return
  endif
  whts(i) = vars%dsigns(i)/dd
end do


dwhts   = sum(whts)
aval    = dot_product(whts,alpha(:,int)) / dwhts
apval   = dot_product(whts,alphader(:,int)) / dwhts

! aval    = sum(whts*alpha(:,int))    / dwhts
! apval    = sum(whts*alphader(:,int)) / dwhts

end subroutine


subroutine freqind_eval3(vars,nints,ab,alpha,alphader,alphader2,x,aval,apval,appval)
implicit double precision (a-h,o-z)
type(freqind_vars_t)             :: vars
double precision                 :: ab(:,:)
double precision                 :: alpha(:,:)
double precision                 :: alphader(:,:)
double precision                 :: alphader2(:,:)
!
!  Evaluate a phase function alpha constructed by one of the routines in this module
!  and its first derivative at a point x on [a,b].
!
!  Input parameters:
!    vars - the structure initialized by freqind_init
!    nints - the number of intervals in the piecewise Chebyshev scheme representing
!      the phase function
!    ab - an (2,nints) array whose jth column gives the endpoints of the
!     jth discretization interval
!    alpha - the values of alpha(x) at the discretization nodes
!    alphader - the values of alpha'(x) at the discretization nodes
!    alphader2 - the values of alpha''(x) at the discretization nodes
!    x - the point at which to evaluate the functions
!
!  Output parameters:
!    aval - the value of alpha(x)
!    apval - the value of alpha'(x)
!    appval - the value of alpha''(x)
!
double precision                 :: whts(vars%k)

eps0 = epsilon(0.0d0)
k    = vars%k

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
do i=1,k
  dd      = xx - vars%xs(i)

  if ( abs(dd) .le. eps0) then
  aval  = alpha(i,int)
  apval = alphader(i,int)
  appval = alphader2(i,int)
  return
  endif

  whts(i) = vars%dsigns(i)/dd
end do

dwhts  = sum(whts)

aval    = dot_product(whts,alpha(:,int))     / dwhts
apval   = dot_product(whts,alphader(:,int))  / dwhts
appval  = dot_product(whts,alphader2(:,int)) / dwhts

! aval    = sum(whts*alpha(:,int))     / dwhts
! apval   = sum(whts*alphader(:,int))  / dwhts
! appval  = sum(whts,alphader2(:,int)) / dwhts

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Local solver routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine freqind_highfreq(vars,ier,eps,maxnewt,maxiters,a,b,k,qs,qders,rs,rders)
implicit double precision (a-h,o-z)
type(freqind_vars_t)         :: vars
double precision             :: qs(k), qders(k)
double complex               :: rs(k), rders(k)
!
!  Solve the equation
!
!    r'(t) + (r(t))^2 + q(t) = 0
!
!  over the interval [a,b] when q is positive and of large magnitude.  This routine
!  uses Newton-Kantorovich iterations and the linearized equations are solved via
!  Jacobi iterations.
!

! double precision, allocatable :: dnorms(:)
! double complex, allocatable   :: amatr(:,:), res(:), delta(:), rder2s(:)
! double complex, allocatable   :: qders(:), qder2s(:)

double precision              :: dnorms(maxiters)
double complex                :: amatr(k,k), res(k), delta(k), bmatr(k,k)
double complex                :: delta1(k)
double complex                :: ima

ier      = 0
eps0     = epsilon(0.0d0)
ima      = (0.0d0,1.0d0)

! use the Liouville-Green phase to initialize
amatr    = 2/(b-a) * vars%adiff
rs       = ima*sqrt(qs) - 0.25*qders/qs
rders    = matmul(amatr,rs)

nextra   = 0
ifexit   = 0


do inewt = 1,maxnewt
  res      = -(rders+rs**2+qs)
  drs      = maxval(abs(rs))

  
  ! Just do two steps of the fixed point iteration
  res    = res / (2*rs)
  delta  = res - matmul(amatr,res)/(2*rs)
  
  ! Allow for more iterations

  ! res    = res/(2*rs)
  ! dres     = maxval(abs(res))
  ! epsjac = max(eps,dres/(100*maxval(abs(rs))))
  ! delta  = res

  ! do iter=1,maxiters
  ! delta1  = res - matmul(amatr,delta)/(2*rs)
  ! diff    = maxval(abs(delta-delta1))
  ! delta   = delta1
  ! if (diff .lt. epsjac) exit
  ! end do

!!!!! SOLVE VIA A QR DECOMPOSITION !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! bmatr = amatr
  ! do i=1,k
  ! bmatr(i,i) = bmatr(i,i)+2*rs(i)
  ! end do
  ! delta=res
  ! call linalg0_solve(k,bmatr,delta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  !  Update the solution
  !
  rs            = rs + delta
  rders         = matmul(amatr,rs)
  dd            = maxval(abs(delta)) / drs
  
  if (dd  .lt. eps) exit

end do

! if (inewt .gt. maxnewt) inewt=inewt-1
! if (dnorms(inewt) .gt. eps) then
! ier = 4
! endif   


end subroutine




subroutine freqind_appell_ivp(vars,a,b,k,ts,qs,qders,rs,rders,ra,rpa)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                :: vars
double precision                    :: ts(k), qs(k), qders(k)
double complex                      :: rs(k), rders(k)
double complex                      :: ra, rpa
!
!  Solve an initial value problem for Appell's equation 
!
!    w'''(t) + 4 q(t) w(t) + 2 q'(t) w(t) = 0,
!
!  where the initial values w(a), w'(a) and w''(a) are specified
!  via the values of a sollution of Riccati's equation.
!

double precision                    :: amatr(k,k), bmatr(k,k)
double precision                    :: rhs(k), sigma(k)
double precision                    :: ws(k), wders(k),wder2s(k)

double complex                      :: ima

real*16                             :: q1,q2,q3,q4,q5,q6,q7

ima    = (0.0d0,1.0d0)

!
!  find the initial values
!
call  freqind_rtow(ra,rpa,wa,wpa,wppa)


!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wa + wpa (t-a) + wppa/2 (t-a)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  a
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        a                         a
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        a
!
!  where
!
!    rhs(t) = - wppa r(t) - ( wpa + wppa(t-a) )  p(t) - (wa + wpa (t-a) + wppa (t-a)^2/2) q(t)
!
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,k
amatr(i,i) = 1.0d0
end do

!
! Handle the p(t) * w'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + 4*qs(i) * vars%aintl2(i,:)*dd2
sigma(i)   = sigma(i) - 4*qs(i)*(wpa + wppa * (ts(i)-a))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + 2*qders(i) * vars%aintl3(i,:)*dd3
sigma(i)   = sigma(i) - 2*qders(i) * (wa + wpa*(ts(i)-a) + wppa * (ts(i)-a)**2/2)
end do

! Solve the system
call freqind_qrsolve(k,amatr,sigma)

! Construct the solutions
wder2s = wppa + dd*matmul(vars%aintl,sigma)
wders  = wpa  + dd*matmul(vars%aintl,wder2s)
ws     = wa   + dd*matmul(vars%aintl,wders)

rs     = ima/ws + 0.5d0*wders/ws
rders  = -ima*wders/ws**2 - wders**2/(2*ws**2) + wder2s/(2*ws)

end subroutine


subroutine freqind_rtow(r,rp,w,wp,wpp)
implicit double precision (a-h,o-z)
double complex           :: r, rp
real*16                  :: q1, q2, q3

w   = 1 / imag(r)

q1  = 2*real(r)*w
wp  = q1

q1  = 2*real(rp)*w
q2  = wp*wp/w
q3  = q1+q2
wpp = q3

end subroutine


subroutine freqind_atow(qval,apval,appval,w,wp,wpp)
implicit double precision (a-h,o-z)

apppval = ( 4*qval*apval**2-4*apval**4+3*appval**2) / (2*apval)
w      = 1.0d0/apval
wp     = -appval/apval**2
wpp    = 2*appval**2/apval**3 - apppval/apval**2

end subroutine



subroutine freqind_appell_ivp2(vars,a,b,k,ts,qs,qders,ws,wders,wder2s,wa,wpa,wppa)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                :: vars
double precision                    :: ts(k), qs(k), qders(k)
double precision                    :: ws(k), wders(k), wder2s(k)
!
!  Solve an initial value problem for Appell's equation 
!
!    w'''(t) + 4 q(t) w(t) + 2 q'(t) w(t) = 0,
!
!  where the initial values w(a), w'(a) and w''(a) are specified directly.
!
double precision                    :: amatr(k,k), bmatr(k,k)
double precision                    :: rhs(k), sigma(k)

double complex                      :: ima


ima    = (0.0d0,1.0d0)

!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wa + wpa (t-a) + wppa/2 (t-a)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  a
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        a                         a
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        a
!
!  where
!
!    rhs(t) = - wppa r(t) - ( wpa + wppa(t-a) )  p(t) - (wa + wpa (t-a) + wppa (t-a)^2/2) q(t)
!
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,k
amatr(i,i) = 1.0d0
end do

!
! Handle the p(t) * w'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + 4*qs(i) * vars%aintl2(i,:)*dd2
sigma(i)   = sigma(i) - 4*qs(i)*(wpa + wppa * (ts(i)-a))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + 2*qders(i) * vars%aintl3(i,:)*dd3
sigma(i)   = sigma(i) - 2*qders(i) * (wa + wpa*(ts(i)-a) + wppa * (ts(i)-a)**2/2)
end do

! Solve the system
call freqind_qrsolve(k,amatr,sigma)

! Construct the solutions
wder2s = wppa + dd*matmul(vars%aintl,sigma)
wders  = wpa  + dd*matmul(vars%aintl,wder2s)
ws     = wa   + dd*matmul(vars%aintl,wders)

end subroutine



subroutine freqind_appell_tvp(vars,a,b,k,ts,qs,qders,rs,rders,rb,rpb)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                :: vars
double precision                    :: ts(k), qs(k), qders(k)
double complex                      :: rs(k), rders(k), rb, rpb
!
!  Solve a termainl value problem for Appell's equation 
!
!    w'''(t) + 4 q(t) w(t) + 2 q'(t) w(t) = 0,
!
!  where the terminal values w(b), w'(b) and w''(b) are specified via a
!  solution of Riccati's equation.
!
double precision                    :: ps0(k), qs0(k)
double precision                    :: ws(k), wders(k), wder2s(k)
double precision                    :: amatr(k,k), bmatr(k,k)
double precision                    :: rhs(k), sigma(k)

double complex                      :: ima

ps0    = 4*qs
qs0    = 2*qders
ima    = (0.0d0,1.0d0)

!
!  find the initial values
!

call  freqind_rtow(rb,rpb,wb,wpb,wppb)


!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wb + wpb (t-b) + wppb/2 (t-b)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  b
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        b                         b
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        b
!
!  where
!
!    rhs(t) = - wppb r(t) - ( wpb + wppb(t-b) )  p(t) - (wb + wpb (t-b) + wppb (t-b)^2/2) q(t)
!
!


!
!  Allocate memory for the procedure and setup some parameters.
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,k
amatr(i,i) = 1.0d0
end do


!
! Handle the p(t) * w'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps0(i) * vars%aintr2(i,:)*dd2
sigma(i)   = sigma(i) - ps0(i)*(wpb + wppb * (ts(i)-b))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs0(i) * vars%aintr3(i,:)*dd3
sigma(i)   = sigma(i) - qs0(i) * (wb + wpb*(ts(i)-b) + wppb * (ts(i)-b)**2/2)
end do

! Solve the system
call freqind_qrsolve(k,amatr,sigma)

! Construct the solutions
wder2s = wppb + dd*matmul(vars%aintr,sigma)
wders  = wpb  + dd*matmul(vars%aintr,wder2s)
ws     = wb   + dd*matmul(vars%aintr,wders)

rs     = ima/ws + 0.5d0*wders/ws
rders  = -ima*wders/ws**2 - wders**2/(2*ws**2) + wder2s/(2*ws)

end subroutine


subroutine freqind_appell_tvp2(vars,a,b,k,ts,qs,qders,ws,wders,wder2s,wb,wpb,wppb)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                :: vars
double precision                    :: ts(k), qs(k), qders(k)
double precision                    :: ws(k), wders(k), wder2s(k)

double precision                    :: ps0(k), qs0(k)
double precision                    :: amatr(k,k), bmatr(k,k)
double precision                    :: rhs(k), sigma(k)

!
!  Solve a termainl value problem for Appell's equation 
!
!    w'''(t) + 4 q(t) w(t) + 2 q'(t) w(t) = 0,
!
!  where the terminal values w(a), w'(a) and w''(a) are specified directly.
!

double complex                      :: ima

ps0    = 4*qs
qs0    = 2*qders
ima    = (0.0d0,1.0d0)


!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wb + wpb (t-b) + wppb/2 (t-b)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  b
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        b                         b
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        b
!
!  where
!
!    rhs(t) = - wppb r(t) - ( wpb + wppb(t-b) )  p(t) - (wb + wpb (t-b) + wppb (t-b)^2/2) q(t)
!
!


!
!  Allocate memory for the procedure and setup some parameters.
!

amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,k
amatr(i,i) = 1.0d0
end do

!
! Handle the p(t) * w'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + ps0(i) * vars%aintr2(i,:)*dd2
sigma(i)   = sigma(i) - ps0(i)*(wpb + wppb * (ts(i)-b))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + qs0(i) * vars%aintr3(i,:)*dd3
sigma(i)   = sigma(i) - qs0(i) * (wb + wpb*(ts(i)-b) + wppb * (ts(i)-b)**2/2)
end do

! Solve the system
call freqind_qrsolve(k,amatr,sigma)


! Construct the solutions
wder2s = wppb + dd*matmul(vars%aintr,sigma)
wders  = wpb  + dd*matmul(vars%aintr,wder2s)
ws     = wb   + dd*matmul(vars%aintr,wders)

end subroutine


subroutine freqind_appell3(vars,a,b,k,ts,qs,qders,ws,wders,wder2s)
implicit double precision (a-h,o-z)
type(freqind_vars_t)                :: vars
double precision                    :: ts(k), qs(k), qders(k)
double precision                    :: ws(k,3), wders(k,3), wder2s(k,3)
!
!  Return three solutions w_1,w_2,w_3 of Appell's differential equation
!
!    w'''(t) + 4 q(t) w(t) + 2 q'(t) w(t) = 0
!
!  on the interval [a,b] such that
!
!    w_1(a) = 1, w_1'(a) = 0, w_1''(a) = 0,
!
!    w_2(a) = 0, w_2'(a) = 1, w_2''(a) = 1   and
!
!    w_3(a) = 0, w_3'(a) = 0, w_3''(a) = 1.
!

double precision                    :: amatr(k,k), wder3s(k,3)
double precision                    :: was(3), wpas(3), wppas(3)

!
!  Use the representation of the solution
! 
!                                                  t
!    w(t) = wa + wpa (t-a) + wppa/2 (t-a)^2 + \int   (t-s)^2/2 sigma(s) ds
!                                                  a
!
!  to solve and initial value problem for (1).  This representation leads to the
!  integral equation
!
!                        t                         t
!    sigma(t) + r(t) \int   sigma(s) ds + p(t) \int  (t-s) sigma(s) ds +
!                        a                         a
!
!                        t
!              q(t)  \int   (t-s)^2/2 sigma(s) ds = rhs(t)
!                        a
!
!  where
!
!    rhs(t) = - wppa r(t) - ( wpa + wppa(t-a) )  p(t) - (wa + wpa (t-a) + wppa (t-a)^2/2) q(t)
!
!

!
!  Construct the discretized integral operator
!
amatr  = 0
sigma  = 0

dd     = (b-a)/2 
dd2    = dd*dd
dd3    = dd*dd2

do i=1,k
amatr(i,i) = 1.0d0
end do

!
! Handle the p(t) * w'(t) term.
!
do i=1,k
amatr(i,:) = amatr(i,:) + 4*qs(i) * vars%aintl2(i,:)*dd2
! sigma(i)   = sigma(i) - 4*qs(i)*(wpa + wppa * (ts(i)-a))
end do

!
!  Handle the q(t) w(t) term
!

do i=1,k
amatr(i,:) = amatr(i,:) + 2*qders(i) * vars%aintl3(i,:)*dd3
! sigma(i)   = sigma(i) - 2*qders(i) * (wa + wpa*(ts(i)-a) + wppa * (ts(i)-a)**2/2)
end do


!
!  Solve three systems involving it
!

was(1)   = 1
wpas(1)  = 0
wppas(1) = 0

was(2)   = 0
wpas(2)  = 1
wppas(2) = 0

was(3)   = 0
wpas(3)  = 0
wppas(3) = 1

do idx=1,3
do i=1,k
wder3s(i,idx)   =             - 4*qs(i)*(wpas(idx) + wppas(idx) * (ts(i)-a))
wder3s(i,idx)   = wder3s(i,idx) - 2*qders(i) * (was(idx) + wpas(idx)*(ts(i)-a) + &
                  wppas(idx) * (ts(i)-a)**2/2)
end do
end do

! Solve three initial value problems
call freqind_qrsolve2(k,amatr,wder3s)
do idx=1,3
wder2s(:,idx) = wppas(idx) + dd*matmul(vars%aintl,wder3s(:,idx))
wders(:,idx)  = wpas(idx)  + dd*matmul(vars%aintl,wder2s(:,idx))
ws(:,idx)     = was(idx)   + dd*matmul(vars%aintl,wders(:,idx))
end do

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Linear algebra routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine freqind_qrsolve(n,a,rhs)
implicit double precision (a-h,o-z)
double precision   :: a(:,:),rhs(:)
!
!  This subroutine uses a version of QR-decomposition to solve the equation
!  A x = b.  
!
!  THE MATRIX A IS DESTROYED BY THIS ROUTINE AND THE INPUT RHS
!  IS OVERWRITTEN.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - a vector of length n speciying the rhs of the system
!
!  Output parameters:
!
!   rhs - upon return, the solution of the linear system


double precision :: aa(2),u(2,2)
! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

!
!  Reduce to upper triangular 
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1)=d1
a(ii,j)=d2
end do
d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix
! 

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

return
end subroutine


subroutine freqind_qrsolve2(n,a,rhs)
implicit double precision (a-h,o-z)
double precision   :: a(:,:),rhs(:,:)
!
!  This subroutine forms a QR decomposition of an (n,n) matrix A
!  using Givens rotations and uses it solve the system of linear
!  equations
!
!    Ax = b
!
!  for several right-hand sides b.
!
!  THE MATRIX A IS DESTROYED BY THIS ROUTINE AND THE INPUT RHS
!  IS OVERWRITTEN.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - an (n,*) array speciying the right-hand sides
!
!  Output parameters:
!   rhs - upon return, the solutions of the linear systems


double precision :: aa(2),u(2,2)
double precision :: dd1(1,size(rhs,2)), dd2(1,size(rhs,2))

! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

!
!  Reduce to upper triangular via Givens rotations
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)
u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1) = d1
a(ii,j)   = d2
end do

dd1(1,:)=u(1,1)*rhs(j-1,:)+u(1,2)*rhs(j,:)
dd2(1,:)=u(2,1)*rhs(j-1,:)+u(2,2)*rhs(j,:)
rhs(j-1,:) = dd1(1,:)
rhs(j,:)   = dd2(1,:)

end do
end do

!
!  Apply the inverse of the triangular matrix
! 

rhs(n,:)=rhs(n,:)/a(n,n)
do i=n-1,1,-1
dd1=0
do j=n,i+1,-1
dd1(1,:)=dd1(1,:)+a(j,i)*rhs(j,:)
end do
rhs(i,:)=(rhs(i,:)-dd1(1,:))/a(i,i)
end do

end subroutine

end module
