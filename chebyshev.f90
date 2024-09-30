!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing Curtis-Clenshaw quadrature rules, and for
!  forming and manipulating univariate Chebyshev expansions.  A univariate Chebyshev
!  expansion is a sum of the form 
!
!              n-1
!      f(x) = \sum   a_i T_i(x),                                                            (1)
!              i=0
!
!  where T_i denotes the Chebyshev polynomial of degree i.  This code provides routines
!  for both real-valued and complex-valued expansions.
!
!  Expansions of the form (1) are represented either via the vector of coefficients
!   
!    ( a_0     )
!    ( a_1     )
!    ( ...     )                                                                            (2)
!    ( a_{n-1} )
!
!  or via the vector 
!
!    ( f(x_1)  )
!    ( f(x_2)  )
!    (  ...    )                                                                            (3)
!    ( f(x_n)  )
!
!  of their values at the "practical" Chebyshev nodes.
!
!  The following routines should be regarded as public:
!      
!    chebyshev_quad - return the nodes and weights of an n-point Curtis-Clenshaw
!      quadrature rule whose nodes are the n-point Chebyshev extrema grid and which
!      integrates polynomials of degrees 0 through n-1 exactly
!
!      ***NOTE*** This routine's running time grows as O(n^2) and a different
!      approach should be used when n is large
!
!    chebyshev_eval - evaluate one or more expansions of the form (1) at a specified
!       point given their coefficient vectors (2)
!
!    chebyshev_evalder - evaluate one or more expansion of the form (1) and their
!      derivatives at a specified point given the vector (2) of coefficients
!
!    chebyshev_interp - evaluate one or more expansions of the form (1) at a specified
!       point given their vectors of values at the Chebyshev extrema grid (3)
!
!    chebyshev_coefsmatrix - return the (n,n) matrix which takes the vector (3) to
!      the vector (2)
!
!    chebyshev_interpmatrix - return an (m,n) matrix which takes the vector of values
!      (3) of a Chebyshev expansion to its values at a user-specified set of points
!      on the interval [-1,1]
!      
!    chebyshev_diffmatrix - return the (n,n) ``spectral differentiation matrix'' which
!      takes the vector (3) of values of an expansion of the form (1) to the vector
!      of values of the derivative of the expansion
!
!    chebyshev_intlmatrix - return the (n,n) spectral integration matrix which
!      takes the vector (3) of values of an expansion of the form (1) to the 
!      the vector of values of 
!
!                        x
!            g(x)  = \int f(x) dt
!                       -1
!
!      at the Curtis-Clenshaw nodes 
!
!    chebyshev_intrmatrix - return the (n,n) spectral integration matrix which
!      takes the vector (3) of values of an expansion of the form (1) to the 
!      the vector of values of 
!
!                        x
!            g(x)  = \int f(x) dt
!                        1
!
!      at the Curtis-Clenshaw nodes 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     chebyshev

interface     chebyshev_interp
module procedure chebyshev_interp1
module procedure chebyshev_interp2
module procedure chebyshev_interp3
module procedure chebyshev_interp4
end interface chebyshev_interp

interface     chebyshev_eval
module procedure chebyshev_eval1
module procedure chebyshev_eval2
module procedure chebyshev_eval3
module procedure chebyshev_eval4
end interface chebyshev_eval

interface     chebyshev_evalder
module procedure chebyshev_evalder1
module procedure chebyshev_evalder2
module procedure chebyshev_evalder3
module procedure chebyshev_evalder4
end interface chebyshev_evalder

contains



subroutine chebyshev_quad(n,xscheb,whtscheb)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)          :: xscheb(:), whtscheb(:)
!
!  Construct a Curtis-Clenshaw quadrature rule whose nodes are the n-point
!  Chebyshev extrema grid and which integrates polynomials of degrees less than
!  or equal to n-1 exactly.
!
!  Input parameters:
!    n - the length of the Curtis-Clenshaw quadrature rule to construct
!
!  Output parameters:
!    xscheb - an array containing the quadrature nodes
!    whtscheb - an array containing the quadrature weights
!
data pi / 3.14159265358979323846264338327950288d0 /
allocate( xscheb(1:n), whtscheb(1:n) )


!
!  Construct the nodes ...
!

do j=1,n
xscheb(n-j+1) = cos( pi * (j-1.0d0) / (n-1.0d0) )
end do

!
!  ... and now the weights
!
do i=1,n
whtscheb(n-i+1) = 2.0d0/(n-1.0d0)
do j=3,n,2
dd          = 1.0d0/(j+0.0d0)-1.0d0/(j-2.0d0)
xx = pi * (j-1.0d0) / (n-1.0d0)
uu = cos((i-1)*xx)*2/(n-1.0d0)
if (j .eq. n) uu=uu/2
whtscheb(n-i+1) = whtscheb(n-i+1) + dd*uu
end do
end do
whtscheb(1) = whtscheb(1)/2
whtscheb(n) = whtscheb(n)/2

end subroutine


subroutine chebyshev_coefsmatrix(n,ucheb)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: ucheb(:,:)
!
!  Return the (n,n) matrix which takes the vector (3) of values of the
!  expansion (1) at the n-point Chebyshev extrema grid to the vector (2)
!  of expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!
!  Output parameters:
!    ucheb - the (n,n) values-to-coefficients matrix
!
data pi / 3.14159265358979323846264338327950288d0 /

allocate(ucheb(n,n))

do j=1,n
jj = n-j+1
xx = pi * (jj-1.0d0) / (n-1.0d0)
do i=1,n
ucheb(i,j) = cos((i-1)*xx)
end do
end do

ucheb      = ucheb * 2.0d0/(n-1.0d0)
ucheb(1,:) = ucheb(1,:)/2
ucheb(n,:) = ucheb(n,:)/2
ucheb(:,1) = ucheb(:,1)/2
ucheb(:,n) = ucheb(:,n)/2

end subroutine


subroutine chebyshev_diffmatrix(n,adiff)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: adiff(:,:)
!
!  Return the (n,n) spectral differentiation matrix which takes the vector (3) of values 
!  of the expansion (1) at the n-point Chebyshev extrema grid to the vector of values
!  of its derivative at the same node.
!
!  Input parameters:
!    n - the number of terms in the expansion
!
!  Output parameters:
!    adiff - the (n,n) values-to-coefficients matrix
!
data pi / 3.14159265358979323846264338327950288d0 /

double precision, allocatable :: u(:,:), vals(:,:), pols(:), ders(:)
allocate(adiff(n,n), u(n,n), vals(n,n), pols(0:n-1), ders(0:n-1) )

call chebyshev_coefsmatrix(n,u)

do i=1,n
ii  = n-i+1
x   = cos( pi * (ii-1.0d0) / (n-1.0d0) )
call chebders(n,x,pols,ders)
do j=1,n
vals(i,j) = ders(j-1)
end do
end do

adiff = matmul(vals,u)

end subroutine


subroutine chebyshev_interpmatrix(n,m,xsout,ainterp)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)  :: ainterp(:,:)
double precision                            :: xsout(:)
!
!  Construct the (m,n) interpolation matrix which takes the vector (3) of values of an 
!  expansion of the form (1) at the Curtis-Clenshaw nodes to its values at a
!  a collection of user-specified nodes.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    m - the number of target nodes
!    xsout - the target nodes
!
!  Output parameters:
!    ainterp - the (m,n) interpolation matrix
!
data pi / 3.14159265358979323846264338327950288d0 /

double precision, allocatable :: whts(:)

eps0 = epsilon(0.0d0)

allocate(ainterp(m,n), whts(n))

ainterp = 0

do i=1,m
x      = xsout(i)
ifdone = 0

dsign = 1.0d0
do j=1,n
jj     = n-j+1
xx      = cos(pi*(jj-1)/(n-1.0d0))
dd      = x - xx

if ( abs(dd) .le. eps0) then
ainterp(i,j) = 1.0d0
ifdone = 1
exit
endif

whts(j) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(n) = whts(n)/2

if (ifdone .eq. 0) then
dsum0 = sum(whts)

do j=1,n
ainterp(i,j) = whts(j)/dsum0
end do
endif

end do

end subroutine


subroutine chebyshev_intlmatrix(n,aint)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) :: aint(:,:)
!
!  Construct the (n,n) "left" spectral integration matrix which takes the values
!  of the an expansion f(x) of the form (1) at the Curtis-Clenshaw nodes
!  nodes to the values of the function
!
!                x
!    g(x) = \int   f(t) dt
!               -1
!
!   at the same nodes.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!
!  Output parameters:
!    aint - the (n+1,n) spectral integration matrix
!
data pi / 3.14159265358979323846264338327950288d0 /
double precision, allocatable :: u(:,:), vals(:,:), pols(:)

allocate(aint(n,n), vals(n,n), u(n,n), pols(1:n+1) )

call chebyshev_coefsmatrix(n,u)

do i=1,n
ii  = n-i+1
x   = cos( pi * (ii-1.0d0) / (n-1.0d0) )
call chebs(n+1,x,pols)

vals(i,1) = (x+1)
vals(i,2) = (x**2-1.0d0)/2

do j=3,n
val1 = (j-1)*pols(j+1)/((j-1.0d0)**2-1.0d0) - x*pols(j)/(j-2.0d0)
val2 = (j-1)*(-1.0d0)**j/((j-1.0d0)**2-1.0d0)+(-1.0d0)**(j-1)/(j-2.0d0)
vals(i,j) = val1-val2
end do
end do

aint = matmul(vals,u)

end subroutine

subroutine chebyshev_intrmatrix(n,aint)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) :: aint(:,:)
!
!  Construct the (n,n) "right" spectral integration matrix which takes the values
!  of the an expansion f(x) of the form (1) at the Curtis-Clenshaw nodes
!  nodes to the values of the function
!
!                x
!    g(x) = \int   f(t) dt
!                1
!
!   at the same nodes.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!
!  Output parameters:
!    aint - the (n+1,n) spectral integration matrix
!
data pi / 3.14159265358979323846264338327950288d0 /
double precision, allocatable :: u(:,:), vals(:,:), pols(:)

allocate(aint(n,n), vals(n,n), u(n,n), pols(1:n+1) )

call chebyshev_coefsmatrix(n,u)

do i=1,n
ii  = n-i+1
x   = cos( pi * (ii-1.0d0) / (n-1.0d0) )
call chebs(n+1,x,pols)

vals(i,1) = (x-1)
vals(i,2) = x**2/2-0.5d0

do j=3,n
val1 = (j-1)*pols(j+1)/((j-1.0d0)**2-1.0d0) - x*pols(j)/(j-2.0d0)
val2 = (j-1.0d0)/((j-1.0d0)**2-1.0d0)-1.0d0/(j-2.0d0)
vals(i,j) = val1-val2
end do
end do

aint = matmul(vals,u)

end subroutine


subroutine chebyshev_interp1(n,xscheb,vals,x,valout)
implicit double precision (a-h,o-z)
double complex             :: vals(:), valout
double precision           :: xscheb(n)
data pi / 3.14159265358979323846264338327950288d0 /
!
!  Use barycentric Lagrange interpolation to evaluate a real-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    vals - the values of (1) at the n Curtis-Clenshaw nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double complex   :: dsum1

eps0  = epsilon(0.0d0)

dsign = 1.0d0
do i=1,n
!ii      = n-i+1
!xx      = cos(pi*(ii-1)/(n-1.0d0))
xx      = xscheb(i)
dd      = x - xx
if ( abs(dd) .le. eps0) then
valout = vals(i)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(n) = whts(n)/2

dsum1 = sum(whts*vals(1:k))
dsum2 = sum(whts)

valout = dsum1/dsum2

end subroutine


subroutine chebyshev_interp2(n,xscheb,vals,x,valout)
implicit double precision (a-h,o-z)
double precision           :: vals(:), valout
double precision           :: xscheb(n)
data pi / 3.14159265358979323846264338327950288d0 /
!
!  Use barycentric Lagrange interpolation to evaluate a real-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    vals - the values of (1) at the n Curtis-Clenshaw nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double precision :: dsum1

eps0  = epsilon(0.0d0)

dsign = 1.0d0
do i=1,n
!ii      = n-i+1
!xx      = cos(pi*(ii-1)/(n-1.0d0))
xx      = xscheb(i)
dd      = x - xx
if ( abs(dd) .le. eps0) then
valout = vals(i)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(n) = whts(n)/2

dsum1 = 0
dsum2 = 0

do i=1,n
dsum1 = dsum1 + whts(i) * vals(i)
dsum2 = dsum2 + whts(i)
end do

valout = dsum1/dsum2

end subroutine

subroutine chebyshev_interp3(n,xscheb,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xscheb(n)
double precision           :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of real-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Curtis-Clenshaw
!      quadrature rule
!    vals - a matrix whose jth column gives the values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
!
double precision              :: whts(n)
!double precision, allocatable :: sums(:)


l = size(vals,2)

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xscheb(i)
if ( abs(dd) .le. eps0) then
valsout(1:l) = vals(i,1:l)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(n) = whts(n)/2

! compute the sums
! do i=1,l
! vals(1:n,i) = vals(1:n,i) * whts(1:n) 
! sums(i)     = sum (vals(:,i))
! vals(1:n,i) = vals(1:n,i) / whts(1:n)
! end do

do j=1,l
valsout(j) = 0
do i=1,n
valsout(j) = valsout(j) + vals(i,j) * whts(i)
end do
end do


! evaluate the formula
whts    = whts
valsout = valsout/sum(whts)

end subroutine


subroutine chebyshev_interp4(n,xscheb,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xscheb(n)
double complex             :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of complex-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Curtis-Clenshaw
!      quadrature rule
!    vals - a matrix whose jth column gives the values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
!
double precision              :: whts(n)
!double precision, allocatable :: sums(:)


l = size(vals,2)

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xscheb(i)
if ( abs(dd) .le. eps0) then
valsout(1:l) = vals(i,1:l)
return
endif
whts(i) = dsign/dd
dsign   = -dsign
end do

whts(1) = whts(1)/2
whts(n) = whts(n)/2

! compute the sums
! do i=1,l
! vals(1:n,i) = vals(1:n,i) * whts(1:n) 
! sums(i)     = sum (vals(:,i))
! vals(1:n,i) = vals(1:n,i) / whts(1:n)
! end do

do j=1,l
valsout(j) = 0
do i=1,n
valsout(j) = valsout(j) + vals(i,j) * whts(i)
end do
end do


! evaluate the formula
whts    = whts
valsout = valsout/sum(whts)

end subroutine


subroutine chebyshev_eval1(n,coefs,x,val)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val
!
!  Evaluate a complex-valued expansion of the form (1) given the vector
!  (2) of its  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)

call chebs(n,x,pols)

val = 0
do i=1,n
val = val + pols(i)*coefs(i)
end do


end subroutine


subroutine chebyshev_eval2(n,coefs,x,val)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:), val
!
!  Evaluate a real-valued expansion of the form (1) given the vector
!  (2) of its  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)

call chebs(n,x,pols)

val = 0
do i=1,n
val = val + pols(i)*coefs(i)
end do

end subroutine


subroutine chebyshev_eval3(n,coefs,x,vals)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:)
!
!  Evaluate a collection of real-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 

double precision :: pols(n+1)

call chebs(n,x,pols)

vals = 0
do i=1,n
vals = vals + pols(i)*coefs(i,:)
end do


end subroutine


subroutine chebyshev_eval4(n,coefs,x,vals)
implicit double precision (a-h,o-z)
double complex                  :: coefs(:,:), vals(:)
!
!  Evaluate a collection of complex-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 

double precision :: pols(n)

call chebs(n,x,pols)

vals = 0
do i=1,n
vals = vals + pols(i)*coefs(i,:)
end do


end subroutine

subroutine chebyshev_evalder1(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val, der
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call chebders(n,x,pols,ders)

val = 0
der = 0

do i=1,n
val = val + pols(i)*coefs(i)
der = der + ders(i)*coefs(i)
end do

end subroutine


subroutine chebyshev_evalder2(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:)
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call chebders(n,x,pols,ders)

val = 0
der = 0

do i=1,n
val = val + pols(i)*coefs(i)
der = der + ders(i)*coefs(i)
end do

end subroutine


subroutine chebyshev_evalder3(n,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call chebders(n,x,pols,ders0)

vals = 0
ders = 0

do i=1,n
vals = vals + pols(i)*coefs(i,:)
ders = ders + ders0(i)*coefs(i,:)
end do

end subroutine


subroutine chebyshev_evalder4(n,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call chebders(n,x,pols,ders0)

vals = 0
ders = 0

do i=1,n
vals = vals + pols(i)*coefs(i,:)
ders = ders + ders0(i)*coefs(i,:)
end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebs(n,x,pols)
implicit double precision (a-h,o-z)
double precision           :: pols(:)
!
!  Return the values of the Chebyshev polynomials of degrees 0 through n-1 at 
!  a specified point.
!
!  Input parameters:
!     n - the number of polynomials to evaluate
!     x - the point at which to evaluate them
!
!  Output parameters:
!     pols - the array of length n containing the values of the first n 
!       polynomials
!

if (x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
end do
return
endif

if (x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
dsign   = -dsign
end do
return
endif

pols(1) = 1.0d0
if (n == 1) return
pols(2) = x
if (n == 2) return


do j=2,n-1
pols(j+1) = 2*x*pols(j)-pols(j-1)
end do

end subroutine


pure subroutine chebders(n,x,pols,ders)
implicit double precision (a-h,o-z)
integer, intent(in)           :: n
double precision, intent(in)  :: x
double precision, intent(out) :: pols(:),ders(:)
!
!  Evaluate the Chebyshev polynomials of degrees 0 through n-1 and their derivatives
!  at a specified point using the standard 3-term recurrence relation.
!
!  Input parameters:
!    n - the number of polynomials to evaluate
!    x - point at which to evaluate the polynomials
!
!  Output parameters:
!   pols - the values of the Chebyshev polynomials
!   ders - the values of their derivatives
!


if (x .eq. 1.0d0) then
pols = 1
do i=1,n
ders(i) = (i-1)**2
end do
return
endif

if (x .eq. -1.0d0) then
pols(1) = 1.0d0
ders(1) = 0.0d0

dsign = -1.0d0
do i=2,n
pols(i) = dsign
ders(i) = -dsign*(i-1)**2
dsign   = -dsign
end do
return
endif


if (n .eq. 0) return

pols(1) = 1.0d0
ders(1) = 0.0d0

if (n .eq. 1) return

pols(2) = x
ders(2) = 1.0d0

do i=2,n-1
pols(i+1) = 2*x*pols(i) - pols(i-1)
ders(i+1) = 2*pols(i) + 2*x*ders(i) - ders(i-1)
end do

! xx = 2*x*pols(n) - pols(n-1)


! !
! !  Compute the derivatives
! !

! do i=2,n-1
! ders(i) = (i-1)*(x*pols(i)-pols(i+1))/(1-x**2)
! end do
! ders(n)   = (n-1)*(x*pols(n)-xx)/(1-x**2)


end subroutine


end module
