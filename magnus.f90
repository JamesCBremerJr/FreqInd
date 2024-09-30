!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for solving initival value problems for second order
!  linear ordinary differential equations of the form
!
!    y''(t) + g(t) y(t) = 0                                                               (2)
!
!  via a 4^th$ order modified Magnus method. 
!
!  The following routines should be regarded as publicly callable:
!
!    magnus_modified - solve a scalar equation of the form y''(t) + g(t)y(t) = 0
!      using a modified Magnus method and an equispaced discretization grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module magnus


interface

subroutine magnus_fun1(t,val,der,pars)
implicit double precision (a-h,o-z)
double precision            :: t
double precision            :: val
double complex, pointer     :: pars(:)
end subroutine

end interface

contains


subroutine magnus_modified(a,b,nsteps,y0,ts,ys,gfun,pars)
implicit double precision (a-h,o-z)
procedure(magnus_fun1)                     :: gfun
double complex, pointer                    :: pars(:)
double complex                             :: y0(2)
double complex, allocatable, intent(out)   :: ys(:,:)
double precision, allocatable, intent(out) :: ts(:)

!
!  A simple implementation of the modified Magnus method for solving a differential 
!  equation of the form (2).  This solver uses an equispaced grid whose density is 
!  determined by the "nsteps" parameter.  
!
!  Input parameters:
!    (a,b) - the interval over which the ODE is to be solved
!    n - the order of the differential equation
!    nsteps - the number of steps to take
!    y0 - the initial value for the solution of the ODE
!    gfun - an external subroutine which supplies the values of the function g
!    pars - a user-supplied array of doubles which is passed to gfun
!
!  Output parameters:
!    ts - an array giving the locations of the discretization grid points
!    ys - an array giving the values of the obtained solutions at the
!     grid points
!  

double complex                             :: y1(2), y2(2)

allocate(ys(2,nsteps+1), ts(nsteps+1))

h       = (b-a)/(nsteps+0.0d0)
ts(1)   = a
ys(:,1) = y0

do i=1,nsteps
t1        = ts(i)
t2        = t1+h
y1        = ys(:,i)

dnu = pars(1)
call magnus_step_modified(t1,y1,t2,y2,gfun,pars)

ts(i+1)   = t2
ys(:,i+1) = y2
end do

end subroutine



subroutine magnus_step_modified(t1,y1,t2,y2,gfun,pars)
implicit double precision (a-h,o-z)
double complex            :: y1(2), y2(2)
double complex, pointer   :: pars(:)
procedure(magnus_fun1)    :: gfun 
!
!  Take a step using a 4th order modified Magnus integrator
!
double complex           :: theta1(2,2), x1(2)
double precision         :: theta2(2,2), omega(2,2)
double precision         :: work(100)
integer                  :: ipiv(100)

n       = 2
h       = t2-t1
t0      = (t1+t2)/2

call gfun(t0,eta,der,pars)
psi = sqrt(eta)*h

theta1(1,1) =  cos(psi)
theta1(1,2) = 1/sqrt(eta)*sin(psi)
theta1(2,1) = -sqrt(eta)*sin(psi)
theta1(2,2) = cos(psi)

call gfun(t1,f0,fp0,pars)
call gfun(t2,fh,fph,pars)
f0=f0-eta
fh=fh-eta

call magnus_filon1(2*sqrt(eta),h,f0,fp0,fh,fph,valsin)
theta2(1,1) =  0.5d0*1/sqrt(eta) * valsin
theta2(2,2) = -0.5d0*1/sqrt(eta) * valsin
call magnus_filon2(sqrt(eta),h,f0,fp0,fh,fph,valcos,valsin)
theta2(1,2) = 1.0d0/eta * valsin
theta2(2,1) = -valcos

call DGPADM(6,2,1.0d0,theta2,2,work,100,ipiv,iexph,ns,iflag )
omega(1,1) = work(iexph)
omega(2,1) = work(iexph+1)
omega(1,2) = work(iexph+2)
omega(2,2) = work(iexph+3)

x1 = matmul(omega,y1)
y2 = matmul(theta1,x1)

end subroutine


subroutine magnus_filon1(omega,h,f0,fp0,fh,fph,valsin)
implicit double precision (a-h,o-z)

valsin = (h*omega*(-0.6d1*fh+0.4d1*fp0*h+0.2d1*fph*h+f0*(6+  &
h**2*omega**2))+h*omega*(0.6d1*f0+0.2d1*fp0*h+0.4d1*fph*h-  &
0.1d1*fh*(6+h**2*omega**2))*Cos(h*omega)+(0.12d2*fh-  &
0.6d1*(2*f0+fp0*h)+fph*h*(-6+  &
h**2*omega**2))*Sin(h*omega))/(h**3*omega**4)

end subroutine

subroutine magnus_filon2(omega,h,f0,fp0,fh,fph,valcos,valsin)
implicit double precision (a-h,o-z)

valsin = (0.6d1*f0*(-3+2*h**4*omega**4+3*Cos(2*h*omega)+  &
3*h*omega*Sin(2*h*omega))+fp0*h*(-9+0.6d1*h**2*omega**2+  &
0.2d1*h**4*omega**4+0.9d1*Cos(2*h*omega)+  &
0.6d1*h*omega*Sin(2*h*omega))+fph*h*(-9-0.2d1*h**4*omega**4+  &
(9-0.6d1*h**2*omega**2)*Cos(0.2d1*h*omega)+  &
0.12d2*h*omega*Sin(2*h*omega))+0.6d1*fh*(3+2*h**4*omega**4-  &
3*Cos(2*h*omega)-h*omega*(3+  &
2*h**2*omega**2)*Sin(2*h*omega)))/(0.48d2*h**3*omega**4)

valcos = -0.20833333333333331d-1*(-0.18d2*f0-0.9d1*fp0*h+  &
0.6d1*fp0*h**3*omega**2-0.12d2*f0*h**4*omega**4-  &
0.2d1*fp0*h**5*omega**4+0.18d2*f0*Cos(2*h*omega)+  &
0.9d1*fp0*h*Cos(2*h*omega)+0.18d2*f0*h*omega*Sin(2*h*omega)+  &
0.6d1*fp0*h**2*omega*Sin(2*h*omega)+fph*h*(-9+  &
0.2d1*h**4*omega**4+(9-  &
0.6d1*h**2*omega**2)*Cos(0.2d1*h*omega)+  &
0.12d2*h*omega*Sin(2*h*omega))-0.6d1*fh*(-3+2*h**4*omega**4+  &
3*Cos(2*h*omega)+h*omega*(3+  &
2*h**2*omega**2)*Sin(2*h*omega)))/(h**3*omega**4)

end subroutine


end module
