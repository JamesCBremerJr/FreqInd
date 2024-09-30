!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains a few utility routines used in our numerical experiments.
!
!  The following routines should be regarded as publicly callable:
!
!    elapsed - return the second elapsed since an arbitrary point in the past;
!      one most systems, the timer this routine relies on has excellent resolution
!
!    plot_functions2 - produce a python file which uses the matplotlib library to
!      generate a PDF file containing a plot of one or more functions
!
!    legeder - evaluate a Legendre polynomial of specified degree and its first 
!      derivative via the well-known 3-term recurrence relations
!
!    legeqder - evaluate the Legendre function of the second kind of a specified
!      degree n and its derivative via the well-known 3-term recurrence relations
!
!    gegpol - evaluate the Gegenbauer polynomial of specified degree n and order
!      using the 3-term recurrence relations
!
!    gegpolw - evaluate the Wronskian-normalized Gegenbauer polynomial of specified
!      degree and order using the 3-term recurrence relation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module experutils

contains


subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
call system_clock(i,irate)

dd = i
dd = dd/irate
t = dd
return
end subroutine

subroutine plot_functions2(scriptname,filename,nfuns,xlabel,ylabel,iflogx,iflogy,x1_,x2_,y1_,y2_, &
   n,xs,ys,legend_loc,legends,linestyles,colors)
implicit double precision (a-h,o-z)
character(len=*)               :: scriptname,filename,legends,legend_loc,linestyles,colors
character(len=*)               :: xlabel,ylabel
double precision               :: xs(n),ys(nfuns,n)
character(len=:), allocatable  :: command, legend, style, color
!
!  Produce a python script which generates a PDF file containing a plot of one or more
!  functions of one variable.
!
!  Input parameters:
!    scriptname - the name of the python script
!    filename - the name of the PDF file to produce
!
!    xlabel - a label for the x-axis (no label will be present if the string is empty)
!    ylabel - a label for the y-axis (no label will be present if the string is empty)
!
!    iflogx - an integer parameter specifying whether the x-axis should be logarithm
!      scale or not
!
!      iflogx = 0    not logarithmic scale
!      iflogx = 1    10^k scale
!      iflogx = 2    2^k scale
!
!    iflogy - an integer parameter specifying whether the y-axis should be logarithm
!      scale or not
!
!      iflogy = 0    not logarithmic scale
!      iflogy = 1    10^k scale
!      iflogy = 2    2^k scale

!    (x1,x2,y1,y2) - extents of the axes ---- if x2 <= x1 then these will be set
!      automatically; likewise if y2 <= y1
!
!    legend_loc - a string specifying the location for the legend --- this can
!       be "upper right", "upper left", "lower right", etc  OR a blank string if 
!       no legend is to be included OR "best" for automatic placement
!
!    legends - a string specifying the legend labels for each function ...
!       each should be terminated by an asterik "*" so "1*2*3*" specifies
!       the label 1 for the first function, 2 for the second, and 3 for the
!       third; this is ignored if legend_loc is blank
!
!    styles - a string specifying styles for each function ... separated as in
!       the legends string ... ignored if empty.  The style strings are passed
!       to matploblib 
!     
!    n - the number of point in the graph of the function to be specified
!    xs - the x-coordinates of the points to plot
!    ys - an (nfuns,n) array whose jth column gives the y-coordinates on
!        the graph of the jth function
!
!  Output parameters:
!    N/A
!


x1 = x1_
x2 = x2_
y1 = y1_
y2 = y2_

if (x2 .le. x1) then
x1 =  1d300
x2 = -1d300
do i=1,n
x1 = min(x1,xs(i))
x2 = max(x2,xs(i))
end do

endif

if (y2 .le. y1) then
y1 =  1d300
y2 = -1d300

do i=1,n
do j=1,nfuns
y1 = min(y1,ys(j,i))
y2 = max(y2,ys(j,i))
end do
end do


if (iflogy .eq. 0) then

if (y1 .gt. 0) then
y1 = y1 * 0.98d0
else
y1 = y1*1.02d0
endif

if (y2 .gt. 0) then
y2 = y2 * 1.02d0
else
y2 = y2 * 0.98d0
endif


endif

endif

iw = 1001
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "import warnings"
write(iw,"(A)") "mpl.use('Agg')"
write(iw,"(A)") "import matplotlib.pyplot as plt"

write(iw,"(A)") 'warnings.filterwarnings("ignore")'
write(iw,"(A,I5,A)") "xs  = np.zeros(",n,")"
write(iw,"(A,I5,A,I5,A)") "ys = np.zeros((",nfuns,",",n,"))"

write(iw,"(A)") "plt.rcParams.update({'font.size': 12, 'figure.autolayout': True})"

do i=1,n
write(iw,"(A,I5,A,ES30.18E3)") "xs[",i-1,"]  = ",xs(i)
end do

do i=1,n
do j=1,nfuns
write(iw,"(A,I5,A,I5,A,ES30.18E3)") "ys[",j-1,",",i-1,"] = ",ys(j,i)
end do
end do

write(iw,"(A)") "fig, ax = plt.subplots()"

if (len(xlabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(xlabel="',xlabel,'")'
endif

if (len(ylabel) .gt. 0) then
write (iw,"(A,A,A)") 'ax.set(ylabel="',ylabel,'")'
endif



!allocate(legend(0),style(0))

idx1 = 1
idx2 = 1

idx3 = 1 
idx4 = 1

idx5 = 1
idx6 = 1

do j=1,nfuns


! find the legend string
if (len(legend_loc) .gt. 0) then
do while (legends(idx2:idx2) .ne. '*') 
idx2 = idx2+1
end do
ll = idx2-idx1
allocate(character(ll) :: legend )
legend(1:ll) = legends(idx1:idx2-1)
idx1 = idx2+1
idx2 = idx1
else
allocate(character(0) :: legend )
endif


! find the style string
if (len(linestyles) .gt. 0) then
do while (linestyles(idx4:idx4) .ne. '*') 
idx4 = idx4+1
end do
ll = idx4-idx3
allocate(character(ll) :: style )
style(1:ll) = linestyles(idx3:idx4-1)
idx3 = idx4+1
idx4 = idx3
else
allocate(character(0) :: style  )
endif


! find the color string
if (len(colors) .gt. 0) then
do while (colors(idx6:idx6) .ne. '*') 
idx6 = idx6+1
end do
ll = idx6-idx5
allocate(character(ll) :: color )
color(1:ll) = colors(idx5:idx6-1)
idx5 = idx6+1
idx6 = idx5
else
allocate(character(0) :: color  )
endif



write(iw,"(A,I5,A,A,A,A,A,A,A)") 'ax.plot(xs,ys[',j-1,',:],linestyle="',style,'",color="',color,'",label="',legend,'")'

deallocate(style)
deallocate(legend)
deallocate(color)
end do


if (len(legend_loc) .gt. 0) then
write(iw,"(A,A,A)") 'plt.legend(loc="',legend_loc,'")'
endif

if (iflogx .eq. 1) then
write(iw,"(A)") 'plt.xscale("log")'
elseif (iflogx .eq. 2) then
write(iw,"(A)") 'ax.set_xscale("log", base=2)'
endif

if (iflogy .eq. 1) then
write(iw,"(A)") 'plt.yscale("log")'
elseif (iflogy .eq. 2) then
write(iw,"(A)") 'ax.set_yscale("log", base=2)'
endif

write(iw,"(A,E24.15,A,E25.15,A)") "plt.xlim([",x1,",",x2,"])"
write(iw,"(A,E24.15,A,E25.15,A)") "plt.ylim([",y1,",",y2,"])"
write(iw,"(A)") 'ax.grid()'
write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'

close(iw)

! allocate(character(7+len(scriptname)) :: command )
! write(command,"(A,A)") "python ",scriptname
! call system(command)

print *,""
print *,"plot_functions2: execute the command `python ",scriptname,"' to produce the "
print *,"                 plot ",filename
print *,""

end subroutine


subroutine legeder(n,x,pol,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre polynomial of degree n and its derivative
!  at the point x.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of P_n(x)
!    der - the value of P_n'(x)
!


if ( x == 1.0d0) then
pol = 1
der = (n)*(n+1.0d0)/2
return
endif

if ( x == -1.0d0) then
pol = (-1.0d0)**n
der = -pol*(n)*(n+1.0d0)/2
return
endif


if (n == 0) then
pol = 1
der = 0
else if (n == 1) then
pol = x
der = 1.0d0
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do
!
pol = p

!
! Compute the derivative using another well-known formula.
!

der=n*(x*p2-p1)/(x**2-1.0d0)
endif

end subroutine


subroutine legeqder(n,x,val,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the Legendre function of the second kind  of degree n and its derivative
!  at the point x.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of Q_n(x)
!    der - the value of Q_n'(x)
!


q1 = 0.5d0*(log(1+x)-log(1-x))
q2 = -1.0d0 + x * 0.5d0*(log(1+x)-log(1-x))

d1 = 0.5d0*(1/(1+x) + 1/(1-x))
d2 = 0.5d0*(2*x/(1-x**2)-log(1-x)+log(1+x))

if (n .eq. 0) then
val = q1
der = d1
return
endif

if (n .eq. 1) then
val = q2
der = d2
return
endif


do j=2,n
q = ((2*j-1)*x*q2-(j-1)*q1)/j
d = ((2*j-1)*q2+(2*j-1)*x*d2-(j-1)*d1)/j

q1 = q2
q2 = q

d1 = d2
d2 = d

end do

val = q
der = d


end subroutine





subroutine gegpol(n,da,x,pol,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the usual Gegenbauer polynomial C_n^a (x) of degree n and 
!  order da at the point x via a recurrence relation.  This routine evaluates
!  the standard, unnormalized polynomial, which is the unique solution of the 
!  Gegenbauer equation which is regular on [-1,1] and such that
! 
!     (a)              Gamma(2a+1)
!    C    (1) =  ----------------------   .
!     n           Gamma(2a) Gamma(n+1)
!  
!  Note that this is not a reasonable normalization when a is large.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    da - the parameter which determines the weight function for
!      the Gegenbauer polynomial
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of C_n^da(x)
!    der - the value of C_n^da'(x)
!

if (n==0) then
pol = 1
der = 0
return
endif

if (n==1) then
pol = 2*da*x
der = 2*da
return
endif

val1 = 1
der1 = 0

val2 = 2*da*x
der2 = 2*da

do i=1,n-1
val  = 2*(i+da)*x*val2 - (i+2*da-1)*val1
val  = val/(i+1.0d0)

der  = 2*(i+da)*val2 + 2*(i+da)*x*der2 -(i+2*da-1)*der1
der  = der/(i+1.0d0)

val1 = val2
val2 = val

der1 = der2
der2 = der

end do

pol = val
end subroutine


subroutine gegpolw(n,da,x,pol,der)
implicit double precision (a-h,o-z)
!
!  Evaluate the "Wronskian normalized" Gegenbauer polynomial of degree n
!  and order da at the point x via a recurrence relation.
!
!  If C_{n,da} is the usual Gegenbauer polynomial of degree n and order a 
!  and D_{n,da} is the usual Gegenbauer function of the second kind of
!  degree n and order a on the cut, then 
!
!    C_{n,da} + i 2/pi D_{n,da}
!
!  has a slowly-varying modulus function on (-1,1).  This routine
!  evaluates the function
!
!     1 / sqrt( W(C_{n,da},2/pi D_{n,da}) ) P_{n,da}
!
!  where W(u,v) is the Wronskian of the pair u,v.  We refer to this
!  as the "Wronskian normalized" Gegenbauer polynomial of order
!  da and degree n.
!
!  This routine is not designed to be efficient and it loses some
!  accuracy as x approaches the singularities at +/- 1 and
!  as da approaches -0.5d0 or becomes excessively large.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    da - the parameter which determines the weight function for
!      the Gegenbauer polynomial
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of C_n^da(x)
!    der - the value of C_n^da'(x)
!
data pi / 3.14159265358979323846264338327950288d0 / 

call gammaratio1(da-0.5d0,val1)
call gammaratio2(da,val2)
call gammaratio3(da,val3)

cc1 = pi**(0.25d0)/sqrt(2.0d0)*sqrt(val1)
cc2 = pi**(0.50d0)*val2
cc3 = -val3
cc4 = -2*(1+da)*cc3


if (n==0) then
pol = cc1
der = 0
return
endif

if (n==1) then
pol = cc2*x
der = cc2
return
endif

if (n==2) then
pol = cc3 + cc4*x**2
der = cc4*2*x
return
endif


val1 = cc2*x
val2 = cc3+cc4*x**2

der1 = cc2
der2 = 2*cc4*x

do i=2,n-1

a    =  (0.2d1*(da+i))/Sqrt((1+i)*(2*da+i))
b    =  (-0.1d1*sqrt((i*(0.1d1+i)*(-0.1d1+0.2d1*da+i))/(0.2d1*da+i)))/(1+i)

pol  = a*x*val2 + b*val1
der  = a*x*der2 + a*val2 + b*der1

val1 = val2
val2 = pol

der1 = der2
der2 = der

end do

end subroutine


subroutine gammaratio1(x,val)
implicit double precision (a-h,o-z)
!
!  Return the value of
!
!     gamma(x+1/2) / gamma(x+1)
!
!  where x > 0.
!

if (x .lt. 100) then
val = exp(log_gamma(x+0.5d0)-log_gamma(x+1.0d0))
else
y   = 1/sqrt(x)
val = y-0.125d0*y**3+0.78125d-2*y**5+0.48828125d-2*y**7-  &
      0.640869140625d-3*y**9-0.1522064208984375d-2*y**11
endif

end subroutine


subroutine gammaratio2(x,val)
implicit double precision (a-h,o-z)

!
!  Return the value of
!
!     2^(-x)  Sqrt(Gamma(1+2x)) / Gamma(x+1/2)   
!
!  for x > 0.
!



if (x .lt. 100) then

dd  = -x*log(2.0d0)+0.5d0*log_gamma(1+2*x) - log_gamma(x+0.5d0)
val = exp(dd)

else

val = -0.878742147871772197957768899945973103d-4/x**0.75d1-  &
0.661216836363539576060570014495473562d-3/x**0.65d1+  &
0.116892039323568394984693775989129872d-3/x**0.55d1+  &
0.8587327721998396243754718825254374d-3/x**0.45d1-  &
0.361571693557827210263356582115973642d-3/x**0.35d1-  &
0.275483195091677874486366919707408489d-2/x**0.25d1+  &
0.440773112146684599178187071531853583d-2/x**0.15d1+  &
0.705236979434695358685099314450965732d-1/Sqrt(x)+  &
0.564189583547756286948079451560772586d0*sqrt(x)
val = sqrt(val)

endif

end subroutine


subroutine gammaratio3(x,val)
implicit double precision (a-h,o-z)

!
!  Return the value of
!
!     2^(x-0.5d0)  Gamma(x+1) / Sqrt( Gamma(2x+2) )
!
!  for x > 0.
!

if (x .lt. 100) then

dd  = log(2.0d0)*(x-0.5d0) + log_gamma(x+1) - 0.5d0*log_gamma(2*x+2)
val = exp(dd)

else


val = 0.536314366705528027569004728149254532d-2/da**0.65d1-  &
0.105426737481095347044550615249621606d-1/da**0.55d1+  &
0.224342417804889762061131275548016347d-1/da**0.45d1-  &
0.454364390490915973794696254274463487d-1/da**0.35d1+  &
0.865455981887458997704183341475168546d-1/da**0.25d1-  &
0.166167548522392127559203201563232361d0/da**0.15d1+  &
0.443113462726379006824541870835286296d0/Sqrt(da)
val = sqrt(val)

endif

end subroutine

end module
