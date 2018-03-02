      program convolute
!
! Gaussian convolution code in FORTRAN 90/95
!

!
! Width is given by ewidth = ...? units??
!
      implicit none

      integer, parameter  	:: dp=8
      integer, parameter 	:: num=8500000
      integer  		 	:: i,j,n,m,npm,new,nskip,jmin,jmax
      real(dp) 		  	:: a,alpha,e,bot,h,hhalpha,rat,pi,sum,x1
      real(dp), allocatable 	:: x(:),y(:),expon(:)
      real(dp), parameter       :: ewidth=5.00e-0_dp
      real(dp), parameter 	:: ymult=1.0_dp,p1=1.0_dp
      real(dp), parameter       :: one=1.0_dp,two=2.0_dp,p25=1.25_dp
      real(dp), parameter 	:: four=4.0_dp,five=5.0_dp,ten=10.0_dp
!
      allocate (x(num))
      allocate (y(num))
      allocate (expon(num))

!     open existing file XPISUM to read data from
!
      open(unit=1,file='XPISUM')
      open(unit=3,file='con')

!      write(6,*)' ymult = ??'
!      read(5,*)ymult

      n=0
   1  read(1,*,end=20)e,x1
      n=n+1
      x(n)=e
      y(n)=x1
      go to 1

   20 continue

		do i=2,n-1
                   bot=max(y(i-1),y(i+1))
                   rat=y(i)/bot
                	if(rat.gt.ymult)then
                           y(i)=(y(i-1)+y(i+1))/two
                	 endif
		enddo

      h=x(2)-x(1)
 !
 !    Convolution value ewidth set in parameter statement
 !
      m=n
      new=int(five/ten * ewidth/h)
      n=n+new

      		if(n.gt.num)then
      		   write(0,*)' error, n=',n,'check dimensions'
      		   stop
     		 endif

      		do i=m+1,n
      		   x(i)=x(i-1)+h
      		   y(i)=y(m)
		enddo

      		 alpha=four * log(two)/(ewidth*ewidth)
     		 pi=four * atan(one)
     		 a=sqrt(alpha/pi)
     		 npm=int(five /ten * ewidth/h)
!
!   use if you want to skip over a number of energy points
!
!     read(5,*)nskip
!     write(6,*)' nskip=??'

      nskip=1
      hhalpha=h*h*alpha
!
      	do i=1,n
           expon(i)=EXP(-real(i*i)*hhalpha)
      	enddo
      		do i=1,n-1,nskip
      		   npm=int(ten / ten * ewidth/h)
      		   jmin=max(1,i-npm)
      	           jmax=min(n,i+npm)
      		   sum=y(i)
     			do j=jmin,i-1
      			   sum=sum+y(j)*expon(i-j)
 			enddo
     				 do j=i+1,jmax
      				    sum=sum+y(j)*expon(j-i)
   				 enddo
      		   sum=sum*a*h
      		   if(x(i).le.x(m))write(3,100)x(i),sum
   		enddo
  100 format(2p,e16.8,1p,e16.6)
      stop
      end program convolute
