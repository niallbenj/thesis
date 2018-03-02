        PROGRAM photoarrange


C


C       Version 1.0      Connor Ballance     May 22 2005


C


C       *make sure that the dimensions of final state photoionisation    *


C       *cross sections are appropriately in stgbf0damp.f / pstgbf0damp.f*


C


C      (1) merges XPISUMXXX files together


C


C          Total photoionisation (including excitation)


C          from a single intial bound term/level


C


C      (2) merges XPIPARXXX files togther when


C          the namelist variables are set, extracts


C          a partial sum (from 1 to all) of final


C          ion state resolved cross sections (Mb)


C


C          see istart,ifinal,ipart below for details


C


C      (3) quick sort algorithm accredited below


C


C      (4) davint (integral package for gaussian convolution) -accredited below


C


        implicit real*8 (a-h,o-z)


        INTEGER :: IONE=1


        INTEGER,allocatable :: POINTER(:)


	INTEGER,allocatable :: ISAT(:),ILAT(:)


        REAL*8,allocatable :: EMESH1(:),PHOTOALL(:,:),photopar(:,:),


     X 	ENAT(:)


        REAL*8,allocatable :: PHOTOPRECON(:),PHOTOCON(:)


        character*1 filec,filed,fileu,filev


        logical bform,dd


	NAMELIST/arrinp/istart,ifinal,ipart,ymult,ewidth,emin,emax


C


C


C


C       istart=first residual ion photoionisation term/level


C               to start summing over


C


C


C       ifinal= last ....


C


C


C       ipart= 1  (ensures the previous 2 commands have effect)


C


C        -------------------------------------------------------


C


C       ie ipart =1, istart=1 & ifinal =1 gives the final state


C                                         photoionisation to the


C                                         groundstate only.


C


C        ------------------------------------------------------


C


C       ymult = value   ( removes numerical spikes, by removing


C                         points that vary from the average background


C                         of those either side by a factor of ymult )


C


C


C        ewidth = if not equal to 0.0d0 then perform gaussian convolution


C                 on the the total (filtered/unfiltered) total photoionisation


C                 cross section   ( ewidth in Ryds )


C


C


C         emin = stop convoluting below this point ... usually the energy


C                    of the final threshold relative to the groundstate


C                    in Ryds


C


C


C        emax = stop convoluting beyond this point ... usually the energy


C                    of the final threshold relative to the groundstate


C                    in Ryds


C


C


        INQUIRE(file='dphoto',exist=dd)


C


C


C


        NPHOTOFINAL=994     !!!!! max number of final state resolved   !!!!


C                                photoionisation cross sections, must


C                                correspond to the compiled value


C				 of pstgbf0damp/stgbf0damp PARAM


C


C


        ipart=0


	istart=1


	ifinal=1


        ymult=-1.0d0


	ewidth=0.0d0


        emax=0.0


        emin=0.0d0


        if(.not.dd)then


        continue


        else


        open(8,file='dphoto',form='formatted',status='unknown')


        read(8,arrinp)


        close(8)            ! for later options


        endif


C


C


C


        open(9,file='XPISUMZ',form='formatted',status='unknown')


C


C


        nproc=9999


	mxe1=0


        K=0


C


C


C


      do iam=0,nproc-1


        ich0 = ichar('0')


        ic = iam/1000


        id = (iam - 1000*ic)/100


        iu = (iam - 1000*ic - 100*id)/10


        iv = (iam - 1000*ic - 100*id -iu*10)


        ic = ic + ich0


        id = id + ich0


        iu = iu + ich0


        iv = iv + ich0


        filec = char(ic)


        filed = char(id)


        fileu = char(iu)


        filev = char(iv)


C


        inquire (file='XPISUM'//filec//filed//fileu//filev,exist=bform)


        if(bform)then


          K=K+1


        open(10+iam,file='XPISUM'//filec//filed//fileu//filev,


     +  status='old')


        write(0,*)'opening xpisum',filec//filed//fileu//filev


        goto 30


        endif


C


        goto 35


 30     continue


            READ(10+iam,*)NSPN,LRGL,NPTY,NBOUND


            READ(10+iam,*)BE,MXE


        MXE1=MXE1+MXE


        rewind(10+iam)


       enddo


C


  35   continue


C


       allocate(EMESH1(mxe1))


       allocate(POINTER(mxe1))


       allocate(PHOTOALL(mxe1,8))


C


C      Re-initialise K   !!!!


C


       nproc=K


       K=0


C


C


c........ suffix for input files xiipsum000,xpisum001 etc


        do iam=0,nproc-1


C


C        print *,'reading file XPITOT',iam


            READ(10+iam,*)NSPN,LRGL,NPTY,NBOUND


            READ(10+iam,*)BE,MXE


C


C


C


C


        DO L=1,MXE


        K=K+1


        READ(10+iam,*)EMESH1(K),PHOTOALL(K,NBOUND)


        POINTER(K)=K


        ENDDO


C


C


        ENDDO


        WRITE(6,*)'TOTAL NUMBER OF POINTS = ',K


C


        call qsort(EMESH1,K,POINTER)


C


C


        WRITE(9,*)NSPN,LRGL,NPTY,NBOUND


        WRITE(9,711)BE,MXE1


  700 FORMAT(1PE14.8,6(1PE11.3)/(14X,6(E11.3)))


  711 FORMAT(E16.6,I6)


C


        do JJ=1,MXE1


	  write(9,700)EMESH1(JJ),PHOTOALL(POINTER(JJ),NBOUND)


	enddo


C


        if(ymult.ne.(-1.0d0))then


          rewind(9)


         READ(9,*)NSPN,LRGL,NPTY,NBOUND


         READ(9,711)BE,MXE1


          do JJ=1,MXE1


          read(9,700)EMESH1(JJ),PHOTOALL(JJ,NBOUND)


          enddo


  720    iflagy=0


         if(mxe1.le.4)stop 'not enough points to filter'


         do JJ=2,MXE1-1


         average=PHOTOALL(JJ+1,NBOUND)+PHOTOALL(JJ-1,NBOUND)/2.0d0


         spikeratio=PHOTOALL(JJ,NBOUND)/average


         if(spikeratio.gt.ymult)then


         grad=PHOTOALL(JJ+1,NBOUND)-PHOTOALL(JJ-1,NBOUND)


         grad=grad/(EMESH1(JJ+1)-EMESH1(JJ-1))


         PHOTOALL(JJ,NBOUND)=PHOTOALL(JJ-1,NBOUND)


     X                       +grad*(EMESH1(JJ)-EMESH1(JJ-1))


         iflagy=1


         endif


         if(iflagy.eq.1)then


         goto 720


         endif


         enddo


         rewind(9)


         WRITE(9,*)NSPN,LRGL,NPTY,NBOUND


         WRITE(9,711)BE,MXE1


          do JJJ=1,MXE1


          write(9,700)EMESH1(JJJ),PHOTOALL(JJJ,NBOUND)


          enddo


        endif


C


C


C        Gaussian convolution section  (ewidth controls use)


C


         if(ewidth.gt.(0.0d0))then      ! start


C


         if(emax.eq.(0.d0))then


         emax=emesh1(mxe1)


         endif


C


         if(emin.eq.(0.d0))then


         emin=emesh1(1)


         endif


C


C


C


C        did i filter photoall or not?  ..... both orderings accounted for


C


         allocate(photoprecon(mxe1),photocon(mxe1))


C


	 do JJJ=1,mxe1


	  if(ymult.ne.(-1.0d0))then


	 photoprecon(JJJ)=photoall(JJJ,NBOUND)


	  else


	 photoprecon(JJJ)=photoall(POINTER(JJJ),NBOUND)


	  endif


         enddo


C


         deallocate(photoall)    !  already written to XPISUMZ .... dumping





C


C


         do JJJ=1,mxe1


	 iflg=0


         photocon(JJJ)=gaussn(JJJ,1,mxe1,mxe1,1,ewidth,iflg,


     1   emesh1,photoprecon,emin,emax)


	 enddo


C


         deallocate(photoprecon)  ! post convolution ..... dumping


C


      open(8,file='XPISUM-CONVOLUTED',form='formatted',status='unknown')


         do JJJ=1,mxe1


	 write(8,*)EMESH1(JJJ),PHOTOCON(JJJ)


	 enddo


C


         deallocate(PHOTOCON)     ! post writing of result .... dumping


	 endif


C


C


C


	close(9)


C


        if(ipart.eq.1)then


	open(9,file='XPIPARZ',form='formatted',status='unknown')


	else


	stop


	endif                            ! endif


C


C       END OF XPISUM arrange


C


         if((istart.ge.1).and.(ifinal.le.NPHOTOFINAL))then


	 continue


	 else


	 write(0,*)'istart or ifinal values inconsistent'


	 write(0,*)'istart=',istart,'ifinal=',ifinal,


     X             'ifinal dimensioned at ', NPHOTOFINAL


         stop


         endif


C


	 allocate(photopar(mxe1,NPHOTOFINAL))


	 JK=0


C


         do 10 iam=0,nproc-1


        ich0 = ichar('0')


        ic = iam/1000


        id = (iam - 1000*ic)/100


        iu = (iam - 1000*ic - 100*id)/10


        iv = (iam - 1000*ic - 100*id -iu*10)


        ic = ic + ich0


        id = id + ich0


        iu = iu + ich0


        iv = iv + ich0


        filec = char(ic)


        filed = char(id)


        fileu = char(iu)


        filev = char(iv)


        open(1000+iam,file='XPIPAR'//filec//filed//fileu//filev,


     +  status='old')


         write(0,*)'reading file', iam


         read(1000+iam,*)NZED,NELC


         TAZ=real(NZED-NELC)


         if(TAZ.eq.0.0d0)TAZ=1.0d0


         TAZ=TAZ**2.0d0


	 read(1000+iam,*)NAST,MXE,IDUM1,IDUM2


	 allocate(ISAT(NAST),ILAT(NAST))


	 read(1000+iam,*)(ISAT(I),ILAT(I),I=1,NAST)


	 allocate(ENAT(NAST))


	 read(1000+iam,*)(ENAT(I),I=1,NAST)


	 read(1000+iam,*)NSPN,LRGL,NPTY,NPIEB


	 read(1000+iam,*)EB,MXXE              ! energy of boundstate


C         write(0,*)EB,MXXE


C


C


         do 20 J=1,mxe


	 read(1000+iam,*)aa


	 backspace(1000+iam)


         ee=0.0d0


         ee=ee+1d-7


	     ee=ee+aa+eb


	     icount=0

C       write(0,*)'-----------------------------'
C       write(0,*)'ee/TAZ:  ',(ee/TAZ)


	   do 25 JJ=1,NAST


	    if((ee/TAZ).gt.ENAT(JJ))then

	    icount=icount+1

C         write(0,*)'test1',icount
C         write(0,*)'ENAT: ', ENAT(JJ)


	    endif


c             write(0,*)'ee =',ee/taz,ENAT(JJ),icount,taz,aa,eb


  25       continue


         itmp=min(icount,NPHOTOFINAL)


         icount=itmp


c         icount=icount-1


c         write(0,*)'icount=',icount


         if(icount.eq.0)continue


  40     iflag=1


         K1=icount/6


C


C


	  if(iflag.eq.1)then


	  JK=JK+1


          icounta=min(6,icount)


	 read(1000+iam,*)eee,(photopar(JK,JJJ),JJJ=1,icounta)


C         write(0,*)J,eee,(photopar(JK,JJJ),JJJ=1,icounta)

      if (iam.eq.0) then

               write(0,*) (JK-(MXXE*iam)), icount
      end if

          endif


C


C


  45      icount=icount-6


  46       if(icount.gt.0)then


          icounta=min(6,icount)


          iflag=iflag+6


     	read(1000+iam,*)(photopar(JK,JJJ),JJJ=iflag,iflag+icounta-1)


C         write(0,*)(photopar(JK,JJJ),JJJ=iflag,iflag+icounta)


          goto 45


           endif


  50      continue


   20      continue


C


          deallocate(ISAT,ILAT,ENAT)


  10      continue


C


        do JJ=1,MXE1


          PARSUM=0.0d0


           do JJJJ=istart,ifinal


          PARSUM=PARSUM+PHOTOPAR(POINTER(JJ),JJJJ)


           enddo


          write(9,700)EMESH1(JJ),PARSUM


        enddo


         deallocate(PHOTOPAR)


C        call system("date")


C


        END


        SUBROUTINE qsort(a, n, t)


!    NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL


!    BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.


!    REAL*8 PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.


        IMPLICIT NONE


        INTEGER, INTENT(IN)    :: n


        REAL*8, INTENT(INOUT)    :: a(n)


        INTEGER, INTENT(INOUT) :: t(n)


!    Local Variables


        INTEGER    :: i, j, k, l, r, s, stackl(15), stackr(15), ww


        REAL*8    :: w, x


        s = 1


        stackl(1) = 1


        stackr(1) = n


!    KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.


 10     CONTINUE


        l = stackl(s)


        r = stackr(s)


        s = s - 1


!    KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.


 20     CONTINUE


        i = l


        j = r


        k = (l+r) / 2


        x = a(k)


!    REPEAT UNTIL I > J.


      DO


      DO


       IF (a(i).LT.x) THEN    ! Search from lower end


      i = i + 1


      CYCLE


       ELSE


       EXIT


       END IF


      END DO


      DO


      IF (x.LT.a(j)) THEN    ! Search from upper end


      j = j - 1


      CYCLE


      ELSE


      EXIT


      END IF


      END DO


      IF (i.LE.j) THEN    ! Swap positions i & j


      w = a(i)


      ww = t(i)


      a(i) = a(j)


      t(i) = t(j)


      a(j) = w


      t(j) = ww


      i = i + 1


      j = j - 1


      IF (i.GT.j) EXIT


      ELSE


      EXIT


      END IF


      END DO


      IF (j-l.GE.r-i) THEN


      IF (l.LT.j) THEN


      s = s + 1


      stackl(s) = l


      stackr(s) = j


      END IF


      l = i


      ELSE


      IF (i.LT.r) THEN


      s = s + 1


      stackl(s) = i


      stackr(s) = r


      END IF


      r = j


      END IF


      IF (l.LT.r) GO TO 20


      IF (s.NE.0) GO TO 10


      RETURN


      END SUBROUTINE qsort


C


      FUNCTION NOPEN(E,ENAT,NAST,IONE,IO)


      implicit real*8(a-h,o-z)


C


      DIMENSION ENAT(NAST)


C


      DO I=1,IO,NAST


       IF(E.LT.ENAT(I)) GOTO 1


      ENDDO


      I=NAST+1


 1      IO=I-1


      NOPEN=(IO*(IO-2*IONE+1))/2


C


      RETURN


      END


c     ******************************************************************


c


c     CPB  2005 - this is a hacked version from Don Griffins  stgsig.f


c                 code.


c


c


      real*8 function gaussn(ie,mxe0,mxe,mxeg,n,ewidth,iflg,


     1 emesh,photoprecon,emin,emax)


      implicit real*8(a-h,o-z)


      dimension photoprecon(mxe)


      dimension vec(mxe),emesh(mxe)


C


C


C


      if(emax.lt.emesh(ie))then


        gaussn=photoprecon(ie)


        return


      endif


C


C


C


       if(emin.gt.emesh(ie))then


        gaussn=photoprecon(ie)


        return


      endif


C


      iv=0


      a=1.6651092/ewidth


      do 1 i=1,mxe


      t=(emesh(ie)-emesh(i))*a


      t=t*t


      if(t.gt.70.0)then


        vec(i)=0.0


      else


        if(iflg.eq.0) then


          vec(i)=photoprecon(i)*exp(-t)


        else


	  write(0,*)'should never get here in gaussn'


	  stop


C          vec(i)=sigma(i)*exp(-t)


        end if


      endif


      if(vec(i).gt.0.0) iv=1


   1  continue


      if(iv.gt.0) go to 50


      gaussn=0.0


      return


   50 ifail=0


      emlo=emesh(1)


      emup=emesh(mxe)


      call davint(emesh,vec,mxe,emlo,emup,ans,ifail)


      if(ifail.gt.1) write(6,100) ifail


 100  format(' ifail=',i2)


      gaussn=ans*a*0.5641895


      return


      end


c


c **********************************************************************


c


c **********************************************************************


c


      subroutine davint (x, y, n, xlo, xup, ans, ierr)


c     ***begin prologue  davint


c     ***purpose  integrate a function tabulated at arbitrarily spaced


c        abscissas using overlapping parabolas.


c     ***from slatec library


c     ***category  h2a1b2


c     ***type      double precision


c     ***keywords  integration, quadrature, tabulated data


c     ***author  jones, r. e., (snla)


c     ***description


c


c     abstract


c         davint integrates a function tabulated at arbitrarily spaced


c         abscissas.  the limits of integration need not coincide


c         with the tabulated abscissas.


c


c         a method of overlapping parabolas fitted to the data is used


c         provided that there are at least 3 abscissas between the


c         limits of integration.  davint also handles two special cases.


c         if the limits of integration are equal, davint returns a


c         result of zero regardless of the number of tabulated values.


c         if there are only two function values, davint uses the


c         trapezoid rule.


c


c     description of parameters


c         the user must dimension all arrays appearing in the call list


c         x(n), y(n)


c


c         input--


c      x    - double precision array of abscissas, which must be in


c             increasing order.


c      y    - double precision array of function values. i.e.,


c                y(i)=func(x(i))


c      n    - the integer number of function values supplied.


c                n .ge. 2 unless xlo = xup.


c      xlo  - double precision lower limit of integration


c      xup  - double precision upper limit of integration.  must have


c              xlo.le.xup


c


c         output--


c      ans  - double precision computed approximate value of integral


c      ierr - a status code


c           --normal code


c                =1 means the requested integration was performed.


c           --abnormal codes


c                =2 means xup was less than xlo.


c                =3 means the number of x(i) between xlo and xup


c                   (inclusive) was less than 3 and neither of the two


c                   special cases described in the abstract occurred.


c                   no integration was performed.


c                =4 means the restriction x(i+1).gt.x(i) was violated.


c                =5 means the number n of function values was .lt. 2.


c                   ans is set to zero if ierr=2,3,4,or 5.


c


c     davint is documented completely in sc-m-69-335


c     original program from *numerical integration* by davis & rabinowitz


c     adaptation and modifications by rondall e jones.


c


c     ***references  r. e. jones, approximate integrator of functions


c       tabulated at arbitrarily spaced abscissas,


c      report sc-m-69-335, sandia laboratories, 1969.


c      690901  date written


c      890831  modified array declarations.  (wrb)


c      890831  revision date from version 3.2


c      891214  prologue converted to version 4.0 format.  (bab)


c      900315  calls to xerror changed to calls to xermsg.  (thj)


c      920501  reformatted the references section.  (wrb)


c     ***end prologue  davint


c


      integer i, ierr, inlft, inrt, istart, istop, n


      real*8 a, ans, b, c, ca, cb, cc, fl, fr, r3, rp5,


     1     slope, sum, syl, syl2, syl3, syu, syu2, syu3, term1, term2,


     2     term3, x, x1, x12, x13, x2, x23, x3, xlo, xup, y


      dimension x(*),y(*)


c     begin block permitting ...exits to 190


c        begin block permitting ...exits to 180


c     ***first executable statement  davint


            ierr = 1


            ans = 0.0d0


            if (xlo .gt. xup) go to 160


               if (xlo .eq. xup) go to 150


                  if (n .ge. 2) go to 10


                     ierr = 5


      write(6,200)


  200 format('** in subroutine davint less than two function values'


     +, ' were supplied **')


c     ...............exit


                     go to 190


   10             continue


                  do 20 i = 2, n


c        ............exit


                     if (x(i) .le. x(i-1)) go to 180


c                 ...exit


                     if (x(i) .gt. xup) go to 30


   20             continue


   30             continue


                  if (n .ge. 3) go to 40


c


c                    special n=2 case


                     slope = (y(2) - y(1))/(x(2) - x(1))


                     fl = y(1) + slope*(xlo - x(1))


                     fr = y(2) + slope*(xup - x(2))


                     ans = 0.5d0*(fl + fr)*(xup - xlo)


c     ...............exit


                     go to 190


   40             continue


                  if (x(n-2) .ge. xlo) go to 50


                     ierr = 3


      write(6,210)


  210 format('** in subroutine davint there were less than three'


     +, ' function values between the limits of integration **')


c     ...............exit


                     go to 190


   50             continue


                  if (x(3) .le. xup) go to 60


      write(6,210)


c     ...............exit


                     go to 190


   60             continue


                  i = 1


   70             if (x(i) .ge. xlo) go to 80


                     i = i + 1


                  go to 70


   80             continue


                  inlft = i


                  i = n


   90             if (x(i) .le. xup) go to 100


                     i = i - 1


                  go to 90


  100             continue


                  inrt = i


                  if ((inrt - inlft) .ge. 2) go to 110


                     ierr = 3


      write(6,210)


c     ...............exit


                     go to 190


  110             continue


                  istart = inlft


                  if (inlft .eq. 1) istart = 2


                  istop = inrt


                  if (inrt .eq. n) istop = n - 1


c


                  r3 = 3.0d0


                  rp5 = 0.5d0


                  sum = 0.0d0


                  syl = xlo


                  syl2 = syl*syl


                  syl3 = syl2*syl


c


                  do 140 i = istart, istop


                     x1 = x(i-1)


                     x2 = x(i)


                     x3 = x(i+1)


                     x12 = x1 - x2


                     x13 = x1 - x3


                     x23 = x2 - x3


                     term1 = y(i-1)/(x12*x13)


                     term2 = -y(i)/(x12*x23)


                     term3 = y(i+1)/(x13*x23)


                     a = term1 + term2 + term3


                     b = -(x2 + x3)*term1 - (x1 + x3)*term2


     1                   - (x1 + x2)*term3


                     c = x2*x3*term1 + x1*x3*term2 + x1*x2*term3


                     if (i .gt. istart) go to 120


                        ca = a


                        cb = b


                        cc = c


                     go to 130


  120                continue


                        ca = 0.5d0*(a + ca)


                        cb = 0.5d0*(b + cb)


                        cc = 0.5d0*(c + cc)


  130                continue


                     syu = x2


                     syu2 = syu*syu


                     syu3 = syu2*syu


                     sum = sum + ca*(syu3 - syl3)/r3


     1                     + cb*rp5*(syu2 - syl2) + cc*(syu - syl)


                     ca = a


                     cb = b


                     cc = c


                     syl = syu


                     syl2 = syu2


                     syl3 = syu3


  140             continue


                  syu = xup


                  ans = sum + ca*(syu**3 - syl3)/r3


     1                  + cb*rp5*(syu**2 - syl2) + cc*(syu - syl)


  150          continue


            go to 170


  160       continue


               ierr = 2


      write(6,220)


  220 format('** in subroutine davint the upper limit of integration'


     +, ' was not greater than the lower limit **')


  170       continue


c     ......exit


            go to 190


  180    continue


         ierr = 4


      write(6,230)


  230 format('** in subroutine davint the abscissas were not strictly'


     +, ' increasing.  must have x(i-1) .lt. x(i) for all i **')


  190 continue


      return


      end


c
