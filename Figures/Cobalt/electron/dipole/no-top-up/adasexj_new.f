c     program adasexj.f  version 1.10                        24 SEP 2013
c
c     calculate effective excitation collision strengths for input
c     into ADAS.  This is the code used for intermediate-coupled
c     level-to-level runs.  This version uses the linear interpolation
c     method of burgess and tully to carry out the effective collision
c     strength integration -- since this is a two-point formula, there
c     are no longer any restriction on the energy mesh.
c
      implicit real*8(a-h,o-z)
      real*4 omegin
c     sun and pentium
      real*4 tarry(2)
      LOGICAL BFORM,BBORN,BSHORT
c
      parameter(nmpts=250000,nstmx=500,ntrmx=nstmx*(nstmx+1)/2,ntmp=20,
     1 nlnmx=100,NMPIN=3000)           !NMPIN SHOULD BE AT LEAST 500
      PARAMETER(MXD1=NMPIN/NTRMX,MXD2=NTRMX/NMPIN,MXD3=MXD1+MXD2,
     X          MXDM=NMPIN*MXD1/MXD3+NTRMX*MXD2/MXD3+1)
C
      PARAMETER (TOLH0=1.D-10)
C
      character iontrm*2,iel*2,term(nstmx)*15,trm*15,mant(ntmp+2)*5,
     1 expon(ntmp+2)*3,name*30,date*30,descrip(nlnmx)*74,lbcd(10)*1,
     2 sbcd(10)*1,isbcd*1,ltbcd*1,XTERM(NSTMX)*31
C
      common/blok1/emesh(NMPIN),ej(NMPIN),omega(MXDM),ccex(6),
     1 enat(nstmx),g(nstmx),fjt(nstmx),omgfit(NMPIN),ue(NMPIN),
     2 emsh(nmpts),omegin(NMPIN,ntrmx),elev(nstmx),
     3 ar(nstmx,nstmx),ups(nstmx,nstmx,ntmp)
      common/blok2/jxtwo(nstmx),lat(nstmx),isat(nstmx),index(nstmx)
      common/const/d0,d1,d2,d3,d4,d30,d100,cnvevk,CONVEV,imskp
C
      dimension tkadas(ntmp),NRD(NMPTS)
C
      namelist/adasex/fipot,iontrm,ielas,irmps,numtmp,mxtmp,irdtmp,
     1 nlevs,iel,iptomg,MXEIN,ILOW,EFITL,EFITH,EFITF,GFMIN,IBORN,ITCC,
     2 irates,imskp
C
      data (tkadas(i),i=1,13)/2.0d+02,5.0d+02,1.0d+03,2.0d+03,5.0d+03,
     11.0d+04,2.0d+04,5.0d+04,1.0d+05,2.0d+05,5.0d+05,1.0d+06,2.0d+06/
      data d0,d1,d2,d3,d4,d30,d100,cnvevk/0.0d+00,1.0d+00,2.0d+00,
     1 3.0d+00,4.0d+00,3.0d+01,1.00d+02,1.16045d+04/
      data (lbcd(i),i=1,10)/"S","P","D","F","G","H","I","K","L","M"/
      data (sbcd(i),i=1,7)/'1','2','3','4','5','6','7'/
c
      open(5,file='adasexj.in')
      open(6,file='adasexj.out')
      open(9,file='adf04')
C
      INQUIRE (FILE='omega',EXIST=BFORM)
      IF(BFORM)THEN
        WRITE(6,*)'*** USING FORMATTED omega FILE'
        open(7,file='omega')
        OPEN(99,STATUS='SCRATCH',FORM='FORMATTED')
      ELSE
        INQUIRE (FILE='omegau',EXIST=BFORM)
        IF(BFORM)THEN
          WRITE(6,*)'*** USING UNFORMATTED omegau FILE'
          open(7,file='omegau',FORM='UNFORMATTED')
          BFORM=.FALSE.
        ELSE
          INQUIRE (FILE='OMEGA',EXIST=BFORM)
          IF(BFORM)THEN
            WRITE(6,*)'*** USING FORMATTED OMEGA FILE'
            open(7,file='OMEGA')
            OPEN(99,STATUS='SCRATCH',FORM='FORMATTED')
          ELSE
            INQUIRE (FILE='OMEGAU',EXIST=BFORM)
            IF(BFORM)THEN
              WRITE(6,*)'*** USING UNFORMATTED OMEGAU FILE'
              open(7,file='OMEGAU',FORM='UNFORMATTED')
              BFORM=.FALSE.
            ELSE
              WRITE(6,*)'***ERROR: NO SUITABLE OMEGA FILE FOUND'
              STOP'***ERROR: NO SUITABLE OMEGA FILE FOUND'
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C
      open(8,status='scratch')
c
c     read input from namelist data.
C     THE FIRST 3 SHOULD BE SPECIFIED FOR A MEANINGFUL ADF04 FILE.
c
c     fipot   -- ionization potential in cm-1
c     iontrm  -- term specification for ion - two characters
c     iel     -- symbol for element - two characters
C
c     nlevs   -- of the levels included in the excitation calculation
c                this is the number for which rates are to be calculated
c                and the number for which level information will be
c                written to adf04
c     ielas   -- the omega file FOR NEUTRALS includes the
c                elastic TRANSITIONS WHILE THE IONIC ONE DOES NOT.
c                TO OVERRIDE THIS SET IELAS .EQ. -1 I.E. IF NO ELASTIC
C                PRESENT (NEUTRALS) OR IF ELASTIC PRESENT (IONS).
C                DEFAULT: 0
c     numtmp  -- number of temperatures at which effective collision
c                strengths should be calculated.
C                NOTE, ADAS CANNOT HANDLE NUMTMP.GT.20.
c                default: 10
c     mxtmp   -- if equal to zero, then do not let the temperatures in
c                the run exceed 1/2 maximum energy in omega file; if
c                not equal to zero, no such limit. But a warning will be
c                printed if the temperature exceeds the maximum energy.
c                default: 1
c     irdtmp  -- if not equal to zero, read in numtmp temperatures in 
c                degree Kelvin following the first input line. these
c                temperature will then be used instead of the default
c                adas values. 
c                default: 0 
C     ILOW    -- IF EQUAL TO ZERO AND NO SUITABLE HIGH ENERGY POINTS FOR
C                THE OMEGA FIT THEN UPSILON ZEROED OUT FOR THE TRANSITION.
C             ***ALSO, WHEN SUITABLE HIGH ENERGY POINTS ARE PRESENT USE 
C                LEAST SQUARES FIT FOR DIPOLE OR QUADRUPOLE (USING INFINITE
C                ENERGY LIMIT) OR EXTRAPOLATION  FOR FORBIDDEN.
C             -- IF EQUAL 1 THEN EVALUATE ONLY THE LOW ENERGY CONTRIBUTION
C                TO UPSILON WHEN NO SUITABLE HIGH ENERGY POINTS PRESENT.
C             ***ALSO, WHEN SUITABLE HIGH ENERGY POINTS PRESENT USE 
C                SIMPLE INTERPOLATION (USING DIPOLE OR QUADRUPOLE INFINITE
C                ENERGY LIMIT) OR EXTRAPOLATION (FORBIDDEN).
C                DEFAULT: 1
C     IBORN   -- IF .GT. 0 THEN BORN LIMITS ARE PRESENT AND SO ZERO LIMITS
C                CORRESPOND TO FORBIDDEN TRANSITIONS, WHICH ARE EXTRAPOLATED
C                AS 1/E**N WITH 1.LE.N.LE.IBORN.
C             -- IF .LE. 0 THEN BORN LIMITS ARE NOT PRESENT AND SO
C                TRANSITIONS WITH ZERO LIMITS ARE ALL TAKEN TO BE ALLOWED
C                AND EXTRAPOLATED AS E**IBORN. (RECOMMEND IBORN=0 HERE.)
C                DEFAULT: 2
C     ITCC    -- IF .NE. 0 THEN ASSUMES INPUT OMEGA FROM AN IC CALCULATION
C                E.G. BREIT-PAULI OR ICFT. ENERGIES ARE TAKEN FROM THE
C                OMEGA FILE.
C             -- IF .EQ. 0 THEN ASSUMES INPUT FROM A PURE JK-COUPLING 
C                CALCULATION WITH STGICF. LEVEL ENERGIES ARE TAKEN FROM
C                ADASEXJ.IN SINCE OMEGA ENERGIES ARE DEGENERATE WITHIN A
C                TERM. 
C                DEFAULT: 1
C
C  N.B. ILOW=1 OVERRIDES ALL FITTING
C
C     GFMIN   -- USE LEAST SQAURES FIT ONLY WHEN GF-VALUE .GT. GFMIN
C                DEFAULT: -1.
c     irmps   -- if not equal to zero then the input files is from an
c                RMPS calculation and the fit to the collision strength
c                will begin at an energy above the last spectroscopic
c                level enat(nlevs) rather than above the last pseudo
c                state [enat(nast)].  however, even though there is
c                no longer any restriction on the energy mesh, for an
c                RMPS calculation, it is recommended that you use a
c                single fine mesh through the resonance region and
c                start the coarse mesh just above the last BOUND
c                pseudo state; the fit to the collision strengths for
c                the dipole transitions will then start at the point
c                where the mesh first changes and that will be above
c                all the resonant contributions from bound states
c                (spectroscopic and pseudo).
c                default:0
C     EFITL   -- ONLY USE POINTS ABOVE EFITL (UNSCALED RYDBERGS) TO DETERMINE
C                THE HIGH ENERGY FITS. DEFAULTS TO ENAT(NAST).
C     EFITH   -- ONLY USE POINTS BELOW EFITH (UNSCALED RYDBERGS) TO DETERMINE
C                THE HIGH ENERGY FITS. DEFAULTS TO ALL POINTS.
C     EFITF   -- AS EFITL BUT ONLY FOR FORBIDDEN TRANSITIONS.
C                .LT.0, NO FITTING.
C                .EQ.0, INTERNALLY SET TO 2/3 MAX SCATTERING ENERGY. 
C                .GE.0, RESET TO MAX(EFITF,EFITL,ENAT(NAST)).
C                DEFAULT: 0
C
c     iptomg  -- if GREATER THAN zero print a comparison of the 
c                calculated and fitted collision strengthS
C             -- IF LESS THAN ZERO PRINT THE ACCURACY OF THE FITS ONLY
C             -- IF EQUAL TO ZERO PRINT THE WORST FIT ONLY
c                default: 0
C     MXEIN   -- MAX NUMBER OF ENEGIES READ-IN AT A TIME FROM OMEGA.
C                THE ENERGIES ABOVE ENAT(NAST) SHOULD BE TREATED IN 
C                THE SAME BATCH.
C                DEFAULT: CHOSEN SO AS TO LEAVE AT LEAST 500 POINTS IN
C                THE LAST BATCH.
C
C     IRATES  -- if not equal to zero, then generate maxwell-averaged 
C                excitation/deexcitation rates ... put into 
C                adf04rij,adf04ji files
C                DEFAULT:0
C
C     IMSKP   -- In the calculation of effective collision strengths,
C                skip iskp mesh points between each mesh point in the 
C                resonance region.  this option can be used to generate 
C                adf04 files with effective collision strengths
C                calculated using different numbers of mesh points.  By
C                comparing such files, one can investigate convergence
C                with the density of mesh points.
C                DEFAULT:1
C        
      convev=1.36058d+01
      convcm=1.097373d+05
      conrad=7.972932d+05
      EINF=1.0D6
      delcon=2.0d-01
      emshtol=1.0d-03
      ratmn=1.0d-04
      conbgn=1.01d+00     
      ielas=0
      numtmp=10
      mxtmp=1
      irmps=0
      irdtmp=0
      iptomg=0
      MXEIN=-1
      ILOW=1
      EFITL=D0
      EFITH=1.0D10
      GFMIN=-1.0D0
      IONTRM='  '
      FIPOT=D0
      IEL='  '
      IBORN=2
      EFITF=D0
      ITCC=1
      irates=0
      imskp=1
C
      read(5,adasex)
C
      if(irates.ne.0)then
        open(10,file='adf04rij')
        open(11,file='adf04rji')
      endif
c
      BBORN=IBORN.GT.0
      IF(NUMTMP.GT.NTMP)THEN
        WRITE(6,535)NUMTMP,NTMP
        NUMTMP=NTMP
      ENDIF
      if(irdtmp.ne.0) then
        read(5,*) (tkadas(i),i=1,numtmp)
        IF(NUMTMP.GT.20)WRITE(6,534)
      ELSE
        IF(NUMTMP.GT.13)THEN
          WRITE(6,536)NUMTMP
          NUMTMP=13
        ENDIF
      end if
c
c     read (CONFIGURATION) term specifications AND j-values.
C     ALSO, ITCC=0 (JK-COUPLING)  energies, and index values in energy order.
c  
c      term(i) --   term specification for the ith level- up to 15
c                   characters.  NOTE the total values of 2S+1 and
c                   and L are determined from the last four characters
c                   within term, and these four characters MUST start
c                   with (, and end with ), and MUST be in the form:
c                   (1S),or (3P), or (5D), etc.
c      fjt(i)  --   j-value for the ith level
C
C    CASE ITCC=0 (JK-COUPLING)
c      elev(i)  --  energy  of the ith level in cm-1.
c      index(i) --  this is the level index used to associate the input
c                   levels with the levels read from the omega file
c                   generated using stgicf.  this is especially useful
c                   in a pure jk coupling run where the levels are
c                   ordered first by term and then by increasing j
c                   value.  if the order of the input levels and those
c                   in the omega file are the same then the indices can
c                   be left blank and the code will automatically set
c                   them equal to the level number read in.
C
C    CASE ITCC=1 (BP/IC)
C      ELEV AND INDEX ARE IGNORED. IF THE CONFIGURATION/TERM ORDER
C      IN ADASEXJ.IN DOES NOT MATCH THAT (IMPLIED) IN OMEGA THEN
C      CUT-AND-PASTE THE CORRECT ORDER INTO ADASEXJ.IN.
c
c check for old or new style term/config list
c
      READ(5,501)TRM,TERM(1),TERM(2)
      BACKSPACE(5)
      BSHORT=.TRUE.
      DO I=LEN(TRM),1,-1
        IF(TRM(I:I).EQ.')')GO TO 3
      ENDDO
C
C ASSUME NEW STYLE, AND NOT BLANK OLD-STYLE...
      READ(5,502)TRM,TERM(1),TERM(2)
      BACKSPACE(5)
      TRM=TERM(2)
      DO I=LEN(TRM),1,-1
        IF(TRM(I:I).EQ.')')GO TO 2
      ENDDO
      DO L=1,NLEVS
        READ(5,542)INDEX(L),TERM(L),ISAT(L),LAT(L),FJT(L),ELEV(L) !SHORT
      ENDDO
      GO TO 15
    2 BSHORT=.FALSE.
      DO L=1,NLEVS
        READ(5,543)INDEX(L),XTERM(L),ISAT(L),LAT(L),FJT(L),ELEV(L) !LONG
      ENDDO
      GO TO 15
c
c old style
    3 lnstr=i                                     !since right justified
      do 10 l=1,nlevs
        read(5,500) term(l),fjt(l),elev(l),index(l)
        write(6,500) term(l),fjt(l),elev(l),index(l)
        if(ITCC.NE.0.OR.index(l).eq.0) index(l)=l
        trm=term(l)
c        do i=len(trm),1,-1
c          if(trm(i:i).eq.')')go to 3
c        enddo
c        stop 'term specification blank'
c    3   lnstr=i
        isbcd=trm(lnstr-2:lnstr-2)
        ltbcd=trm(lnstr-1:lnstr-1)
        do i=1,7
          if(sbcd(i).eq.isbcd) then
            isat(l)=i
            go to 6
          end if
        enddo
        stop 'input of terms invalid'
    6   do i=1,10
        if(lbcd(i).eq.ltbcd) then
          lat(l)=i-1
          go to 10
        end if
        enddo
        stop 'input of terms invalid'
   10 enddo    
c
c     read information describing data
c
c      name    -- name of person generating the excitation data - up to
c                 30 characters
c      date    -- date - up to 30 characters
c      descrip -- description of data - up to 72 characters per line
c                 end description with a period in the first column of
c                 the last line.  limited to nlnmx lines, including the
c                 period on the last line.
c
   15 read(5,510) name
      read(5,510) date
      i=1
   20 read(5,530) descrip(i)
      if(descrip(i).eq.'.') go to 30 
      i=i+1
      if(i.gt.nlnmx) then
        write(6,540) nlnmx
        stop 'INCREASE THE PARAMETER NLNMX (OR SHORTEN DESCRIPTION)'
      end if
      go to 20
   30 nmlines=i-1      
c
c     read basic data from file omega
c
      IF(BFORM)read(7,*) nzed,nelc
      IF(.NOT.BFORM)read(7) nzed,nelc
      taz=nzed-nelc
      iq=int(taz)
      iq1=iq+1
      fiq1=iq1
      if(irdtmp.eq.0) then
        do it=1,numtmp
          tkadas(it)=tkadas(it)*(fiq1**2)
        enddo
      end if
      taz=taz*taz
      if(taz.eq.d0) taz=d1
C
      IF(BFORM)read(7,*) nast,mxe,NOMWRT
      IF(.NOT.BFORM)read(7) nast,mxe,NOMWRT
      if(nast.gt.nstmx) then
        write(6,550) nast
        stop 'increase the parameter nstmx'
      end if
      if(nlevs.gt.nast) then
        write(6,555) nast,nlevs
        stop 'NO. OF CC LEVELS DOES NOT MATCH THOSE IN ADASEXJ.IN'
      end if
      if(mxe.gt.nmpts) then
        write(6,560) mxe
        stop 'increase the parameter nmpts'
      end if
C
      IF(BFORM)read(7,*) (idum,jxtwo(i),i=1,nast)
      IF(.NOT.BFORM)read(7) (idum,jxtwo(i),i=1,nast)
      IF(BFORM)read(7,*) (enat(i),i=1,nast)
      IF(.NOT.BFORM)read(7) (enat(i),i=1,nast)
C
      IF(IELAS.LT.0)THEN
        IF(NZED.EQ.NELC)THEN
          NELC=-NELC
        ELSE
          NELC=NZED
        ENDIF
      ENDIF
C
      IF(NZED.EQ.NELC)THEN
        NOMT=(NAST*(NAST+1))/2
        IONE=0
      ELSE
        NOMT=(NAST*(NAST-1))/2
        IONE=1
      ENDIF
      IF(NOMWRT.NE.0)NOMT=IABS(NOMWRT)
C
      IF(NOMWRT.GE.0)THEN
        DO I=1,MXE
          NRD(I)=NOMT
        ENDDO
      ELSE
        IF(BFORM)READ(7,*)E0                 !CHECK FIRST E
        IF(.NOT.BFORM)READ(7)E0                 !CHECK FIRST E
        BACKSPACE(7)
        I=1
        IF(E0.GE.ENAT(NAST))GO TO 1    !ALL OPEN - EARLY EXIT
        IF(BFORM)THEN
          I0=1
          NO=NOPEN(E0,ENAT,NAST,IONE,I0)
          NO=MIN(NO,NOMT)
          READ(7,*)E0,(DUM,N=1,NO)
          READ(7,*)E                   !=E0+EINCR
          EINCR=E-E0
          E0=E0-EINCR
          EINCH=EINCR*0.5D0
          NREC=2+(NO-1)/6
          DO N=1,NREC
            BACKSPACE(7)
          ENDDO
        ENDIF
        I0=1
        DO I=1,MXE
          IF(BFORM)THEN
            E=E0+I*EINCR
            WRITE(99,600)E
            BACKSPACE(99)
            READ(99,600)E
            BACKSPACE(99)
          ELSE
            READ(7)E
          ENDIF
          IF(E.GE.ENAT(NAST))GO TO 1
          NRD(I)=NOPEN(E,ENAT,NAST,IONE,I0)
          NRD(I)=MIN(NRD(I),NOMT)
        ENDDO
        I=MXE+1
   1    I0=I
        DO I=I0,MXE
          NRD(I)=NOMT
        ENDDO
        DO N=1,NOMT
          DO I=1,NMPIN
            OMEGIN(I,N)=D0
          ENDDO
        ENDDO
        IF(.NOT.BFORM)THEN
          REWIND(7)
          NREC=4
          DO N=1,NREC
            READ(7)
          ENDDO
        ENDIF
      ENDIF
c
c     set the lowest energy for the least squares fit to the highest
c     energy threshold
c
      if(irmps.eq.0) then
        EFITL=MAX(enat(nast)*taz,EFITL)*convev
      else
        EFITL=MAX(enat(nlevs)*taz,EFITL)*convev
      end if
c
c     print level information in the ADAS format
c
      write(9,570) iel,iq,nzed,iq1,fipot,iontrm
      if(irates.ne.0)then
        write(10,570) iel,iq,nzed,iq1,fipot,iontrm
        write(11,570) iel,iq,nzed,iq1,fipot,iontrm
      endif
      do l=1,nlevs
        g(l)=d2*fjt(l)+d1
        IWJ=NINT(G(L))
        lx=index(l)
        IF(IWJ.NE.IABS(JXTWO(lx)))THEN
          WRITE(6,575)L,IWJ,lx,IABS(JXTWO(lx))
          STOP 'MIS-MATCH OF TARGET LEVELS: ADASEXJ.IN AND OMEGA'
        ENDIF
        IF(ITCC.NE.0)ELEV(L)=ENAT(L)*TAZ*CONVCM
        IF(BSHORT)THEN
          write(9,580) l,term(l),isat(l),lat(l),fjt(l),elev(l)
          if(irates.ne.0)then
            write(10,580) l,term(l),isat(l),lat(l),fjt(l),elev(l)
            write(11,580) l,term(l),isat(l),lat(l),fjt(l),elev(l)
          endif
        ELSE
          WRITE(9,581) L,XTERM(L),ISAT(L),LAT(L),FJT(L),ELEV(L)
          IF(IRATES.NE.0)THEN
            WRITE(10,581) L,XTERM(L),ISAT(L),LAT(L),FJT(L),ELEV(L)
            WRITE(11,581) L,XTERM(L),ISAT(L),LAT(L),FJT(L),ELEV(L)
          ENDIF
        ENDIF
      enddo
      write(9,590)
      if(irates.ne.0)then
        write(10,590)
        write(11,590)
      endif
C
C     SET-UP BATCHES OF ENERGIES TO BE READ SO AS NOT TO EXCEED
C     THE NMPIN DIMENSION OF OMEGIN. (MXE.LE.NMPTS STILL REQUIRED)
C
      IF(MXE.LE.NMPIN)THEN
        NBATCH=1
        MBATCH=MXE
      ELSE
        IF(MXEIN.GT.NMPIN)THEN
          WRITE(6,*)'USER INPUT ERROR: REDUCE MXEIN TO ',NMPIN
          STOP 'REDUCE MXEIN'
        ENDIF
        III=1
        IF(BFORM)III=0
        IF(MXEIN.LE.0)MXEIN=NMPIN
        MBATCH=MXEIN-III
        IF(ILOW.EQ.0.OR.IBORN.GT.0)THEN
   55     NBATCH=1+(MXE-1-III)/MBATCH
          MNB=MBATCH*(NBATCH-1)
          MHIGH=MXE-MNB
          IF(MHIGH.GT.0.AND.MHIGH.LT.500)THEN
            MBATCH=MBATCH-MBATCH/5
            IF(NBATCH.GT.100.OR.MBATCH.LT.500)THEN
              WRITE(6,*)' BATCH ENERGY SET-UP OUT OF CONTROL'
              WRITE (6,*)'TRY OVERRIDING WITH INPUT MXEIN'
              STOP ' BATCH ENERGY SET-UP OUT OF CONTROL'
            ENDIF
            GO TO 55
          ENDIF
        ENDIF
      ENDIF
C
C     INITIALIZE UPSILON
C
      DO K=1,NUMTMP
        DO J=1,NLEVS
          DO I=1,NLEVS
            UPS(I,J,K)=D0
          ENDDO
        ENDDO
      ENDDO
c
c     read electron-impact excitation collision strengths
c
      IE1=0
   60 IE0=IE1+1
      IE1=IE1+MBATCH
      IE1=MIN(IE1,MXE)
      IOM=0
c      write(77,*)'*****',ie0,ie1,mbatch
      IF(BFORM)THEN
        do i=IE0,IE1
          IOM=IOM+1
          read(7,600) emsh(i),(omegin(IOM,n),n=1,NRD(I))
c        write(77,601)i,emsh(i)
c 601    format(i6,e18.8)
          IF(NOMWRT.LT.0)THEN
            ET=E0+I*EINCR
            IF(ABS(EMSH(I)-ET).GT.EINCH)THEN
              IF(NRD(I).LT.NOMT)THEN
              WRITE(6,*)'*** MIS-MATCH DURING READ OF FORMATTED OMEGAS'
              WRITE(6,*)I,NRD(I),ET,EMSH(I)
              STOP'*** MIS-MATCH DURING READ OF FORMATTED OMEGAS'
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE
        DO I=IE0,IE1
          IOM=IOM+1
          read(7) emsh(i),(omega(n),n=1,NRD(I))           !OMEGA*8
c        write(77,601)i,emsh(i)
          DO N=1,NRD(I)
            omegin(IOM,n)=OMEGA(N)                          !OMEGIN*4
          ENDDO
        ENDDO
      ENDIF
C
C TEST ENERGY MESH FOR SEQUENTIAL POINTS
C
      TOLH=TOLH0*TAZ*CONVEV
      DO I=IE0+1,IE1
        H=EMSH(I)-EMSH(I-1)
        IF(H.LT.-TOLH0)THEN
          WRITE(6,*)'*** MESH ERROR:',I-1,EMSH(I-1),I,EMSH(I)
          STOP '*** MESH ERROR'
        ENDIF
      ENDDO
C
      IFIN=IE1
      IF(IE1.EQ.MXE)THEN
        IF(ILOW.EQ.0.AND.EMSH(IE0).GT.ENAT(NAST))THEN
          WRITE(6,*)'NEED ALL ENERGIES ABOVE ENAT(NAST) IN SAME BATCH'
          WRITE(6,*)'TRY SETTING MXEIN = ',2*MXEIN/3
          STOP 'NEED ALL ENERGIES ABOVE ENAT(NAST) IN THE SAME BATCH'
        ENDIF
        if(emsh(mxe).LT.EINF.AND.EMSH(MXE).GE.1.5*ENAT(NAST)) then
          write(6,610) 
        stop 'ERROR: INFINITE ENERGY INFO MISSING'
        end if
C
        ifin=mxe-1
c
c     test energy mesh to see if a higher value of efit can be used
c     where the mesh spacing changes the first time
c
        hold=emsh(IE0+1)-emsh(IE0)
        do i=IE0+2,ifin
          hnew=emsh(i)-emsh(i-1)
          hdif=abs(hnew-hold)
          havg=(hnew+hold)/d2
          ratio=hdif/havg
          if(ratio.gt.emshtol) then
            efitnew=emsh(i-1)*taz*convev
            if(efitnew.gt.EFITL) EFITL=efitnew
            go to 75       
          end if
        enddo
   75   EFITH=EFITH*CONVEV
        IF(EFITH.LT.1.E5)WRITE(6,614)EFITL,EFITH
        IF(EFITF.EQ.D0)EFITF=EMSH(IFIN)*TAZ*2/3
        EFITF=EFITF*CONVEV
        IF(EFITF.GE.D0)EFITF=MAX(EFITF,EFITL)
c
c     check to see that the maximum temperature does
c     not exceed half the maximum energy in the omega.  
c     if it does, then if mxtmp.eq.0
c     reduce the maximum temperature
c     else warn if 
c
        emxtmp=emsh(ifin)*taz*convev*cnvevk/d2
        do it=1,numtmp
          if(emxtmp.lt.tkadas(it)) go to 615
        end do
        go to 617
  615   if(mxtmp.eq.0) then                             !reduce
          write(6,616) tkadas(numtmp),tkadas(it-1)
          numtmp=it-1
        else                                         !just warn
          if(emxtmp.lt.tkadas(numtmp)/d2)
     x       write(6,618) tkadas(numtmp),emxtmp*d2
        endif
c
c     write temperatures for ADAS
c
  617   write(8,630) (tkadas(it),it=1,numtmp)
        rewind (8)
        read(8,640) (mant(it),expon(it),it=1,numtmp)
        itype=3
        write(9,650) fiq1,itype,(mant(it),expon(it),it=1,numtmp)
        if(irates.ne.0)then
          itype=0                                  !not yet defined
          write(10,650) fiq1,itype,(mant(it),expon(it),it=1,numtmp)
          write(11,650) fiq1,itype,(mant(it),expon(it),it=1,numtmp)
        endif
      ENDIF
c
      li0=0
      lf0=0
      avgpd0=d0
c
c *** begin loops over levels   ***
c
      do 150 li=1,nlevs-1
      loi=index(li)
      ei=enat(loi)
      do 140 lf=li+1,nlevs
      lof=index(lf)
c      if(enat(lof).gt.enat(loi)) then
      if(itcc.ne.0) then
        ef=enat(lof)
        ei=enat(loi)
        emshmn=ef
      else
        ef=elev(lf)/(convcm*taz)
        ei=elev(li)/(convcm*taz)
        emshmn=ef
        if(enat(loi).gt.ef) emshmn=enat(loi)
      end if
      ethrsh=ef-ei
      eavg=(ef+ei)/d2
      ratio=ethrsh/eavg
      fdeg=d0
      if(ratio.lt.ratmn) then
        ethrsh=ratmn*eavg
        eadd=conbgn*ethrsh
        fdeg=d1
      end if
      ethrsh=ethrsh*taz*convev
      if(loi.lt.lof) then
        ILI=loi
        ILF=lof
      else
        ILI=lof
        ILF=loi
      end if
      IF(NOMWRT.LT.0)NTR=((ILF-1)*(ILF-2*IONE))/2+ILI
      IF(NOMWRT.GT.0)NTR=ILF+NAST*(ILI-1)-(ILI*(ILI-1+2*IONE))/2
      IF(NTR.GT.NOMT)GO TO 140
c
C     SET-UP INFINITE ENERGY LIMIT IF LAST OMEGA IS NON-ZERO
C     AND SET TRANSITION TYPE.
c
      IF(IE1.EQ.MXE)THEN
        idip=0
        C1=D0
        c5=d0
        ar(li,lf)=d0
        ar(lf,li)=d0
        omgasym=omegin(mxe-IE0+1,ntr)
        if(omgasym.ne.d0) then
          IF(BBORN.AND.OMGASYM.LT.D0.OR..NOT.BBORN)THEN
            ISPN=-1
            idip=1
            c5=ABS(omgasym)/log(EMSH(MXE)*taz)
            dele=(elev(lf)-elev(li))*convev/convcm
            ar(li,lf)=conrad*(dele**3)*c5/g(lf)
            AR(LF,LI)=-C5
          ELSE
            ISPN=0
            C1=OMGASYM
            AR(LF,LI)=C1
          ENDIF
        ELSE
          ISPN=IABS(IBORN)
C TO ALLOW FIT TO FORBIDDEN, WITH ZERO LIMIT POINT, UNCOMMENT NEXT LINE
          IF(BBORN.AND.EFITF.GT.0.)C1=1.D-20
        end if
      ENDIF
c
c     SET LIMIT POINTS FOR UPEX2 INTERPOLATION.
c
      CCEX(1)=C1 
      CCEX(5)=C5
c
c     begin loop over mesh energy
c
      ie=0
      ieft=0
      IOM=0
      EFITL0=EFITL
      IF(ISPN.GT.0.AND.EFITF.GT.0.)EFITL0=EFITF
      emshmn=emshmn+fdeg*eadd
      do 80 i=IE0,ifin
        IOM=IOM+1
        if(emsh(i).le.emshmn) go to 80
        ie=ie+1
        emesh(ie)=(emsh(i)-ei)*taz*convev
        ej(ie)=(emsh(i)-ef)*taz*convev
        emshev=emsh(i)*taz*convev
        omega(ie)=omegin(IOM,ntr)
        if(emshev.gt.EFITL0.AND.EMSHEV.LT.EFITH) then
          ieft=ieft+1
          omgfit(ieft)=omega(ie)
          ue(ieft)=emesh(ie)/ethrsh 
        end if
   80 enddo
      mxetot=ie
C
      IF(IE1.EQ.MXE)THEN
        if(IEFT.EQ.0) then
          IF(ILOW.EQ.1)THEN
            WRITE(6,*)'***NO HIGH ENERGY POINTS, UPSILON INCOMPLETE'
     X,               ' FOR TRANSITION',LI,'  TO',LF
            GO TO 120
          ENDIF
          WRITE(6,*)'***NO HIGH ENERGY POINTS, UPSILON SET TO ZERO'
     X,             ' FOR TRANSITION',LI,'  TO',LF
          ar(li,lf)=d0
          ar(lF,lI)=d0
          do it=1,numtmp
            ups(li,lf,it)=d0
          enddo
          go to 140
        end if
        mxe2=ieft
c
c     set coefficient for INTER/extrapolation of collision strengths
c     above the last mesh point.  
c
        casym=d0
        UASYM=UE(MXE2)
        ALF=D0
        if(idip.eq.0) THEN
          if(ispn.EQ.0) then
            IF(ILOW.EQ.1.OR.C1.EQ.D0)casym=omgfit(mxe2)
          elseIF(ISPN.GT.0)THEN
            IF(ILOW.EQ.1.OR.C1.EQ.D0)THEN
      if(irmps.ge.0)then
              IF(MXE2.LT.10)THEN
                WRITE(6,*)'***ERROR, INSUFFICIENT POINTS IN LAST BATCH,'
     X                    ,'OR TOO FEW HIGH ENERGY POINTS: MXE2=',MXE2
                STOP '***ERROR, INSUFFICIENT HIGH ENERGY POINTS'
              ENDIF
              MXE1=MXE2-MXE2/3
c              TI=LOG(OMGFIT(MXE1)/OMGFIT(MXE2))/LOG(UE(MXE2)/UE(MXE1))
              TJ=LOG(OMGFIT(MXE1)/OMGFIT(MXE2))
     X          /LOG((UE(MXE2)-D1)/(UE(MXE1)-D1))
              ISPN0=NINT(TJ)
c              WRITE(6,*)LI,LF,ISPN0,TI,TJ
              IF(ISPN0.NE.ISPN)THEN
                IF(ISPN0.LT.ISPN)THEN
                  ISPN=ISPN-1
                  T=ISPN
                  ALF=MAX(T,TJ)
                ENDIF
                IF(ISPN0.GT.ISPN)THEN
                  ISPN=ISPN+1
                  T=ISPN
                  ALF=MIN(T,TJ)
                ENDIF
              ELSE
                ALF=TJ
              ENDIF
      else
        alf=ispn
      endif
              casym=omgfit(mxe2)*(ue(mxe2)-D1)**ALF
            ENDIF
          ELSE
            STOP 'DIPOLE TRANSITION, SHOULD NOT BE HERE...'
          end if
        ELSE
          IF(ILOW.EQ.1.OR.C5.LT.GFMIN)
     X       CASYM=OMGFIT(MXE2)/LOG(UE(MXE2)-D1+EXP(D1))
        ENDIF
c
c     skip the least-squares fit on omegas if CASYM.NE.0
c
        if(CASYM.NE.D0) go to 120
C
C
        if(mxe2.lt.6) then
          write(6,660) mxe2,li,lf
          stop 'INSUFFICIENT POINTS FOR OMEGFIT'
        end if
C
C     FIT OMEGA
C
        call omegfit(omgfit,ue,mxe2,ccex,idip,C1,c5,ifail)
C
        if(ifail.ne.0) THEN
          WRITE(6,*)'FAILURE IN OMEGFIT FOR TRANSITION',LI,'  TO',LF
          stop 'FAILURE IN OMEGFIT'
        ENDIF
c
c     print results on cross section fit
c 
        if(iptomg.GT.0) then
          write(6,670) li,lf,(i,ccex(i),i=1,6)
          write(6,680)
        end if
        avgpd=d0
        do ie=1,mxe2
          e=ue(ie)*ethrsh
          omgcal=ccex(1)+ccex(2)/ue(ie)+ccex(3)/ue(ie)**2
     1          +ccex(4)/ue(ie)**3+ccex(5)*log(ue(ie))
     2          +ccex(6)*log(ue(ie))/ue(ie)
          pdiff=d100*abs((omgfit(ie)-omgcal)/(omgfit(ie)+omgcal))
          avgpd=avgpd+pdiff
          if(iptomg.GT.0) write(6,690) e,omgcal,omgfit(ie),pdiff
        enddo
        fmxe2=mxe2
        avgpd=avgpd/fmxe2
        if(avgpd.gt.avgpd0)then
          avgpd0=avgpd
          li0=li
          lf0=lf
        endif
        if(iptomg.LT.0) write(6,700) li,lf,avgpd
c
c     calculate and print the extrapolated values of omega
c
        IF(IPTOMG.GT.0)THEN
          write(6,710) avgpd
          u=ue(mxe2)
          delu=delcon*u
          do iu=1,10
            u=u+delu
            e=u*ethrsh
            omgcal=ccex(1)+ccex(2)/u+ccex(3)/u**2+ccex(4)/u**3
     1            +ccex(5)*log(u)+ccex(6)*log(u)/u
            write(6,720) e,omgcal
            if(omgcal.lt.d0) THEN
              WRITE(6,721)LI,LF,E
              GO TO 120
            ENDIF
            IF(IDIP.EQ.0)THEN
              IF(ISPN.GT.0)THEN
                OMXTRP=OMGFIT(MXE2)*(UE(MXE2)-D1)**ALF
                OMXTRP=OMXTRP/(U-D1)**ALF
              ELSE
                OMXTRP=OMGFIT(MXE2)
              ENDIF
              IF(OMGCAL.GT.2.0D0*OMXTRP)THEN
                WRITE(6,722)IU,LI,LF
                GO TO 120
              ENDIF
            ELSE
              OMXTRP=OMGFIT(MXE2)/LOG(UE(MXE2)-D1+EXP(D1))
              OMXTRP=OMXTRP*LOG(U-D1+EXP(D1))
              IF(OMGCAL.GT.1.5D0*OMXTRP)THEN
                WRITE(6,723)IU,LI,LF
                GO TO 120
              ENDIF
            ENDIF
          enddo
        ENDIF
c
      ENDIF
c
c     calculate and print effective collision strengths for excitation
c     as a function of temperature
c
  120 iunit=1
      do it=1,numtmp
        temp=tkadas(it)
        ups1=upex1(ej,omega,1,mxetot,temp,iunit,TOLH)
        ups(li,lf,it)=ups(li,lf,it)+ups1
      ENDDO
      IF(IE1.EQ.MXE.AND.IEFT.GT.0)THEN
        emin=emesh(mxetot)
        IFAIL=0
        do it=1,numtmp
          IP=IT/NUMTMP
          IP=IP*IPTOMG
          temp=tkadas(it)
          ups2=upex2(emin,ethrsh,ccex,idip,temp,iunit,UASYM,casym,ispn
     X              ,ALF,LI,LF,IFAIL,IP)
          ups(li,lf,it)=ups(li,lf,it)+ups2
        ENDDO
      ENDIF
C
  140 continue
  150 continue
C
C *** END LOOP OVER LEVELS ***
C
      IF(IE1.LT.MXE)THEN
        IF(.NOT.BFORM)THEN        !MINOR REFINEMENT
          BACKSPACE(7)
          IE1=IE1-1
        ENDIF
        GO TO 60
      ENDIF
C
C     END OF ENERGIES
c
      IF(AVGPD0.GT.0.)THEN
        write(6,725)
        write(6,700)li0,lf0,avgpd0
      ENDIF
c
      do li=1,nlevs-1
        do lf=li+1,nlevs
          rewind(8)
          ar(li,lf)=max(ar(li,lf),1.0d-30)
          write(8,730) ar(li,lf),(ups(li,lf,it),it=1,numtmp),AR(LF,LI)
          rewind(8)
          read(8,740) (mant(i),expon(i),i=1,numtmp+2)
          write(9,750) lf,li,(mant(i),expon(i),i=1,numtmp+2)    !eff col
C
          if(irates.ne.0)then
            do iiii=1,numtmp
              ups(li,lf,iiii)=2.1716E-8*ups(li,lf,iiii)/real(JXTWO(loi))
              EDELIJ=(ef-ei)*TAZ
              tempkel=tkadas(iiii)/cnvevk
              tempkel=tempkel/13.6058
              tmpa=exp(-(EDELIJ/tempkel))
              tmpb=sqrt((1.0d0/tempkel))
              ups(li,lf,iiii)=ups(li,lf,iiii)*tmpa*tmpb
            enddo
c            write(10,*) lf,li,ar(li,lf),(ups(li,lf,it),it=1,numtmp)  ! q_i-j 
C
            rewind(8)
            write(8,730) ar(li,lf),(ups(li,lf,it),it=1,numtmp),AR(LF,LI)
            rewind(8)
            read(8,740) (mant(i),expon(i),i=1,numtmp+2)
            write(10,750) lf,li,(mant(i),expon(i),i=1,numtmp+2)    !eff col
C
            do iiii=1,numtmp
              tempkel=tkadas(iiii)/cnvevk
              tempkel=tempkel/13.6058
              tmpa=exp((EDELIJ/tempkel))
              ups(li,lf,iiii)=ups(li,lf,iiii)*
     X                   (real(JXTWO(loi))/real(JXTWO(lof)))*tmpa
            enddo
c            write(11,*) lf,li,ar(li,lf),(ups(li,lf,it),it=1,numtmp)  ! q_j-i
C
            rewind(8)
            write(8,730) ar(li,lf),(ups(li,lf,it),it=1,numtmp),AR(LF,LI)
            rewind(8)
            read(8,740) (mant(i),expon(i),i=1,numtmp+2)
            write(11,750) lf,li,(mant(i),expon(i),i=1,numtmp+2)    !eff col
	  endif 
        ENDDO
      ENDDO
c
c     print description of data
c
      write(9,760)
      write(9,770) 
      if(nmlines.gt.0) then
        do i=1,nmlines
          write(9,780) descrip(i)
        enddo
      end if
      write(9,790) name,date
      if(irates.ne.0)then
        write(10,760)
        write(10,770)
        write(11,760)
        write(11,770)
        if(nmlines.gt.0) then
          do i=1,nmlines
            write(10,780) descrip(i)
            write(11,780) descrip(i)
          enddo
        end if
        write(10,790) name,date
        write(11,790) name,date
      endif
c
c     sun and pentium
      dum=dtime(tarry)
      time=tarry(1)
c
c     cray
c     call second(time)
c
      TIME=TIME/60.0
      write(6,800) time
C
C FORMATS
C
  500 format(a15,f5.1,f11.0,i5)
  501 FORMAT(3A15)
  502 FORMAT(5X,3A15)
  510 format(a30)
  530 format(a72)
  534 FORMAT(/5X,'***ATTENTION, THE NUMBER OF TEMPERATURES EXCEEDS'
     X ,' THE MAX (20) THAT CAN BE HANDLED BY ADAS'/5X,'***THE ADF04'
     X ,' SHOULD NOT BE PUT INTO THE ADAS DATABASE.'/)
  535 FORMAT(/5X,'***ATTENTION, THE NUMBER OF TEMPERATURES HAS BEEN'
     X ,' REDUCED ON DIMENSION GROUNDS FROM',I4,' TO',I4/)
  536 FORMAT(/5X,'***ATTENTION, THE NUMBER OF TEMPERATURES HAS BEEN'
     X ,' REDUCED TO THE ADAS DEFAULT I.E. FROM',I4,' TO  13'/)
  540 format(/5x,
     1 '********************** FATAL INPUT ERROR! *******************'/ 
     2 5x,'the description is longer than ',i2,' lines. increase the '/
     3 5x,'value of nlnmx in the code or shorten the description'/5x,
     4 '************************************************************'/)
  542 FORMAT(I5,1X,A15,4X,I1,1X,I1,1X,F4.1,1X,F21.4)
  543 FORMAT(I5,1X,A31,4X,I1,1X,I1,1X,F4.1,1X,F21.4)
  550 format(/5x,
     1 '****************** FATAL OMEGA INPUT ERROR! ****************'/ 
     2 5x,'the number of states in the close-coupling calculation is '/
     2 5x,'larger than the current dimension.  increase the parameter'/
     3 5x,' nstmx to ',i3/5x,
     4 '************************************************************'/)
  555 format(/5x,
     1 '********************** FATAL INPUT ERROR! *******************'/
     2 5x,'the number of states in the close-coupling calculation is',
     3 ' equal to ',i5/5x,'but the number of levels in the file',
     4 ' adasexj.in is equal to ',i5)
  560 format(/5x,
     1 '****************** FATAL OMEGA INPUT ERROR! ****************'/ 
     2 5x,'the number of energy points in the mesh is larger than the'/
     2 5x,' current dimension. increase the parameter nmpts to ',i6/5x,
     4 '************************************************************'/)
  570 format(a2,'+',i2,2i10,f15.0,'(',a2,')')
  575 FORMAT('MIS-MATCH OF TARGET LEVELS: ADASEXJ.IN AND OMEGA'/4I5)
cold  580 format(i5,3x,a15,' (',i1,')',i1,'(',f4.1,')',f17.0)
  580 FORMAT(I5,1X,A15,2X,' (',I1,')',I1,'(',F4.1,')',F21.4)
  581 FORMAT(I5,1X,A31,2X,' (',I1,')',I1,'(',F4.1,')',F21.4)
  590 format(3x,'-1')
  600 format(1Pe14.8,6(e11.3)/(14x,6(e11.3)))
  610 format(/5x,
     1 '****************** FATAL OMEGA INPUT ERROR! ****************'/ 
     2 5x,'the omega file must contain omegas for dipole transitions',
     3 5x,' calculated from a long range dipole approximation',5x,
     4 '************************************************************'/)
  614 FORMAT(' FITTING CARRIED-OUT IN THE ENERGY RANGE (EV): ',2F10.3/)
  616 format(/5x,'****************************** WARNING!',
     1       '******************************',/
     2       5x,'The maximum temperature was decreased from ',1pe9.2,
     3       ' K to ',e9.2,' K'/
     4       5x,'In light of the maximum energy contained in the omega',
     5       ' file, a higher',/
     6       5x,'temperature could yield inaccurate effective collision'
     7       ,' strengths',/
     9       5x,'*****************************************************',
     *       '*****************'/) 
  618 format(/5x,'****************************** WARNING!',
     1       '******************************',/
     2       5x,'The maximum temperature should be decreased from ',
     3       1pe9.2,' K to ',e9.2,' K'/
     4       5x,'In light of the maximum energy contained in the omega',
     5       ' file, a higher',/
     6       5x,'temperature could yield inaccurate effective collision'
     7       ,' strengths',/
     9       5x,'*****************************************************',
     *       '*****************'/) 
  630 format(20(1pe9.2))
  640 format(20(a5,1x,a3))
  650 format(f5.2,4x,i1,6x,20(a5,a3))
  660 format(5x/
     1 '*********************** FATAL ERROR! ***********************'/
     2 5x,i1,' is not enough points to perform a least squares fit', 
     3 ' for'/5x,'the transition from level number ',i3,' to level ',
     4 'number ',i3,5x/
     5 '************************************************************')
  670 format(/'Fitting Coefficients For Transition ',i5,' To ',
     1   i5/3(3x,'c(',i1,')= ',e11.4)/
     2   3(3x,'c(',i1,')= ',e13.6))
  680 format(/'       Energy          Omega from fit   ',
     1   '  Calculated Omega   Percent Difference')
  690 format(6x,f8.2,10x,f12.6,8x,f12.6,10x,f8.4)
  700 format(/3x,'For the transition from level',i4,' to level',i4,
     1   ' the average percent difference'/3x,'between the calculated',
     2   ' and fitted omega is equal to ',f8.4)
  710 format(66x,'________',/66x,f8.4/)
  720 format(6x,f8.2,10x,f12.6)
  721 FORMAT('***EXTRAPOLATED OMEGA HAS GONE NEGATIVE FOR TRANSTION: '
     X,2I5,' AT ENERGY ',F10.3)
  722  FORMAT('NON-DIPOLE EXTRAPOLATION MAYBE DIVERGING AT IU=',I2,
     X      '  FOR TRANSITION: ', 2I5)
  723  FORMAT('DIPOLE EXTRAPOLATION MAYBE DIVERGING AT IU=',I2,
     X       '  FOR TRANSITION: ',  2I5)           
  725 format(/3x,'The ** worst ** fit is given by:')
  730 format(22(1pe9.2))
  740 format(22(a5,1x,a3))
  750 format(2i4,22(a5,a3))
  760 format(2x,'-1'/2x,'-1  -1')
  770 format('C',79('-')/'C')
  780 format('C',5x,a74)
  790 format('C'/'C',5x,a30,4x,a30/'C',79('-'))
  800 format(/'cpu time =',f9.3, ' min')
      end
c
c     ******************************************************************
c
      subroutine omegfit(omg,u,nfit,c,idip,C1,c5,ifail)
c
c     subroutine to perform least squares fit of omega to e in threshold
c     units.  
c     if dipole use omega=c1+c2/u+c3/u**2+c4/u**3+c5*log(u)+c6*log(u)/u,
c     with c5 determined from a long range dipole approximation for
c     omega.  if non-dipole, use omega=c1+c2/u+c3/u**2+c4/u**3
C     WHERE C1 MAY BE DETERMINED FROM BORN LIMIT.
C
      implicit real*8(a-h,o-z)
c
      common/const/d0,d1,d2,d3,d4,d30,d100,cnvevk,CONVEV,imskp
      dimension omg(nfit),u(nfit),c(6)
      dimension cfit(5),fna(5,5),fnb(5),p(5)
c
      do i=1,5
        cfit(i)=d0
        fnb(i)=d0
      enddo
      do j=1,5
        do i=1,5
          fna(i,j)=d0
        enddo
      enddo
c
      do n=1,nfit
c
        fx1=d1/u(n)
        fx2=fx1*fx1
        fx3=fx1*fx2
        if(idip.eq.0) then
          IF(C1.EQ.D0)THEN
            isz=4
            fna(1,1)=fna(1,1)+d1
            fna(1,2)=fna(1,2)+fx1
            fna(1,3)=fna(1,3)+fx2
            fna(1,4)=fna(1,4)+fx3
            fna(2,2)=fna(2,2)+fx1*fx1
            fna(2,3)=fna(2,3)+fx1*fx2
            fna(2,4)=fna(2,4)+fx1*fx3
            fna(3,3)=fna(3,3)+fx2*fx2
            fna(3,4)=fna(3,4)+fx2*fx3
            fna(4,4)=fna(4,4)+fx3*fx3
            fy=omg(n)
            fnb(1)=fnb(1)+fy
            fnb(2)=fnb(2)+fy*fx1
            fnb(3)=fnb(3)+fy*fx2
            fnb(4)=fnb(4)+fy*fx3
          ELSE
            isz=3
            fna(1,1)=fna(1,1)+fx1*fx1
            fna(1,2)=fna(1,2)+fx1*fx2
            fna(1,3)=fna(1,3)+fx1*fx3
            fna(2,2)=fna(2,2)+fx2*fx2
            fna(2,3)=fna(2,3)+fx2*fx3
            fna(3,3)=fna(3,3)+fx3*fx3
            fy=omg(n)-C1
            fnb(1)=fnb(1)+fy*fx1
            fnb(2)=fnb(2)+fy*fx2
            fnb(3)=fnb(3)+fy*fx3
          ENDIF
        else               !ASSUME C5 NON-ZERO
          isz=5
          fx4=log(u(n))*fx1
          fna(1,1)=fna(1,1)+d1
          fna(1,2)=fna(1,2)+fx1
          fna(1,3)=fna(1,3)+fx2
          fna(1,4)=fna(1,4)+fx3
          fna(1,5)=fna(1,5)+fx4
          fna(2,2)=fna(2,2)+fx1*fx1
          fna(2,3)=fna(2,3)+fx1*fx2
          fna(2,4)=fna(2,4)+fx1*fx3
          fna(2,5)=fna(2,5)+fx1*fx4
          fna(3,3)=fna(3,3)+fx2*fx2
          fna(3,4)=fna(3,4)+fx2*fx3
          fna(3,5)=fna(3,5)+fx2*fx4
          fna(4,4)=fna(4,4)+fx3*fx3
          fna(4,5)=fna(4,5)+fx3*fx4
          fna(5,5)=fna(5,5)+fx4*fx4
          fy=omg(n)-c5*log(u(n))
          fnb(1)=fnb(1)+fy
          fnb(2)=fnb(2)+fy*fx1
          fnb(3)=fnb(3)+fy*fx2
          fnb(4)=fnb(4)+fy*fx3
          fnb(5)=fnb(5)+fy*fx4
        end if
c
      enddo
c
c     call f04ase(fna,5,fnb,isz,cfit,w1,w2,ifail)
      call choldc(fna,isz,5,p,ifail)
c
      if(ifail.ne.0) then
        write(6,220)
        return
      end if
c
      call cholsl(fna,isz,5,p,fnb,cfit)
c
      IF(ISZ.EQ.3)THEN
        c(1)=C1
        c(2)=cfit(1)
        c(3)=cfit(2)
        c(4)=cfit(3)
        C(5)=D0
        C(6)=D0
      ELSEIF(ISZ.EQ.4)THEN
        c(1)=cfit(1)
        c(2)=cfit(2)
        c(3)=cfit(3)
        c(4)=cfit(4)
        C(5)=D0
        C(6)=D0
      ELSE
        c(1)=cfit(1)
        c(2)=cfit(2)
        c(3)=cfit(3)
        c(4)=cfit(4)
        c(5)=c5
        c(6)=cfit(5)
      ENDIF
c
      return
c
  220   format(/5x,
     1 '*********************** FATAL ERROR! ***********************'/
     2 5x,'the cholesky routine used to solve the normal equations',
     3 ' for'/5x,' the cross section fit failed'/5x,
     4 '************************************************************'/)
      end
c
c     ******************************************************************
c
      real*8 function upex1(e,omg,iei,ief,temp,iunit,TOLH)
c
c     calculate the contribution to the effective collision strength
c     from that portion of the collision strength with energies
c     between e(iei) and e(ief) using the burgess and tully
c     method of linear interpolation on the collision strength
c     (Atron. Astrophys. 254 436-452 (1992) and Burgess et al.
c      J. Phys. B 30 33-57 (1997).
c
      implicit real*8 (a-h,o-z)
      common/const/d0,d1,d2,d3,d4,d30,d100,cnvevk,CONVEV,imskp
      dimension e(*),omg(*)
      ypar=d30
c     
c     if temperature is in degrees Kelvin convert to eV for calculation
c
      tmp=temp
      if(iunit.ne.0) tmp=temp/cnvevk
      imn=iei
      imx=ief-imskp
      sum=d0
      do 10 i=imn,imx,imskp
      h=e(i+imskp)-e(i)
      IF(ABS(H).LT.TOLH)GO TO 10
      del=h/tmp
      expdel=exp(-del)
      y=e(i)/tmp
      if(y.ge.ypar) go to 10
      sum=sum+(omg(i)-omg(i+imskp)*expdel-(omg(i)-omg(i+imskp))*
     1 ((d1-expdel)/del))*exp(-y)
ctest sum=sum+(omg(i)+omg(i+imskp)*expdel)*exp(-y)*del*0.5         !trapezoidal
   10 continue
      upex1=sum
      return
      end
c
c     ******************************************************************
c
      real*8 function upex2(emin,ethrsh,cc,idip,temp,iunit,UASYM,
     1 casym,ispn,ALF,LI,LF,IFAIL,IPTOMG)
c
c     calculate the contribution to the effective collision strength
c     from the portion of the collision strength with energies greater
c     than emin from the least-squares fitting coefficients
c
      implicit real*8 (a-h,o-z)
      common/const/d0,d1,d2,d3,d4,d30,d100,cnvevk,CONVEV,imskp
      dimension cc(6),ESTORE(1000),OMGSTR(1000)
      fnmsh=9.99d+02
      ypar=d30
c     
c     if temperature is in degrees Kelvin convert to eV for calculation
c
      tmp=temp
      if(iunit.ne.0) tmp=temp/cnvevk
      eji=emin-ethrsh
      x1=log(eji)
      u1=emin/ethrsh
      y1=eji/tmp
      if(y1.ge.ypar) then
        upex2=d0
        return
      end if
      ejf=ypar*tmp
      h=log(ejf/eji)/fnmsh
      ymx=ypar
      umx=(ejf+ethrsh)/ethrsh
C
      if(idip.eq.0) then
        IF(CASYM.EQ.D0)THEN
          omg1=cc(1)+cc(2)/u1+cc(3)/u1**2+cc(4)/u1**3
          omgmx=cc(1)+cc(2)/umx+cc(3)/umx**2+cc(4)/umx**3
        ELSE
          if(ispn.LE.0) then
            IF(CC(1).EQ.D0)THEN
              omg1=casym
              omgmx=casym
            ELSE
              XASYM=UASYM/(UASYM+D1)
              XR1=U1/(U1+D1)
              OMG1=CASYM*(D1-XR1)/(D1-XASYM)
     X            +CC(1)*(XASYM-XR1)/(XASYM-D1)
              XRMX=UMX/(UMX+D1)
              OMGMX=CASYM*(D1-XRMX)/(D1-XASYM)
     X             +CC(1)*(XASYM-XRMX)/(XASYM-D1)
            ENDIF
          else
            omg1=casym/(u1-D1)**ALF
            omgmx=casym/(umx-D1)**ALF
          end if
        ENDIF
      else
        IF(CASYM.EQ.D0)THEN
          omg1=cc(1)+cc(2)/u1+cc(3)/u1**2+cc(4)/u1**3+cc(5)*log(u1)+
     1         cc(6)*log(u1)/u1
          omgmx=cc(1)+cc(2)/umx+cc(3)/umx**2+cc(4)/umx**3+cc(5)*log(umx)
     1         +cc(6)*log(umx)/umx
        ELSE              !ASSUME CC(5).NE.0
          XASYM=LOG(UASYM-D1+EXP(D1))
          XASYM=D1-D1/XASYM
          XR1=LOG(U1-D1+EXP(D1))
          XR1=D1-D1/XR1
          OMG1=CASYM*(D1-XR1)/(D1-XASYM)+CC(5)*(XASYM-XR1)/(XASYM-D1)
          OMG1=OMG1*LOG(U1-D1+EXP(D1))
          XRMX=LOG(UMX-D1+EXP(D1))
          XRMX=D1-D1/XRMX
          OMGMX=CASYM*(D1-XRMX)/(D1-XASYM)+CC(5)*(XASYM-XRMX)/(XASYM-D1)
          OMGMX=OMGMX*LOG(UMX-D1+EXP(D1))
        ENDIF
      end if
C
      OMGSTR(1)=OMG1
      ESTORE(1)=EMIN/CONVEV
      sum=exp(-y1)*eji*omg1
      if(omgmx.gt.D0) sum=sum+exp(-ymx)*ejf*omgmx
      fac=d2
      x=x1
C
      do j=2,1000
      x=x+h
      ej=exp(x)
      y=ej/tmp
      if(fac.eq.d2) then
        fac=d4
      else 
        fac=d2
      end if
      u=(ej+ethrsh)/ethrsh
      if(idip.eq.0) then
        IF(CASYM.EQ.D0)THEN
          omg=cc(1)+cc(2)/u+cc(3)/u**2+cc(4)/u**3
        ELSE
          if(ispn.LE.0) then
            IF(CC(1).EQ.D0)THEN
              omg=casym
            ELSE
              XR=U/(U+D1)
              OMG=CASYM*(D1-XR)/(D1-XASYM)+CC(1)*(XASYM-XR)/(XASYM-D1)
            ENDIF
          else
            omg=casym/(u-D1)**ALF
          end if
        ENDIF
      else
        IF(CASYM.EQ.D0)THEN
          omg=cc(1)+cc(2)/u+cc(3)/u**2+cc(4)/u**3+cc(5)*log(u)+
     1        cc(6)*log(u)/u
        ELSE       !ASSUME CC(5).NE.0
          XR=LOG(U-D1+EXP(D1))
          XR=D1-D1/XR
          OMG=CASYM*(D1-XR)/(D1-XASYM)+CC(5)*(XASYM-XR)/(XASYM-D1)
          OMG=OMG*LOG(U-D1+EXP(D1))
        ENDIF
      end if
C
      ESTORE(J)=U*ETHRSH/CONVEV
      OMGSTR(J)=OMG
c
c     test for problems with fit
c
      if(omg.lt.D0) then
        IFAIL=IFAIL+1
        IF(IFAIL.EQ.1)THEN
          e=u*ethrsh
          write(6,100) e,tmp,li,lf
          T=ETHRSH/CONVEV
          WRITE(10,110)LI,LF,T
          DO K=1,J
            WRITE(10,120)D0,ESTORE(K),ESTORE(K)-T,OMGSTR(K)
          ENDDO
        ENDIF
        go to 20
      end if
      sum=sum+fac*exp(-y)*ej*omg
C
      ENDDO
   20 upex2=sum*h/(d3*tmp)
C
C  DETAILED XTRAP OUTPUT, FORMATTED FOR BURG PROG.
C
      IF(IPTOMG.LT.0)THEN
        T=ETHRSH/CONVEV
        WRITE(33,110)LI,LF,T
        DO K=1,1000,10
          WRITE(33,120)D0,ESTORE(K),ESTORE(K)-T,OMGSTR(K)
        ENDDO
      ENDIF
C
      return
C
  100 format(/5x,
     1 '************************* WARNING! *************************'/
     2 5x,'omega went negative at an energy of ',1pe12.6,' eV for a'/
     3 5x,'temperature of ',e10.4,' eV in the rate calculation.  the'/
     4 5x,'rate calculation was terminated at this energy for'
     5 ,' transition:',2I5/5x,
     6 '******************** NO MORE WARNINGS! *********************'/)
  110  FORMAT('#',I3,' -',I4,' TRANSITION, ENERGY=',F12.6)
  120  FORMAT(1X,4(1PE14.6))
      end
c
c     ******************************************************************
c
      subroutine choldc(a,n,np,p,ifail)
      integer n,np
      real*8 a(np,np),p(n),sum
      integer i,j,k,ifail
c
c     Given a positive definite symmetric matrix a(1:n,1:n), with
c     physical dimension np, this routine constructs its Cholesky
c     decomposition A=L.L^T.  On input, only the upper triangle of A
c     needs to be given; it is not modified.  the Cholesky factor L
c     is returned in the lower triangle of A, except for its diagonal
c     elements which are returned in p(1:n)
c
      ifail=0   
      do 13 i=1,n
      do 12 j=i,n
      sum=a(i,j)
      do 11 k=i-1,1,-1
      sum=sum-a(i,k)*a(j,k)
   11 continue
      if(i.eq.j) then
        if(sum.le.0) then
          write(6,*) 'choldec failed:sum= ',sum
          write(6,*) 'i=',i,' j=',j
          ifail=1
          return
        endif
        p(i)=sqrt(sum)
      else
        a(j,i)=sum/p(i)
      endif
   12 continue
   13 continue
      return
      end
c
c     **************************************************************** 
c
      subroutine cholsl(a,n,np,p,b,x)
      integer n,np
      real*8 a(np,np),b(n),p(n),x(n),sum
      integer i,k
c 
c     Solves the set of n linear equations Ax=b, where A is a positive
c     definite symmetric matrix with physical dimension np.  A and p
c     are input as the output of the routine choldc.  Only the lower
c     triangle of A is accessed.  b(1:n) is input as the right-hand
c     side vector.  The solution vector is returned in x(1:n).  A, n,
c     np, and p are not modified and can be left in place for
c     successive calls with different right-hand sides b.  b is not
c     modified unless you identify b and x in the calling sequence,
c     which is allowed.
c
      do 12 i=1,n
      sum=b(i)
      do 11 k=i-1,1,-1
      sum=sum-a(i,k)*x(k)
   11 continue
      x(i)=sum/p(i)
   12 continue
      do 14 i=n,1,-1
      sum=x(i)
      do 13 k=i+1,n
      sum=sum-a(k,i)*x(k)
   13 continue
      x(i)=sum/p(i)
   14 continue
      return
      end
C
C************************************************
C
      INTEGER FUNCTION NOPEN(E,ENAT,NAST,IONE,I0)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ENAT(NAST)
C
      DO I=I0,NAST
        IF(E.LT.ENAT(I))GO TO 1
      ENDDO
      I=NAST+1
   1  I0=I-1
      NOPEN=(I0*(I0-2*IONE+1))/2
C
      RETURN
      END
