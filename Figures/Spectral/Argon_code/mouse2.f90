                           program RadT
!
  implicit real*8(A-H,O-Z)
!
! Crucial information for code to run
  NAMELIST /COLLINP/ Nodens,notemp,sDens,fDens,IC,grad, &
           worder,petol,mantmp,wrngup,wrngdn,pops,corapp,lte
! Optional Line ratio stuff
  NAMELIST/LIRA/uplev,uplevl,lolev,lolevl,uppec,lopec,ppec,line, &
                & lvart,lvard,diaglvl
!
  integer :: i, ii, j, jj, jjj, k, kk, kkk, order, Ne, prnt, pops, corapp, &
           & IC, grad, icount, worder, dum, dumm, einf, jcount, &
           & levels, notemp, Z, tripi, tripj, Nodens, sDens, fDens, &
           & uplev, lolev, kcount, petol, lopec, uppec, dumrr, &
           & ppec, mantmp, wrngup, wrngdn, line, uplevl, lolevl, lte, &
           & collprint, cup, clo, lvart, lvard, lvari, lvarj, tk, Nee, &
           & diaglvl
  real*8 :: hold, IH, dumr, sDensn, start, finish, start1, &
          & finish1, start2, finish2, start3, finish3, petolr, asu(2), &
          & const, boltz, factor, plank, sos, Me, zero, newDen, corappv, &
          & ltev, DE1, DE2, totAup, totAlo, S1, S2, S3
  real*8, allocatable :: aij(:,:), rate(:,:,:), temp(:), holde(:), ED(:), &
                       & XJ(:), weight(:), en(:), err(:), lratio(:,:), &
                       & COLL(:,:,:,:), LHS(:,:,:), pop(:,:,:), wlength(:,:), &
                       & wlenor(:), pec(:,:,:,:), pecor(:,:,:)
  integer, allocatable :: ind(:), mult(:), bEll(:), icoarr(:)
  character*13, allocatable :: config(:)
  character*3 :: charin,chr
!
  parameter(MXEN=700)         ! NUMBER OF TARGET STATES
  parameter(EINF=0)           ! = 0 for no infinite energies
                              ! = 1 otherwise
!
  open(15, file='mouse',status='unknown', form='formatted')
  open(16, file='mout',status='unknown', form='formatted')
  open(404, file='const_dens',status='unknown', form='formatted')
  open(405, file='const_temp',status='unknown', form='formatted')
  open(406, file='populations',status='unknown', form='formatted')
!
  call cpu_time(start)
  call cpu_time(start1)
!
  worder = 0
  corapp = 0
    pops = 0
   petol = 0
   lopec = 1
   uppec = 2
    ppec = 0
  mantmp = 1
  wrngup = 0
  wrngdn = 0
    line = 0
    grad = 0
     clo = 0
     cup = 0
!
  read(15,collinp)
  read(15,lira)
!
  if ((Nodens.lt.2).or.(Nodens.gt.999)) then
      print *, 'WARNING, Nodens in RATinp MUST BE 2 < Nodens < 999'
  endif
!
! Set up the mesh of Densities
!
  allocate(ED(Nodens))
!
  newDen = real(sDens)
!
  if (Nodens.eq.2) then
     ED(1) = 10.0d0**(real(sDens))
     ED(2) = 10.0d0**(real(fDens))
  else
        if ((real(fDens) - real(sDens).lt.1.1)) then
           Me = ((real(fDens) - real(sDens)))/(real(Nodens) - 1.0)
        else
           Me = ((real(fDens) - real(sDens)))/(real(Nodens) - 1.0)
        endif
     ED(1) = 10.0d0**(real(sDens))
        do i=2,Nodens
           newDen = newDen + Me
           ED(i) = 10.0d0**(real(newDen))
           print *, 'Density: ', i, ED(i)
        enddo
  endif
!
  close(15)
!
! Input constants
!
  const = 2.1716d-8     ! 2pia_0 const. for rates (q_ij)
     IH = 13.605698     ! Ionization of Hydrogen
  boltz = 8.61733247d-5 ! Boltzmann constant in eV
  plank = 4.1357d-15    ! Plancks constant to convert to Ang
    sos = 3.0d8         ! Speed of sound
!
!   open(10, file ='input1', status = 'old')
  open(10, file ='input2', status = 'old')
  read(10,1000) Z
!
  levels = 0
!
! Calculate the number of levels in adf04 without having to specify
!
  do i=1,MXEN
     read(10,*) dum
        if(dum.eq.(-1))then
           go to 200
        else
           continue
        endif
     levels = levels + 1
  enddo
!
! Reset reading of input back to the first level
!
  200 continue
     do i=1,(levels+1)
        backspace(10)
     enddo
!
  allocate(ind(levels))
    allocate(en(levels))
      allocate(mult(levels))
        allocate(bEll(levels))
          allocate(XJ(levels))
            allocate(weight(levels))
            allocate(config(levels))
          allocate(aij(levels,levels))
        allocate(rate(levels,levels,notemp))
      allocate(temp(notemp))
    allocate(holde(notemp+einf))
  allocate(err(levels))
!
  write(16,7000)
  write(16,7001) Z, notemp, Nodens, sDens, fDens, Me
  write(16,7002)
  write(16,7003) levels
!
! Calculate rates from Upsilons in ADF04
!
  do i=1,levels
!
! WARNING:
! This may need adjusted depending on the output of ADF04!
!
! read(10,1001) ind(i), config(i), mult(i), bEll(i), XJ(i), en(i)
     read(10,2001) ind(i), config(i), mult(i), bEll(i), XJ(i), en(i)
!
        if (i.eq.levels) then
           if (ind(i).eq.levels) then
              print *, 'Read of levels correct, total in ADF04: ', levels
           else
              print *, 'Read of crucial input data INCORRECT.'
           endif
        endif
!
! cm-1 -> eV conversion
     en(i) = en(i)*0.000123984
     write(16,7004), ind(i), config(i), en(i)
! Statistical weightings..
! LSpi / Jpi :
        if (IC.eq.0) then
           weight(i) = (mult(i))*(2*bEll(i)+1)
        else
           weight(i) = (2*XJ(i)+1)
        endif

        if ((i.eq.1).or.(i.eq.levels)) then
          print *, 'Check the line is read correctly:'
          print *, 'Index: ', ind(i), 'Config: ', config(i), &
          & 'Stat weight: ', weight(i), 'Energy: ', (en(i)/0.000123984)
        endif
  enddo
!
  read(10,*) dumm
!
! Just give exact number of temperatures during input and double check
!
  read(10,*) dumr, dumrr, (temp(i), i=1,(notemp))
  print *, 'Check these temperatures are correct:', (temp(i), i=1,(notemp))
!
! I use tripi/j, hold and holde for holding variables
! to compute the ex/de-excitation rates and A-values.
!
  if (grad.eq.1) then
     open(5,file='adf04rad',status='unknown', form='formatted')
  endif

   aij = 0.0
  rate = 0.0
!
  order = ((levels**2-levels)/2)
!
  do i=1,order
     read(10,2003) tripj, tripi, hold, (holde(j),  j=1,notemp+einf)
!     print *, tripj, tripi, hold, (holde(j),  j=1,notemp+einf)
!
! If grad = 1 then read A-values from a GRASPRAD ADF04RAD file
!
        if (grad.eq.1) then
           read(5,*) dumrr, dumm, hold
if ((dumrr.eq.6).and.(dumm.eq.2)) then
  print *, hold
  print *, tripj, tripi
  print *, dumrr, dumm
endif


              if (hold.lt.1e-20) then   ! If A-value is too low.
                 hold = 0.0
              endif
        endif
!
    aij(tripj,tripi) = hold
    aij(tripi,tripj) = hold
!
! Excitation rate, depends on the temperature but not density.
!
     do k=1,notemp
        rate(tripi,tripj,k) = (const/weight(tripi))* &
           & sqrt((IH)/(boltz*temp(k)))* &
           & exp((en(tripi)-en(tripj))/(boltz*temp(k)))*holde(k)
!
! De-excitation rate.
!
        rate(tripj,tripi,k) = (weight(tripi)/weight(tripj))* &
           & exp(-(en(tripi)-en(tripj))/(boltz*temp(k)))* &
           & rate(tripi,tripj,k)
     enddo
!
  enddo

  DE1 = en(uplev)-en(uplevl)
  DE2 = en(lolev)-en(lolevl)
!
! Important formats for ADF04 file
! These may need changed depending how
!
  1000 format(13X,I2)
  1001 format(2X,I3,3X,A13,4X,I1,1X,I1,1X,F4.1,3X,F9.1)
  1003 format(2X,I2,2X,I2,25F8.5)
!
  2001 format(2X,I3,3X,A13,4X,I1,1X,I1,1X,F4.1,4X,F12.4)
  2003 format(1X,I3,1X,I3,25F8.5)
  2010 format(1X,I3,1X,I3,1X,16F8.6)
!
     if (allocated(mult)) deallocate(mult)
     if (allocated(bEll)) deallocate(bEll)
     if (allocated(XJ)) deallocate(XJ)
!     if (allocated(weight)) deallocate(weight)
!
! Reading routine complete..
! I've got all A-values stored in aij(i,j)
! and all rates in rate(i,j,k)
! #############################################################
! #############################################################
!
! Next step is to fill the Collision matrix COLL..
! #############################################################
! #############################################################
!
  edens = sDensn
  allocate(COLL(levels,levels,notemp,Nodens))
  COLL = 0.0d0
!
print *, 'Begin filling collision matrix..'
!
  do Ne=1,Nodens
     edens = edens + Me
     print *, Ne, edens
        do k=mantmp,notemp
           do i=1,levels
              do j=1,levels
                 if (i.eq.j) then
                    do kk=1,j
                       COLL(i,j,k,Ne) = COLL(i,j,k,Ne) - aij(j,kk)
                    enddo
                    do kkk=1,levels
                       if (kkk.eq.j) then
                          go to 99
                       else
                          COLL(i,j,k,Ne) = COLL(i,j,k,Ne) - &
                          & ((ED(Ne))*(rate(j,kkk,k)))
                       endif
                          99 continue
                    enddo
                 else if (i.gt.j) then
                    COLL(i,j,k,Ne) = COLL(i,j,k,Ne) + &
                    & ((ED(Ne))*(rate(j,i,k)))
                 else if (i.lt.j) then
                    COLL(i,j,k,Ne) = COLL(i,j,k,Ne) + aij(i,j) + &
                    & ((ED(Ne))*(rate(j,i,k)))
                 endif
              enddo
           enddo
        enddo
  enddo
!

print *, 'COLL HERE: ', COLL(1,1,7,30)

collprint = 1

if (collprint.eq.1) then
     open(1000,file='level_ind',status='unknown', form='formatted')

icount = 0
jcount = 0

     do ii=1,3
!     do Ne=1,1
       icount = icount + 1
       jcount = jcount + 1
!       print *, 'ED(29) = ', ED(29)
  do i=1,levels
            if (ii.eq.1) then
            write(1000,*) i
            endif

  write(chr,'(I3)') icount
  chr = adjustl(chr)

     if (ii.eq.1) then
       k = 4
     elseif (ii.eq.2) then
       k = 6
     elseif (ii.eq.3) then
       k = 7
     endif

open(1001+jcount,file='COLL'//trim(chr),status='unknown', form='formatted')
    write(1001+jcount,*) (abs(COLL(i,j,k,30)), j=1,levels)
    if (ii.eq.3) then
      print *, 'COLL HERE: ', COLL(1,1,7,30)
    endif
!enddo
enddo

enddo
endif


  call cpu_time(finish1)
  print '("Basic Read = ",f6.3," seconds.")',finish1-start1
  call cpu_time(start2)
!
! Gaussian Elimination
! #############################################################
! #############################################################
!
  allocate(LHS((levels),notemp,Nodens))
!
  LHS = 0.0d0
!
! I know at this point I put in LHS(1,k,kk) but obviously with the
! dN/dt .ne. 0 then I can't write the LHS for the first one as -C11N1
! So I just won't use it later in the Gauss elimination.
!
  do k=mantmp,notemp
     do kk=1,Nodens
        do kkk=1,(levels)
           LHS(kkk,k,kk) = -COLL((kkk),1,k,kk)
        enddo
     enddo
  enddo
!
  print *, 'Begin Guassian elimination..'
  print *, 'Density      Temperature'
  print *, '________________________'
!
  do Ne=1,Nodens
     print *,'  ', Ne
        do kk=mantmp,notemp
           print *, '           ', kk
              do i=2,(levels-1)
                 do j=i+1,levels
                    factor=COLL(j,i,kk,Ne)/COLL(i,i,kk,Ne)
                    COLL(j,i,kk,Ne)=0.0d0
                    LHS(j,kk,Ne) = LHS(j,kk,Ne) - factor*LHS(i,kk,Ne)
                       do k=i+1,levels
                          COLL(j,k,kk,Ne) = COLL(j,k,kk,Ne) - &
                          & factor*COLL(i,k,kk,Ne)
                       enddo
                 enddo
              enddo
        enddo
  enddo
!
! Check matrix is diagonal...
! Add up all the lower triangular terms and see if they
! are < 0 or some tolerance.. (TOL)
!
! Clearly the first (i,j) has a value since it is diagonal.
! But we also reduce the matrix to (n-1) x (n-1), so the second i
! is the new diagonal, so the first 0 will occur from i = 3. We
! run this down to (levels).
!
  zero = 0.0
  icount = 0
!
  do k=mantmp,notemp
     do l=1,Nodens
        do i=3,levels
           do j=2,(i-1)
              zero = zero + COLL(i,j,k,l)
           enddo
        enddo
        !
        if (zero.gt.0.0000000001d0) then
           print *, 'WARNING, GAUSSIAN ELIMINATION MAY HAVE FAILED'
           icount = icount + 1
        endif
        !
     enddo
  enddo
!
  if (icount.eq.0) then
     print *, 'Guassian elimination successful'
  endif
!
  call cpu_time(finish2)
  print '("Gaussian elimination = ",f6.3," seconds.")',finish2-start2
!
! If matrix is small can visualise the diagonalised output by uncommenting
! this below.
!
!  do k=2,levels
!     print *, '|', (COLL(k,kk,4,1),kk=2,levels), '||', LHS(k,4,1), '|'
!  enddo
!
  write(16,7005) levels
  jcount = 0
!
  do k=2,levels
     if (mod(k,8).eq.0) then
        write(16,7006) (LHS(kk+(jcount*8),mantmp,1), kk=1,8)
        jcount = jcount + 1
     else if (k.eq.levels) then
        write(16,7006) (LHS(levels-mod(levels,8)+kk,mantmp,1), &
        & kk=1,(mod(levels,8)))
     else
        continue
     endif
  enddo
!
! Determine populations relative to the ground state, this is
! taking the diagonal COLL matrix and using back substitution
! to work out each of the N_i values
!
  icount = 0
  allocate(pop(levels,notemp,Nodens))
  pop = 0.0
!
  do Ne=1,Nodens
     do kk=mantmp,notemp
        do i=levels,2,-1
           pop(i,kk,Ne) = (LHS(i,kk,Ne))/(COLL(i,i,kk,Ne))
              if (i.ne.levels) then
                 do j=levels,(i+1),-1
                    pop(i,kk,Ne) = pop(i,kk,Ne) - &
                    & ((COLL(i,j,kk,Ne)*pop(j,kk,Ne))/(COLL(i,i,kk,Ne)))
                 enddo
              endif
        enddo
     enddo
  enddo
!
! Before line ratios, plot the population as a function
! of density for a given temperature, this is useful
! to see what lines may be good as a density diagnostic
!
  if (pops.eq.1) then
     do i=1,levels
        write(charin,'(I3)') i
        charin = adjustl(charin)
        open (410 + i, file='population_'//trim(charin), &
        & status='unknown', form='formatted')
            write(410 + i,7011) (temp(j), j=mantmp,notemp)
           do Ne=1,Nodens
              write(410 + i,7012) ED(Ne), (pop(i,k,Ne), k=mantmp,notemp)
!               write(410 + i,*) (pop(i,k,Ne), k=mantmp,notemp)
           enddo
              close(410+i)
     enddo
  endif
!
!
! Loop over densities and temperatures now in Population arrays
! to obtain a specific line ratio defined by lolev and uplev
!
!
   allocate(lratio(Nodens,notemp))
!
!
  if ((lvart.eq.1).or.(lvard.eq.1)) then
    continue
    if (mod(Nodens,2).eq.1) then
       Nee = (Nodens+1)/2
     else
       Nee = (Nodens/2)
     endif
     if (mod((notemp-mantmp+1),2).eq.1) then
        tk = (notemp+1)/2
      else
        tk = (notemp/2)
      endif

  else
    go to 2211
  endif
  open(2223,file='lrat_var',status='unknown',form='formatted')
  if (lvard.eq.1) then
  write(2223,*) 'Constant density..'
  else
    write(2223,*) 'Constant temperature..'
  endif
!
!

do ii=1,diaglvl
do jj=1,diaglvl
  if (aij(ii,jj).lt.0.000000001d0) then
    go to 2210
  elseif (ii.lt.jj) then
    go to 2210
  endif

write(2223,*) ' ------------'
write(2223,*) '|', ii, jj
write(2223,*) ' ------------'

do i=1,diaglvl
  do j=1,diaglvl
               if (aij(i,j).lt.0.000000001d0) then
                 go to 2222
               elseif (i.lt.j) then
                 go to 2222
               endif

    if (lvard.eq.1) then
    do kk=mantmp,notemp
       lratio(Nee,kk) = (pop(ii,kk,Nee)*aij(ii,jj))/ &
       & (pop(i,kk,Nee)*aij(i,j))*(en(ii)-en(jj))/(en(i)-en(j))
!               print *, lratio(Ne,kk)*(DE1/DE2)
    enddo
  elseif (lvart.eq.1) then
    do Ne=1,Nodens
       lratio(Ne,tk) = (pop(ii,tk,Ne)*aij(ii,jj))/ &
       & (pop(i,tk,Ne)*aij(i,j))*(en(ii)-en(jj))/(en(i)-en(j))
     enddo
   endif

if (lvard.eq.1) then
if ((abs(lratio(Nee,mantmp)).lt.(0.1)).or. &
& (abs(lratio(Nee,mantmp)).gt.(90.0)).or. &
& (abs(lratio(Nee,tk)).lt.(0.1)).or.(abs(lratio(Nee,tk)).gt.(1.0))) then
      go to 2222
  else
     print *, i, j, ii, jj, lratio(Nee,1)
write(2223,*) i, j
  endif
elseif (lvart.eq.1) then
  if ((abs(lratio(Nee,tk)-lratio(Nee,tk)).gt.1.0d1).or. &
  & (abs(lratio(Nee,tk)).lt.(0.1)).or.(abs(lratio(Nee,tk)).gt.(1.0d1))) then
        go to 2222
    else
  write(2223,*) i, j
    endif
endif

          2222 continue
  enddo
enddo

 2210 continue

enddo
enddo

2211 continue
!
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  edens = sDensn
  write(16,7007) uplev, lolev
!
!
! Keeping Density constant:
!
  prnt = 0
!
  do Ne=1,Nodens
     if (prnt.eq.1) then
        write(16,7008) Ne, ED(Ne)
           do kk=mantmp,notemp
              lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
              & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
           enddo
     elseif (prnt.eq.0) then
        if (Ne.eq.1) then
           write(16,*) '-> Increasing density'
           write(16,*)  '      Temperature          ', (ED(jj), jj=1,Nodens)
        else
           continue
        endif
!
           do kk=mantmp,notemp
!             print *, 'get here', aij(uplev,uplevl)
!             print *, 'get here', lolev, lolevl, aij(lolev,lolevl)
              lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
              & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
           enddo
!
     endif
        if (Ne.eq.Nodens) then
           do kk=mantmp,notemp
              write(16,*) temp(kk), (lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
              write(404,*) temp(kk), (lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
           enddo
        endif
  enddo
!
! Keeping Temperature constant:
!
  do i=1,2
     write(16,*)
  enddo
!
  do kk=mantmp,notemp
     if (prnt.eq.1) then
        write(16,7010) kk, temp(kk)
           do Ne=1,Nodens
              lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
              & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
              write(16,*) ED(Ne), lratio(Ne,kk)*(DE1/DE2)
           enddo
     elseif (prnt.eq.0) then
        do Ne=1,Nodens
           lratio(Ne,kk) = (pop(uplev,kk,Ne)*aij(uplev,uplevl))/ &
           & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
        enddo
           if (kk.eq.mantmp) then
              write(16,*) '-> Increasing temperature'
              write(16,*)  '      Density          ', &
              & (temp(jj), jj=mantmp,notemp)
           else
              continue
           endif
           !
           if (kk.eq.notemp) then
              do Ne=1,Nodens
          write(16,*) ED(Ne), (lratio(Ne,ii)*(DE1/DE2), ii=mantmp,notemp)
          write(405,*) ED(Ne), (lratio(Ne,ii)*(DE1/DE2), ii=mantmp,notemp)
              enddo
           endif
          !
     endif
  enddo

!
! Corronal Approx and LTE:
! DENSITY LIMITS, LOOP OVER T
!

if (lte.eq.1) then
write(16,*) 'LTE values:'

do k=mantmp,notemp
ltev = (weight(uplev)/weight(lolev))* &
       & exp(-(en(uplev)-en(lolev))/(boltz*temp(k)))
       write(16,*) k, 10.0e8, &
                 & ltev*((aij(uplev,uplevl))/(aij(lolev,lolevl)))* &
                 & (DE1/DE2)
enddo
endif

if (corapp.eq.1) then
!
totAup = 0.0
totAlo = 0.0
!
  write(16,*) 'Coronal Approx values:'

  do i=1,(uplev-1)
    totAup = totAup + aij(uplev,i)
  enddo
  do i=1,(lolev-1)
    totAlo = totAlo +  aij(lolev,i)
  enddo

  print *, 'tots A lo:', totAup
  print *, 'tots A:', totAlo

do k=mantmp,notemp

corappv = (rate(1,uplev,k)/rate(1,lolev,k))*(totAlo/totAup)* &
         & (aij(uplev,uplevl)/aij(lolev,lolevl))
         write(16,*) k, 10.0e2, &
         & (corappv*(DE1/DE2))
enddo


endif


!


!
! Have to generate all the pecs (= pop * A-value)
!
! lambda = c*h / (Ej - Ei) where lambda will be in m
! 1m = 10^10 Ang
!
! Need to define new arrays to save for printing
!
  call cpu_time(start3)
!
  petolr = real(10.0**(real(petol)))
!
  allocate(wlength(levels,levels-1))
    allocate(pec(levels,levels-1,notemp,Nodens))
      allocate(wlenor(order))
    allocate(icoarr(order))
  allocate(pecor(order,notemp,Nodens))
!
  open(18, file='lineID',status='unknown', form='formatted')
  write(18,*) 'Upper state Lower State  Upper config &
  &    Lower config         Upper energy              Lower energy  &
  &              wavelength'
!
  open(19, file='PEC_plot',status='unknown', form='formatted')
!
  do kk=mantmp,notemp
     if (ppec.eq.0) then
        write(19,*) 'Temperature: ', kk
        write(19,*) '--------------------'
        write(19,*)  'Wavelength                 Density     Value  &
                    &         PEC'
     endif
        do Ne=1,Nodens
           do j=2,levels ! This double loop covers 'order' repititions.
              do i=1,j-1
                 if ((kk.eq.mantmp).and.(Ne.eq.1)) then
!
! Only need to compute this once for wavelengths
!
                    wlength(j,i) = (sos*plank)/(abs(en(j)-en(i)))
                    wlength(j,i) = wlength(j,i)*1d10
                    write(18,*) j, i, config(j), ' -> ', config(i), &
                                & en(j), en(i), wlength(j,i)
!
! Wavelength between each state j and i.
!
                 endif
                    pec(j,i,kk,Ne) = pop(j,kk,Ne)*aij(j,i)
!
! a 5 nested if statement. Checks the following
!               1) wavelength ordering is off
!               2) ppec is for all PEC's.
!               3) Only PEC's above a certain tolerance
!               4) Only wavlengths in a certain range
                 !
                 if (worder.eq.0) then
                    if (ppec.eq.0) then
                       if ((pec(j,i,kk,Ne)).gt.petolr) then
                          if ((wrngup.gt.0).or.(wrgndn.gt.0)) then
                             if ((wlength(j,i).gt.wrngdn).and. &
                                & (wlength(j,i).lt.wrngup)) then
                                write(19,*) wlength(j,i), Ne, ED(Ne), &
                                            & pec(j,i,kk,Ne)
                             endif
                          endif
                       endif
                    endif
                 endif
                 !
              enddo
           enddo
        enddo
  enddo
!
  edens = sDensn
!
! ppec = 1 singles out a unique j -> i transition and loops over
! both temperature and density again.
!
  if (ppec.eq.1) then
     write(19,*) 'ppec = ', ppec
     write(19,*) 'Looking at PEC for transition', uppec, '->', lopec
        do kk=mantmp,notemp
           write(19,*) 'Temp ', temp(kk)
              do Ne=1,Nodens
                 write(19,*) ED(Ne), pec(uppec,lopec,kk,Ne)
              enddo
        enddo
           write(19,*) 'Now loop with a constant density'
        do Ne=1,Nodens
           write(19,*) 'Dens ', ED(Ne)
              do kk=mantmp,notemp
                 write(19,*) temp(kk), pec(uppec,lopec,kk,Ne)
              enddo
        enddo
  endif
!
! This is all the wavelengths computed, but not in energy order
! Need to sort by wavelength..
!
! This takes a long time, so only use if absolutely necessary!
!
  if (worder.eq.1) then
     print *, 'Wavelength ordering requested, with worder = ', worder
     continue
  else
     print *, 'No wavelength ordering requested, with worder = ', worder
     go to 801
  endif
!
  icount = 1
  kcount = 1
  jcount = 0
  wlenor = 0.0
!
  do kk=mantmp,notemp
     print *, 'Starting Temperature ', kk
        do Ne=1,Nodens
           print *, 'Starting Density ', Ne
              do jj=2,levels
                 do ii=1,jj-1
                    do j=2,levels
                       do i=1,j-1
                          if ((jj.eq.j).and.(ii.eq.i)) then
                             go to 800
                          else
                             if (wlength(jj,ii).gt.wlength(j,i)) then
                                icount = icount + 1
                             endif
                          endif
                             800 continue
                       enddo
                    enddo
!
                       icoarr(kcount) = icount
!
! Only the pec are dependent on temperature and density.
! So only adjust the wavelengths for one of the temperatures
! and densities.
!
                    if (kcount.gt.1) then
                       do jjj=1,(kcount-1)
                          if (icount.eq.icoarr(jjj)) then
                             jcount = jcount + 1
                          endif
                       enddo
                       !
                       pecor(icount+jcount,kk,Ne) = pec(jj,ii,kk,Ne)
                       !
                       if ((kk.eq.1).and.(Ne.eq.1)) then
                          wlenor(icount+jcount) = wlength(jj,ii)
                       endif
                    else
                       pecor(icount,kk,Ne) = pec(jj,ii,kk,Ne)
                    !
                          if ((kk.eq.1).and.(Ne.eq.1)) then
                             wlenor(icount) = wlength(jj,ii)
                          endif
                    endif
!
! Reinitiliaze counters
!
  icount = 1
  kcount = kcount + 1
  jcount = 0
!
                 enddo
              enddo
           print *, '                 ', Ne
        enddo
     print *, '                     ', kk
  enddo
!
! Now begin to print out line ratios into output unit
!

  do kk=mantmp,notemp
     write(19,*) 'Temperature: ', kk
        do Ne=1,Nodens
           write(19,*) 'Density: ', Ne
              do k=1,order
                 if (pecor(k,kk,Ne).gt.petolr) then
                    write(19,*) k, wlenor(k), pecor(k,kk,Ne)
                 endif
              enddo
        enddo
  enddo
call cpu_time(finish3)
print '("Rearranging Wavelength Order = ",f6.2," minutes.")', &
       & (finish3-start3)/60.0
!
! This marker skips the whole wavelength ordering. Probably
! recommend as this is probably the lengthiest portion to compute.
!
  801 continue
!
  call cpu_time(finish)
  print '("Total Time = ",f6.2," minutes.")',(finish-start)/60.0
!
  if (line.eq.1) then
     !
     open(950, file ='lratio' , status = 'unknown',form='formatted')
     open(951, file ='dens' , status = 'unknown',form='formatted')
     open(952, file ='temp' , status = 'unknown',form='formatted')
     !
     write(950,*) 'Electron Density | Electron Temperature | Line Ratio'
     write(950,*) '----------------------------------------------------'
!
! These vectors are needed for 3D plotting in Matlab for example
!
        do Ne=1,Nodens
           write(951,*) ED(Ne)
        enddo
        do kk=mantmp,notemp
           write(952,*) temp(kk)
        enddo
        do Ne=1,Nodens
           write(950,*) (lratio(Ne,kk),kk=mantmp,notemp)
        enddo
  endif
 print *, 'last thing'
  write(16,*) 'Print off some SXB ratio stuff for comparison!'
  S1 = 1.1d-11
  S2 = 1.3d-8
  S3 = 3.3d-8

  do Ne=1,Nodens
  write(16,*) ED(Ne), (S1*ED(Ne))/(aij(3,2)*pop(3,4,Ne)), &
& (S2*ED(Ne))/(aij(3,2)*pop(3,6,Ne)), (S3*ED(Ne))/(aij(3,2)*pop(3,7,Ne))
  enddo

!
! Format specifiers for output file
!
  7000 format(' --- RAdiative Transfer Output ---',/,/,'  - Basic Read of &
            & Input - ', /, &
            & 8X, '_____________________________________________')
  7001 format(8X, '| Atomic Number = ', 22X, I3, ' |', / ,&
            & 8X, '| Number of Temperatures = ', 13X, I3, ' |', / ,&
            & 8X, '| Number of Densities = ', 15X, I4, ' |', / ,&
            & 8X, '| Starting at 10^', I1, 'K and finishing at 10^', I2, 'K', &
            & 1X, '|', /, 8X,'| with a mesh increment of ', F16.4, ' |', &
            & /, 8X, '_____________________________________________', /)
  7002 format('RESERVED FOR SPECIFIC LINES NAMELIST LINE.')
  7003 format('  - Basic Level Information -  ', /, /,  &
            & 'Total number of levels in ADF04 file: ', I3, /, /, &
            & '_____________________________________________', /, &
            & '| IND.|   CONFIGURATION     |  ENERGY (eV)  |', /, &
            & '_____________________________________________')
  7004 format('| ',I3,' | '3X,A13,3X,' | ',F12.4,'  |')
  7005 format('_____________________________________________', /, /, &
            & 'The new left hand side values from N_i = 1 - ', I3, &
            & ' for a particular temperature and density (normally 1 and 1)',/)
  7006 format(8D12.4)
  7007 format(/, 'Looking at the specific line ratio of', I3, '->1/', I3, &
            & '->1 over temperatures and densities.')
  7008 format(/, '-------------- Density ...', 4X, I3, '(',E10.4,')')
  7009 format(2X,E10.4,4X,F20.16)
  7010 format(/, '-------------- Temperature ...', 4X, I3, '(',E10.4,')')
!
  7011 format(3X,'Density',1X,40ES11.4)
  7012 format(999ES11.4)
!
end program
