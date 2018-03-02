                           program RadT
!
  implicit real*8(A-H,O-Z)
!
! Crucial information for code to run
  NAMELIST/COLLINP/Nodens,notemp,sDens,fDens,IC,popt,grad, &
          & worder,petol,mantmp,wrngup,wrngdn,pops,corapp,lte
! Optional Line ratio stuff
  NAMELIST/LIRA/uplev,uplevl,lolev,lolevl,uppec,lopec,ppec, &
          & line
!
  integer :: i,ii,j,jj,jjj,k,kk,kkk,order,Ne,prnt,pops, &
          &  corapp,popt,IC,grad,icount,worder,dum,dumm, &
          &  einf,jcount,levels,notemp,Z,tripi,tripj, &
          &  Nodens,sDens,fDens,uplev,lolev,kcount,petol, &
          &  lopec,uppec,dumrr,ppec,mantmp,wrngup,wrngdn, &
          &  line,uplevl,lolevl,lte,collprint,cup,clo,tk,Nee

  real*8 :: hold,IH,dumr,start,finish,start1,finish1, &
          &  start2,finish2,start3,finish3,petolr, &
          &  const,boltz,factor,plank,sos,Me,zero, &
          &  newDen,corappv,ltev,DE1,DE2,totAup,totAlo

  real*8, allocatable :: aij(:,:),rate(:,:,:),temp(:), &
          &  holde(:),ED(:),XJ(:),weight(:),en(:),err(:), &
          &  lratio(:,:),COLL(:,:,:,:),LHS(:,:,:), &
          &  pop(:,:,:),wlength(:,:),wlenor(:), &
          &  pec(:,:,:,:),pecor(:,:,:)

  integer, allocatable :: ind(:),mult(:),bEll(:),icoarr(:)
  character*13, allocatable :: config(:)
  character*3 :: charin,chr
!
  parameter(MXEN=700)         ! NUMBER OF TARGET STATES
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
    popt = 0
     clo = 0
     cup = 0
    einf = 0
!
  read(15,collinp)
  read(15,lira)
!
  if ((Nodens.lt.2).or.(Nodens.gt.999)) then
      print *, 'WARNING, Nodens in RATinp MUST BE 2 &
               & < Nodens < 999'
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
           Me = ((real(fDens) - real(sDens)))/ &
                & (real(Nodens) - 1.0)
        else
           Me = ((real(fDens) - real(sDens)))/ &
                & (real(Nodens) - 1.0)
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
  open(10, file ='adf04', status = 'old')
  read(10,1000) Z
!
  levels = 0
!
! Calculate the number of levels in adf04 without having
! to specify
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
  read(10,2001) ind(i), config(i), mult(i), bEll(i), XJ(i), &
             &  en(i)
!
        if (i.eq.levels) then
           if (ind(i).eq.levels) then
             print *, 'Read of levels correct, total in &
                    &  ADF04: ', levels
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
  enddo
!
  read(10,*) dumm
!
! Just give exact number of temperatures
! during input and double check
!
  read(10,*) dumr, dumrr, (temp(i), i=1,(notemp))
  print *, 'Check these temperatures are correct:', &
        &  (temp(i), i=1,(notemp))
!
! I use tripi/j, hold and holde for holding variables
! to compute the ex/de-excitation rates and A-values.
!
  if (grad.eq.1) then
  open(5,file='adf04rad',status='unknown', form='formatted')
  endif
!
   aij = 0.0
  rate = 0.0
!
  order = ((levels**2-levels)/2)
!
  do i=1,order
     read(10,2003) tripj, tripi, hold, (holde(j), &
                &  j=1,notemp+einf)
!
! If grad = 1 then read A-values
! from a GRASPRAD ADF04RAD file
!
        if (grad.eq.1) then
           read(5,*) dumrr, dumm, hold
if ((dumrr.eq.6).and.(dumm.eq.2)) then
  print *, hold
  print *, tripj, tripi
  print *, dumrr, dumm
endif


              if (hold.lt.1e-20) then
                 hold = 0.0
              endif
        endif
!
    aij(tripj,tripi) = hold
    aij(tripi,tripj) = hold
!
! Excitation rate, depends on the
! temperature but not density.
!
     do k=1,notemp
      rate(tripi,tripj,k) = (const/weight(tripi))* &
           & sqrt((IH)/(boltz*temp(k)))* &
           & exp((en(tripi)-en(tripj))/ &
           & (boltz*temp(k)))*holde(k)
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
  2001 format(2X,I3,3X,A13,4X,I1,1X,I1,1X,F4.1,10X,F12.4)
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
! ########################################################
! ########################################################
!
! Next step is to fill the Collision matrix COLL..
! ########################################################
! ########################################################
!
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
                       COLL(i,j,k,Ne) = COLL(i,j,k,Ne) - &
                                      & aij(j,kk)
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
                    COLL(i,j,k,Ne) = COLL(i,j,k,Ne) + &
                    &  aij(i,j) + ((ED(Ne))*(rate(j,i,k)))
                 endif
              enddo
           enddo
        enddo
  enddo
!
  call cpu_time(finish1)
  print '("Basic Read = ",f6.3," seconds.")',finish1-start1
  call cpu_time(start2)
!
! Gaussian Elimination
! #########################################################
! #########################################################
!
  allocate(LHS((levels),notemp,Nodens))
!
  LHS = 0.0d0
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
                    LHS(j,kk,Ne) = LHS(j,kk,Ne) - &
                                 &  factor*LHS(i,kk,Ne)
                       do k=i+1,levels
                          COLL(j,k,kk,Ne) = COLL(j,k,kk,Ne) &
                                &  - factor*COLL(i,k,kk,Ne)
                       enddo
                 enddo
              enddo
        enddo
  enddo
!
  call cpu_time(finish2)
  print '("Gaussian elimination = ",f6.3," seconds.")', &
        &  finish2-start2
!
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
                    &  ((COLL(i,j,kk,Ne)*pop(j,kk,Ne))/ &
                    &  (COLL(i,i,kk,Ne)))
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
           do Ne=1,Nodens
            write(410 + i,*) (pop(i,k,Ne), k=mantmp,notemp)
           enddo
            close(410+i)
     enddo
  endif
!
!
! Loop over densities and temperatures now in
! Population arrays to obtain a specific line
! ratio defined by lolev and uplev
!
   allocate(lratio(Nodens,notemp))
!
! Keeping Density constant:
!
  prnt = 0
!
  do Ne=1,Nodens
     if (prnt.eq.1) then
        write(16,7008) Ne, ED(Ne)
           do kk=mantmp,notemp
              lratio(Ne,kk) = (pop(uplev,kk,Ne)* &
              &  aij(uplev,uplevl))/ &
              &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
           enddo
     elseif (prnt.eq.0) then
        if (Ne.eq.1) then
           write(16,*) '-> Increasing density'
           write(16,*)  '      Temperature          ', &
                    &  (ED(jj), jj=1,Nodens)
        else
           continue
        endif
!
           do kk=mantmp,notemp
              lratio(Ne,kk) = (pop(uplev,kk,Ne)* &
              & aij(uplev,uplevl))/ &
              & (pop(lolev,kk,Ne)*aij(lolev,lolevl))
           enddo
!
     endif
        if (Ne.eq.Nodens) then
           do kk=mantmp,notemp
              write(16,*) temp(kk), &
              & (lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
              write(404,*) temp(kk), &
              & (lratio(ii,kk)*(DE1/DE2), ii=1,Nodens)
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
              lratio(Ne,kk) = (pop(uplev,kk,Ne)* &
              &  aij(uplev,uplevl))/ &
              &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
              write(16,*) ED(Ne), lratio(Ne,kk)*(DE1/DE2)
           enddo
     elseif (prnt.eq.0) then
        do Ne=1,Nodens
           lratio(Ne,kk) = (pop(uplev,kk,Ne)*  &
           &  aij(uplev,uplevl))/ &
           &  (pop(lolev,kk,Ne)*aij(lolev,lolevl))
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
           write(16,*) ED(Ne), (lratio(Ne,ii)*(DE1/DE2), &
           &  ii=mantmp,notemp)
          write(405,*) ED(Ne), (lratio(Ne,ii)*(DE1/DE2), &
           &  ii=mantmp,notemp)
              enddo
           endif
          !
     endif
  enddo
!
! Coronal Approx and LTE:
!
  if (lte.eq.1) then
     write(16,*) 'LTE values:'

    do k=mantmp,notemp
       ltev = (weight(uplev)/weight(lolev))* &
          & exp(-(en(uplev)-en(lolev))/(boltz*temp(k)))
         write(16,*) k, 10.0e8, ltev*((aij(uplev,uplevl))/ &
                          &  (aij(lolev,lolevl)))*(DE1/DE2)
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
!
    do i=1,(lolev-1)
      totAlo = totAlo +  aij(lolev,i)
    enddo

    print *, 'tots A lo:', totAup
    print *, 'tots A:', totAlo

     do k=mantmp,notemp
        corappv = (rate(1,uplev,k)/rate(1,lolev,k))* &
         &  (totAlo/totAup)* &
         &  (aij(uplev,uplevl)/aij(lolev,lolevl))
         write(16,*) k, 10.0e2,(corappv*(DE1/DE2))
     enddo
  endif
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
  open(18,file='lineID',status='unknown',form='formatted')
  write(18,*) 'Upper state Lower State  Upper config &
  &    Lower config         Upper energy              &
  & Lower energy                wavelength'
!
  open(19,file='PEC_plot',status='unknown',form='formatted')
!
  do kk=mantmp,notemp
     if (ppec.eq.0) then
        write(19,*) 'Temperature: ', kk
        write(19,*) '--------------------'
        write(19,*)  'Wavelength                 Density  &
                  &   Value           PEC'
     endif
        do Ne=1,Nodens
           do j=2,levels
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
! ppec = 1 singles out a unique j -> i transition and
! loops over both temperature and density again.
!
  if (ppec.eq.1) then
     write(19,*) 'ppec = ', ppec
     write(19,*) 'Looking at PEC for transition', &
              &   uppec, '->', lopec
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
! This is all the wavelengths computed, but not
! in energy order. Need to sort by wavelength..
!
  if (worder.eq.1) then
     print *, 'Wavelength ordering requested, &
           &   with worder = ', worder
     continue
   else
     print *, 'No wavelength ordering requested, &
           &   with worder = ', worder
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
! The PEC are only dependent on temperature and
! density. So only adjust the wavelengths for
! one of the temperatures and densities.
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
                    !
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
!
  call cpu_time(finish3)
  print '("Rearranging Wavelength Order = ",f6.2," &
       & minutes.")', (finish3-start3)/60.0
!
  801 continue
!
  call cpu_time(finish)
  print '("Total Time = ",f6.2," minutes.")', &
       & (finish-start)/60.0
!
  if (line.eq.1) then
     !
  open(950,file='lratio',status='unknown',form='formatted')
  open(951,file='dens',status='unknown',form='formatted')
  open(952,file='temp',status='unknown',form='formatted')
     !
  write(950,*) 'Electron Density | Electron Temperature &
             & | Line Ratio'
  write(950,*) '-----------------------------------------&
             &-----------'
!
! Can use these vectors for 3D
! plotting in Matlab for example
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
!
! Format specifiers for output file
!
  7000 format(' --- RAdiative Transfer Output ---',/,/,' &
   &  - Basic Read of Input - ', /, &
   & 8X, '_____________________________________________')
  7001 format(8X, '| Atomic Number = ', 22X, I3, ' |', / ,&
   & 8X, '| Number of Temperatures = ', 13X, I3, ' |', / ,&
   & 8X, '| Number of Densities = ', 15X, I4, ' |', / ,&
   & 8X, '| Starting at 10^', I1, 'K and finishing at 10^',&
   & I2, 'K',1X, '|', /, 8X,'| with a mesh increment of ',&
   & F16.4, ' |',/,8X, &
   & '_____________________________________________', /)
  7002 format('RESERVED FOR SPECIFIC LINES NAMELIST LINE.')
  7003 format('  - Basic Level Information -  ', /, /,&
   & 'Total number of levels in ADF04 file: ', I3, /, /,&
   & '_____________________________________________', /,&
   & '| IND.|   CONFIGURATION     |  ENERGY (eV)  |', /,&
   & '_____________________________________________')
  7004 format('| ',I3,' | '3X,A13,3X,' | ',F12.4,'  |')
  7005 format('__________________________________________',&
   & /, /, 'The new left hand side values from N_i = 1 - ',&
   & I3,' for a particular temperature and density &
   & (normally 1 and 1)',/)
  7006 format(8D12.4)
  7007 format(/, 'Looking at the specific line ratio of', I3,&
   & '->1/', I3,'->1 over temperatures and densities.')
  7008 format(/, '-------------- Density ...', 4X, I3, &
   & '(',E10.4,')')
  7009 format(2X,E10.4,4X,F20.16)
  7010 format(/, '-------------- Temperature ...', 4X, I3,&
   & '(',E10.4,')')
!
end program
