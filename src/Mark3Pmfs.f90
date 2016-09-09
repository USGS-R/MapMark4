
! Mark3Pmf
!
! Description: Calculates a probability distribution for the number of deposits, using the
! methodology in Root and others (1992).
!
! Arguments:
! ND()      Lists the thresholds for the calculated probability distribution (see Details).
! sizeND    Size of vector ND
! INOD      Number of specified thresholds. The only allowed values are 3, 5, 7, and 9.
! XX()      Lists the calculated probabilities.
! sizeXX    Size of vector XX.
! status    Indicates whether any of the input varibables are incorrect (see Details).
!
! Details:
! The configurations for ND depends on INOD:
!
!                                      Description
!    Element    INOD = 3                   INOD = 5                   INOD = 7                   INOD = 9
!    ND(1)      not used                   not used                   not used                   not used
!    ND(2)      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9
!    ND(3)      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5
!    ND(4)      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1
!    ND(5)      not used                   P( N >= ND(5) ) = 0.05     P( N >= ND(5) ) = 0.05     P( N >= ND(5) ) = 0.05
!    ND(6)      not used                   P( N >= ND(6) ) = 0.01     P( N >= ND(6) ) = 0.01     P( N >= ND(6) ) = 0.01
!    ND(7)      not used                   not used                   P( N >= ND(7) ) = 0.005    P( N >= ND(7) ) = 0.005
!    ND(8)      not used                   not used                   P( N >= ND(8) ) = 0.001    P( N >= ND(8) ) = 0.001
!    ND(9)      not used                   not used                   not used                   P( N >= ND(9) ) = 0.0005
!    ND(10)     not used                   not used                   not used                   P( N >= ND(10) ) = 0.0001
!
! The expression P( N >= ND(2) ) = 0.9 means the probability that the number of the deposits (N)
! will be at least ND(2) is 0.9. That is, there is a 0.9 probability (90% chance) of finding ND(2) or more
! deposits. The elements of ND must be nondecreasing. For example,
! if INOD is 3, then ND(2) <= ND(3) <= ND(4). Although it seems that the elements of ND should
! be strictly increasing, this is not required by the algorithm. The size of ND must be at least
! INOD+1. Unfortunately, ND is modified within the subroutine.
!
! If ND(2) is specified as 0, then, according to the previous definition, P( N >= 0 ) = 0.9. Of course,
! this is wrong because P( N >= 0 ) = 1. To address this problem, the algorithm is implemented such that
! P( N >= 1 ) < 0.9. (This doesn't make sense to me. KJE)
!
! The configuration of XX is:
!
!                                         Probability of having
!    Element              (INOD = 3)       (INOD = 5)       (INOD = 7)       (INOD = 9)
!    XX(1)                0 deposits       0 deposits       0 deposits       0 deposits
!    XX(2)                1 deposit        1 deposit        1 deposit        1 deposit
!    XX(3)                2 deposits       2 deposits       2 deposits       2 deposits
!    ...
!    XX( ND(INOD+1)+1 )   ND(4) deposits   ND(6) deposits   ND(8) deposits   ND(10) deposits
!
! The size of XX must be at least ND(INOD+1)+1. (One is added because the first element of XX is
! the probability of zero deposits.
!
! If status equals 1, then INOD is not 3, 5, 7, or 9. If status equals 10, then the elements of ND are not
! nondecreasing. If status equals 11, then both of these errors occured.
!
! This code was extracted from the Mark3B program (Root and others, 1998). In that program
! are subroutines DISTRIB3, DISTRIB5, and DISTRIB9. The parts of these subroutines that calculate the
! probabilitydistribution are practically identical, so they were replaced by this one subroutine. There were
! only three modifications: (1) explicit declaration of the varibles, (2) checks for errors
! in the input variables, (3) replaced variables DISTP3/DISTP5/DISTP9 with a single variable DISTP,
! and (4) revision of the last executable statement (so that it works with all three values of INOD).
!
! All variables (except four) are identical to the corresponding variables in program Mark3B, to make it
! easy to compare the codes. The exceptions are new variables: sizeND, sizeXX, status, and error, which are
! needed to check for errors.
!
! Dates: 26 Jan 2012, Karl J. Ellefsen: Adapted for use in a dll
!
! References:
! Root, D.H., Menzie, W.D., and Scott, W.A., 1992, Computer Monte Carlo Simulation in Quantitative
! Resource Estimation: Nonrenewable resources, v. 1, no. 2, p. 125-138.
!
! Root, D.H., Scott, W.A., Jr., Schruben, P.G., 1998, MARK3B Resource assessment program for Macintosh:
! U.S. Geological Survey, Open-file report 98-356.

      SUBROUTINE Mark3Pmf(ND,sizeND,INOD,XX,sizeXX,status)
    implicit none
    integer, intent(in) :: sizeXX, sizeND,INOD
    integer, intent(in out), dimension(sizeND) :: ND
    integer, intent(out) :: status
    real(kind=8), intent(out), dimension(sizeXX) :: XX

    real(kind=8), dimension(10) :: D
    real(kind=8), dimension(9) :: DISTP
    real(kind=8) :: DUM
    integer :: I, L, IDUM, error

!   DISTP lists the probabilities associated with the thresholds in ND.
!
!   When INOD = 3, the first 3 elements of DISTP are used.
!   When INOD = 5, the first 5 elements of DISTP are used.
!   When INOD = 7, the first 7 elements of DISTP are used.
!   When INOD = 9, all 9 elements of DISTP are used.
    DISTP = [ 0.90, 0.50, 0.10, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001 ]

    status = 0
    if( INOD /= 3 .AND. INOD /= 5 .AND. INOD /= 7 .AND. INOD /= 9 ) status = 1

    ! the elements of ND must be non-decreasing
    error = 0
    do i = 2, inod
        if( nd(i) > nd(i+1) ) error = 10
    end do
    if( error == 10 ) status = status + 10

    if( status /= 0 ) return

20    ND(1) = 0
      DO 1 I = 1,INOD+1
        ND(I) = ND(I) + 1
        D(I) = ND(I)
1     CONTINUE

      DO 2 L = 1, ND(INOD+1)
        XX(L) = 0.
2     CONTINUE

      IDUM = ND(2) - ND(1)
      IF (IDUM.EQ.0) XX(ND(1)) = 1.-DISTP(1)
      IF (IDUM.EQ.0) GO TO 12
      IF (IDUM.GT.0) DUM = (1.-DISTP(1))/(2*(D(2)-D(1))+1)
      DO 3 L = ND(1), ND(2)
        IF (L.EQ.ND(2)) XX(L) = XX(L)+DUM
        IF (L.LT.ND(2)) XX(L) = XX(L) + 2.*DUM
3     CONTINUE

12    DO 4 I = 2,INOD
        IDUM = ND(I+1)-ND(I)
        IF (IDUM.EQ.0) DUM = DISTP(I-1)-DISTP(I)
        IF (IDUM.EQ.1) DUM = (DISTP(I-1)-DISTP(I))/2.
        IF (IDUM.GT.1) DUM = (DISTP(I-1)-DISTP(I))/(2.*(IDUM-1.)+2.)
        DO 5 L = ND(I), ND(I + 1)
            IF (L.EQ.ND(I).OR.L.EQ.ND(I+1)) XX(L) = XX(L) + DUM
            IF (L.GT.ND(I).AND.L.LT.ND(I+1)) XX(L) = XX(L)+2.*DUM
5       CONTINUE
4     CONTINUE

      XX(ND(INOD+1)) = XX(ND(INOD+1)) + DISTP(INOD)

      RETURN
      END

! Mark3RevisedPmf
!
! Description: Calculates a probability distribution for the number of deposits, using an
! algorithm developed and implemented by Jeffery D. Phillips.
!
! Arguments:
! ND()      Lists the thresholds for the calculated a probability distribution (see Details).
! sizeND    Size of vector ND
! INOD      Number of specified thresholds. The only allowed values are 3, 5, 7, and 9.
! maxd      Maximum number of deposits for which the probability is non-zero.
! XX()      Lists the calculated probabilities.
! sizeXX    Size of vector XX.
! status    Indicates whether any of the input varibables are incorrect (see Details).
!
! Details:
! The algorithm is based on the following three criteria:
! 1. The calculated probabilities should be as close as possible to the specified probabalities. (See the
! table below.
! 2. The maximum number of deposits will be maxd.
! 3. The probabilities between the specified quantilies will be constant. This criterion requires
! qualification, which will be presented with an example: Assume that INOD is 3, ND(3) is 5 and ND(4) is 9.
! Then, when the number of deposits is 6, 7, or 8, the associated probabilities will be equal (and hence
! constant).
!
! The configurations for ND depends on INOD:
!
!                                      Description
!    Element    INOD = 3                   INOD = 5                   INOD = 7                   INOD = 9
!    ND(1)      not used                   not used                   not used                   not used
!    ND(2)      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9      P( N >= ND(2) ) = 0.9
!    ND(3)      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5      P( N >= ND(3) ) = 0.5
!    ND(4)      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1      P( N >= ND(4) ) = 0.1
!    ND(5)      not used                   P( N >= ND(5) ) = 0.05     P( N >= ND(5) ) = 0.05     P( N >= ND(5) ) = 0.05
!    ND(6)      not used                   P( N >= ND(6) ) = 0.01     P( N >= ND(6) ) = 0.01     P( N >= ND(6) ) = 0.01
!    ND(7)      not used                   not used                   P( N >= ND(7) ) = 0.005    P( N >= ND(7) ) = 0.005
!    ND(8)      not used                   not used                   P( N >= ND(8) ) = 0.001    P( N >= ND(8) ) = 0.001
!    ND(9)      not used                   not used                   not used                   P( N >= ND(9) ) = 0.0005
!    ND(10)     not used                   not used                   not used                   P( N >= ND(10) ) = 0.0001
!
! The expression P( N >= ND(2) ) = 0.9 means the probability that the number of the deposits (N)
! will be at least ND(2) is 0.9. That is, there is a 0.9 probability (90% chance) of finding ND(2) or more
! deposits. The elements of ND must be nondecreasing. For example,
! if INOD is 3, then ND(2) <= ND(3) <= ND(4). Although it seems that the elements of ND should
! be strictly increasing, this is not required by the algorithm. The size of ND must be at least
! INOD+1.
!
! If ND(2) is specified as 0, then, according to the previous definition, P( N >= 0 ) = 0.9. Of course,
! this is wrong because P( N >= 0 ) = 1. To address this problem, the algorithm is implemented such that
! P( N >= 1 ) < 0.9.
!
! The configuration of XX is:
!
!    Element              Probability of having
!    XX(1)                0 deposits
!    XX(2)                1 deposit
!    XX(3)                2 deposits
!    ...
!    XX( maxd+1 )         maxd deposits
!
! maxd must be greater than or equal to ND(INOD+1). The size of XX must be at least maxd+1. (One is added
! because the first element of XX is the probability of zero deposits.
!
! If status equals 1, then INOD is not 3, 5, or 9. (Error 1)
! If status equals 10, then the elements of ND are not nondecreasing. (Error 2)
! If status equals 100, then maxd is less than ND(INOD+1). (Error 3)
! If status equals 11, then Errors 1 and 2 occured.
! If status equals 101, then Errors 1 and 3 occured.
! If status equals 110, then Errors 2 and 3 occured.
! If status equals 111, then Errors 1, 2 and 3 occured.
!
! Dates: 2011, Jeffrey D. Phillips, written.
!        30 Jan 2012, Karl J. Ellefsen: Adapted for use in a dll
!
! References:
! Currently none.

    subroutine Mark3RevisedPmf( ND, sizeND, INOD, maxd, XX, sizeXX, status )

    implicit none
    integer, intent(in) :: sizeXX, sizeND, INOD, maxd
    integer, intent(in), dimension(sizeND) :: ND
    integer, intent(out) :: status
    real(kind=8), intent(out), dimension(sizeXX) :: XX

    real(kind=8) :: frac1
    integer :: I, error

    status = 0
    if( INOD /= 3 .AND. INOD /= 5 .AND. INOD /= 7 .AND. INOD /= 9 ) status = 1

    ! the elements of ND must be non-decreasing
    error = 0
    do i = 2, inod
        if( nd(i) > nd(i+1) ) error = 10
    end do
    if( error == 10 ) status = status + 10

    if( maxd < ND(INOD+1) ) status = status + 100

    if( status /= 0 ) return


      do i=1,sizeXX
        xx(i)=0.
      enddo
!
      if(inod.eq.3) then
        frac1=0.1
      else if(inod.eq.5) then
        frac1=0.01
      else if(inod.eq.7) then
        frac1=0.001
      else if(inod.eq.9) then
        frac1=0.0001
      endif
! distribute probabilities at the high end
      if(maxd.gt.nd(inod+1)) then
        do i=nd(inod+1)+1,maxd+1
          xx(i)=frac1/float(maxd+1-nd(inod+1))
        enddo
      else
        xx(nd(inod+1)+1)=frac1
      endif
      if(inod.eq.9) then
! distribute the next .04% of deposits
        if(nd(9).gt.0) then
          if(nd(10).gt.nd(9)) then
            do i=nd(9)+1,nd(10)
            xx(i)=.0004/float(nd(10)-nd(9))
            enddo
          else
            xx(nd(10)+1)=xx(nd(10)+1)+0.0004
          endif
! distribute the next .05% of deposits
          if(nd(8).gt.0) then
            if(nd(9).gt.nd(8)) then
              do i=nd(8)+1,nd(9)
              xx(i)=.0005/float(nd(9)-nd(8))
              enddo
            else
              xx(nd(9)+1)=xx(nd(9)+1)+0.0005
            endif
            endif
        endif
      endif
      if(inod.ge.7) then
! distribute the next .4% of deposits
        if(nd(7).gt.0) then
          if(nd(8).gt.nd(7)) then
            do i=nd(7)+1,nd(8)
            xx(i)=.004/float(nd(8)-nd(7))
            enddo
          else
            xx(nd(8)+1)=xx(nd(8)+1)+0.004
          endif
! distribute the next .5% of deposits
          if(nd(6).gt.0) then
            if(nd(7).gt.nd(6)) then
              do i=nd(6)+1,nd(7)
              xx(i)=.005/float(nd(7)-nd(6))
              enddo
            else
              xx(nd(7)+1)=xx(nd(7)+1)+0.005
            endif
          endif
        endif
      endif
      if(inod.ge.5) then
! distribute the next 4% of deposits
        if(nd(5).gt.0) then
          if(nd(6).gt.nd(5)) then
            do i=nd(5)+1,nd(6)
            xx(i)=.04/float(nd(6)-nd(5))
            enddo
          else
            xx(nd(6)+1)=xx(nd(6)+1)+0.04
          endif
! distribute the next 5% of deposits
          if(nd(4).gt.0) then
            if(nd(5).gt.nd(4)) then
              do i=nd(4)+1,nd(5)
              xx(i)=.05/float(nd(5)-nd(4))
              enddo
            else
              xx(nd(5)+1)=xx(nd(5)+1)+0.05
            endif
          endif
        endif
      endif
      if(inod.ge.3) then
! distribute the 40% of deposits above the median
        if(nd(3).gt.0.) then
          if(nd(4).gt.nd(3)) then
            do i=nd(3)+1,nd(4)
            xx(i)=.4/float(nd(4)-nd(3))
            enddo
          else
            xx(nd(4)+1)=xx(nd(4)+1)+0.4
          endif
! distribute the 40% of deposits below the median
          if(nd(2).gt.0) then
            if(nd(3).gt.nd(2)) then
              do i=nd(2)+1,nd(3)
              xx(i)=.4/float(nd(3)-nd(2))
              enddo
            else
              xx(nd(3)+1)=xx(nd(3)+1)+0.4
            endif
          endif
        endif
      endif
! distribute the 10%, 50%, or 90% of deposits at the low end
      if(nd(2).gt.0) then
        do i=1,nd(2)
        xx(i)=0.1/float(nd(2))
        enddo
      else if(nd(3).gt.0) then
        do i=1,nd(3)
        xx(i)=0.5/float(nd(3))
        enddo
      else if(nd(4).gt.0) then
        do i=1,nd(4)
        xx(i)=0.9/float(nd(4))
        enddo
      endif

      RETURN
      END
