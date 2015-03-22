! --------------------------------------------------------------------------
! erlasso.f90: the GCD algorithm for expectile regression.
! --------------------------------------------------------------------------
! call dblepr('x', -1, x, nobs * nvars)
! call intpr()
! USAGE:
! 
! call erlasso(tau, nobs, nvars, x, y, pf, pf2, dfmax, pmax, nlam, &
! & flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, &
! & nbeta, alam, npass, jerr)
! 
! INPUT ARGUMENTS:
!
!    tau = the asymmetry coefficient in the expectile regression model.
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of covariates, dimension N * p; row -> observation.
!    y(nobs) = response variable.
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                pf2(j) = 0 => jth variable unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    beta(pmax, nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS:
! YUWEN GU 04/18/2014
! ------------------------------------------------------------------------------------- !
SUBROUTINE erlasso(w, tau, nobs, nvars, x, y, pf, pf2, dfmax, pmax, &
& nlam, flmin, ulam, eps, maxit, nalam, beta, ibeta, nbeta, theta, &
& itheta, ntheta, alam, npass, jerr)
! ------------------------------------------------------------------------------------- !
      IMPLICIT NONE
! ----- passed variables ----- !
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: itheta (pmax)
      INTEGER :: ntheta (nlam)
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: tau
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: theta (pmax, nlam)
      DOUBLE PRECISION :: alam (nlam)
! ----- local declarations ----- !
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! ---- allocate variables ----- !

      

      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju) !!! check sample covariance ?= 0


      IF (maxval(ju) <= 0) THEN
        jerr = 7777
        RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
        jerr = 10000
        RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
        jerr = 10000
        RETURN
      END IF
      pf = Max(0.0D0, pf)
      pf2 = Max(0.0D0, pf2)


      DO j=1, nvars                                  
        IF(ju(j)==1) THEN
          maj(j) = dot_product(x(:,j), x(:,j))/nobs
        END IF
      END DO


! ------- call erlassopath -------- !
      CALL erlassopath (w, tau, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
      & pmax, nlam, flmin, ulam, eps, maxit, nalam, beta, ibeta, nbeta, &
      & theta, itheta, ntheta, alam, npass, jerr)
      IF (jerr > 0) RETURN !!! check error after calling function
      DEALLOCATE (ju, maj)
      RETURN
END SUBROUTINE erlasso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ---------------- erlassopath ---------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE erlassopath (w, tau, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, beta, ibeta, nbeta, &
& theta, itheta, ntheta, alam, npass, jerr)
! ---------------------------------------------------------------------------------- !
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: big = 9.9D30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0D-06
      INTEGER, PARAMETER :: mnlam = 6
! ----- passed variables ----- !
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      INTEGER :: itheta (pmax)
      INTEGER :: ntheta (nlam)
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: w
      DOUBLE PRECISION :: tau
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: theta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
      DOUBLE PRECISION :: flmin
! ----- local declarations ----- !
      DOUBLE PRECISION :: bigm
      DOUBLE PRECISION :: d
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: oldth
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: dl (nobs)
      INTEGER :: mnl
      INTEGER :: i
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ctr
      INTEGER :: ierr
      INTEGER :: nib
      INTEGER :: nith
      INTEGER :: me
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: th
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldtheta
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r1
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r2
      INTEGER, DIMENSION (:), ALLOCATABLE :: mmb
      INTEGER, DIMENSION (:), ALLOCATABLE :: mmth
! ----- allocate variables ----- !
      ALLOCATE (b(1:nvars), STAT=jerr)
      ALLOCATE (oldbeta(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (th(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (oldtheta(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mmb(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mmth(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r1(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r2(1:nobs), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
! ----- some initial setup ----- !
      r1 = y
      r2 = y
      b = 0.0D0
      oldbeta = 0.0D0
      th = 0.0D0
      oldtheta = 0.0D0
      ibeta = 0
      itheta = 0
      mmb = 0
      mmth = 0
      npass = 0
      nib = npass
      nith = npass
      mnl = Min(mnlam, nlam)
      bigm = 2.0D0 * Max((1.0D0 - tau), tau)
      IF (flmin < 1.0D0) THEN
        flmin = Max(mfl, flmin)
        alf = flmin ** (1.0D0/(nlam-1.0D0))
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! ----------------- lambda loop ------------------- !
          ! ------------ This is the outmost loop ------------ !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO l = 1, nlam !!! begin lambda loop (the outmost loop)
      ! --------- computing lambda ------------------- !

        IF (flmin >= 1.0D0) THEN
          al = ulam (l)
        ELSE
          IF (l > 2) THEN
            al = al * alf
          ELSE IF (l == 1) THEN
            al = big
          ELSE IF (l == 2) THEN
            al = 0.0D0
            DO i = 1, nobs
              IF (r2(i) < 0.0D0) THEN
                dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
              ELSE
                dl(i) = 2.0D0 * tau * r2(i)
              END IF
            END DO
            DO j = 1, nvars
              IF (ju(j) /= 0) THEN
                IF (pf(j) > 0.0D0 .AND. pf2(j) > 0.0D0) THEN
                  u = dot_product(dl, x(:, j))
                  al = Max(al, Abs(u) / pf2(j), Abs(w * &
                       dot_product(r1, x(:, j)) + u) / pf(j))
                END IF
              END IF
            END DO
            al = al * alf / nobs
          END IF
        END IF
        ctr = 0

        ! ------------------ outer loop -------------------- !
        DO !!! begin outer loop
          IF (nib > 0) oldbeta(ibeta(1:nib)) = b(ibeta(1:nib))
          IF (nith > 0) oldtheta(itheta(i:nith)) = th(itheta(1:nith))

        ! ----------------- middle loop -------------------- !

          DO !!! begin middle loop
            npass = npass + 1

            call intpr('passes', -1, npass, 1)

            dif = 0.0D0
            DO k = 1, nvars
              IF (ju(k) /= 0) THEN
                oldb = b(k)
                u = 0.0D0
                DO i = 1, nobs
                  IF (r2(i) < 0.0D0) THEN
                    dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
                  ELSE
                    dl(i) = 2.0D0 * tau * r2(i)
                  END IF
                  u = u + (w * r1(i) + dl(i)) * x (i, k)
                END DO
                u = maj(k) * b(k) * (w + bigm) + u / nobs
                v = al * pf(k)
                v = Abs(u) - v
                IF (v > 0.0D0) THEN
                  b(k) = sign(v, u) / (maj(k) * (w + bigm))
                ELSE
                  b(k) = 0.0D0
                END IF
                d = b(k) - oldb
                IF (ABS(d) > 0.0D0) THEN
                  dif = MAX(dif, bigm * d**2)
                  r1 = r1 - x(:, k) * d
                  r2 = r2 - x(:, k) * d
                  IF (mmb(k) == 0) THEN
                    nib = nib + 1
                    IF (nib > pmax) EXIT
                    mmb(k) = nib
                    ibeta(nib) = k !indicate which coefficient is non-zero
                  END IF
                END IF
              END IF
            END DO

            DO k = 1, nvars
              IF (ju(k) /= 0) THEN
                oldth = th(k)
                u = 0.0D0
                DO i = 1, nobs
                  IF (r2(i) < 0.0D0) THEN
                    dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
                  ELSE
                    dl(i) = 2.0D0 * tau * r2(i)
                  END IF
                  u = u + dl(i) * x(i, k)
                END DO
                u = maj(k) * th(k) * bigm + u / nobs
                v = al * pf2(k)
                v = Abs(u) - v
                IF (v > 0.0D0) THEN
                  th(k) = sign(v, u) / (maj(k) * bigm)
                ELSE
                  th(k) = 0.0D0
                END IF
                d = th(k) - oldth
                IF (ABS(d) > 0.0D0) THEN
                  dif = MAX(dif, bigm * d**2)
                  r2 = r2 - x(:, k) * d
                  IF (mmth(k) == 0) THEN
                    nith = nith + 1
                    IF (nith > pmax) EXIT
                    mmth(k) = nith
                    itheta(nith) = k !indicate which coefficient is non-zero
                  END IF
                END IF
              END IF
            END DO
            IF (dif < eps) EXIT
            IF (nib > pmax) EXIT
            IF (nith > pmax) EXIT
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ! -----------------inner loop ------------------- !
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             DO !!! begin inner loop
!               npass = npass + 1
!               dif = 0.0D0
!               DO j = 1, nib
!                 k = ibeta(j)
!                 oldb = b(k)
!                 u = 0.0D0
!                 DO i = 1, nobs
!                   IF (r2(i) < 0.0D0) THEN
!                     dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
!                   ELSE
!                     dl(i) = 2.0D0 * tau * r2(i)
!                   END IF
!                   u = u + (w * r1(i) + dl(i)) * x(i, k)
!                 END DO
!                 u = maj(k) * b(k) * (w + bigm) + u / nobs
!                 v = al * pf(k)
!                 v = Abs(u) - v
!                 IF (v > 0.0D0) THEN
!                   b(k) = sign(v, u) / (maj(k) * (w + bigm))
!                 ELSE
!                   b(k) = 0.0D0
!                 END IF
!                 d = b(k) - oldb
!                 IF (ABS(d) > 0.0D0) THEN
!                   dif = MAX(dif, bigm * d**2)
!                   r1 = r1 - x(:, k) * d
!                   r2 = r2 - x(:, k) * d
!                 END IF
!               END DO      

!               DO j = 1, nith
!                 k = itheta(j)
!                 oldth = th(k)
!                 u = 0.0D0
!                 DO i = 1, nobs
!                   IF (r2(i) < 0.0D0) THEN
!                     dl (i) = 2.0D0 * (1.0D0 - tau) * r2(i)
!                   ELSE
!                     dl (i) = 2.0D0 * tau * r2(i)
!                   END IF
!                   u = u + dl(i) * x(i, k)
!                 END DO
!                 u = maj(k) * th(k) * bigm + u / nobs
!                 v = al * pf2(k)
!                 v = Abs(u) - v
!                 IF (v > 0.0D0) THEN
!                   th(k) = sign(v, u) / (maj(k) * bigm)
!                 ELSE
!                   th(k) = 0.0D0
!                 END IF
!                 d = th(k) - oldth
!                 IF (ABS(d) > 0.0D0) THEN
!                   dif = MAX(dif, bigm * d**2)
!                   r2 = r2 - x(:, k) * d
!                 END IF
!               END DO 
!               IF (dif < eps) EXIT
!             END DO !!! end inner loop
          END DO !!! end middle loop
          IF (nib > pmax) EXIT
          IF (nith > pmax) EXIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! -------------- final check ---------------- !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          vrg = 1
          DO j = 1, nib
            IF ((b(ibeta(j))-oldbeta(ibeta(j)))**2 >= eps) THEN
              vrg = 0
              EXIT
            END IF
          END DO
          DO j = 1, nith
            IF ((th(itheta(j))-oldtheta(itheta(j)))**2 >= eps) THEN
              vrg = 0
              EXIT
            END IF
          END DO
          IF (vrg == 1) EXIT
          ctr = ctr + 1
          IF (ctr > maxit) THEN
            jerr = - l
            RETURN
          END IF
        END DO !!! end outer loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ----------- final update & save results ------------ !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (nib > pmax) THEN
          jerr = - 10000 - l
          EXIT
        END IF
        IF (nith > pmax) THEN
          jerr = - 10000 - l
          EXIT
        END IF
        IF (nib > 0) beta(1:nib, l) = b(ibeta(1:nib))
        nbeta(l) = nib
        IF (nith > 0) theta(1:nith, l) = th(itheta(1:nith))
        ntheta(l) = nith
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        me = count(ABS(beta(1:nib, l)) > 0.0D0)
        IF (me > dfmax) EXIT
        me = count(ABS(theta(1:nith, l)) > 0.0D0)
        IF (me > dfmax) EXIT
      END DO !!! end lambda loop (the outmost loop)
      DEALLOCATE(b, oldbeta, th, oldtheta, r1, r2, mmb, mmth)
      RETURN
END SUBROUTINE erlassopath
