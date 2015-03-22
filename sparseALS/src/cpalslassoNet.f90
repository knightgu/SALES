! -------------------------------------------------------------------------------
! cpalslasso.f90: the coordinate descent algorithm for the couple ALS regression.
! -------------------------------------------------------------------------------
!
! USAGE:
! 
! call cpalslasso(tau, nobs, nvars, x, y, pf, pf2, dfmax, pmax, nlam, &
! & flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, &
! & nbeta, alam, npass, jerr)
! 
! INPUT ARGUMENTS:
!
!    tau = the asymmetry coefficient in the expectile regression model.
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs,nvars) = matrix of covariates, dimension N * p; row -> observation.
!    y(nobs) = response variable.
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                pf2(j) = 0 => jth variable unpenalized
!    dfmax = limit the MAXimum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the MAXimum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => MINimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    maxit = MAXimum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    beta(pmax,nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => MAXval(vp) <= 0.0
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
! YUWEN GU 03/13/2015
!
! ------------------------------------------------------------------------------------- !
SUBROUTINE cpalslasso(w, tau, nobs, nvars, x, y, jd, pf, pf2, dfmax, pmax, &
& nlam, flmin, ulam, eps, isd, maxit, nalam, b0,beta, ibeta, nbeta, &
& t0, theta, itheta, ntheta, alam, npass, jerr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs,nvars
  INTEGER :: dfmax,pmax
  INTEGER :: jd(*)
  INTEGER :: nlam,nalam
  INTEGER :: isd,npass,maxit,jerr
  INTEGER :: ibeta(pmax)
  INTEGER :: nbeta(nlam)
  INTEGER :: itheta(pmax)
  INTEGER :: ntheta(nlam)
  DOUBLE PRECISION :: flmin
  DOUBLE PRECISION :: eps
  DOUBLE PRECISION :: w
  DOUBLE PRECISION :: tau
  DOUBLE PRECISION :: x(nobs, nvars)
  DOUBLE PRECISION :: y(nobs)
  DOUBLE PRECISION :: pf(nvars)
  DOUBLE PRECISION :: pf2(nvars)
  DOUBLE PRECISION :: ulam(nlam)
  DOUBLE PRECISION :: beta(pmax,nlam),b0(nlam)
  DOUBLE PRECISION :: theta(pmax,nlam),t0(nlam)
  DOUBLE PRECISION :: alam(nlam)
  ! --------- LOCAL DECLARATIONS ---------- !
  INTEGER :: j
  INTEGER :: l
  INTEGER :: nk,nkk
  INTEGER :: ierr
  INTEGER, DIMENSION(:), ALLOCATABLE :: ju
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean,xnorm,maj
  ! --------- ALLOCATE VARIABLES ---------- !
  ALLOCATE(ju(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(xmean(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(maj(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(xnorm(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  IF (jerr /= 0) RETURN
  CALL chkvars(nobs, nvars, x, ju)
  IF (jd(1) > 0) ju(jd(2:(jd(1) + 1))) = 0 ! EXCLUDED VARIABLES
  IF (MAXVAL(ju) <= 0) THEN
    jerr = 7777
    RETURN
  END IF
  IF (MAXVAL(pf) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  IF (MAXVAL(pf2) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  pf = MAX(0.0D0, pf)
  pf2 = MAX(0.0D0, pf2)
  CALL standard(nobs, nvars, x, ju, isd, xmean, xnorm, maj)
  ! ------- CALL cpalslassopath -------- !
  CALL cpalslassopath(w, tau, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
  & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, nbeta, &
  & t0, theta, itheta, ntheta, alam, npass, jerr)
  IF (jerr > 0) RETURN ! CHECK ERROR AFTER CALLING FUNCTION
  ! ----------- TRANSFORM BETA BACK TO THE ORIGINAL SCALE ----------- !
  DO l = 1, nalam
    nk = nbeta(l)
    nkk = ntheta(l)
    IF (isd == 1) THEN
      DO j = 1, nk
        beta(j,l) = beta(j,l)/xnorm(ibeta(j))
      END DO
      DO j = 1, nkk
        theta(j,l) = theta(j,l)/xnorm(itheta(j))
      END DO
    END IF
    b0(l) = b0(l) - DOT_PRODUCT(beta(1:nk,l),xmean(ibeta(1:nk)))
    t0(l) = t0(l) - DOT_PRODUCT(theta(1:nkk,l),xmean(itheta(1:nkk)))
  END DO
  DEALLOCATE(ju,xmean,xnorm,maj)
  RETURN
END SUBROUTINE cpalslasso

! ---------------------------------- cpalslassopath -------------------------------- !
SUBROUTINE cpalslassopath(w, tau, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, nbeta, &
& t0, theta, itheta, ntheta, alam, npass, jerr)

  IMPLICIT NONE     
  ! --------------- INPUT VARIABLES --------------- !
  INTEGER :: nobs,nvars,dfmax,pmax
  INTEGER :: nlam,nalam
  INTEGER :: maxit,npass,jerr
  INTEGER :: ju(nvars)
  INTEGER :: ibeta(pmax),nbeta(nlam)
  INTEGER :: itheta(pmax),ntheta(nlam)
  DOUBLE PRECISION :: eps,w,tau
  DOUBLE PRECISION :: x(nobs,nvars),y(nobs)
  DOUBLE PRECISION :: pf(nvars),pf2(nvars)
  DOUBLE PRECISION :: beta(pmax,nlam)
  DOUBLE PRECISION :: theta(pmax,nlam)
  DOUBLE PRECISION :: ulam(nlam)
  DOUBLE PRECISION :: alam(nlam)
  DOUBLE PRECISION :: maj(nvars)
  DOUBLE PRECISION :: flmin
  ! ------------- LOCAL DECLARATIONS -------------- !
  INTEGER, PARAMETER :: mnlam = 6
  DOUBLE PRECISION, PARAMETER :: big = 9.9D30,mfl = 1.0D-06
  DOUBLE PRECISION :: bigm
  DOUBLE PRECISION :: d
  DOUBLE PRECISION :: dif
  DOUBLE PRECISION :: oldb,oldth
  DOUBLE PRECISION :: u
  DOUBLE PRECISION :: v
  DOUBLE PRECISION :: al
  DOUBLE PRECISION :: alf
  DOUBLE PRECISION :: dl(nobs)
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
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: b,th
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oldbeta,oldtheta
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r1,r2
  INTEGER, DIMENSION (:), ALLOCATABLE :: mmb,mmth
  ! ----------------- ALLOCATE VARIABLES ------------------ !
  ALLOCATE (b(0:nvars), STAT=jerr)
  ALLOCATE (oldbeta(0:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE (th(0:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE (oldtheta(0:nvars), STAT=ierr)
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
  ! -------------------- INITIALIZATION ------------------- !
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
  mnl = MIN(mnlam, nlam)
  bigm = 2.0D0 * MAX((1.0D0 - tau), tau)
  IF (flmin < 1.0D0) THEN
    flmin = MAX(mfl, flmin)
    alf = flmin ** (1.0D0/(nlam-1.0D0))
  END IF
  ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
  DO l = 1, nlam
    ! --------------- COMPUTE LAMBDA ------------------- !
    IF (flmin >= 1.0D0) THEN
      al = ulam(l)
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
            IF (pf(j) > 0.0D0) THEN
              u = DOT_PRODUCT(dl, x(:,j))
              al = MAX(al,ABS(w*DOT_PRODUCT(r1,x(:,j))+u)/pf(j))
            END IF
            IF (pf2(j) > 0.0D0) THEN
              IF (pf(j) <= 0.0D0) u = DOT_PRODUCT(dl, x(:,j))
              al = MAX(al, ABS(u)/pf2(j))
            END IF
          END IF
        END DO
        al = al * alf/nobs
      END IF
    END IF
    ctr = 0
    ! ------------------ OUTER LOOP -------------------- !
    DO
      oldtheta(0) = b(0)
      IF (nib > 0) oldbeta(ibeta(1:nib)) = b(ibeta(1:nib))
      oldtheta(0) = th(0)
      IF (nith > 0) oldtheta(itheta(1:nith)) = th(itheta(1:nith))
      ! ----------------- MIDDLE LOOP -------------------- !
      DO
        npass = npass + 1
        dif = 0.0D0
        DO k = 1, nvars
          IF (ju(k) /= 0) THEN
            oldb = b(k)
            DO ! BEGIN PROXIMAL GRADIENT DESCENT
              u = 0.0D0
              DO i = 1, nobs
                IF (r2(i) < 0.0D0) THEN
                  dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
                ELSE
                  dl(i) = 2.0D0 * tau * r2(i)
                END IF
                u = u + (w * r1(i) + dl(i)) * x(i,k)
              END DO
              u = (w + bigm) * maj(k) * b(k) + u/nobs
              v = ABS(u) - al * pf(k)
              IF (v > 0.0D0) THEN
                b(k) = SIGN(v, u)/((w + bigm) * maj(k))
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
            END DO ! END PROXIMAL GRADIENT DESCENT
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
              th(k) = SIGN(v, u) / (maj(k) * bigm)
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
!                   b(k) = SIGN(v, u) / (maj(k) * (w + bigm))
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
!                   th(k) = SIGN(v, u) / (maj(k) * bigm)
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
END SUBROUTINE cpalslassopath
