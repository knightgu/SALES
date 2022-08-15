! ------------------------------------------------------------------------------
! cpalslassoNet.f90: COORDINATE DESCENT ALGORITHM FOR THE COUPLE ALS REGRESSION.
! ------------------------------------------------------------------------------
!
! USAGE:
!
! CALL cpalslassoNET(w, tau, lam2, nobs, nvars, x, y, jd, pfmean, pfscale,&
!      & pf2mean, pf2scale, dfmax, pmax, nlam, flmin, ulam, eps, isd, intr,&
!      & maxit, nalam, b0, beta, ibeta, nbeta, t0, theta, itheta, ntheta, alam,&
!      & npass, jerr)
!
! INPUT ARGUMENTS:
!
!    tau = the asymmetry coefficient in the coupled sparse ALS model.
!    w = weight assigned to the first part of the loss function
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs,nvars) = matrix of covariates, dimension N * p; row -> observation.
!    y(nobs) = response variable.
!    pfmean(nvars) = relative L1 penalties for each mean coefficient
!                pfmean(j) = 0 => jth mean coefficient unpenalized
!    pfscale(nvars) = relative L1 penalties for each scale coefficient
!                pfscale(j) = 0 => jth scale coefficient unpenalized
!    pf2mean(nvars) = relative L2 penalties for each mean coefficient
!                pfmean(j) = 0 => jth mean coefficient unpenalized
!    pf2scale(nvars) = relative L2 penalties for each scale coefficient
!                pfscale(j) = 0 => jth scale coefficient unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero.
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
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original scales.
!    intr = intercept flag:
!           intr = 0 => intercepts are always set to be zero
!           intr = 1 => intercepts are calculated
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda
!    value.
!            (suggested values, maxit = 100000)
!
! OUTPUT:
!
!    nalam = actual number of lambda values (solutions)
!    beta(pmax,nalam) = compressed mean coefficient values for each solution
!    ibeta(pmax) = pointers to compressed mean coefficients
!    nbeta(nalam) = number of compressed mean coefficients for each solution
!    theta(pmax,nalam) = compressed scale coefficient values for each solution
!    itheta(pmax) = pointers to compressed scale coefficients
!    ntheta(nalam) = number of compressed scale coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                jerr < 7777 => memory allocation error
!                jerr = 7777 => all used predictors have zero variance
!                jerr = 10000 => maxval(vp) <= 0.0
!           jerr < 0 => non fatal error - partial output:
!                Solutions for larger lambdas (1:(k-1)) returned.
!                jerr = -k => convergence for kth lambda value not reached
!                  after maxit (see above) iterations.
!                jerr = -10000-k => number of nonzero coefficients along path
!                  exceeds pmax (see above) at kth lambda value.
!
! LICENSE: GNU GPL (version 2 or later)
!
! AUTHORS:
!
! YUWEN GU (yuwen.gu@uconn.edu), HUI ZOU (zouxx019@umn.edu)
!   DEPARTMENT OF STATISTICS, UNIVERSITY OF CONNECTICUT
!   SCHOOL OF STATISTICS, UNIVERSITY OF MINNESOTA
!
! REFERENCES:
!
! Gu, Y., and Zou, H. (2016).
! High-dimensional generalizations of asymmetric least squares regression and
! their applications. The Annals of Statistics, 44(6), 2661–2694.
!
!----------------------------------------------------------------------------- !
SUBROUTINE cpalslassoNET(w, tau, lam2, nobs, nvars, x, y, jd, pfmean, pfscale,&
     & pf2mean, pf2scale, dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit&
     &, nalam, b0, beta, ibeta, nbeta, t0, theta, itheta, ntheta, alam, npass,&
     & jerr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs, nvars, dfmax, pmax, jd(*), nlam, nalam
  INTEGER :: isd, npass, maxit, jerr, intr
  INTEGER :: ibeta(pmax), nbeta(nlam), itheta(pmax), ntheta(nlam)
  DOUBLE PRECISION :: w, tau, lam2, flmin, eps, ulam(nlam), alam(nlam)
  DOUBLE PRECISION :: x(nobs, nvars), y(nobs)
  DOUBLE PRECISION :: pfmean(nvars), pfscale(nvars)
  DOUBLE PRECISION :: pf2mean(nvars), pf2scale(nvars)
  DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam), theta(pmax, nlam), t0(nlam)
  ! --------- LOCAL DECLARATIONS ---------- !
  INTEGER :: j, l, nk, nkk, ierr
  INTEGER,  DIMENSION(:),  ALLOCATABLE :: ju
  DOUBLE PRECISION,  DIMENSION(:),  ALLOCATABLE :: xmean, xnorm, maj
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
  IF (MAXVAL(pfmean) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  IF (MAXVAL(pfscale) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  IF (MAXVAL(pf2mean) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  IF (MAXVAL(pf2scale) <= 0.0D0) THEN
    jerr = 10000
    RETURN
  END IF
  pfmean = MAX(0.0D0, pfmean)
  pfscale = MAX(0.0D0, pfscale)
  pf2mean = MAX(0.0D0, pf2mean)
  pf2scale = MAX(0.0D0, pf2scale)
  CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
  ! ------------------------ CALL cpalslassopath -------------------------- !
  CALL cpalslassoNETpath(w, tau, lam2, maj, nobs, nvars, x, y, ju, pfmean,&
       & pfscale, pf2mean, pf2scale ,dfmax, pmax, nlam, flmin, ulam, eps, maxit&
       &, nalam, b0, beta, ibeta, nbeta, t0, theta, itheta, ntheta, alam, npass&
       &, jerr, intr)
  IF (jerr > 0) RETURN ! CHECK ERROR AFTER CALLING FUNCTION
  ! ----------- TRANSFORM BETA BACK TO THE ORIGINAL SCALE ----------- !
  DO l = 1, nalam
    nk = nbeta(l)
    nkk = ntheta(l)
    IF (isd == 1) THEN
      DO j = 1, nk
        beta(j, l) = beta(j, l) / xnorm(ibeta(j))
      END DO
      DO j = 1, nkk
        theta(j, l) = theta(j, l) / xnorm(itheta(j))
      END DO
    END IF
    b0(l) = b0(l) - DOT_PRODUCT(beta(1:nk, l), xmean(ibeta(1:nk)))
    t0(l) = t0(l) - DOT_PRODUCT(theta(1:nkk, l), xmean(itheta(1:nkk)))
  END DO
  DEALLOCATE(ju, xmean, xnorm, maj)
  RETURN
END SUBROUTINE cpalslassoNET

! ----------------------------- cpalslassoNETpath --------------------------- !
SUBROUTINE cpalslassoNETpath(w, tau, lam2, maj, nobs, nvars, x, y, ju, pfmean,&
     & pfscale, pf2mean, pf2scale ,dfmax, pmax, nlam, flmin, ulam, eps, maxit,&
     & nalam, b0, beta, ibeta, nbeta, t0, theta, itheta, ntheta, alam, npass,&
     & jerr, intr)

  IMPLICIT NONE
  ! --------------- INPUT VARIABLES --------------- !
  INTEGER :: nobs, nvars, dfmax, pmax, nlam, nalam, maxit, npass, jerr, intr
  INTEGER :: ju(nvars), ibeta(pmax), nbeta(nlam), itheta(pmax), ntheta(nlam)
  DOUBLE PRECISION :: eps, w, tau, lam2, flmin
  DOUBLE PRECISION :: x(nobs, nvars), y(nobs), maj(nvars)
  DOUBLE PRECISION :: pfmean(nvars), pfscale(nvars), ulam(nlam), alam(nlam)
  DOUBLE PRECISION :: pf2mean(nvars), pf2scale(nvars)
  DOUBLE PRECISION :: beta(pmax, nlam), b0(nlam), theta(pmax, nlam), t0(nlam)
  ! ------------- LOCAL DECLARATIONS -------------- !
  INTEGER,  PARAMETER :: mnlam = 6
  INTEGER :: mnl, i, k, j, l, vrg, ctr, ierr, nib, nith, me
  INTEGER,  DIMENSION (:),  ALLOCATABLE :: mmb, mmth
  DOUBLE PRECISION,  PARAMETER :: big = 9.9D30, mfl = 1.0D-06
  DOUBLE PRECISION :: tmp, bigm, d, dif, oldb, oldth, al, dl(nobs)
  DOUBLE PRECISION :: alf = 1.0D0, u = 0.0D0, v
  DOUBLE PRECISION,  DIMENSION(:),  ALLOCATABLE :: b, th, oldbeta, oldtheta
  DOUBLE PRECISION,  DIMENSION(:),  ALLOCATABLE :: r1, r2
  ! ----------------- ALLOCATE VARIABLES ------------------ !
  ALLOCATE (b(0:nvars), STAT = jerr)
  ALLOCATE (oldbeta(0:nvars), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (th(0:nvars), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (oldtheta(0:nvars), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (mmb(1:nvars), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (mmth(1:nvars), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (r1(1:nobs), STAT = ierr)
  jerr = jerr + ierr
  ALLOCATE (r2(1:nobs), STAT = ierr)
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
    alf = flmin ** (1.0D0 / (nlam - 1.0D0))
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
            IF (pfmean(j) > 0.0D0) THEN
              u = DOT_PRODUCT(dl, x(:, j))
              al = MAX(al, ABS(w * DOT_PRODUCT(r1, x(:, j)) + u) / pfmean(j))
            END IF
            IF (pfscale(j) > 0.0D0) THEN
              IF (pfmean(j) <= 0.0D0) u = DOT_PRODUCT(dl, x(:, j))
              al = MAX(al, ABS(u) / pfscale(j))
            END IF
          END IF
        END DO
        al = al * alf / nobs
      END IF
    END IF
    ctr = 0
    ! ------------------ OUTER LOOP -------------------- !
    DO
      IF (intr == 1) THEN
        oldbeta(0) = b(0)
        oldtheta(0) = th(0)
      END IF
      IF (nib > 0) oldbeta(ibeta(1:nib)) = b(ibeta(1:nib))
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
                u = u + (w * r1(i) + dl(i)) * x(i, k)
              END DO
              u = (w + bigm) * maj(k) * b(k) + u / nobs
              v = ABS(u) - al * pfmean(k)
              IF (v > 0.0D0) THEN
                tmp = SIGN(v, u) / ((w + bigm) * maj(k) + lam2 * pf2mean(k))
              ELSE
                tmp = 0.0D0
              END IF
              d = tmp - b(k)
              IF (bigm * d**2 < eps) EXIT
              b(k) = tmp
              r1 = r1 - x(:, k) * d
              r2 = r2 - x(:, k) * d
            END DO ! END PROXIMAL GRADIENT DESCENT
            d = b(k) - oldb
            IF (ABS(d) > 0.0D0) THEN
              dif = MAX(dif, bigm * d**2)
              IF (mmb(k) == 0) THEN
                nib = nib + 1
                IF (nib > pmax) EXIT
                mmb(k) = nib
                ibeta(nib) = k ! RECORD ACTIVE VARIABLES
              END IF
            END IF
          END IF
        END DO
        DO k = 1, nvars
          IF (ju(k) /= 0) THEN
            oldth = th(k)
            DO ! BEGIN PROXIMAL GRADIENT DESCENT
              u = 0.0D0
              DO i = 1, nobs
                IF (r2(i) < 0.0D0) THEN
                  dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
                ELSE
                  dl(i) = 2.0D0 * tau * r2(i)
                END IF
                u = u + dl(i) * x(i,k)
              END DO
              u = bigm * maj(k) * th(k) + u/nobs
              v = ABS(u) - al * pfscale(k)
              IF (v > 0.0D0) THEN
                tmp = SIGN(v, u) / (bigm * maj(k) + lam2 * pf2scale(k))
              ELSE
                tmp = 0.0D0
              END IF
              d = tmp - th(k)
              IF (bigm * d**2 < eps) EXIT
              th(k) = tmp
              r2 = r2 - x(:, k) * d
            END DO ! END PROXIMAL GRADIENT DESCENT
            d = th(k) - oldth
            IF (ABS(d) > 0.0D0) THEN
              dif = MAX(dif, bigm * d**2)
              IF (mmth(k) == 0) THEN
                nith = nith + 1
                IF (nith > pmax) EXIT
                mmth(k) = nith
                itheta(nith) = k ! RECORD ACTIVE VARIABLES
              END IF
            END IF
          END IF
        END DO
        IF (nib > pmax) EXIT
        IF (nith > pmax) EXIT
        IF (intr == 1) THEN
          oldb = b(0)
          oldth = th(0)
          DO ! BEGIN NEWTON-RAPHSON
            DO i = 1, nobs
              IF (r2(i) < 0.0D0) THEN
                dl(i) = (1.0D0 - tau) * r2(i)
              ELSE
                dl(i) = tau * r2(i)
              END IF
            END DO
            d = SUM(r1) / nobs
            v = SUM(dl) / (nobs * tau + (1.0D0 - 2.0D0 * tau) * COUNT(r2 <&
                 & 0.0D0)) - d
            IF (bigm * d**2 < eps .AND. bigm * v**2 < eps) EXIT
            b(0) = b(0) + d
            th(0) = th(0) + v
            r1 = r1 - d
            r2 = r2 - d - v
          END DO ! END NEWTON-RAPHSON
          d = b(0) - oldb
          IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
          d = th(0) - oldth
          IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
        END IF
        IF (dif < eps) EXIT
        ! ! ---------- INNER LOOP (ACTIVE SET ACCELERATION) ---------- !
        ! DO
        !   npass = npass + 1
        !   dif = 0.0D0
        !   DO j = 1, nib
        !     k = ibeta(j)
        !     oldb = b(k)
        !     DO ! BEGIN PROXIMAL GRADIENT DESCENT
        !       u = 0.0D0
        !       DO i = 1, nobs
        !         IF (r2(i) < 0.0D0) THEN
        !           dl(i) = 2.0D0 * (1.0D0 - tau) * r2(i)
        !         ELSE
        !           dl(i) = 2.0D0 * tau * r2(i)
        !         END IF
        !         u = u + (w * r1(i) + dl(i)) * x(i,k)
        !       END DO
        !       u = (w + bigm) * maj(k) * b(k) + u/nobs
        !       v = ABS(u) - al * pfmean(k)
        !       IF (v > 0.0D0) THEN
        !         tmp = SIGN(v, u) / ((w + bigm) * maj(k) + lam2 * pf2mean(k))
        !       ELSE
        !         tmp = 0.0D0
        !       END IF
        !       d = tmp - b(k)
        !       IF (bigm * d**2 < eps) EXIT
        !       b(k) = tmp
        !       r1 = r1 - x(:,k) * d
        !       r2 = r2 - x(:,k) * d
        !     END DO ! END PROXIMAL GRADIENT DESCENT
        !     d = b(k) - oldb
        !     IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
        !   END DO
        !   DO j = 1, nith
        !     k = itheta(j)
        !     oldth = th(k)
        !     u = 0.0D0
        !     DO ! BEGIN PROXIMAL GRADIENT DESCENT
        !       DO i = 1, nobs
        !         IF (r2(i) < 0.0D0) THEN
        !           dl (i) = 2.0D0 * (1.0D0 - tau) * r2(i)
        !         ELSE
        !           dl (i) = 2.0D0 * tau * r2(i)
        !         END IF
        !         u = u + dl(i) * x(i,k)
        !       END DO
        !       u = maj(k) * th(k) * bigm + u / nobs
        !       v = ABS(u) - al * pfscale(k)
        !       IF (v > 0.0D0) THEN
        !         tmp = SIGN(v, u) / (bigm * maj(k) + lam2 * pf2scale(k))
        !       ELSE
        !         tmp = 0.0D0
        !       END IF
        !       d = tmp - th(k)
        !       IF (bigm * d**2 < eps) EXIT
        !       th(k) = tmp
        !       r2 = r2 - x(:,k) * d
        !     END DO ! END PROXIMAL GRADIENT DESCENT
        !     d = th(k) - oldth
        !     IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
        !   END DO
        !   IF (intr == 1) THEN
        !     oldb = b(0)
        !     oldth = th(0)
        !     DO ! BEGIN NEWTON-RAPHSON
        !       DO i = 1, nobs
        !         IF (r2(i) < 0.0D0) THEN
        !           dl(i) = (1.0D0 - tau) * r2(i)
        !         ELSE
        !           dl(i) = tau * r2(i)
        !         END IF
        !       END DO
        !       d = SUM(r1) / nobs
        !       v = SUM(dl) / (nobs * tau + (1.0D0 - 2.0D0 * tau) *&
        !            & COUNT(r2 < 0.0D0)) - d
        !       IF (bigm * d**2 < eps .AND. bigm * v**2 < eps) EXIT
        !       b(0) = b(0) + d
        !       th(0) = th(0) + v
        !       r1 = r1 - d
        !       r2 = r2 - d - v
        !     END DO ! END NEWTON-RAPHSON
        !     d = b(0) - oldb
        !     IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
        !     d = th(0) - oldth
        !     IF (ABS(d) > 0.0D0) dif = MAX(dif, bigm * d**2)
        !   END IF
        !   IF (dif < eps) EXIT
        ! END DO ! -----> END INNER LOOP
      END DO ! ----> END MIDDLE LOOP
      IF (nib > pmax) EXIT
      IF (nith > pmax) EXIT
      ! -------------- FINAL CHECK ---------------- !
      vrg = 1
      IF (intr == 1) THEN
        IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
        IF ((th(0) - oldtheta(0))**2 >= eps) vrg = 0
      END IF
      DO j = 1, nib
        IF ((b(ibeta(j)) - oldbeta(ibeta(j)))**2 >= eps) THEN
          vrg = 0
          EXIT
        END IF
      END DO
      DO j = 1, nith
        IF ((th(itheta(j)) - oldtheta(itheta(j)))**2 >= eps) THEN
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
    END DO ! -------> END OUTER LOOP
    ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
    IF (nib > pmax) THEN
      jerr = - 10000 - l
      EXIT
    END IF
    IF (nith > pmax) THEN
      jerr = - 10000 - l
      EXIT
    END IF
    IF (nib > 0) beta(1:nib, l) = b(ibeta(1:nib))
    b0(l) = b(0)
    nbeta(l) = nib
    IF (nith > 0) theta(1:nith, l) = th(itheta(1:nith))
    t0(l) = th(0)
    ntheta(l) = nith
    alam(l) = al
    nalam = l
    IF (l < mnl) CYCLE
    IF (flmin >= 1.0D0) CYCLE
    me = COUNT(ABS(beta(1:nib, l)) > 0.0D0)
    IF (me > dfmax) EXIT
    me = COUNT(ABS(theta(1:nith, l)) > 0.0D0)
    IF (me > dfmax) EXIT
  END DO ! ----------> END LAMBDA LOOP
  DEALLOCATE(b, oldbeta, th, oldtheta, r1, r2, mmb, mmth)
  RETURN
END SUBROUTINE cpalslassoNETpath
