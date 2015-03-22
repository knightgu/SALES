! --------------------------------------------------------------------------
! alslassoNET.f90: the coordinate descent algorithm for the ALS regression.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! call alslassoNET(tau, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, &
! & pmax, nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, &
! & nbeta, alam, npass, jerr)
! 
! INPUT ARGUMENTS:
!
!    tau = the parameter in the ALS regression
!    lam2 = regularization parameter for the L2 penalty
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; 
!                     each row is an observation
!    y(nobs) = response variable
!    jd(jd(1) + 1) = predictor variable deletion flag
!                    jd(1) = 0  => use all variables
!                    jd(1) != 0 => do not use variables jd(2), ..., jd(jd(1) + 1)
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                 pf2(j) = 0 => jth variable unpenalized
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
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nalam) = intercept values for each solution
!    beta(pmax,nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr = 0 => no error
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
!    Yuwen Gu (guxxx192@umn.edu), Yi Yang (yiyang@umn.edu), and Hui Zou (hzou@stat.umn.edu), 
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H. (2012). An Efficient Algorithm for Computing The HHSVM and Its 
!      Generalizations. Journal of Computational and Graphical Statistics, 22, 396-415.

! ------------------------------------------------------------------------------------ !
SUBROUTINE alslassoNET(tau, lam2, nobs, nvars, x, y, jd, pf, pf2, dfmax, pmax, &
& nlam, flmin, ulam, eps, isd, maxit, nalam, b0, beta, ibeta, nbeta, &
& alam, npass, jerr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: nobs,nvars,dfmax,pmax,nlam,nalam,isd,npass,jerr,maxit
  INTEGER :: jd(*)
  INTEGER :: ibeta(pmax),nbeta(nlam)
  DOUBLE PRECISION :: lam2,flmin,eps,tau
  DOUBLE PRECISION :: x(nobs,nvars),y(nobs)
  DOUBLE PRECISION :: pf(nvars),pf2(nvars)
  DOUBLE PRECISION :: beta(pmax,nlam),b0(nlam)
  DOUBLE PRECISION :: ulam(nlam),alam(nlam)
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER :: j,l,nk,ierr
  INTEGER, DIMENSION(:), ALLOCATABLE :: ju
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xmean,xnorm,maj
  ! -------- ALLOCATE VARIABLES -------- !
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
  ! -------------------- CALL alslassoNET --------------------- !
  CALL alslassoNETpath(tau, lam2, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
  & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, &
  & nbeta, alam, npass, jerr)
  IF (jerr > 0) RETURN ! CHECK ERROR AFTER CALLING FUNCTION
  ! ----------- TRANSFORM BETA BACK TO THE ORIGINAL SCALE ----------- !
  DO l = 1, nalam
    nk = nbeta(l)
    IF (isd == 1) THEN
      DO j = 1, nk
        beta(j,l) = beta(j,l)/xnorm(ibeta(j))
      END DO
    END IF
    b0(l) = b0(l) - DOT_PRODUCT(beta(1:nk,l),xmean(ibeta(1:nk)))
  END DO
  DEALLOCATE(ju,xmean,xnorm,maj)
  RETURN
END SUBROUTINE alslassoNET


! --------------------------------- alslassoNETpath --------------------------------- !
SUBROUTINE alslassoNETpath(tau, lam2, maj, nobs, nvars, x, y, ju, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, nbeta, alam, &
& npass, jerr)

  IMPLICIT NONE
  ! -------- INPUT VARIABLES -------- !
  INTEGER :: mnl,nobs,nvars,dfmax,pmax,nlam,maxit,nalam,npass,jerr
  INTEGER :: ju(nvars),m(pmax),nbeta(nlam)
  DOUBLE PRECISION :: lam2,eps,tau
  DOUBLE PRECISION :: x(nobs,nvars),y(nobs),maj(nvars)
  DOUBLE PRECISION :: pf(nvars),pf2(nvars)
  DOUBLE PRECISION :: beta(pmax,nlam),b0(nlam)
  DOUBLE PRECISION :: ulam(nlam),alam(nlam)
  ! -------- LOCAL DECLARATIONS -------- !
  INTEGER, PARAMETER :: mnlam = 6
  DOUBLE PRECISION, PARAMETER :: big = 9.9D30, mfl = 1.0D-06
  DOUBLE PRECISION :: tmp,bigm,d,dif,oldb,u,v,al,alf,flmin,dl(nobs)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: b,oldbeta,r
  INTEGER :: i,k,j,l,vrg,ctr,ierr,ni,me
  INTEGER, DIMENSION(:), ALLOCATABLE :: mm
  ! -------- ALLOCATE VARIABLES -------- !
  ALLOCATE(b(0:nvars), STAT=jerr)
  ALLOCATE(oldbeta(0:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(mm(1:nvars), STAT=ierr)
  jerr = jerr + ierr
  ALLOCATE(r(1:nobs), STAT=ierr)
  jerr = jerr + ierr
  IF (jerr /= 0) RETURN
  ! ---------------- INITIALIZATION ---------------- !
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  m = 0
  mm = 0
  npass = 0
  ni = npass
  mnl = MIN(mnlam,nlam)
  bigm = 2.0D0 * MAX((1.0D0-tau),tau)
  maj = bigm * maj
  IF (flmin < 1.0D0) THEN
    flmin = MAX(mfl,flmin)
    alf = flmin ** (1.0D0/(nlam-1.0D0))
  END IF
  ! ----------------- LAMBDA LOOP (OUTMOST LOOP) ------------------- !
  DO l = 1, nlam
  ! ----------------- COMPUTE LAMBDA ------------------- !
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
          IF (r(i) < 0.0D0) THEN
            dl(i) = 2.0D0 * (1.0D0 - tau) * r(i)
          ELSE
            dl(i) = 2.0D0 * tau * r(i)
          END IF
        END DO
        DO j = 1, nvars
          IF (ju(j) /= 0) THEN
            IF (pf(j) > 0.0D0) THEN
              u = DOT_PRODUCT(dl,x(:,j))
              al = MAX(al, ABS(u)/pf(j))
            END IF
          END IF
        END DO
        al = al * alf/nobs
      END IF
    END IF
    ctr = 0
    ! ------------------ OUTER LOOP -------------------- !
    DO
      oldbeta(0) = b(0)
      IF (ni > 0) oldbeta(m(1:ni)) = b(m(1:ni))
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
                IF (r(i) < 0.0D0) THEN
                  dl(i) = 2.0D0 * (1.0D0 - tau) * r(i)
                ELSE
                  dl(i) = 2.0D0 * tau * r(i)
                END IF
                u = u + dl(i) * x(i,k)
              END DO
              u = maj(k) * b(k) + u/nobs
              v = al * pf(k)
              v = ABS(u) - v
              IF (v > 0.0D0) THEN
                tmp = SIGN(v,u)/(maj(k) + pf2(k) * lam2)
              ELSE
                tmp = 0.0D0
              END IF
              d = tmp - b(k)
              IF (bigm * d**2 < eps) EXIT
              b(k) = tmp
              r = r - x(:,k) * d
            END DO ! END PROXIMAL GRADIENT DESCENT
            d = b(k) - oldb
            IF (ABS(d) > 0.0D0) THEN
              dif = MAX(dif, bigm * d**2)
              IF (mm(k) == 0) THEN
                ni = ni + 1
                IF (ni > pmax) EXIT
                mm(k) = ni
                m(ni) = k ! RECORD ACTIVE VARIABLES
              END IF
            END IF
          END IF
        END DO
        IF (ni > pmax) EXIT
        oldb = b(0)
        DO ! BEGIN GRADIENT DESCENT
          DO i = 1, nobs
            IF (r(i) < 0.0D0) THEN
              dl(i) = 2.0D0 * (1.0D0 - tau) * r(i)
            ELSE
              dl(i) = 2.0D0 * tau * r(i)
            END IF
          END DO
          d = SUM(dl)/(nobs * bigm)
          IF (bigm * d**2 < eps) EXIT
          b(0) = b(0) + d
          r = r - d
        END DO ! END GRADIENT DESCENT
        d = b(0) - oldb
        IF (ABS(d) > 0.0D0) THEN
          dif = MAX(dif, bigm * d**2)
        END IF
        IF (dif < eps) EXIT
        ! ----------------- INNER LOOP ------------------- !
        DO
          npass = npass + 1
          dif = 0.0D0
          DO j = 1, ni
            k = m(j)
            oldb = b(k)
            DO ! BEGIN PROXIMAL GRADIENT DESCENT
              u = 0.0D0
              DO i = 1, nobs
                IF (r(i) < 0.0D0) THEN
                  dl(i) = 2.0D0 * (1.0D0 - tau) * r(i)
                ELSE
                  dl(i) = 2.0D0 * tau * r(i)
                END IF
                u = u + dl(i) * x(i,k)
              END DO
              u = maj(k) * b(k) + u/nobs
              v = al * pf(k)
              v = ABS(u) - v
              IF (v > 0.0D0) THEN
                tmp = SIGN(v,u)/(maj(k) + pf2(k) * lam2)
              ELSE
                tmp = 0.0D0
              END IF
              d = tmp - b(k)
              IF (bigm * d**2 < eps) EXIT
              b(k) = tmp
              r = r - x(:,k) * d
            END DO ! END PROXIMAL GRADIENT DESCENT
            d = b(k) - oldb
            IF (ABS(d) > 0.0D0) THEN
              dif = MAX(dif, bigm * d**2)
            END IF
          END DO
          oldb = b(0)
          DO ! BEGIN GRADIENT DESCENT     
            DO i = 1, nobs
              IF (r(i) < 0.0D0) THEN
                dl(i) = 2.0D0 * (1.0D0 - tau) * r(i)
              ELSE
                dl(i) = 2.0D0 * tau * r(i)
              END IF
            END DO
            d = SUM(dl)/(nobs * bigm)
            IF (bigm * d**2 < eps) EXIT
            b(0) = b(0) + d
            r = r - d
          END DO ! END GRADIENT DESCENT
          d = b(0) - oldb
          IF (ABS(d) > 0.0D0) THEN
              dif = MAX(dif, bigm * d**2)
          END IF
          IF (dif < eps) EXIT
        END DO ! ----------> END INNER LOOP
      END DO ! ----------> END MIDDLE LOOP
      IF (ni > pmax) EXIT
      ! -------------- FINAL CHECK ---------------- !
      vrg = 1
      IF ((b(0) - oldbeta(0))**2 >= eps) vrg = 0
      DO j = 1, ni
        IF ((b(m(j)) - oldbeta(m(j)))**2 >= eps) THEN
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
    IF (ni > pmax) THEN
      jerr = - 10000 - l
      EXIT
    END IF
    IF (ni > 0) beta(1:ni,l) = b(m(1:ni))
    nbeta(l) = ni
    b0(l) = b(0)
    alam(l) = al
    nalam = l
    IF (l < mnl) CYCLE
    IF (flmin >= 1.0D0) CYCLE
    me = COUNT(ABS(beta(1:ni,l)) > 0.0D0)
    IF (me > dfmax) EXIT
  END DO ! -------> END LAMBDA LOOP
  DEALLOCATE(b,oldbeta,r,mm)
  RETURN
END SUBROUTINE alslassoNETpath