SUBROUTINE TScomp(W, V, nx, ny, mx, my, GL, nGL, BDs, TS)
IMPLICIT NONE
! W: (nx, mx) matrix
! V: (ny, my) matrix
! Wdiff: (nx, mWdiff) matrix of differences
! Vdiff: (ny, mVdiff) matrix of differences
! GH: (nGH, 2) matrix with abscissa in the 1st col, weight in the 2nd col for Gauss-Hermite quadrature
! GL: (nGl, 2) matrix with abscissa in the 1st col, weight in the 2nd col for Gauss-Legendre quadrature
! VARs: nVARs vector with variances for normal weight function
! BDs: nBDs vector with bounds for uniform weight function
! TS: nVARs + nBDs vector of results
INTEGER :: nx, ny, mx, my, nGL, N_x, N_y
real(8) :: W(nx, mx), V(ny, my)
real(8) :: GL(nGL, 2), BDs(2)
real(8) :: avgx(nx), avgy(ny), TS(2)

INTEGER :: i, j, k
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: Wcos, Vcos, Wsin, Vsin, denomx, denomy, x, wei, TStmp

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

!! This is a little bit faster than below	
avgx(:) = SUM(W, DIM = 2)/DBLE(mx)
avgy(:) = SUM(V, DIM = 2)/DBLE(my)

!	DO i = 1, nx
!		avgx(i) = SUM(W(i, :))/DBLE(mx)
!	END DO

!	DO i = 1, ny
!		avgy(i) = SUM(V(i, :))/DBLE(my)
!	END DO

k = 0
DO i = 1, (mx-1)
	DO j = (i+1), mx
		k = k + 1
		Wdiff(:, k) = W(:, j) - W(:, i)
	END DO
END DO

k = 0
DO i = 1, (my-1)
	DO j = (i+1), my
		k = k + 1
		Vdiff(:, k) = V(:, j) - V(:, i)
	END DO
END DO


TS(1) = 0
TS(2) = 0
DO j = 1, nGL
	x = (BDs(1)+BDs(2))/DBLE(2) + (BDs(2)-BDs(1))*GL(j,1)/DBLE(2)
	wei = GL(j,2)
	denomx = sqrt(abs(SUM(cos(x*Wdiff/DBLE(mx)))/DBLE(N_x)))**(DBLE(mx))
	denomy = sqrt(abs(SUM(cos(x*Vdiff/DBLE(my)))/DBLE(N_y)))**(DBLE(my))
	Wcos = (SUM(cos(x*avgx(:)))/DBLE(nx))/denomx
	Vcos = (SUM(cos(x*avgy(:)))/DBLE(ny))/denomy
	Wsin = (SUM(sin(x*avgx(:)))/DBLE(nx))/denomx
	Vsin = (SUM(sin(x*avgy(:)))/DBLE(ny))/denomy

	TStmp = wei*((Wcos - Vcos)**2 + (Wsin - Vsin)**2)
	TS(1) = TS(1) + TStmp
	TS(2) = TS(2) + TStmp*EXP(-x**2/DBLE(2))
END DO
TS(1) = (BDs(2)-BDs(1))*DBLE(nx)*TS(1)/DBLE(2)
TS(2) = (BDs(2)-BDs(1))*DBLE(nx)*TS(2)/DBLE(2)


RETURN
END SUBROUTINE TScomp


SUBROUTINE AMISE(W, V, nx, ny, mx, my, GL, nGL, Wseq, Vseq, nseq,&
				Wobj, Vobj, varx, vary, Wdiffvar, Vdiffvar)
IMPLICIT NONE
! W: (nx, mx) matrix
! V: (ny, my) matrix
! Wdiff: (nx, mWdiff) matrix of differences
! Vdiff: (ny, mVdiff) matrix of differences
! GL: (nGl, 2) matrix with abscissa in the 1st col, weight in the 2nd col for Gauss-Legendre quadrature
! Wobj: evaluated values corresponding to Wseq (for dataset W)
! Vobj: evaluated values corresponding to Wseq (for dataset V)
INTEGER :: nx, ny, mx, my, nGL, N_x, N_y, nseq
real(8) :: W(nx, mx), V(ny, my), GL(nGL, 2), Wseq(nseq), Vseq(nseq)
real(8) :: Wobj(nseq), Vobj(nseq)

INTEGER :: i, j, k
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: denomx, denomy, x, wei, Wvar, Vvar, uxvar, uyvar
real(8) :: varx, vary, Wdiffvar, Vdiffvar
REAL(8), PARAMETER :: Pi = 3.141593

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

k = 0
DO i = 1, (mx-1)
	DO j = (i+1), mx
		k = k + 1
		Wdiff(:, k) = W(:, j) - W(:, i)
	END DO
END DO

k = 0
DO i = 1, (my-1)
	DO j = (i+1), my
		k = k + 1
		Vdiff(:, k) = V(:, j) - V(:, i)
	END DO
END DO

Wvar = (SUM(W**2) - SUM(W)**2/DBLE(nx*mx))/DBLE(nx*mx-1)
Vvar = (SUM(V**2) - SUM(V)**2/DBLE(ny*my))/DBLE(ny*my-1)
uxvar = SUM((SUM(W**2, DIM = 2) - (SUM(W, DIM = 2)**2/DBLE(mx)))/DBLE(mx-1))/DBLE(nx)
uyvar = SUM((SUM(V**2, DIM = 2) - (SUM(V, DIM = 2)**2/DBLE(my)))/DBLE(my-1))/DBLE(ny)

varx = Wvar - uxvar
vary = Vvar - uyvar

Wdiffvar = (SUM(Wdiff**2) - SUM(Wdiff)**2/DBLE(N_x))/DBLE(N_x-1)
Vdiffvar = (SUM(Vdiff**2) - SUM(Vdiff)**2/DBLE(N_y))/DBLE(N_y-1)


DO i = 1, nseq
	Wobj(i) = 0
	Vobj(i) = 0
	DO j = 1, nGL/2
		x = GL(j, 1)
		wei = GL(j, 2)
		denomx = sqrt(abs(SUM(cos(x*Wdiff/(DBLE(mx)*Wseq(i))))/DBLE(N_x)))**(DBLE(mx))
		denomy = sqrt(abs(SUM(cos(x*Vdiff/(DBLE(my)*Vseq(i))))/DBLE(N_y)))**(DBLE(my))
		Wobj(i) = Wobj(i) + wei*(1 - (1 - x**2)**3/denomx)**2/x**2
		Vobj(i) = Vobj(i) + wei*(1 - (1 - x**2)**3/denomy)**2/x**2			
	END DO
	Wobj(i) = Wseq(i)*Wobj(i)/(Pi*nx) + 9*Wseq(i)**4/(4*sqrt(Pi)*varx**(1.5))
	Vobj(i) = Vseq(i)*Vobj(i)/(Pi*nx) + 9*Vseq(i)**4/(4*sqrt(Pi)*vary**(1.5))
END DO


RETURN
END SUBROUTINE AMISE


SUBROUTINE Fhat(W, V, nx, ny, mx, my, GL, nGL, zseq, nseq, bwx, bwy, Fx, Fy)
IMPLICIT NONE
! W: (nx, mx) matrix
! V: (ny, my) matrix
! Wdiff: (nx, mWdiff) matrix of differences
! Vdiff: (ny, mVdiff) matrix of differences
! GL: (nGl, 2) matrix with abscissa in the 1st col, weight in the 2nd col for Gauss-Legendre quadrature
! Wobj: evaluated values corresponding to Wseq (for dataset W)
! Vobj: evaluated values corresponding to Wseq (for dataset V)
! Fx: evaluated true signal distribution funciton estimates for dataset W at zseq
! Fy: evaluated true signal distribution funciton estimates for dataset W at zseq
INTEGER :: nx, ny, mx, my, nGL, N_x, N_y, nseq
real(8) :: W(nx, mx), V(ny, my), GL(nGL, 2), zseq(nseq), bwx, bwy
real(8) :: Fx(nseq), Fy(nseq)

INTEGER :: i, j, k
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: denomx, denomy, x, wei
REAL(8), PARAMETER :: Pi = 3.141593

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

k = 0
DO i = 1, (mx-1)
	DO j = (i+1), mx
		k = k + 1
		Wdiff(:, k) = W(:, j) - W(:, i)
	END DO
END DO

k = 0
DO i = 1, (my-1)
	DO j = (i+1), my
		k = k + 1
		Vdiff(:, k) = V(:, j) - V(:, i)
	END DO
END DO

DO i = 1, nseq
	Fx(i) = 0
	Fy(i) = 0
	DO j = 1, nGL/2
		x = GL(j, 1)
		wei = GL(j, 2)
		denomx = sqrt(abs(SUM(cos(x*Wdiff/(DBLE(mx)*bwx)))/DBLE(N_x)))**(DBLE(mx))
		denomy = sqrt(abs(SUM(cos(x*Vdiff/(DBLE(my)*bwy)))/DBLE(N_y)))**(DBLE(my))
		Fx(i) = Fx(i) + wei*(1-x**2)**3*SUM(sin(x*(zseq(i) - (SUM(W, DIM = 2)/&
			DBLE(mx)))/bwx))/(denomx*x)
		Fy(i) = Fy(i) + wei*(1-x**2)**3*SUM(sin(x*(zseq(i) - (SUM(V, DIM = 2)/&
			DBLE(my)))/bwy))/(denomy*x)
	END DO
	Fx(i) = Fx(i)/(Pi*nx) + 0.5
	Fy(i) = Fy(i)/(Pi*ny) + 0.5
	IF (Fx(i) > 1) THEN 
		Fx(i) = 1 
	END IF
	IF (Fx(i) < 0) THEN 
		Fx(i) = 0 
	END IF
	IF (Fy(i) > 1) THEN 
		Fy(i) = 1 
	END IF
	IF (Fy(i) < 0) THEN 
		Fy(i) = 0 
	END IF		
END DO	



RETURN
END SUBROUTINE Fhat


SUBROUTINE Uhat(W, V, nx, ny, mx, my, GL, nGL, errseq, nseq, Ux, Uy)
IMPLICIT NONE
! W: (nx, mx) matrix
! V: (ny, my) matrix
! Wdiff: (nx, mWdiff) matrix of differences
! Vdiff: (ny, mVdiff) matrix of differences
! GL: (nGl, 2) matrix with abscissa in the 1st col, weight in the 2nd col for Gauss-Legendre quadrature
! Wobj: evaluated values corresponding to Wseq (for dataset W)
! Vobj: evaluated values corresponding to Wseq (for dataset V)
! Ux: evaluated error distribution funciton estimates for dataset W at zseq
! Uy: evaluated error distribution funciton estimates for dataset W at zseq
! Wdiffvar: variance of N_x differences of W
! Vdiffvar: variance of N_y differences of V
INTEGER :: nx, ny, mx, my, nGL, N_x, N_y, nseq
real(8) :: W(nx, mx), V(ny, my), GL(nGL, 2), errseq(nseq)
real(8) :: Ux(nseq), Uy(nseq)

INTEGER :: i, j, k
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: hoptW, hoptV, denomx, denomy, x, wei
REAL(8), PARAMETER :: Pi = 3.141593

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

k = 0
DO i = 1, (mx-1)
	DO j = (i+1), mx
		k = k + 1
		Wdiff(:, k) = W(:, j) - W(:, i)
	END DO
END DO

k = 0
DO i = 1, (my-1)
	DO j = (i+1), my
		k = k + 1
		Vdiff(:, k) = V(:, j) - V(:, i)
	END DO
END DO

hoptW = (1.06/EXP(1.0))*((SUM(Wdiff**2) - SUM(Wdiff)**2/DBLE(N_x))/DBLE(N_x-1))*N_x**(-0.2)
hoptV = (1.06/EXP(1.0))*((SUM(Vdiff**2) - SUM(Vdiff)**2/DBLE(N_y))/DBLE(N_y-1))*N_y**(-0.2)

DO i = 1, nseq
	Ux(i) = 0
	Uy(i) = 0
	DO j = 1, nGL
		x = GL(j, 1) + 1
		wei = GL(j, 2)
		denomx = sqrt(abs(SUM(cos(x*Wdiff/(2*hoptW)))/DBLE(N_x)))
		denomy = sqrt(abs(SUM(cos(x*Vdiff/(2*hoptV)))/DBLE(N_y)))
		Ux(i) = Ux(i) + wei*sin(errseq(i)*x/(2*hoptW))*denomx*(1 - x**2/4)**(1.5)/x
		Uy(i) = Uy(i) + wei*sin(errseq(i)*x/(2*hoptV))*denomy*(1 - x**2/4)**(1.5)/x
	END DO
	Ux(i) = Ux(i)/Pi + 0.5
	Uy(i) = Uy(i)/Pi + 0.5
	IF (Ux(i) > 1) THEN 
		Ux(i) = 1 
	END IF
	IF (Ux(i) < 0) THEN 
		Ux(i) = 0 
	END IF
	IF (Uy(i) > 1) THEN 
		Uy(i) = 1 
	END IF
	IF (Uy(i) < 0) THEN 
		Uy(i) = 0 
	END IF	
END DO

RETURN
END SUBROUTINE Uhat

SUBROUTINE BootGen(nx, ny, mx, my, nboot, F, Ux, Uy, zseq, errseq,& 
				nseq, Xb, Yb, Xerr, Yerr)
IMPLICIT NONE

INTEGER :: nx, ny, mx, my, nboot, nseq
real(8) :: F(nseq), Ux(nseq), Uy(nseq), zseq(nseq), errseq(nseq)
real(8) :: Xb(nx*nboot), Yb(ny*nboot), Xerr(nx*mx*nboot), Yerr(ny*my*nboot)

INTEGER :: i, idx
real(8) :: u

! ----- variables for portable seed setting -----
INTEGER :: i_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
INTEGER, DIMENSION(1:8) :: dt_seed
! ----- end of variables for seed setting -----

! ----- Set up random seed portably -----
CALL RANDOM_SEED(size=i_seed)
ALLOCATE(a_seed(1:i_seed))
CALL RANDOM_SEED(get=a_seed)
CALL DATE_AND_TIME(values=dt_seed)
a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
CALL RANDOM_SEED(put=a_seed)
DEALLOCATE(a_seed)
! ----- Done setting up random seed -----

!	CALL RANDOM_NUMBER(u)
!	idx = MINLOC(F, DIM = 1, MASK = (F >= u))
!	tmp = zseq(idx)

DO i = 1, nx*nboot
	CALL RANDOM_NUMBER(u)
	idx = MINLOC(F, DIM = 1, MASK = (F >= u))
	Xb(i) = zseq(idx)
END DO	

DO i = 1, ny*nboot
	CALL RANDOM_NUMBER(u)
	idx = MINLOC(F, DIM = 1, MASK = (F >= u))
	Yb(i) = zseq(idx)
END DO	

DO i = 1, nx*mx*nboot
	CALL RANDOM_NUMBER(u)
	idx = MINLOC(F, DIM = 1, MASK = (Ux >= u))
	Xerr(i) = errseq(idx) 
END DO	
	
DO i = 1, ny*my*nboot
	CALL RANDOM_NUMBER(u)
	idx = MINLOC(F, DIM = 1, MASK = (Uy >= u))
	Yerr(i) = errseq(idx) 
END DO	

RETURN
END SUBROUTINE BootGen



SUBROUTINE BootDist(nx, ny, mx, my, nboot, Xb, Yb, Xerr, Yerr, &
				GL, nGL, BDs, Bdist)
IMPLICIT NONE

INTEGER :: nx, ny, mx, my, nboot, nGH, nGL, nVARs
real(8) :: Xb(nx*nboot), Yb(ny*nboot), Xerr(nx*mx*nboot), Yerr(ny*my*nboot)
real(8) :: Bdist(nboot,2)
real(8) :: GL(nGL, 2), BDs(2)

real(8) :: Wtmp(nx, mx), Vtmp(ny, my), TS(2)
INTEGER :: i, j

!	DO i = 1, mx
!		Wtmp(:,i) = Xb + RESHAPE(Xerr(((i-1)*(nx*nboot)+1):(nx*nboot*i)), (/nx*nboot/))
!	END DO
!	DO i = 1, my
!		Vtmp(:,i) = Yb + RESHAPE(Yerr(((i-1)*(ny*nboot)+1):(ny*nboot*i)), (/ny*nboot/) )
!	END DO

!	DO i = 1, nboot
!		CALL TScomp(Wtmp(((i-1)*nx+1):(i*nx), :), Vtmp(((i-1)*ny+1):(i*nx), :),&
!			nx, ny, mx, my, GH, GL, nGH, nGL, VARs, BDs, nVARs, nBDs, TS)
!		Bdist(i, :) = TS
!	END DO 

DO i = 1, nboot
	DO j = 1, mx
		Wtmp(:,j) = Xb(((i-1)*nx + 1):(i*nx))
	END DO
	Wtmp = Wtmp + RESHAPE(Xerr(((i-1)*(nx*mx)+1):(i*nx*mx)), (/nx, mx/))
	DO j = 1, my
		Vtmp(:,j) = Yb(((i-1)*ny + 1):(i*ny))
	END DO
	Vtmp = Vtmp + RESHAPE(Yerr(((i-1)*(ny*my)+1):(i*ny*my)), (/ny, my/))
	CALL TScomp(Wtmp, Vtmp, nx, ny, mx, my, GL, nGL, BDs, TS)
	Bdist(i,:) = TS(:)
END DO

RETURN
END SUBROUTINE BootDist	

SUBROUTINE ABhat(W, V, nx, ny, mx, my, points, nt0, awhat, bwhat, avhat, bvhat)
IMPLICIT NONE

INTEGER :: nx, ny, nGL, nt0, N_x, N_y, mx, my
real(8) :: W(nx, mx), V(ny, my)
real(8) :: points(nt0)
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: avgx(nx), avgy(ny)

real(8) :: awhat(nt0), bwhat(nt0), avhat(nt0), bvhat(nt0)
INTEGER :: i, j, k
real(8) :: denomx, denomy

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

avgx(:) = SUM(W, DIM = 2)/DBLE(mx)
avgy(:) = SUM(V, DIM = 2)/DBLE(my)
	
k = 0
DO i = 1, (mx-1)
	DO j = (i+1), mx
		k = k + 1
		Wdiff(:, k) = W(:, j) - W(:, i)
	END DO
END DO

k = 0
DO i = 1, (my-1)
	DO j = (i+1), my
		k = k + 1
		Vdiff(:, k) = V(:, j) - V(:, i)
	END DO
END DO

DO i = 1, nt0
	denomx = sqrt(abs(SUM(cos(points(i)*Wdiff/DBLE(mx)))/DBLE(N_x)))**(DBLE(mx))
	denomy = sqrt(abs(SUM(cos(points(i)*Vdiff/DBLE(my)))/DBLE(N_y)))**(DBLE(my))
	awhat(i) = (SUM(cos(points(i)*avgx(:)))/DBLE(nx))/denomx
	bwhat(i) = (SUM(sin(points(i)*avgx(:)))/DBLE(nx))/denomx
	avhat(i) = (SUM(cos(points(i)*avgy(:)))/DBLE(ny))/denomy
	bvhat(i) = (SUM(sin(points(i)*avgy(:)))/DBLE(ny))/denomy
END DO

RETURN
END SUBROUTINE ABhat


SUBROUTINE ABhat1(W, V, nx, ny, mx, my, points, nt0, Bidx, nboot, &
				awhat, bwhat, avhat, bvhat)
IMPLICIT NONE

INTEGER :: nx, ny, nGL, nt0, N_x, N_y, mx, my, nboot, Bidx(nboot*nx)
real(8) :: W(nx, mx), V(ny, my), Wtmp(nx, mx), Vtmp(ny, my)
real(8) :: points(nt0)
real(8) :: Wdiff(nx, mx*(mx-1)/2), Vdiff(ny, my*(my-1)/2)
real(8) :: avgx(nx), avgy(ny)

real(8) :: awhat(nt0, nboot), bwhat(nt0, nboot), avhat(nt0, nboot), bvhat(nt0, nboot)
real(8) :: awhat_tmp(nt0), bwhat_tmp(nt0), avhat_tmp(nt0), bvhat_tmp(nt0)
INTEGER :: i, j, k, l
real(8) :: denomx, denomy

N_x = nx*mx*(mx-1)/2
N_y = ny*my*(my-1)/2

DO l = 1, nboot

	Wtmp(:, :) = W(Bidx(((l-1)*nx+1):(l*nx)), :)
	Vtmp(:, :) = V(Bidx(((l-1)*nx+1):(l*nx)), :)

	avgx(:) = SUM(Wtmp, DIM = 2)/DBLE(mx)
	avgy(:) = SUM(Vtmp, DIM = 2)/DBLE(my)

	k = 0
	DO i = 1, (mx-1)
		DO j = (i+1), mx
			k = k + 1
			Wdiff(:, k) = Wtmp(:, j) - Wtmp(:, i)
		END DO
	END DO

	k = 0
	DO i = 1, (my-1)
		DO j = (i+1), my
			k = k + 1
			Vdiff(:, k) = Vtmp(:, j) - Vtmp(:, i)
		END DO
	END DO

	DO i = 1, nt0
		denomx = sqrt(abs(SUM(cos(points(i)*Wdiff/DBLE(mx)))/DBLE(N_x)))**(DBLE(mx))
		denomy = sqrt(abs(SUM(cos(points(i)*Vdiff/DBLE(my)))/DBLE(N_y)))**(DBLE(my))
		awhat_tmp(i) = (SUM(cos(points(i)*avgx(:)))/DBLE(nx))/denomx
		bwhat_tmp(i) = (SUM(sin(points(i)*avgx(:)))/DBLE(nx))/denomx
		avhat_tmp(i) = (SUM(cos(points(i)*avgy(:)))/DBLE(ny))/denomy
		bvhat_tmp(i) = (SUM(sin(points(i)*avgy(:)))/DBLE(ny))/denomy
	END DO
	
	awhat(:, l) = awhat_tmp(:)
	bwhat(:, l) = bwhat_tmp(:)
	avhat(:, l) = avhat_tmp(:)
	bvhat(:, l) = bvhat_tmp(:)
	
END DO

RETURN
END SUBROUTINE ABhat1


