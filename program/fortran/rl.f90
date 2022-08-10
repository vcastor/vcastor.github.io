PROGRAM linearregression

!***********************************************************************
!------      This program compute a linear regression and the     ------
!------                 uncertainties with n-points               ------
!------                                                           ------
!------                       ~VCastor 2022                       ------
!***********************************************************************

IMPLICIT NONE
INTEGER                                 :: i, npoints
REAL(KIND=8)                            :: m, b, sx, sx2, sy, sy2, sxy
REAL(KIND=8)                            :: symxb2, ub, um, r
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: x, y

!=======================================================================
OPEN(108,FILE="data.dat")
  READ(108,*) npoints
  ALLOCATE(x(npoints), y(npoints))
  DO i = 1, npoints
    READ(108,*) x(i), y(i)
  ENDDO
CLOSE(108)
!=======================================================================

sx = 0.d0; sy = 0.d0; sx2 = 0.d0; sy2 = 0.d0; sxy = 0.d0

DO i = 1, npoints
  sx  = sx  + x(i)
  sy  = sy  + y(i)
  sx2 = sx2 + x(i)*x(i)
  sy2 = sy2 + y(i)*y(i)
  sxy = sxy + x(i)*y(i)
ENDDO

m = (npoints*sxy - sx*sy)/(npoints*sx2 - sx*sx)
b = ((sy*sx2) - (sx*sxy))/(npoints*sx2 - (sx*sx))
r = (sxy - sx*sy/npoints)/DSQRT((sx2 - sx*sx/npoints)*&
                                                  (sy2 - sy*sy/npoints))
!UNCERTAINTIES
symxb2 = 0.d0
DO i = 1, npoints
  symxb2 = symxb2 + (y(i) -m*x(i) -b)*(y(i) -m*x(i) -b)
ENDDO

um = DSQRT(npoints*symxb2/((npoints-2)*(npoints*sx2 -sx*sx)))
ub = DSQRT(sx2*symxb2/((npoints-2)*(npoints*sx2 -sx*sx)))

!=======================================================================
OPEN(109,FILE="data.csv")
  DO i = 1, npoints
    WRITE(109,FMT='(F14.8,A1,F14.8)') x(i), ',', y(i)
  ENDDO
CLOSE(109)

WRITE(*,FMT='(A6,F15.8)')"b    =", b
WRITE(*,FMT='(A6,F15.8)')"m    =", m
WRITE(*,FMT='(A6,F15.8)')"U(m) =", um
WRITE(*,FMT='(A6,F15.8)')"U(b) =", ub
WRITE(*,FMT='(A6,F15.8)')"r    =", r
WRITE(*,FMT='(A6,F15.8)')"r^2  =", r*r
!***********************************************************************
ENDPROGRAM linearregression
