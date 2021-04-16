PROGRAM polynomial

!-----------------------------------------------------------------------------
!---------   This program compute the polynomial regresion and the   ---------
!---------                 r^2 for n-pair of points.                 ---------
!---------                                                           ---------
!---------                       ~VCastor 2020                       ---------
!-----------------------------------------------------------------------------

IMPLICIT NONE
INTEGER :: j, k, n, m
REAL :: rr
REAL, DIMENSION(:,:), ALLOCATABLE :: MAT, D 
CHARACTER(LEN=60) :: f, str


!-----------------------------------------------------------------------------
!---- The data will be read

CALL Reading(D,n)

!-----------------------------------------------------------------------------
!---- Can we compute?

WRITE(*,FMT='(A)',ADVANCE='NO') "Which degree for the polynomial regresion &
                                 you want? "
READ(*,*) m

IF ( m >= n ) THEN
    WRITE(*,*) '**************************************'
    WRITE(*,FMT='(A26,I1,A35,I1,A14)') "Is not possible compute a ",m,"th &
                                        degree polynomial equation &
                                        with ", n, " pairs numbers"
    WRITE(*,*) "exist more than one equation that fulfill that criterion"
    WRITE(*,*) '**************************************'
    STOP 'Try with lower polynomial degree'
ENDIF

!-----------------------------------------------------------------------------
!---- Compute the values of the matrix and get the matrix 

CALL Setting(MAT,D,n,m)

!-----------------------------------------------------------------------------
!---- Diagonalizate the matrix

CALL Solve(MAT,m)

!-----------------------------------------------------------------------------
!---- Write the coeficientes, (too fancy UuUr)

f='(A3,'//TRIM(str(m))//'(F6.3,A2,I1,A3),F6.3)'
WRITE(*,*) '**************************************'
WRITE(*,FMT=f) 'y= ', (MAT(j,m+2), 'x^',j-1, '+ ', j=m+1,2,-1), MAT(1,m+2)
WRITE(*,*) '**************************************'

!-----------------------------------------------------------------------------
!---- Also compute and write R^2

CALL R2(D,n,m,rr)
WRITE(*,FMT='(A11,F6.4)') 'With r^2 = ', rr
WRITE(*,*) '**************************************'

!------------------------------------------------------------------------------
!---- The subroutines that we need:
CONTAINS
    SUBROUTINE Reading(D,n)
    IMPLICIT NONE
    INTEGER :: i, j, n
    REAL, DIMENSION(:,:), ALLOCATABLE :: D
    CHARACTER(LEN=50) :: archivo

    WRITE(*,*) "This program needs a n-pair of points &
                to compute a polynomial regresion, "
    WRITE(*,*) "The file data should have in the first line &
                how many pairs we will have, then the pairs"
    WRITE(*,FMT='(A)',ADVANCE='NO') "File name of the data: "
    READ(*,*) archivo
    OPEN(11,FILE=archivo) 
    READ(11,*) n                  !The 1st line: the number of pairs we have
    ALLOCATE(D(n,2))

    DO i = 1, n                   !Read every pair as a row in a Matrix D
           READ(11,*) (D(i,j), j=1,2)
    ENDDO
    ENDSUBROUTINE Reading
!------------------------------------------------------------------------------
    SUBROUTINE Setting(MAT,D,n,m)
    IMPLICIT NONE
    INTEGER :: i, j, n, m
    REAL :: ss, cc
    REAL, DIMENSION(n,2) :: D
    REAL, DIMENSION(:,:), ALLOCATABLE :: MAT
    REAL, DIMENSION(:), ALLOCATABLE :: S, C

    ALLOCATE(MAT(m+1,m+2))        !The matrix will have the dimension m+1 & m+2

    !---- Compute the necesary values
    ALLOCATE(S(2*m))              !These vector will have all sums possibles for x

    DO i = 1, 2*m
        ss = 0.                   !Reset the summ every k-loop
        DO k = 1, n
            ss = ss + D(k,1)**i
        ENDDO
        S(i) = ss
    ENDDO                         !Now we have the values of \sum x_i^j
                                  !in the vector, that was for no compute
                                  !more than one time the same sum

    !----- Recurrent values will be in the correct positions of the matrix
    DO k = 1, 2*m                 !Size of S
        DO i = 1, m+1             !Part of MAT which use x_i^j
            DO j = 1, m+1         !Part of MAT which use x_i^j
                IF ( (i+j) == (k+2) ) THEN
                    MAT(i,j) = S(k)
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    MAT(1,1) = n                  !The first value

    !---- The values which use y_i
    ALLOCATE(C(m+1))
    DO i = 1, m+1
        cc = 0.                   !Restart the sum for the j-loop
        DO j = 1, n
            !           y    *   x   ** i-1
            cc = cc + (D(j,2)*(D(j,1)**(i-1)))
        ENDDO
        C(i) = cc
    ENDDO

    !---- The matrix will be the contatenation of the both before parts
    DO i = 1, m+1
        MAT(i,m+2) = C(i)
    ENDDO
    ENDSUBROUTINE Setting
!------------------------------------------------------------------------------
    SUBROUTINE Solve(MAT,m)
    IMPLICIT NONE
    INTEGER :: i, j, k, m
    REAL :: aux, aux2
    REAL, DIMENSION(:), ALLOCATABLE :: V
    REAL, DIMENSION(m+1,m+2) :: MAT

    ALLOCATE(V(m+2))              !axuliar to get the numbers of the matrix rows

    !---- First up diagonal
    DO i = 1, m                   !for all rows except the last one
        aux = MAT(i,i)            !diagonal elemnt
        DO j = 1, m+2             !all elements in the row
            V(j) = MAT(i,j)/aux   !divided by the diagonal elemnt
        ENDDO
        DO k = i+1, m+1           !only with the rows under itself
            aux2 = MAT(k,i)
            DO j = 1, m+2         !all elements in the row
                MAT(k,j) = MAT(k,j) - V(j)*aux2
            ENDDO
        ENDDO
    ENDDO

    !---- 2nd zeros in up diagonal
    DO i = 1, m                   !for all rows except the last one
        DO k = i+1, m+1           !only with the rows under itself
            aux = MAT(k,k)
            DO j = 1, m+2         !all elements in the row
                V(j) = MAT(k,j)/aux
            ENDDO
            aux2 = MAT(i,k)
            DO j = 1, m+2         !all elements in the row
                MAT(i,j) = MAT(i,j) - V(j)*aux2
            ENDDO
        ENDDO
    ENDDO

    !---- 3rd diagonal elements equals to one
    DO i = 1, m+1 
        aux = MAT(i,i)            !diagonal elemnt
        DO j = 1, m+2 
            MAT(i,j) = MAT(i,j)/aux
        ENDDO
    ENDDO
    ENDSUBROUTINE Solve
!------------------------------------------------------------------------------
    SUBROUTINE R2(D,n,m,rr)
    IMPLICIT NONE
    INTEGER :: i, j, n, m
    REAL :: tot, res, ss, rr, ytes
    REAL, DIMENSION(n) :: pred
    REAL, DIMENSION(n,2) :: D

    ss = 0.
    DO i = 1, n
        ss = ss + D(i,2)
    ENDDO
    ytes = ss/REAL(n)              !average of y data values

    pred(:) = 0.                   !predicted by the regresion equation
    DO i = 1, n                    
        DO j = m+1, 2, -1          !differents coefficients of equation
            !                   coeficient*   x  ** j-1
            pred(i) = pred(i) + MAT(j,m+2)*D(i,1)**(j-1)
        ENDDO
        pred(i) = pred(i) + MAT(1,m+2)     !plus constant coefficient
    ENDDO

    tot = 0. ; res = 0.                    !start the sum
    DO i = 1, n
        tot = tot + (D(i,2) - ytes)**2     !squared diff between y & average
        res = res + (D(i,2) - pred(i))**2  !squared diff between y & predicted
    ENDDO

    rr = 1 - res/tot                !R^2 is easy so compute
    ENDSUBROUTINE R2
ENDPROGRAM polynomial
!------------------------------------------------------------------------------
!---- For fancy style, transform a number to charecter
CHARACTER(LEN=60) FUNCTION str(q)
    INTEGER, INTENT(IN) :: q
    WRITE (str,*) q
    str = ADJUSTL(str)
ENDFUNCTION str 
