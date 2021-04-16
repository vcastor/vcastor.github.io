PROGRAM huckel

!--------------------------------------------------------------------
!----       This program use the Jacobi algorithm to compute     ----
!----       the eigenvalues and eigenvectors in a subroutine     ----
!----                to the Hückel approximation                 ----
!----                                                            ----
!----                       ~VCastor 2020                        ----
!--------------------------------------------------------------------

IMPLICIT NONE
INTEGER :: i, j, k, l, atoms, bonds, pi_e
REAL*8, PARAMETER :: e = 1.0e-08
REAL*8 :: alpha, beta, energy_s
INTEGER, DIMENSION(:,:), ALLOCATABLE :: D
REAL*8, DIMENSION(:,:), ALLOCATABLE :: H, X, Info
REAL*8, DIMENSION(:), ALLOCATABLE :: Gama, V
CHARACTER(LEN=40) :: f
CHARACTER(LEN=60) :: archivo, title, outputfile

!--------------------------------------------------------------------
!---- Data

WRITE(*,*) 'This program compute the MO energy and &
            coefficients, using Hückel approximation'
WRITE(*,*) 'The data file would be written as the templates &
            given in the zip file with these code'
WRITE(*,FMT='(A)',ADVANCE='NO') 'Where is the data? '
READ(*,*) archivo

OPEN(11,FILE=archivo)
READ(11,*) alpha, beta
READ(11,*) title
READ(11,*) atoms, bonds, pi_e
ALLOCATE(D(bonds,2))
DO i = 1, bonds
    READ(11,*) D(i,1), D(i,2)
ENDDO

WRITE(*,*) 
WRITE(*,FMT='(A11,*(A))') 'Job title: ', title
WRITE(*,*)

!--------------------------------------------------------------------
!---- Let's starts with the Hamiltonian

ALLOCATE(H(atoms,atoms),X(atoms,atoms))
H(:,:) = 0.                    !Set zeros

DO i = 1, atoms
    H(i,i) = alpha             !diagonal
ENDDO

DO i = 1, bonds
    H(D(i,1),D(i,2)) = beta    !adjacents
    H(D(i,2),D(i,1)) = beta
ENDDO

WRITE(*,*) ' '
WRITE(*,*) '****************************************'
WRITE(*,*) '**********The hamiltonian is: **********'
DO i = 1, atoms
    WRITE(*,FMT='(*(F7.1))') (H(i,j), j=1,atoms)
ENDDO

!--------------------------------------------------------------------
!---- Calling the subroutine, where H is the initial
!---- Matrix, X will have the eigenvectors, e is the
!---- tolerance to stop, and atoms is the dimension of
!---- the matrix

CALL jacobi(H,X,e,atoms)

!---- Sort the resoutls, the lowest energy is the last row

ALLOCATE(Info(atoms,atoms+1))  !to have: (E1|phi1,ph2...)...
DO i = 1, atoms
    Info(i,1) = H(i,i)         !energy values
    DO j = 2, atoms+1
        Info(i,j) = X(j-1,i)   !column vectors to row vectors
    ENDDO
ENDDO

ALLOCATE(V(atoms+1))           !auxiliar to have the rows of the matrix
DO i = 1, atoms
    DO j = i+1, atoms
        IF (Info(i,1) < Info(j,1)) THEN  !larger value at the top
            V(:) = Info(i,:)             
            Info(i,:) = Info(j,:)        !flip the rows
            Info(j,:) = V(:)
        ENDIF
    ENDDO
ENDDO

!---- Energy of the pi system, 

energy_s = 0.
DO i = 1, pi_e/2               !two electrons per orbital
    energy_s = energy_s + Info(atoms+1-i,1)
ENDDO
energy_s = 2*energy_s          !two electrons per orbital

IF (.NOT.(MOD(pi_e,2).EQ.0)) THEN     !if we have a odd e-number
    energy_s = energy_s + Info(atoms-pi_e/2,1)
ENDIF

!---- Let's have the eigenvalues as a sum of alpha & beta
ALLOCATE(Gama(atoms))
DO i = 1, atoms
    Gama(i) = (Info(i,1)-alpha)/beta
ENDDO

!---- Write reoults

outputfile=TRIM(title)//'.out'
OPEN(11,FILE=outputfile)

WRITE(11,*) ' '
WRITE(11,*) '****************************************'
WRITE(11,*) '************Eigenvalues (eV)************'
WRITE(11,*) ' '
f = '(A2,I2,A10,F5.3,A8,F8.4)'
DO i = 1, atoms
    IF (Gama(i) < 0.0) THEN    !to print + & - values without doubles + or -
        WRITE(11,FMT=f) 'E_',atoms+1-i,'= alpha -', -Gama(i), 'beta = ', Info(i,1)
    ELSE
        WRITE(11,FMT=f) 'E_',atoms+1-i,'= alpha +', Gama(i), 'beta = ', Info(i,1)
    ENDIF
ENDDO

WRITE(11,*) ' '
WRITE(11,*) '****************************************'
WRITE(11,*) '**************Eigenvectors**************'
WRITE(11,*) ' '
!to have table with values of phi_n for Psi_n
WRITE(11,FMT='(A11,*(A4,I2,A2))') '*******    ', ('phi_',j,'  ', j=1,atoms)
DO i = 1, atoms
    WRITE(11,FMT='(A4,I2,A3,*(F8.4))') 'Psi_',atoms+1-i,'= ',(Info(i,j+1), j=1,atoms)
ENDDO

WRITE(11,*) ' '
WRITE(11,*) '****************************************'
WRITE(11,FMT='(A29,F9.4,A3)') ' The energy of pi system is: ', energy_s, ' eV'
WRITE(11,*) '****************************************'
WRITE(11,*) ' '

WRITE(*,*) 'The resoults are in: ', outputfile
ENDPROGRAM huckel

!--------------------------------------------------------------------
!---- Jacobi algorithm
!---- X will have the eigenvectors as column vectors 
!---- H will have the eigenvalues in the diagonal
!--------------------------------------------------------------------
SUBROUTINE jacobi(H,X,e,n)
IMPLICIT NONE
INTEGER :: i, j, k, n
REAL*8 :: e, b2, bar, ttheta, coeff, c, s, cs, sc, t, aux
REAL*8, DIMENSION(n,n) :: H, X

!---- Starts with the Identity matrix
X = 0.
DO i = 1, n
    X(i,i) = 1.
ENDDO

b2 = 0.             !sum for no-diagonal elements squared
DO i = 1, n
    DO j = 1, n
        IF ( i /= j ) THEN
            b2 = b2 + H(i,j)**2
        ENDIF
    ENDDO
ENDDO

IF (b2 <= e) RETURN  !e says if it's ok, i.e., when no-
                     !diagonal is like 0. Diagonalizaded

bar = 0.5*b2/FLOAT(n*n)   !average for off-diagonal /2

DO WHILE (b2 > e)         !until converge
    DO i = 1, n-1
        DO j = i+1, n
            IF (H(j,i)**2 <= bar) CYCLE  !no with small elements
            b2 = b2 - 2.0*H(j,i)**2
            bar = 0.5*b2/FLOAT(n*n)
            !Compute with the given matrix
            ttheta = (H(j,j)-H(i,i))/(2.0*H(j,i))
            coeff = 0.5*ttheta/SQRT(1. + ttheta**2)
            s = SQRT(MAX(0.5 + coeff, 0.))                 !sin
            c = SQRT(MAX(0.5 - coeff, 0.))                 !cos

            DO k = 1, n   !Recalculate rows
                cs =  c*H(i,k) + s*H(j,k)           ! cos + sin
                sc = -s*H(i,k) + c*H(j,k)           !-sin + cos
                H(i,k) = cs
                H(j,k) = sc
            ENDDO

            DO k = 1, n   !The new matrix H, and the eigenvectors
                cs =  c*H(k,i) + s*H(k,j)
                sc = -s*H(k,i) + c*H(k,j)
                H(k,i) = cs
                H(k,j) = sc
                cs =  c*X(k,i) + s*X(k,j)
                sc = -s*X(k,i) + c*X(k,j)
                X(k,i) = cs
                X(k,j) = sc
            ENDDO
        ENDDO
    ENDDO
ENDDO

DO i = 1, n            !Normalize if we are not
    aux = NORM2(X(i,:))
    X(i,:) = X(i,:)/aux
ENDDO

RETURN
ENDSUBROUTINE
