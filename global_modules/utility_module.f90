MODULE utility_module

    USE parameter_module, ONLY: DP, ZERO, ONE, TWO, HALF, ERROR

    IMPLICIT NONE

    CONTAINS
            
        SUBROUTINE PrintMatrix(realMat, intMat, unit)
            !Purpose:
            ! Prints the Matrix A row by row in the specified unit 
            ! 

            REAL(DP), OPTIONAL, INTENT(IN) :: realMat(:,:)
            INTEGER,  OPTIONAL, INTENT(IN) :: intMat(:,:)
            INTEGER,  OPTIONAL, INTENT(IN) :: unit

            INTEGER :: i, j
            INTEGER :: Ncol 

            CHARACTER(len=12) :: fmt

            IF (PRESENT(realMat)) THEN 

                Ncol = SIZE(realMat,2)
                WRITE(fmt, '(A, I0, A)') '(5X', Ncol, 'ES15.7)'
                IF (.NOT. PRESENT(unit)) THEN
                    DO i = 1, SIZE(realMat,1)
                        WRITE(*, fmt) (realMat(i,j), j = 1,Ncol)
                    END DO
                END IF
                
                IF (PRESENT(unit)) THEN
                    DO i = 1, SIZE(realMat,1)
                        WRITE(unit, fmt) (realMat(i,j), j = 1,Ncol)
                    END DO
                END IF

            ELSE IF (PRESENT(intMat)) THEN 
                
                Ncol = SIZE(intMat,2)
                WRITE(fmt, '(A, I0, A)') '(5X', Ncol, 'I10)'

                IF (.NOT. PRESENT(unit)) THEN
                    DO i = 1, SIZE(intMat,1)
                        WRITE(*, fmt) (intMat(i,j), j = 1,Ncol)
                    END DO
                END IF
                
                IF (PRESENT(unit)) THEN
                    DO i = 1, SIZE(intMat,1)
                        WRITE(unit, fmt) (intMat(i,j), j = 1,Ncol)
                    END DO
                END IF

            ELSE
                RETURN
            END IF

        END SUBROUTINE PrintMatrix

        INTEGER FUNCTION newunit() RESULT(n)

            LOGICAL :: inuse
            INTEGER :: nmin = 101
            INTEGER :: nmax = 500

            DO n = nmin, nmax

                INQUIRE(UNIT=n, OPENED=inuse)

                IF (.NOT. inuse) THEN
                    RETURN
                END IF
            END DO

        END FUNCTION

        PURE FUNCTION KDelta(i, j) 
            !Purpose
            ! Kronecker Delta function
            
            REAL(DP) :: KDelta

            INTEGER, INTENT(IN)   :: i
            INTEGER, INTENT(IN)   :: j

            IF (i .EQ. j) KDelta = ONE

            IF (i .NE. j) KDelta = ZERO

        END FUNCTION 

        PURE FUNCTION DoubleContraction42(A, B) RESULT(C)
            !Purpose: 
            ! Define the double contraction product of two tensors of order 4 and 2 
            ! respectively
            !
            !Inputs:
            !   A : Tensor of order 4
            !   B : Tensor of order 2
            !
            !Outputs:
            !   C : Tensor of order 2

            REAL(DP),  INTENT(IN)  :: A(:,:,:,:)
            REAL(DP),  INTENT(IN)  :: B(:,:)
        
            REAL(DP), ALLOCATABLE :: C(:,:)

            ! Loop varibales
            INTEGER :: i, j, k, l

            ! Allocating the size of C
            ALLOCATE(C(SIZE(A,1), SIZE(A,2)))

            C(:,:) = ZERO

            DO i = 1, SIZE(A,1)
             DO j = 1, SIZE(A,2)
              DO k = 1, SIZE(A,3)
               DO l = 1, SIZE(A,4)
                C(i,j) = C(i,j) + A(i,j,k,l) * B(k,l)
               END DO
              END DO
             END DO
            END DO

            DEALLOCATE(C)

        END FUNCTION DoubleContraction42

        PURE FUNCTION DoubleContraction24(A, B) RESULT(C)
            !Purpose: 
            ! Define the double contraction product of two tensors of order 2 and 4 
            ! respectively
            !
            !Inputs:
            !   A : Tensor of order 2
            !   B : Tensor of order 4
            !
            !Outputs:
            !   C : Tensor of order 2

            REAL(DP),  INTENT(IN)  :: A(:,:)
            REAL(DP),  INTENT(IN)  :: B(:,:,:,:)
        
            REAL(DP), ALLOCATABLE :: C(:,:)

            ! Loop varibales
            INTEGER :: i, j, k, l

            ! Allocating the size of C
            ALLOCATE(C(SIZE(A,1), SIZE(A,2)))

            C(:,:) = ZERO

            DO i = 1, SIZE(B,1)
             DO j = 1, SIZE(B,2)
              DO k = 1, SIZE(B,3)
               DO l = 1, SIZE(B,4)
                C(k,l) = C(k,l) + A(i,j) * B(i,j,k,l)
               END DO
              END DO
             END DO
            END DO

            DEALLOCATE(C)

        END FUNCTION DoubleContraction24

        
        PURE FUNCTION DoubleContraction22(A, B) RESULT(prod)
            !Purpose: 
            ! Define the double contraction product of two tensors of order 2 
            !
            !Inputs:
            !   A : Tensor of order 2
            !   B : Tensor of order 2
            !
            !Outputs:
            !   C : scalar productive
            REAL(DP),  INTENT(IN) :: A(:,:)
            REAL(DP),  INTENT(IN) :: B(:,:)
            
            REAL(DP) :: prod

            ! Loop varibales
            INTEGER :: i, j

            ! Allocating the size of C
            prod = ZERO

            DO i = 1, SIZE(A,1)
                DO j = 1, SIZE(A,2)
                    prod = prod + A(i,j) * B(i,j)
                END DO
            END DO
        END FUNCTION DoubleContraction22

        PURE FUNCTION TensorProduct2(A, B) RESULT(C)
            !Purpose:
            ! Tensor product of a 2nd order tensors
            !Note : Make sure that the dimensions of A and B are same

            REAL(DP), INTENT(IN) :: A(:,:)
            REAL(DP), INTENT(IN) :: B(:,:)
            
            REAL(DP), ALLOCATABLE :: C(:,:,:,:)

            ! Loop variables
            INTEGER :: i, j, k, l, dim 

            dim = SIZE(A,1)

            ALLOCATE(C(dim, dim, dim, dim))

            DO i = 1, dim
             DO j = 1, dim
              DO k = 1, dim
               DO l = 1, dim
                C(i,j,k,l) = A(i,j) * B(k,l)
               END DO
              END DO
             END DO
            END DO

            DEALLOCATE(C)

        END FUNCTION TensorProduct2 

        PURE FUNCTION Inverse2(A, det) RESULT(Ainv)
            !Purpose:
            ! Find the inverse of a 2,2 matrix 

            REAL(DP), INTENT(IN) :: A(2,2)
            REAL(DP), INTENT(IN) :: det

            ! Output
            REAL(DP) :: Ainv(2,2)

            Ainv(1,1) = (1/det) * A(2,2)
            Ainv(2,2) = (1/det) * A(1,1)
            Ainv(1,2) = - (1/det) * A(1,2)
            Ainv(2,1) = - (1/det) * A(2,1)

        END FUNCTION Inverse2

        PURE FUNCTION Inverse3(A, det) RESULT(Ainv)

            REAL(DP), INTENT(IN) :: A(3,3)
            REAL(DP), INTENT(IN) :: det

            REAL(DP) :: Ainv(3,3)

            Ainv = ZERO            

            Ainv(1,1) =  (ONE/det) * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
            Ainv(2,1) = -(ONE/det) * (A(2,1)*A(3,3) - A(2,3)*A(3,2))
            Ainv(3,1) =  (ONE/det) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            Ainv(1,2) = -(ONE/det) * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
            Ainv(2,2) =  (ONE/det) * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
            Ainv(3,2) = -(ONE/det) * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
            Ainv(1,3) =  (ONE/det) * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
            Ainv(1,3) = -(ONE/det) * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
            Ainv(1,3) =  (ONE/det) * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

        END FUNCTION Inverse3 

        PURE FUNCTION checkSymm(A) RESULT(Issym)

            REAL(DP), INTENT(IN) :: A(:,:)

            LOGICAL :: Issym

            INTEGER :: i, j

            DO i = 1, SIZE(A,1)
                DO j = 1, SIZE(A,2)
                    IF (ABS(A(i,j) - A(j,i)) .LE. ERROR) THEN 
                        Issym = .TRUE.
                    ELSE
                        Issym = .FALSE.
                        RETURN
                    END IF
                END DO
            END DO
            
        END FUNCTION

        PURE FUNCTION checkEq(A, B) RESULT(IsEq)

            REAL(DP), INTENT(IN) :: A(:,:)
            REAL(DP), INTENT(IN) :: B(:,:)

            LOGICAL :: IsEq

            REAL(DP) :: Diff(SIZE(A,1), SIZE(A,2))

            INTEGER :: i, j

            IF ((SIZE(A,1)) .EQ. SIZE(B,1)) THEN
                IsEq = .FALSE. 
                RETURN
            END IF

            IF ((SIZE(A,2)) .EQ. SIZE(B,2)) THEN
                IsEq = .FALSE. 
                RETURN
            END IF

            Diff = A - B

            DO i = 1, SIZE(A,1)
                DO j = 1, SIZE(A,2)
                    IF (Diff(i,j) .GT. ERROR) THEN
                        IsEq = .FALSE.
                        RETURN
                    END IF
                END DO
            END DO

            IsEq = .TRUE.

        END FUNCTION checkEq 

        PURE SUBROUTINE Exp2ss(A, expA, fstat)  
            !Pupose: 
            ! Exponential of a 2x2 skew symmetric matrix
            REAL(DP),  INTENT(IN) :: A(2,2)
            REAL(DP), INTENT(OUT) :: expA(2,2)
            INTEGER,  INTENT(OUT) :: fstat

            REAL(DP) :: th

            fstat = 0

            ! Check if the matrix is skew symmetric
            IF ((A(1,1) > ZERO + ERROR) .OR. (A(1,1) < ZERO - ERROR)) fstat=1 
            IF ((A(2,2) > ZERO + ERROR) .OR. (A(1,1) < ZERO - ERROR)) fstat=1
            IF ((A(1,2) > -A(2,1) + ERROR) .OR. (A(1,2) < -A(2,1) - ERROR)) fstat=1

            th = A(1,2)

            ! Calculating the exponential of the matrix
            expA(1,1) = COS(th)
            expA(2,2) = COS(th)
            expA(1,2) = SIN(th)
            expA(2,1) = -SIN(th)

        END SUBROUTINE Exp2ss

        PURE FUNCTION Voigtmap(dim) RESULT(map)
            
            CHARACTER(len=2), INTENT(IN) :: dim
            
            INTEGER, ALLOCATABLE ::  map(:,:)

            IF (dim .EQ. '2D') THEN
                ALLOCATE(map(3,2))
                map(1,1) = 1
                map(1,2) = 1
                map(2,1) = 2
                map(2,2) = 2
                map(3,1) = 1
                map(3,2) = 2
                
            ELSE IF (dim .EQ.'3D') THEN
                ALLOCATE(map(6,2))
                map(1,:) = (/1,1/)
                map(2,:) = (/2,2/)
                map(3,:) = (/3,3/)
                map(4,:) = (/2,3/)
                map(5,:) = (/1,3/)
                map(6,:) = (/1,2/)

            ELSE
                ! Incorrect dimension
                ALLOCATE(map(1,1))
                map = 0
            END IF
        END FUNCTION Voigtmap

        PURE FUNCTION Voigt4(C) RESULT(Cvoigt)
            !Purpose:
            ! Use Voigt kinetic rule to reduce a 4th order tensor to a 2D matrix
        
            REAL(DP), INTENT(IN) :: C(:,:,:,:)
            REAL(DP), ALLOCATABLE :: Cvoigt(:,:)

            INTEGER :: dim
            INTEGER :: i, j
            INTEGER, ALLOCATABLE :: map(:,:)

            dim = SIZE(C, 1)

            IF (dim .EQ. 2) THEN
                ALLOCATE(map(3,2))
                ALLOCATE(Cvoigt(3,3))
                map = Voigtmap('2D')
                DO i = 1,3
                    DO j = 1,3
                        Cvoigt(i,j) = C(map(i,1), map(i,2), map(j,1), map(j,2))
                    END DO
                END DO

            ELSE IF (dim .EQ. 3) THEN
                ALLOCATE(map(6,2))
                ALLOCATE(Cvoigt(6,6))
                map = Voigtmap('3D')
                DO i = 1,6
                    DO j = 1,6
                        Cvoigt(i,j) = C(map(i,1), map(i,2), map(j,1), map(j,2))
                    END DO
                END DO
            END IF

            DEALLOCATE(Cvoigt)
            DEALLOCATE(map)

        END FUNCTION Voigt4

        PURE FUNCTION Voigt2(T, state) RESULT(vT)
            ! Purpose:
            !   Calculates Voigt form of a 2nd order tensor
            !
            REAL(DP), INTENT(IN) :: T(:,:)
            INTEGER,  INTENT(IN) :: state
            

            REAL(DP) :: vT(SIZE(T,1) * (SIZE(T,1) + 1) / 2)
            REAL(DP), ALLOCATABLE :: vec(:)

            INTEGER :: dim, i
            INTEGER, ALLOCATABLE :: map(:,:)

            dim = SIZE(T,1)

            IF (dim .EQ. 2) THEN
                ALLOCATE(map(3,2))
                ALLOCATE(vec(3))
                map = Voigtmap('2D')
                DO i = 1,3
                    vec(i) = T(map(i,1), map(i,2))
                END DO
                IF (state .EQ. 1) THEN
                    vT = vec
                ELSE IF (state .EQ. 2) THEN
                    vT(1:2) = vec(1:2)
                    vT(3) = TWO * vec(3)
                END IF
            
            ELSE IF (dim .EQ. 3) THEN
                ALLOCATE(map(6,2))
                ALLOCATE(vec(6))
                map = Voigtmap('3D')
                DO i = 1,6
                    vec(i) = T(map(i,1), map(i,2))
                END DO
                IF (state .EQ. 1) THEN
                    vT = vec
                ELSE IF (state .EQ. 2) THEN
                    vT(1:3) = vec(1:3)
                    vT(4:6) = TWO * vec(4:6)
                END IF
            END IF 

            DEALLOCATE(map)
            DEALLOCATE(vec)

        END FUNCTION Voigt2

        PURE FUNCTION Voigt2Tensor(vT, state) RESULT(T)
            ! Purpose:
            !   Converts 2nd order Voigt vector to tensor
            REAL(DP), INTENT(IN) :: vT(:)
            INTEGER,  INTENT(IN) :: state
            
            REAL(DP), ALLOCATABLE :: T(:,:)

            INTEGER :: dim
            INTEGER :: i
            INTEGER, ALLOCATABLE :: map(:,:)
            
            dim = SIZE(vT)

            IF (dim == 3) THEN
                ALLOCATE(T(2,2))
                T = ZERO
                ALLOCATE(map(3,2))
                map = Voigtmap('2D')
                
                DO i = 1,2
                    T(map(i,1), map(i,2)) = vT(i)
                END DO

                IF (state .EQ. 1) THEN
                    T(map(3,1), map(3,2)) = vT(3)
                    T(map(3,2), map(3,1)) = vT(3)
                ELSE IF (state .EQ. 2) THEN
                    T(map(3,1), map(3,2)) = HALF * vT(3)
                    T(map(3,2), map(3,1)) = HALF * vT(3)
                END IF

            ELSE IF (dim == 6) THEN
                ALLOCATE(map(6,2))
                ALLOCATE(T(3,3))
                T = ZERO
                map = Voigtmap('3D')

                DO i = 1,3
                    T(map(i,1), map(i,2)) = vT(i)
                END DO

                DO i = 4,6
                    IF (state .EQ. 1) THEN
                        T(map(i,1), map(i,2)) = vT(i)
                        T(map(i,2), map(i,1)) = vT(i)
                    ELSE IF (state .EQ. 2) THEN
                        T(map(i,1), map(i,2)) = HALF * vT(i)
                        T(map(i,2), map(i,1)) = HALF * vT(i)
                    END IF
                END DO
            END IF

            DEALLOCATE(map)
            
        END FUNCTION Voigt2Tensor 

        PURE FUNCTION Identity(n) RESULT(Id)
            ! Purpose:
            !   Creates an Identity matrix of size n
            INTEGER, INTENT(IN) :: n

            REAL(DP) :: Id(n,n)

            INTEGER :: i

            Id = ZERO

            DO i = 1, n
                Id(i,i) = ONE
            END DO

        END FUNCTION Identity

        PURE FUNCTION Tensor2Output(dim, T) RESULT(oT)
            INTEGER,  INTENT(IN) :: dim
            REAL(DP), INTENT(IN) :: T(dim, dim)

            REAL(DP) :: oT(6)

            oT(:) = ZERO

            IF (dim .EQ. 2) THEN
                oT(1) = T(1,1)
                oT(2) = T(2,2)
                oT(6) = T(1,2)

            ELSE IF (dim .EQ. 3) THEN
                oT(1) = T(1,1)
                oT(2) = T(2,2)
                oT(3) = T(3,3)
                oT(4) = T(2,3)
                oT(5) = T(1,3)
                oT(6) = T(1,2)
            END IF

        END FUNCTION Tensor2Output

        SUBROUTINE inv(A, Ainv, stat)
            ! Calculates the inverse of an n,n matrix
            REAL(DP),  INTENT(IN) :: A(:,:)
            REAL(DP), INTENT(OUT) :: Ainv(SIZE(A,1), SIZE(A,1))
            INTEGER,  INTENT(OUT) :: stat

            REAL(DP) :: det

            INTEGER :: n, info

            INTEGER :: ipiv(SIZE(A,1))

            EXTERNAL dgetrf
            EXTERNAL dgetri

            stat = 0
            det = ZERO

            Ainv = A
            n = SIZE(A,1)

            IF (n .EQ. 2) THEN
                det = A(1,1) * A(2,2) -  A(1,2) * A(2,1)
                IF (det .EQ. ZERO) THEN
                    info = 1
                    RETURN
                END IF
                Ainv = Inverse2(A, det)
            ELSE IF (n .EQ. 3) THEN
                det =   A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)     &
                    &  - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)    &
                    &  + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
                IF (det .EQ. ZERO) THEN
                    info = 1
                    RETURN
                END IF
                Ainv = Inverse3(A, det)
            ELSE 
                ! Compute the LU decomposition  of a general m by n matrix
                !   using partial pivoting with row exchange method
                CALL dgetrf(n, n, Ainv, n, ipiv, info)
                
                IF (info .NE. 0) THEN
                    stat = 1
                    RETURN
                END IF

                ! Compute the inverse of the matrix using LU factorization 

                CALL dgetri(n, Ainv, n, ipiv, n, info)

                IF (info .NE. 0) THEN
                    stat = 2
                    RETURN
                END IF
            END IF
            
        END SUBROUTINE inv

        PURE FUNCTION Identity4s(dim) RESULT(Id4s)
            !Purpose:
            ! The function calculates the 4th order symmetric identity matrix 
            ! in Voigt notation given by:
            !       I_ijkl = 1/2 (delta_ik * delta_jl + delta_il * delta_jk)
            !

            INTEGER, INTENT(IN) :: dim
            
            REAL(DP) :: Id4s(dim * (dim + 1)/2, dim * (dim + 1)/2)

            INTEGER, ALLOCATABLE :: map(:,:)
            INTEGER :: i, j

            Id4s = ZERO

            IF (dim .NE. 2 .AND. dim .NE. 3) RETURN

            IF (dim .EQ. 2) THEN
                ALLOCATE(map(3,2))
                map = Voigtmap('2D')
            ELSE IF (dim .EQ. 3) THEN
                ALLOCATE(map(6,2))
                map = Voigtmap('3D')
            END IF

            DO i = 1, SIZE(Id4s, 1)
                DO j = 1, SIZE(Id4s, 2)

                    Id4s(i,j) = HALF * ( KDelta(map(i,1), map(j,1))) * &
                    &           KDelta(map(i,2), map(j,2)) +           &
                    &           KDelta(map(i,1), map(j,2)) *           &
                    &           KDelta(map(i,2), map(j,1))       
                    
                END DO
            END DO

            DEALLOCATE(map)
        
        END FUNCTION Identity4s

        PURE FUNCTION IdoxId(dim) RESULT(mat)
            !Purpose:
            ! The function calculates the tensor product of two identity matrix
            ! in voigt notation:
            !       I x I = delta_ij * delta_kl
            !

            INTEGER, INTENT(IN) :: dim

            REAL(DP) :: mat(dim * (dim + 1)/2, dim * (dim + 1)/2)

            INTEGER, ALLOCATABLE :: map(:,:)
            INTEGER :: i, j

            mat = ZERO

            IF (dim .NE. 2 .AND. dim .NE. 3) RETURN

            IF (dim .EQ. 2) THEN
                ALLOCATE(map(3,2))
                map = Voigtmap('2D')
            ELSE IF (dim .EQ. 3) THEN
                ALLOCATE(map(6,2))
                map = Voigtmap('3D')
            END IF

            DO i = 1, SIZE(mat, 1)
                DO j = 1, SIZE(mat, 2)

                    mat(i,j) = KDelta(map(i,1), map(i,2)) * &   
                    &          KDelta(map(j,1), map(j,2))
                    
                END DO
            END DO
            
            DEALLOCATE(map)

        END FUNCTION IdoxId

        PURE FUNCTION TensorCrossProduct(A, B) RESULT(C)
            !Purpose:
            ! The function takes in two 2D tensor in matrix notaion and return the
            ! tensor product in a Voigt notation  

            REAL(DP), INTENT(IN) :: A(:,:)
            REAL(DP), INTENT(IN) :: B(:,:)

            REAL(DP) :: C(SIZE(A,1)*(SIZE(A,1) + 1)/2, SIZE(B,1)*(SIZE(B,1) + 1)/2)

            REAL(DP) :: Av(SIZE(A,1)*(SIZE(A,1) + 1)/2, 1)
            REAL(DP) :: Bv(1, SIZE(B,1)*(SIZE(B,1) + 1)/2)

            Av(:,1) = Voigt2(A, 1)
            Bv(1,:) = Voigt2(B, 1)

            C = MATMUL(Av, Bv)

        END FUNCTION TensorCrossProduct

END MODULE utility_module
