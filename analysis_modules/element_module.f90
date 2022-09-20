
MODULE element_module 
    !Purpose:
    ! The module stores the ESFEM element object
    
    USE EdgeSFEM,         ONLY: calc_SF_Derivative_FE, calc_SF_Derivative_IE
    USE parameter_module, ONLY: DP, ONE, NDIM, HALF, ZERO, NDIM_STRESS, &
                            &   NDIM_STRESS_Voigt, NDIM_Voigt
    USE material_module,  ONLY: IG_point
    ! USE material_module,  ONLY: IG_point
    
    IMPLICIT NONE 

    TYPE, PUBLIC :: ESFEM_elem
        PRIVATE
        INTEGER               :: Nnodes 
        REAL(DP), ALLOCATABLE :: u(:)
        REAL(DP), ALLOCATABLE :: du(:) 
        REAL(DP), ALLOCATABLE :: dNdX0(:,:)      
        REAL(DP), ALLOCATABLE :: dNdx(:,:)      
        REAL(DP)              :: F_n(NDIM,NDIM)
        REAL(DP)              :: F_np1(NDIM,NDIM)
        REAL(DP)              :: Q_np1(NDIM_STRESS,NDIM_STRESS)
        REAL(DP)              :: dtD(NDIM_STRESS_Voigt)
        REAL(DP)              :: vol
        REAL(DP)              :: j
        REAL(DP)              :: rot
        TYPE(IG_point)        :: igp
    END TYPE

    INTERFACE update
        MODULE PROCEDURE update_ESFEM_elem
    END INTERFACE update
    
    INTERFACE display
        MODULE PROCEDURE display_ESFEM_elem
    END INTERFACE display
    
    INTERFACE extract
        MODULE PROCEDURE extract_ESFEM_elem
    END INTERFACE extract

    INTERFACE integrate
       MODULE PROCEDURE integrate_ESFEM_elem
    END INTERFACE integrate

    PUBLIC :: update, extract, display, integrate

    CONTAINS

        SUBROUTINE update_ESFEM_elem(elem, NodesCoord, u, du, igp, info)
            !Purpose:
            ! Updates the element 

            USE utility_module, ONLY: Identity, inv, Exp2ss, Voigt2

            TYPE(ESFEM_elem),         INTENT(INOUT) :: elem
            REAL(DP), DIMENSION(:,:), INTENT(IN) :: NodesCoord
            REAL(DP),   DIMENSION(:), INTENT(IN) :: u
            REAL(DP),   DIMENSION(:), INTENT(IN) :: du
            TYPE(IG_point),           INTENT(IN) :: igp
            INTEGER,                 INTENT(OUT) :: info
            
            ! Element variables 
            ! Nnodes : Number of nodes associated to the Smoothed region
            ! du     : Increment in displacement stored in 1D array 
            !            [du_1x, du_1y, du_2x, du_2y .....]
            ! du2    : Increment in displacement stored in 2D array with separated x and y column
            !               [[du_1x, du_1y]
            !                [du_2x, du_2y]...]
            ! dNdX0  : Derivative of the shape function w.r.t. material coordinate frame, 
            !           2D array (Nnodes, NDIM)
            ! dNdx   : Derivative of the shape function w.r.t. current coordinate frame
            !           2D array (Nnodes, NDIM)
            ! F      : Deformation gradient 
            ! dtD    : Incremental rate of deformation
            ! vol    : Volume of the element
            ! jdet   : Deformation Jacobian for the element 

            INTEGER               :: Nnodes
            REAL(DP), ALLOCATABLE :: du2(:,:)
            REAL(DP), ALLOCATABLE :: dNdX0(:,:)
            REAL(DP), ALLOCATABLE :: dNdx(:,:)
            REAL(DP)              :: F_n(NDIM, NDIM)
            REAL(DP)              :: H(NDIM, NDIM)
            REAL(DP)              :: invF_n(NDIM, NDIM)
            REAL(DP)              :: F_np1(NDIM,NDIM)
            REAL(DP)              :: invF_np1(NDIM,NDIM)
            REAL(DP)              :: delF(NDIM,NDIM)
            REAL(DP)              :: invdelF(NDIM,NDIM)
            REAL(DP)              :: del0du(NDIM,NDIM)
            REAL(DP)              :: delndu(NDIM,NDIM)
            REAL(DP)              :: dtW(NDIM,NDIM)
            REAL(DP)              :: dtD(NDIM_STRESS,NDIM_STRESS)
            REAL(DP)              :: Q_np1(NDIM_STRESS,NDIM_STRESS)
            REAL(DP)              :: Vsk
            REAL(DP)              :: jdet
        
            INTEGER               :: i  ! Loop variable for the loop
            REAL(DP)              :: Id(2,2) ! Variable for identity matrix

            Nnodes = SIZE(NodesCoord, 1)

            ALLOCATE(dNdx(Nnodes, 2))
            ALLOCATE(dNdX0(Nnodes, 2))
            ALLOCATE(du2(Nnodes, 2))

            F_n       = ZERO
            invF_n    = ZERO
            F_np1     = ZERO
            invF_np1  = ZERO
            delF      = ZERO
            invdelF   = ZERO
            del0du    = ZERO
            delndu    = ZERO
            jdet      = ZERO
            dtW       = ZERO
            dtD       = ZERO
            Q_np1     = Identity(NDIM_STRESS)

            IF (Nnodes == 3) THEN
                ! For free edge
                CALL calc_SF_Derivative_FE(NodesCoord, dNdX0, Vsk)

            ELSE IF (Nnodes == 4) THEN
                ! For inside edge 
                CALL calc_SF_Derivative_IE(NodesCoord, dNdX0, Vsk)
            END IF
            
            DO i = 1, Nnodes
                H(1,1) = H(1,1) + dNdX0(i,1) * u(2*i-1)
                H(1,2) = H(1,2) + dNdX0(i,2) * u(2*i-1)
                H(2,1) = H(2,1) + dNdX0(i,1) * u(2*i)
                H(2,2) = H(2,2) + dNdX0(i,2) * u(2*i)
                du2(i,1) = du(2*i-1)
                du2(i,2) = du(2*i)
            END DO

            ! Identity matrix
            Id = Identity(NDIM)

            F_n = Id + H
            
            jdet = F_n(1,1)*F_n(2,2) - F_n(1,2)*F_n(2,1)

            CALL inv(F_n, invF_n, info)
            IF (info .NE. 0) RETURN

            dNdx = MATMUL(dNdX0, invF_n)

            ! Calculating derivative of incremental displacement : 
            ! In current domain 
            delndu = MATMUL(TRANSPOSE(dNdx), du2)
            ! In the material domain
            del0du = MATMUL(TRANSPOSE(dNdX0), du2)

            ! Calculating increment in the deformation gradient 
            delF = Id + TRANSPOSE(delndu)

            ! Deformation gradient in the next step
            F_np1 = MATMUL(delF, F_n)

            ! Inverse of deformation gradient for the next step
            CALL inv(F_np1, invF_np1, info)
            IF (info .NE. 0) RETURN

            ! Inverse of delF
            CALL inv(delF, invdelF, info)
            IF (info .NE. 0) RETURN

            ! Calculating effective spin tensor 
            dtW = HALF * (MATMUL(TRANSPOSE(del0du), invF_np1) - &
                    &    MATMUL(TRANSPOSE(invF_np1), del0du))

            CALL Exp2ss(dtW, Q_np1(1:NDIM,1:NDIM), info)

            dtD(1:NDIM, 1:NDIM) = HALF * (Id - &
                            &     MATMUL(TRANSPOSE(invdelF), invdelF))

            IF (info .NE. 0) RETURN

            IF (.NOT. ALLOCATED(elem%u))    ALLOCATE(elem%u(SIZE(u)))
            IF (.NOT. ALLOCATED(elem%du))   ALLOCATE(elem%du(SIZE(du)))
            IF (.NOT. ALLOCATED(elem%dNdX0)) THEN
                ALLOCATE(elem%dNdX0(SIZE(dNdX0,1), SIZE(dNdX0,2)))
            END IF
            IF (.NOT. ALLOCATED(elem%dNdx)) THEN
                ALLOCATE(elem%dNdx(SIZE(dNdx,1), SIZE(dNdx,2)))
            END IF

            elem%Nnodes = Nnodes
            elem%u      = u
            elem%du     = du
            elem%dNdX0  = dNdX0
            elem%dNdx   = dNdx
            elem%F_n    = F_n
            elem%F_np1  = F_np1
            elem%dtD    = Voigt2(dtD,2)
            elem%Q_np1  = Q_np1
            elem%vol    = Vsk
            elem%j      = jdet
            elem%rot    = dtW(1,2)
            elem%igp    = igp

        END SUBROUTINE update_ESFEM_elem

        PURE SUBROUTINE extract_ESFEM_elem(elem, Nnodes, u, du, dNdX0, dNdx, F_n, F_np1, &
                                        &  Q_np1, dtD, vol, jdet, rot,igp)
            !Purpose to extract the ESFEM element 

            TYPE(ESFEM_elem),          INTENT(IN) :: elem
            INTEGER,        OPTIONAL, INTENT(OUT) :: Nnodes
            REAL(DP),       OPTIONAL, INTENT(OUT) :: u(:)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: du(:)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: dNdX0(:,:)        
            REAL(DP),       OPTIONAL, INTENT(OUT) :: dNdX(:,:)        
            REAL(DP),       OPTIONAL, INTENT(OUT) :: F_n(NDIM,NDIM)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: F_np1(NDIM,NDIM)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: Q_np1(NDIM_STRESS,NDIM_STRESS)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: dtD(NDIM_STRESS_Voigt)
            REAL(DP),       OPTIONAL, INTENT(OUT) :: vol
            REAL(DP),       OPTIONAL, INTENT(OUT) :: jdet
            REAL(DP),       OPTIONAL, INTENT(OUT) :: rot
            TYPE(IG_point), OPTIONAL, INTENT(OUT) :: igp

            IF (PRESENT(Nnodes)) Nnodes = elem%Nnodes

            IF (PRESENT(u)) THEN
                IF (ALLOCATED(elem%u)) u = elem%u
            END IF

            IF (PRESENT(du)) THEN
                IF (ALLOCATED(elem%du)) du = elem%du
            END IF

            IF (PRESENT(dNdX0)) THEN
                IF (ALLOCATED(elem%dNdX0)) dNdX0 = elem%dNdX0
            END IF

            IF (PRESENT(dNdX0)) THEN
                IF (ALLOCATED(elem%dNdX)) dNdx = elem%dNdx
            END IF

            IF (PRESENT(F_n))   F_n   = elem%F_n
            IF (PRESENT(F_np1)) F_np1 = elem%F_np1
            IF (PRESENT(dtD))   dtD   = elem%dtD
            IF (PRESENT(Q_np1)) Q_np1 = elem%Q_np1
            IF (PRESENT(vol))   vol   = elem%vol
            IF (PRESENT(jdet))  jdet  = elem%j
            IF (PRESENT(rot))   rot   = elem%rot
            IF (PRESENT(igp))   igp   = elem%igp
            
        END SUBROUTINE extract_ESFEM_elem

        SUBROUTINE integrate_ESFEM_elem(elem, K, fvec, sf)
            !Purpose:
            ! Integrates the ESFEM element  
            USE material_module,    ONLY: integrate, extract
            USE stress_module,      ONLY: Stress, extract, display
            USE utility_module,     ONLY: Voigt2, PrintMatrix
        
            TYPE(ESFEM_elem), INTENT(INOUT) :: elem
            REAL(DP),           INTENT(OUT) :: K(:,:)
            REAL(DP),           INTENT(OUT) :: fvec(:)
            CHARACTER(len=1),    INTENT(IN) :: sf
            !Internal Variables

            TYPE(IG_point)           :: igp
            TYPE(Stress)             :: sigma
            REAL(DP), ALLOCATABLE    :: Bgeo(:,:)
            REAL(DP), ALLOCATABLE    :: B(:,:)
            REAL(DP), ALLOCATABLE    :: dNdx(:,:) 
            REAL(DP), ALLOCATABLE    :: Kgeo(:,:)
            REAL(DP), ALLOCATABLE    :: Kmat(:,:)
            REAL(DP)                 :: S(NDIM_STRESS, NDIM_STRESS)
            REAL(DP)                 :: SI(2*NDIM, 2*NDIM)
            REAL(DP)                 :: Cep(NDIM_Voigt, NDIM_Voigt)
            
            REAL(DP) :: jdet, vol              
            INTEGER  :: Nnodes
            INTEGER  :: i

            Nnodes=elem%Nnodes
            igp = elem%igp
            jdet = elem%j
            vol = elem%vol

            ALLOCATE(dNdx(Nnodes, 2))
            !Allocating B matrices size
            ALLOCATE(Bgeo(4, NDIM * Nnodes))
            ALLOCATE(B(3, NDIM * Nnodes))
            ! Allocating K matrices
            ALLOCATE(Kgeo(Nnodes*NDIM, Nnodes*NDIM))
            ALLOCATE(Kmat(Nnodes*NDIM, Nnodes*NDIM))

            dNdx = elem%dNdx

            CALL integrate(igp, Cep, elem%Q_np1, elem%dtD, elem%j, &
                            &   algo='I', sf=sf, cj2t='Y')


            CALL extract(igp, sigma=sigma)
            CALL extract(sigma, stress_2d=S)

            ! Calculating the B matrix
            Bgeo(:,:) = ZERO

            DO i = 1, Nnodes
                
                Bgeo(1, 2*i-1) = dNdx(i,1)
                Bgeo(2, 2*i-1) = dNdx(i,2)

                Bgeo(3, 2*i) = dNdx(i,1)
                Bgeo(4, 2*i) = dNdx(i,2)

            END DO

            B(:,:) = ZERO

            DO i = 1, Nnodes
                
                B(1, 2*i-1) = dNdx(i,1)
                B(3, 2*i-1) = dNdx(i,2)

                B(2, 2*i) = dNdx(i,2)
                B(3, 2*i) = dNdx(i,1)

            END DO

            SI = ZERO
            SI(1:NDIM, 1:NDIM) = S(1:NDIM, 1:NDIM)
            SI(NDIM+1:2*NDIM, NDIM+1:2*NDIM) = S(1:NDIM, 1:NDIM)

            IF (sf .EQ. 'K') THEN
                Kmat = MATMUL(MATMUL(TRANSPOSE(B), Cep), B) * vol * jdet 
                Kgeo = MATMUL(MATMUL(TRANSPOSE(Bgeo), SI), Bgeo) * vol
                fvec = MATMUL(TRANSPOSE(B), Voigt2(S(1:NDIM,1:NDIM),1)) * vol
                
            ELSE IF (sf .EQ. 'C') THEN
                Kmat = MATMUL(MATMUL(TRANSPOSE(B), Cep), B) * vol * jdet 
                Kgeo = MATMUL(MATMUL(TRANSPOSE(Bgeo), SI), Bgeo) * vol * jdet
                fvec = MATMUL(TRANSPOSE(B), Voigt2(S(1:NDIM,1:NDIM),1)) * vol * jdet
            END IF

            K = Kmat + Kgeo
            elem%igp = igp
            
        END SUBROUTINE integrate_ESFEM_elem 

        SUBROUTINE display_ESFEM_elem(elem, elemno, inc, unit)
            ! Purpose:
            ! Displays the element object for debugging

            TYPE(ESFEM_elem),  INTENT(IN) :: elem
            INTEGER,           INTENT(IN) :: elemno
            INTEGER,           INTENT(IN) :: inc
            INTEGER, OPTIONAL, INTENT(IN) :: unit

            INTEGER :: i, j

            IF (PRESENT(unit)) THEN

                WRITE(unit,*) ''
                WRITE(unit,*) '****************ELEMENT DISPLAY INFORMATION****************'
                WRITE(unit,*) '==========================================================='
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A,I5,A,I5)') 'Element Number : ', elemno,      & 
                                    &    'for increment number : ', inc        
                        
                WRITE(unit,'(2X,A,I2)') 'Number of nodes in the element is : ', &
                                    &   elem%Nnodes
                WRITE(unit,'(1X, A)')''
                WRITE(unit,'(2X,A)') 'Displacement of nodes at the end of last inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(unit, '(5X,A,I2,A,ES10.3,A,ES10.3)')                          &
                                    &  'Displacement for node', i, ' in X dir is ',     &
                                    &  elem%u(2*i-1), ' and in Y dir is ', elem%u(2*i)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Displacement increment of nodes for current inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(unit,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Increment for node', i,     &
                        & ' in X dir is ', elem%du(2*i-1), ' and in Y dir is ', elem%du(2*i)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Derivative of shape function w.r.t. material reference frame : '
                DO i = 1, elem%Nnodes
                    WRITE(unit,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Derivative for Shape function', i, & 
                            &   ' wrt X is ', elem%dNdX0(i,1), ' and in Y dir is ', elem%dNdX0(i,2)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Derivative of shape function w.r.t. last converged inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(unit,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Derivative for Shape function', i, &
                        & ' wrt x is ', elem%dNdx(i,1), ' and in y dir is ', elem%dNdx(i,2)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Deformation gradient of last converged inc.: '
                DO i = 1,NDIM
                    WRITE(unit, '(5X, 2ES15.7)') (elem%F_n(i,j), j=1,NDIM)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Deformation gradient after the end of current inc.: '
                DO i = 1,NDIM
                    WRITE(unit, '(5X, 2ES15.7)') (elem%F_np1(i,j), j=1,NDIM)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Incremental rate of deformation : '
                DO i = 1,NDIM_STRESS
                    WRITE(unit, '(5X, 2ES15.7)') elem%dtD(i)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(2X, A)') 'Rotation tensor : '
                DO i = 1,NDIM_STRESS
                    WRITE(unit, '(5X, 2ES15.7)') (elem%Q_np1(i,j), j=1,NDIM)
                END DO
                WRITE(unit,'(1X, A)')''

                WRITE(unit,'(1X, A, ES10.3)')'The volume of Smoothed domain : ', elem%vol

                WRITE(unit,*) 'The Deformation Jacobian for the element : ', elem%j
                WRITE(unit,*) ''

                WRITE(unit,*) 'Rotation angle of the element : ', elem%rot
                WRITE(unit,*) ''

                WRITE(unit,*) '==========================================================='

            ELSE
                WRITE(*,*) ''
                WRITE(*,*) '****************ELEMENT DISPLAY INFORMATION****************'
                WRITE(*,*) '==========================================================='
                WRITE(*,*) ''

                WRITE(*,'(2X,A,I5,A,I5)') 'Element Number : ', elemno,      & 
                                    &    'for increment number : ', inc      
                        
                WRITE(*,'(2X,A,I2)') 'Number of nodes in the element is : ', &
                                    &   elem%Nnodes
                WRITE(*,'(1X, A)')''
                WRITE(*,'(2X,A)') 'Displacement of nodes at the end of last inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(*, '(5X,A,I2,A,ES10.3,A,ES10.3)')                          &
                                    &  'Displacement for node', i, ' in X dir is ',     &
                                    &  elem%u(2*i-1), ' and in Y dir is ', elem%u(2*i)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Displacement increment of nodes for current inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(*,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Increment for node', i,     &
                        & ' in X dir is ', elem%du(2*i-1), ' and in Y dir is ', elem%du(2*i)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Derivative of shape function w.r.t. material reference frame : '
                DO i = 1, elem%Nnodes
                    WRITE(*,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Derivative for Shape function', i, & 
                            &   ' wrt X is ', elem%dNdX0(i,1), ' and in Y dir is ', elem%dNdX0(i,2)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Derivative of shape function w.r.t. last converged inc. : '
                DO i = 1, elem%Nnodes
                    WRITE(*,'(5X,A,I2,A,ES10.3,A,ES10.3)') 'Derivative for Shape function', i, &
                        & ' wrt x is ', elem%dNdx(i,1), ' and in y dir is ', elem%dNdx(i,2)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Deformation gradient of last converged inc.: '
                DO i = 1,NDIM
                    WRITE(*, '(5X, 2ES15.7)') (elem%F_n(i,j), j=1,NDIM)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Deformation gradient after the end of current inc.: '
                DO i = 1,NDIM
                    WRITE(*, '(5X, 2ES15.7)') (elem%F_np1(i,j), j=1,NDIM)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Incremental rate of deformation : '
                DO i = 1,NDIM_STRESS
                    WRITE(*, '(5X, 2ES15.7)') elem%dtD(i)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(2X, A)') 'Rotation tensor : '
                DO i = 1,NDIM_STRESS
                    WRITE(*, '(5X, 2ES15.7)') (elem%Q_np1(i,j), j=1,NDIM)
                END DO
                WRITE(*,'(1X, A)')''

                WRITE(*,'(1X, A, ES10.3)')'The volume of Smoothed domain : ', elem%vol

                WRITE(*,*) 'The Deformation Jacobian for the element : ', elem%j
                WRITE(*,*) ''

                WRITE(*,*) 'Rotation angle of the element : ', elem%rot
                WRITE(*,*) ''

                WRITE(*,*) '==========================================================='

            END IF
        END SUBROUTINE display_ESFEM_elem
        
END MODULE element_module
