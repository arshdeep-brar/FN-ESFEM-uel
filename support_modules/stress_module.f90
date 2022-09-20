MODULE stress_module
    
    USE parameter_module,   ONLY: DP, ZERO, HALF, ONE_THIRD, ONE, TWENTYSEVEN, &
                            &     TWO, NDIM, NDIM_STRESS, NDIM_STRESS_Voigt,   &
                            &     ERROR, PI, THREE, NINE
    USE utility_module,     ONLY: Voigt2, Voigt2Tensor, Identity

    IMPLICIT NONE    
    
    PRIVATE

    TYPE, PUBLIC ::Stress
        PRIVATE
        REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)
        REAL(DP) :: stress_v(NDIM_STRESS_Voigt)
        REAL(DP) :: stress_VM
        REAL(DP) :: stress_mean
        REAL(DP) :: stress_deviatoric(NDIM_STRESS, NDIM_STRESS)
        REAL(DP) :: stress_triaxiality
        REAL(DP) :: lode_angle
        REAL(DP) :: lode_angle_parameter
    END TYPE

    INTERFACE update
        MODULE PROCEDURE update_stress
    END INTERFACE update

    INTERFACE extract
        MODULE PROCEDURE extract_stress
    END INTERFACE extract

    INTERFACE display
        MODULE PROCEDURE display_stress
    END INTERFACE display

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE real_to_stress
        MODULE PROCEDURE array2d_to_stress
        MODULE PROCEDURE array_to_stress
    END INTERFACE
    
    INTERFACE OPERATOR(+)
        MODULE PROCEDURE add_with_array2d
        MODULE PROCEDURE add_with_array
        MODULE PROCEDURE add_with_stress
    END INTERFACE

    INTERFACE OPERATOR(-)
        MODULE PROCEDURE sub_with_array2d
        MODULE PROCEDURE sub_with_array
        MODULE PROCEDURE sub_with_stress
    END INTERFACE
    
    PUBLIC :: ASSIGNMENT(=), update, extract, display,  &
            & OPERATOR(+), OPERATOR(-)

    CONTAINS

        SUBROUTINE update_stress(this, stress_2d, stress_v)
            USE utility_module, ONLY: PrintMatrix

            TYPE(Stress),    INTENT(INOUT) :: this
            REAL(DP), OPTIONAL, INTENT(IN) :: stress_2d(NDIM_STRESS, NDIM_STRESS)
            REAL(DP), OPTIONAL, INTENT(IN) :: stress_v(NDIM_STRESS_Voigt)

            REAL(DP) :: stress_VM2, stress_VM
            REAL(DP) :: stress_mean
            REAL(DP) :: J3, ratio
            REAL(DP) :: stress_principal(3)
            REAL(DP) :: stress_2dL(NDIM_STRESS, NDIM_STRESS)
            REAL(DP) :: stress_dev(NDIM_STRESS, NDIM_STRESS)

            INTEGER :: i

            ! Lapack subroutine to calculate eigen value of a 3x3 matrix
            EXTERNAL dsyev     
            
            stress_mean = ZERO

            IF (PRESENT(stress_2d)) THEN
                
                this%stress_2d = stress_2d
                this%stress_v = Voigt2(stress_2d, 1)

                DO i = 1, NDIM_STRESS
                    stress_mean = stress_mean + ONE_THIRD * stress_2d(i,i)
                END DO

                IF (NDIM_STRESS .EQ. 2) THEN
                    
                    CALL get_principal_stress(stress_2d, stress_principal)

                    stress_VM2 = HALF * ((stress_principal(1) - stress_principal(2))**2 + &
                                    &    (stress_principal(2) - stress_principal(3))**2 + &
                                    &    (stress_principal(3) - stress_principal(1))**2)
                    
                    J3 = (TWENTYSEVEN/TWO) * (stress_principal(1) - stress_mean) * &
                                         &   (stress_principal(2) - stress_mean) * & 
                                         &   (stress_principal(3) - stress_mean)
                    
                    stress_dev = stress_2d - stress_mean * Identity(NDIM_STRESS)

                ELSE IF (NDIM_STRESS .EQ. 3) THEN
                    
                    stress_dev = stress_2d - stress_mean * Identity(NDIM_STRESS)

                    stress_VM2 = (THREE / TWO) * SUM(stress_dev * stress_dev)

                    J3 = (NINE/TWO) * SUM(MATMUL(stress_dev, stress_dev) * &
                                        &  stress_dev)

                END IF

                stress_VM = SQRT(stress_VM2)

                this%stress_VM          = stress_VM
                this%stress_mean        = stress_mean
                this%stress_deviatoric  = stress_dev

                IF (stress_VM .NE. ZERO) THEN
                
                    ratio = J3 / stress_VM**3

                    IF (ABS(ratio - ONE) < ERROR) ratio = ONE
                    IF (ABS(ratio + ONE) < ERROR) ratio = -ONE

                    this%stress_triaxiality   = stress_mean / stress_VM
                    this%lode_angle           = ONE_THIRD * DACOS(ratio)
                    this%lode_angle_parameter = ONE - (TWO/PI) * DACOS(ratio)

                ELSE
                    this%stress_triaxiality   = ZERO
                    this%lode_angle           = ZERO
                    this%lode_angle_parameter = ZERO
                END IF

            ELSE IF (PRESENT(stress_v)) THEN
    
                stress_2dL = Voigt2Tensor(stress_v, 1)
                this%stress_v = stress_v
                
                this%stress_2d = stress_2dL

                DO i = 1, NDIM_STRESS
                    stress_mean = stress_mean + ONE_THIRD * stress_2dL(i,i)
                END DO

                IF (NDIM_STRESS .EQ. 2) THEN
                    
                    CALL get_principal_stress(stress_2dL, stress_principal)

                    stress_VM2 = HALF * ((stress_principal(1) - stress_principal(2))**2 + &
                                    &    (stress_principal(2) - stress_principal(3))**2 + &
                                    &    (stress_principal(3) - stress_principal(1))**2)
                    
                    J3 = (TWENTYSEVEN/TWO) * (stress_principal(1) - stress_mean) * &
                                         &   (stress_principal(2) - stress_mean) * & 
                                         &   (stress_principal(3) - stress_mean)

                    stress_dev = stress_2dL - stress_mean * Identity(NDIM_STRESS)
                
                ELSE IF (NDIM_STRESS .EQ. 3) THEN

                    stress_dev = stress_2dL - stress_mean * Identity(NDIM_STRESS)

                    stress_VM2 = (THREE / TWO) * SUM(stress_dev * stress_dev)

                    J3 = (NINE/TWO) * SUM(MATMUL(stress_dev, stress_dev) * &
                                        &  stress_dev)
                
                END IF

                stress_VM = SQRT(stress_VM2)

                this%stress_VM          = stress_VM
                this%stress_mean        = stress_mean
                this%stress_deviatoric  = stress_dev

                IF (stress_VM .NE. ZERO) THEN
                     
                    ratio = J3 / (stress_VM**3)

                    IF (ABS(ratio - ONE) < ERROR) ratio = ONE
                    IF (ABS(ratio + ONE) < ERROR) ratio = -ONE

                    this%stress_triaxiality   = stress_mean / stress_VM
                    this%lode_angle           = ONE_THIRD * DACOS(ratio)
                    this%lode_angle_parameter = ONE - (TWO/PI) * DACOS(ratio)

                ELSE
                    this%stress_triaxiality   = ZERO
                    this%lode_angle           = ZERO
                    this%lode_angle_parameter = ZERO
                END IF

            ELSE
                this%stress_2d            = ZERO
                this%stress_v             = ZERO
                this%stress_VM            = ZERO
                this%stress_mean          = ZERO
                this%stress_deviatoric    = ZERO
                this%stress_triaxiality   = ZERO
                this%lode_angle           = ZERO
                this%lode_angle_parameter = ZERO
            END IF

        END SUBROUTINE update_stress

        SUBROUTINE real_to_stress(this, realnum)
            
            TYPE(Stress), INTENT(OUT) :: this
            REAL(DP),      INTENT(IN) :: realnum

            this%stress_2d            = realnum
            this%stress_v             = realnum
            this%stress_VM            = realnum
            this%stress_mean          = realnum
            this%stress_deviatoric    = realnum
            this%stress_triaxiality   = realnum
            this%lode_angle           = realnum
            this%lode_angle_parameter = realnum

        END SUBROUTINE real_to_stress

        SUBROUTINE array2d_to_stress(this, mat)
            TYPE(Stress), INTENT(OUT) :: this
            REAL(DP),      INTENT(IN) :: mat(NDIM_STRESS, NDIM_STRESS)

            CALL update(this, stress_2d=mat)
            
        END SUBROUTINE array2d_to_stress

        SUBROUTINE array_to_stress(this, vstress)
            TYPE(Stress), INTENT(OUT) :: this
            REAL(DP),      INTENT(IN) :: vstress(NDIM_STRESS_Voigt)

            CALL update(this, stress_v=vstress)
        
        END SUBROUTINE array_to_stress
        
        FUNCTION add_with_array2d(stress1, mat) RESULT(stress2)
            
            TYPE(Stress), INTENT(IN) :: stress1 
            REAL(DP),     INTENT(IN) :: mat(NDIM_STRESS, NDIM_STRESS)

            TYPE(Stress) :: stress2

            REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)

            stress_2d = stress1%stress_2d + mat

            CALL update(stress2, stress_2d=stress_2d)

        END FUNCTION add_with_array2d

        FUNCTION add_with_array(stress1, vstress) RESULT(stress2)
            
            TYPE(Stress), INTENT(IN) :: stress1 
            REAL(DP),     INTENT(IN) :: vstress(NDIM_STRESS_Voigt)

            TYPE(Stress) :: stress2

            REAL(DP) :: stress_v(NDIM_STRESS_Voigt)

            stress_v = stress1%stress_v + vstress

            CALL update(stress2, stress_v=stress_v)

        END FUNCTION add_with_array

        FUNCTION add_with_stress(stress1, stress2) RESULT(stress3)

            TYPE(Stress), INTENT(IN) :: stress1
            TYPE(Stress), INTENT(IN) :: stress2

            TYPE(Stress) :: stress3

            REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)

            stress_2d = stress1%stress_2d + stress2%stress_2d

            CALL update(stress3, stress_2d=stress_2d)

        END FUNCTION add_with_stress

        FUNCTION sub_with_array2d(stress1, mat) RESULT(stress2)
            
            TYPE(Stress), INTENT(IN) :: stress1 
            REAL(DP),     INTENT(IN) :: mat(NDIM_STRESS, NDIM_STRESS)

            TYPE(Stress) :: stress2

            REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)

            stress_2d = stress1%stress_2d - mat

            CALL update(stress2, stress_2d=stress_2d)

        END FUNCTION sub_with_array2d

        FUNCTION sub_with_array(stress1, vstress) RESULT(stress2)
            
            TYPE(Stress), INTENT(IN) :: stress1 
            REAL(DP),     INTENT(IN) :: vstress(NDIM_STRESS_Voigt)

            TYPE(Stress) :: stress2

            REAL(DP) :: stress_v(NDIM_STRESS_Voigt)

            stress_v = stress1%stress_v - vstress

            CALL update(stress2, stress_v=stress_v)

        END FUNCTION sub_with_array

        FUNCTION sub_with_stress(stress1, stress2) RESULT(stress3)

            TYPE(Stress), INTENT(IN) :: stress1
            TYPE(Stress), INTENT(IN) :: stress2

            TYPE(Stress) :: stress3

            REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)

            stress_2d = stress1%stress_2d - stress2%stress_2d

            CALL update(stress3, stress_2d=stress_2d)

        END FUNCTION sub_with_stress

        SUBROUTINE extract_stress(this, stress_2d, stress_v, stress_principal, stress_VM,  &
                                &      stress_mean, stress_deviatoric, stress_triaxiality, &
                                &      lode_angle, lode_angle_parameter, info)
            
            TYPE(Stress),        INTENT(IN) :: this
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_2d(NDIM_STRESS, NDIM_STRESS)
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_v(NDIM_STRESS_Voigt)
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_principal(3)
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_VM
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_mean
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_deviatoric(NDIM_STRESS, NDIM_STRESS)
            REAL(DP), OPTIONAL, INTENT(OUT) :: stress_triaxiality
            REAL(DP), OPTIONAL, INTENT(OUT) :: lode_angle
            REAL(DP), OPTIONAL, INTENT(OUT) :: lode_angle_parameter
            INTEGER,  OPTIONAL, INTENT(OUT) :: info

            INTEGER :: info_stat

            info_stat = 0

            IF (PRESENT(stress_2d)) stress_2d = this%stress_2d
            IF (PRESENT(stress_v)) stress_v = this%stress_v
            IF (PRESENT(stress_principal)) THEN
                CALL get_principal_stress(this%stress_2d, stress_principal, info=info_stat)
            END IF
            IF (PRESENT(stress_VM)) stress_VM = this%stress_VM
            IF (PRESENT(stress_mean)) stress_mean = this%stress_mean
            IF (PRESENT(stress_deviatoric)) stress_deviatoric = this%stress_deviatoric
            IF (PRESENT(stress_triaxiality)) stress_triaxiality = this%stress_triaxiality
            IF (PRESENT(lode_angle)) lode_angle = this%lode_angle
            IF (PRESENT(lode_angle_parameter)) lode_angle_parameter = &
                                                            &   this%lode_angle_parameter
            IF (PRESENT(info)) info = info_stat
    
        END SUBROUTINE extract_stress

        SUBROUTINE get_principal_stress(this, stress_principal, info)
            
            REAL(DP),            INTENT(IN) :: this(NDIM_STRESS, NDIM_STRESS)
            REAL(DP),           INTENT(OUT) :: stress_principal(3)
            INTEGER,  OPTIONAL, INTENT(OUT) :: info
 
            REAL(DP) :: stress_2d(NDIM_STRESS, NDIM_STRESS)
            REAL(DP) :: stress_1p2, stress_1m2
            INTEGER  :: lwork, info_stat

            REAL(DP), ALLOCATABLE :: work(:)

            lwork = 3*NDIM_STRESS - 1
            ALLOCATE(work(lwork))
            
            info_stat = 0

            stress_principal = ZERO

            stress_2d = this

            stress_1p2 = HALF * (stress_2d(1,1) + stress_2d(2,2))
            stress_1m2 = HALF * (stress_2d(1,1) - stress_2d(2,2))

            IF (NDIM_STRESS .EQ. 2) THEN
                
                stress_principal(1) = stress_1p2 - SQRT(stress_1m2**2 + &
                                            &   stress_2d(1,2)**2) 
                stress_principal(2) = stress_1p2 + SQRT(stress_1m2**2 + &
                                            &   stress_2d(1,2)**2)
                stress_principal(3) = ZERO

            ELSE IF (NDIM_STRESS .EQ. 3) THEN
                CALL dsyev('N', 'U', NDIM_STRESS, stress_2d, NDIM_STRESS, &
                            &   stress_principal, work, lwork, info_stat)
            END IF

            IF (PRESENT(info)) info = info_stat


        END SUBROUTINE get_principal_stress

        SUBROUTINE display_stress(this, unit)

            USE utility_module, ONLY: PrintMatrix

            TYPE(Stress),      INTENT(IN) :: this
            INTEGER, OPTIONAL, INTENT(IN) :: unit

            INTEGER :: i

            REAL(DP) :: stress_principal(3)
            INTEGER  :: info

            CALL get_principal_stress(this%stress_2d, stress_principal, info=info)

            IF (PRESENT(unit)) THEN

                WRITE(unit,*) ''
                WRITE(unit,*) '***************STRESS DISPLAY INFORMATION***************'
                WRITE(unit,*) '========================================================'
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Stress 2D : '
                CALL PrintMatrix(realMat=this%stress_2d, unit=unit)
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Stress Voigt : '
                DO i = 1, NDIM_STRESS_Voigt
                    WRITE(unit,'(5X,3ES15.7)') this%stress_v(i)
                END DO
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Principal Stresses : '
                DO i = 1, 3
                    WRITE(unit,'(5X,3ES15.7)') stress_principal(i)
                END DO
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Von mises stress/Effective stress  : '
                WRITE(unit,'(5X,3ES15.7)') this%stress_VM
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Mean Stress : '
                WRITE(unit,'(5X,3ES15.7)') this%stress_mean
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Deviatoric stress : '
                CALL PrintMatrix(realMat=this%stress_deviatoric, unit=unit)
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Stress Triaxiality : '
                WRITE(unit,'(5X,3ES15.7)') this%stress_triaxiality
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Lode angle : '
                WRITE(unit,'(5X,3ES15.7)') this%lode_angle
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A)') 'Lode angle parameter : '
                WRITE(unit,'(5X,3ES15.7)') this%lode_angle_parameter
                WRITE(unit,*) ''

                WRITE(unit,'(2X,A,I2)') 'INFO STATUS : ', info 
                WRITE(unit,*) ''

                WRITE(unit,*) '========================================================'

            ELSE
                WRITE(*,*) ''
                WRITE(*,*) '***************STRESS DISPLAY INFORMATION***************'
                WRITE(*,*) '========================================================'
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Stress 2D : '
                CALL PrintMatrix(realMat=this%stress_2d)
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Stress Voigt : '
                DO i = 1, NDIM_STRESS_Voigt
                    WRITE(*,'(5X,3ES15.7)') this%stress_v(i)
                END DO
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Principal Stresses : '
                DO i = 1, 3
                    WRITE(*,'(5X,3ES15.7)') stress_principal(i)
                END DO
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Von mises stress/Effective stress  : '
                WRITE(*,'(5X,3ES15.7)') this%stress_VM
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Mean Stress : '
                WRITE(*,'(5X,3ES15.7)') this%stress_mean
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Deviatoric stress : '
                CALL PrintMatrix(realMat=this%stress_deviatoric)
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Stress Triaxiality : '
                WRITE(*,'(5X,3ES15.7)') this%stress_triaxiality
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Lode angle : '
                WRITE(*,'(5X,3ES15.7)') this%lode_angle
                WRITE(*,*) ''

                WRITE(*,'(2X,A)') 'Lode angle parameter : '
                WRITE(*,'(5X,3ES15.7)') this%lode_angle_parameter
                WRITE(*,*) ''

                WRITE(*,'(2X,A,I2)') 'INFO STATUS : ', info 
                WRITE(*,*) ''

                WRITE(*,*) '========================================================'

            END IF

        END SUBROUTINE display_stress
    
END MODULE stress_module