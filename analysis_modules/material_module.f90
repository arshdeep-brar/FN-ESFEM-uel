MODULE material_module
    !Pupose:
    ! The module contains the object material containing information of the sdv
    !and material constants for plasticity model

    USE parameter_module,   ONLY: DP, ZERO, ONE, PI, ONE_SIXTH, THREE,         &
                            &     TWO, ONE_THIRD, HALF, TWO, TWENTYSEVEN,      &
                            &     SIX, NDIM, NDIM_STRESS, ERROR, MSGLENGTH,    &
                            &     NDIM_Voigt, NDIM_STRESS_Voigt, ERROR_SMALL,  &
                            &     EXPLICIT_MARGIN
    USE stress_module,      ONLY: Stress, update, extract
    USE utility_module,     ONLY: Identity
    USE log_module,         ONLY: Logs, add_log_message

    IMPLICIT NONE

    PRIVATE
    
    TYPE, PUBLIC :: IG_point
        PRIVATE
        TYPE(Stress) :: sigma
        REAL(DP)     :: strain(NDIM_STRESS_Voigt)
        REAL(DP)     :: Ep(NDIM_STRESS_Voigt)
        REAL(DP)     :: trailStress(NDIM_STRESS_Voigt)
        REAL(DP)     :: eps
    END TYPE IG_point

    INTERFACE update
        MODULE PROCEDURE update_ig_point
    END INTERFACE update

    INTERFACE extract
        MODULE PROCEDURE extract_ig_point
    END INTERFACE extract

    INTERFACE integrate 
        MODULE PROCEDURE integrate_material
    END INTERFACE integrate

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE sdv_to_IGpoint
        MODULE PROCEDURE IGpoint_to_sdv
    END INTERFACE 

    TYPE(Logs), PUBLIC, SAVE, ALLOCATABLE :: material_logs
    
    PUBLIC :: ASSIGNMENT(=), update, extract, integrate, &
            & start_material_logging, end_material_logging  

    CONTAINS
    
        SUBROUTINE update_ig_point(this, sigma, sigma_v, strain, Ep, &
                                &   trialStress, eps)
            !Purpose:
            ! Updates material integration point required for the integration 
            ! at a material point
            !
            TYPE(IG_point),      INTENT(INOUT) :: this
            TYPE(Stress), OPTIONAL, INTENT(IN) :: sigma
            REAL(DP),     OPTIONAL, INTENT(IN) :: sigma_v(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(IN) :: eps
            REAL(DP),     OPTIONAL, INTENT(IN) :: strain(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(IN) :: trialStress(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(IN) :: Ep(NDIM_STRESS_Voigt)
            
            IF (PRESENT(sigma)) this%sigma = sigma

            IF (PRESENT(sigma_v)) THEN
                CALL update(this%sigma, stress_v=sigma_v)
            END IF

            IF (PRESENT(Ep))          this%Ep          = Ep
            IF (PRESENT(strain))      this%strain      = strain
            IF (PRESENT(trialStress)) this%trailStress = trialStress
            IF (PRESENT(eps))         this%eps         = eps
            
        END SUBROUTINE update_ig_point

        PURE SUBROUTINE extract_ig_point(this, sigma, strain, Ep,   &   
                                    &    trialStress, eps)
            !Purpose:
            ! Extarct integration point variable 

            TYPE(IG_point),          INTENT(IN) :: this
            TYPE(Stress), OPTIONAL, INTENT(OUT) :: sigma
            REAL(DP),     OPTIONAL, INTENT(OUT) :: strain(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(OUT) :: Ep(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(OUT) :: trialStress(NDIM_STRESS_Voigt)
            REAL(DP),     OPTIONAL, INTENT(OUT) :: eps

            IF (PRESENT(sigma))       sigma       = this%sigma
            IF (PRESENT(strain))      strain      = this%strain
            IF (PRESENT(Ep))          Ep          = this%Ep
            IF (PRESENT(trialStress)) trialStress = this%trailStress
            IF (PRESENT(eps))         eps         = this%eps

        END SUBROUTINE

        SUBROUTINE sdv_to_IGpoint(this, arr)

            TYPE(IG_point), INTENT(OUT) :: this
            REAL(DP),        INTENT(IN) :: arr(3*NDIM_STRESS_Voigt + 1)

            CALL update(this, sigma_v=arr(1:NDIM_STRESS_Voigt),                &
                    &     strain=arr(NDIM_STRESS_Voigt+1:2*NDIM_STRESS_Voigt), &
                    &     Ep=arr(2*NDIM_STRESS_Voigt+1:3*NDIM_STRESS_Voigt),   &
                    &     eps=arr(3*NDIM_STRESS_Voigt+1))

        END SUBROUTINE sdv_to_IGpoint

        SUBROUTINE IGpoint_to_sdv(arr, this)

            REAL(DP),       INTENT(OUT) :: arr(3*NDIM_STRESS_Voigt + 1)
            TYPE(IG_point),  INTENT(IN) :: this

            REAL(DP) :: stress_v(NDIM_STRESS_Voigt)

            CALL extract(this%sigma, stress_v=stress_v)

            arr(1:NDIM_STRESS_Voigt)                       = stress_v
            arr(NDIM_STRESS_Voigt+1:2*NDIM_STRESS_Voigt)   = this%strain
            arr(2*NDIM_STRESS_Voigt+1:3*NDIM_STRESS_Voigt) = this%Ep
            arr(3*NDIM_STRESS_Voigt+1)                     = this%eps

        END SUBROUTINE IGpoint_to_sdv

        SUBROUTINE start_material_logging(level, elem, inc)
            USE log_module,     ONLY: add_log_message, set_log_elem, &
                                &     set_log_kinc, set_log_level

            INTEGER, INTENT(IN) :: level, elem, inc

            ALLOCATE(material_logs)
            CALL set_log_level(material_logs, level)
            CALL set_log_elem(material_logs, elem)
            CALL set_log_kinc(material_logs, inc)

            CALL add_log_message(material_logs, info=       &
                    &   'Initiating log record for material module')            

        END SUBROUTINE start_material_logging

        SUBROUTINE end_material_logging
            USE log_module,     ONLY: add_log_message, write_logs

            CALL add_log_message(material_logs, info=       &
            &   'Ending log record for material module')
            CALL write_logs(material_logs)

            DEALLOCATE(material_logs)
        END SUBROUTINE end_material_logging
        
        SUBROUTINE integrate_material(igp, Cep, Q_np1, dtD, detF, algo, sf, cj2t)
            !Purpose: 
            ! Performs explicit integration for plasticity model
            !
            !Inputs:
            !   igp   : Integration point object 
            !   Q_np1 : Rotation tensor for stress
            !   dtD   : Incremental rate of deformation
            !   algo  : Alogorithm type used 
            !           ('E' - Explicit, 'I' - Implicit) 
            !   sf    : Stress formulation used
            !           ('K' - Kirchoff stress formulation, 
            !            'C' - Cauchy Stress formulation  ) 
            !   cj2t  : Flag to convert Jaumann moduli to Truesdell
            !           ('Y' - Convert, 'N' - don't convert)
            !Output:
            !   Cep : Elastoplastic modulus 
            USE parameter_module,   ONLY: PLASTICITY_FLAG
            USE stress_module,      ONLY: Stress, ASSIGNMENT(=)
            USE utility_module,     ONLY: Voigt2 
            USE matlib_module,      ONLY: get_yield_function, get_elasticity_moduli

            IMPLICIT NONE 

            TYPE(IG_point), INTENT(INOUT) :: igp
            REAL(DP),         INTENT(OUT) :: Cep(NDIM_Voigt,NDIM_Voigt) 
            REAL(DP),          INTENT(IN) :: Q_np1(NDIM_STRESS, NDIM_STRESS)
            REAL(DP),          INTENT(IN) :: dtD(NDIM_STRESS_Voigt)
            REAL(DP),          INTENT(IN) :: detF
            CHARACTER(len=1),  INTENT(IN) :: algo
            CHARACTER(len=1),  INTENT(IN) :: sf
            CHARACTER(len=1),  INTENT(IN) :: cj2t
  
            !Kinetic variables:
            ! stress_n         : Stress onject at the start of increment
            ! stress_np1       : Stress object at the end of integration
            ! stress_np1_trail : Trial stress object
            ! Ep               : Copy of last converged pastic strain 
            ! strain           : Copy of last converged strain
            ! eps              : Copy of the last converged effective plastic strain
            ! stress2d_n       : Copy of the stress in array 2d
            ! stressV_np1      : Copy of the stress in Voigt form 
            TYPE(Stress) :: stress_n
            TYPE(Stress) :: stress_np1_trial
            REAL(DP)     :: Ep(NDIM_STRESS_Voigt)
            REAL(DP)     :: strain(NDIM_STRESS_Voigt)
            REAL(DP)     :: eps
            REAL(DP)     :: stress2d(NDIM_STRESS,NDIM_STRESS)
            REAL(DP)     :: stressV(NDIM_STRESS_Voigt)

            !Tangent Modulus: 
            ! CelJ   : Elastic tangent moduli for Jaumann rate
            ! CepJ   : Elastoplastic moduli for Jaumann rate
            ! CepJ2d : 2-Dimensional copy of CepJ
            ! CepT   : Elastoplastic moduli for Truesdell rate 
            REAL(DP)     :: CelJ(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP)     :: CepJ(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP)     :: CepJ2d(NDIM_Voigt, NDIM_Voigt)
            REAL(DP)     :: CepT(NDIM_Voigt, NDIM_Voigt)
            
            !Yield function variables:
            !   yield_func : Value of the yield function obained 
            REAL(DP)     :: yield_func

            CHARACTER(len=MSGLENGTH) :: logmsg
            
            ! ---------------------------------------------------------------------- !
            !                         PREAMBLE ENDS HERE 
            ! ---------------------------------------------------------------------- !

            ! Initiating varibles
            stress_n         = ZERO
            stress_np1_trial = ZERO
            Ep               = ZERO
            strain           = ZERO
            eps              = ZERO
            stress2d         = ZERO
            stressV          = ZERO

            CelJ    = ZERO
            CepJ    = ZERO
            CepJ2d  = ZERO
            CepT    = ZERO

            yield_func  = ZERO

            logmsg = ''

            ! Updating variable from integration point object
            stress_n = igp%sigma
            strain   = igp%strain
            eps      = igp%eps
            Ep       = igp%Ep

            IF (ALLOCATED(material_logs)) THEN
                CALL add_log_message(material_logs, &
                        info='Material integration begins')
            END IF

            ! Extracting the copy of 2d stress and storing it at
            !  locat varibale stress2d_n
            CALL extract(stress_n, stress_2d=stress2d)

            CelJ = get_elasticity_moduli(NDIM_STRESS)
            ! Obtaining the traial stress using incremental objectivity
            stress_np1_trial = Voigt2(MATMUL(MATMUL(Q_np1, stress2d), &
                           &          TRANSPOSE(Q_np1)),1) +            &
                           &   MATMUL(CelJ, dtD)

            CALL extract(stress_np1_trial, stress_v=stressV)

            igp%trailStress = stressV

            ! Obtaining yield function value of the trial stress if the analysis
            !  is plastic
            IF (PLASTICITY_FLAG .EQ. 1) THEN
                CALL get_yield_function(yield_func, stress_np1_trial, eps)
                WRITE(logmsg,'(A,ES15.7)')  &
                    &   'Plastic analysis, value of yield function : ',  &
                    &   yield_func 
            END IF

            IF (ALLOCATED(material_logs)) THEN
                CALL add_log_message(material_logs, info=logmsg)
            END IF

            ! Checking if the yield function value of trail stress is 
            !   the yield surface:
            !       if yield_fucn < 0 => only elastic loading
            !       else ==> plastic loading

            IF (yield_func .LE. ZERO) THEN
                ! No plastic loading 
                igp%sigma  = stress_np1_trial
                igp%Ep     = Ep
                igp%eps    = eps
                igp%strain = strain + dtD ! Check again
                IF (ALLOCATED(material_logs)) THEN
                    CALL add_log_message(material_logs, &
                        &   info='Elastic loading')
                END IF
                IF (PLASTICITY_FLAG .EQ. 1) THEN
                    CepJ = CelJ
                ELSE 
                    Cep = get_elasticity_moduli(NDIM)
                    RETURN
                END IF


            ELSE 
                ! Plastic Loading
                IF (ALLOCATED(material_logs)) THEN
                    CALL add_log_message(material_logs, &
                        &   info='Plastic loading')
                END IF
                IF (algo .EQ. 'E') THEN
                    IF (ALLOCATED(material_logs)) THEN
                        CALL add_log_message(material_logs, &
                            &   info='Explicit integration')
                    END IF
                    CALL explicit_integration(igp, CepJ, stress_np1_trial, dtD, CelJ)
                ELSE IF (algo .EQ. 'I') THEN
                    IF (ALLOCATED(material_logs)) THEN
                        CALL add_log_message(material_logs, &
                            &   info='Implicit integration')
                    END IF
                    CALL implicit_integration(igp, CepJ, stress_np1_trial, dtD, CelJ)
                END IF

            END IF

            CepJ2d(1:NDIM, 1:NDIM)         = CepJ(1:NDIM, 1:NDIM)
            CepJ2d(NDIM_Voigt, 1:NDIM)     = CepJ(NDIM_STRESS_Voigt, 1:NDIM)
            CepJ2d(1:NDIM, NDIM_Voigt)     = CepJ(1:NDIM, NDIM_STRESS_Voigt)
            CepJ2d(NDIM_Voigt, NDIM_Voigt) = CepJ(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)

            IF (cj2t .EQ. 'Y') THEN 
                CALL extract(igp%sigma, stress_v=stressV)
                IF (sf .EQ. 'C') THEN
                    CepT = Jaumann_to_Tressdell(CepJ2d, stressV, 1, detF=detF)
                ELSE IF (sf .EQ. 'K') THEN
                    CepT = Jaumann_to_Tressdell(CepJ2d, stressV, 2, detF=detF)
                END IF
            ELSE 
                CepT = CepJ2d
            END IF
            
            Cep = CepT
            ! CALL add_log_message(elemlogs, info='Material Integration Successful!\n')

        END SUBROUTINE integrate_material

        SUBROUTINE explicit_integration(igp, Cep, stress_trial, dtD, Cel)
            !Purpose:
            ! Explicit integration for material
            USE matlib_module,      ONLY: gets_yield_derivative
            USE stress_module,      ONLY: OPERATOR(-) 
            
            TYPE(IG_point), INTENT(INOUT) :: igp
            REAL(DP),         INTENT(OUT) :: Cep(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            TYPE(Stress),      INTENT(IN) :: stress_trial
            REAL(DP),          INTENT(IN) :: dtD(NDIM_STRESS_Voigt)
            REAL(DP),          INTENT(IN) :: Cel(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)

            ! Internal variable
            REAL(DP) :: f_stress(NDIM_STRESS_Voigt)
            REAL(DP) :: f_q
            REAL(DP) :: f_stressT(1,NDIM_STRESS_Voigt)
            REAL(DP) :: Celr(NDIM_STRESS_Voigt,1)
            REAL(DP) :: rCel(1,NDIM_STRESS_Voigt)
            REAL(DP) :: dtDp(NDIM_STRESS_Voigt)
            REAL(DP) :: dellambda
            REAL(DP) :: eps
            REAL(DP) :: Ep(NDIM_STRESS_Voigt)
            REAL(DP) :: strain(NDIM_STRESS_Voigt)
            REAL(DP) :: nom
            REAL(DP) :: denom

            eps      = igp%eps
            strain   = igp%strain
            Ep       = igp%Ep

            CALL gets_yield_derivative(f_q, f_stress, igp%sigma, igp%eps)

            f_stressT(1,:) = f_stress
            
            nom = SUM(f_stress * MATMUL(Cel, dtD))

            denom = - f_q + SUM(f_stress * MATMUL(Cel, f_stress))

            dellambda = nom/denom

            rCel = MATMUL(f_stressT, Cel)
            Celr(:,1) = MATMUL(Cel, f_stress)

            Cep = Cel - (ONE/denom) * MATMUL(Celr, rCel)

            dtDp = dellambda * f_stress
            
            igp%sigma  = stress_trial - MATMUL(Cel, dtDp)
            igp%eps    = eps + dellambda
            igp%Ep     = Ep + dtDp
            igp%strain = strain + dtD
        
        END SUBROUTINE explicit_integration

        SUBROUTINE semi_implicit_integration(igp, Cep, stress_trial, dtD, Cel)
            
            USE matlib_module,      ONLY: gets_yield_derivative, get_yield_function 
            USE stress_module,      ONLY: OPERATOR(-)

            TYPE(IG_point), INTENT(INOUT) :: igp
            REAL(DP),         INTENT(OUT) :: Cep(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)
            TYPE(Stress),      INTENT(IN) :: stress_trial
            REAL(DP),          INTENT(IN) :: dtD(NDIM_STRESS_Voigt)
            REAL(DP),          INTENT(IN) :: Cel(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)

            TYPE(Stress) :: stress_k
            REAL(DP)     :: eps_k
            REAL(DP)     :: Ep_k(NDIM_STRESS_Voigt)
            REAL(DP)     :: strain_k(NDIM_STRESS_Voigt)

            REAL(DP) :: Amatrix_k(NDIM_STRESS_Voigt+1, NDIM_STRESS_Voigt+1)
            REAL(DP) :: lambda_k
            REAL(DP) :: tildaf_k(NDIM_STRESS_Voigt+1)
            REAL(DP) :: tildar_k(NDIM_STRESS_Voigt+1)
            REAL(DP) :: f_k
            REAL(DP) :: f_stress(NDIM_STRESS_Voigt)
            REAL(DP) :: dtDp(NDIM_STRESS_Voigt)
            REAL(DP) :: f_q
            REAL(DP) :: dellambda

            REAL(DP) :: rn(NDIM_STRESS_Voigt)
            REAL(DP) :: f_stressT(1, NDIM_STRESS_Voigt)
            
            REAL(DP) :: Celrn(NDIM_STRESS_Voigt,1)
            REAL(DP) :: fsigCel(1, NDIM_STRESS_Voigt)
            REAL(DP) :: denom

            INTEGER  :: numitr 
            
            stress_k = stress_trial
            eps_k    = igp%eps
            strain_k = igp%strain
            Ep_k     = igp%Ep
            
            CALL gets_yield_derivative(f_q, f_stress, stress_k, eps_k)
            
            rn = f_stress

            tildar_k(1:NDIM_STRESS_Voigt) = f_stress
            tildar_k(NDIM_STRESS_Voigt+1) = ONE

            lambda_k = ZERO

            dellambda = ZERO

            Amatrix_k = ZERO
            Amatrix_k(1:NDIM_STRESS_Voigt, 1:NDIM_STRESS_Voigt) = Cel
            Amatrix_k(NDIM_STRESS_Voigt+1, NDIM_STRESS_Voigt+1) = -ONE

            numitr = 0

            DO 

                tildaf_k(1:NDIM_STRESS_Voigt) = f_stress
                tildaf_k(NDIM_STRESS_Voigt+1) = f_q
                
                CALL get_yield_function(f_k, stress_k, eps_k)

                lambda_k = f_k/(SUM(tildaf_k * MATMUL(Amatrix_k, tildar_k)))

                dellambda = dellambda + lambda_k
                dtDp      = dellambda * rn
                stress_k  = stress_trial - MATMUL(Cel, dtDp)
                eps_k     = eps_k + lambda_k
                Ep_k      = Ep_k + dtDp
                strain_k  = strain_k + dtD

                numitr = numitr + 1

                CALL gets_yield_derivative(f_q, f_stress, stress_k, eps_k)

                IF (ABS(f_k) .LT. ERROR) EXIT
                
            END DO
            ! WRITE(logmsg, '(A, I5)') 'Number of interation required for convergence: ', numitr
            ! CALL add_log_message(elemlogs, info=logmsg)
            
            Celrn(:,1) = MATMUL(Cel, rn)
            f_stressT(1,:) = f_stress 
            fsigCel =  MATMUL(f_stressT, Cel)

            denom = -f_q + SUM(f_stress * MATMUL(Cel, rn))

            Cep = Cel - (ONE/denom) * MATMUL(Celrn, fsigCel)

            igp%sigma  = stress_k
            igp%eps    = eps_k
            igp%strain = strain_k
            igp%Ep     = Ep_k

        END SUBROUTINE semi_implicit_integration

        SUBROUTINE implicit_integration(igp, Cep, stress_trial, dtD, Cel)

            USE matlib_module,      ONLY: gets_yield_derivative, get_yield_function, &
                                    &     plasticty_type 
            USE stress_module,      ONLY: OPERATOR(-)

            TYPE(IG_point), INTENT(INOUT) :: igp
            REAL(DP),         INTENT(OUT) :: Cep(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)
            TYPE(Stress),      INTENT(IN) :: stress_trial
            REAL(DP),          INTENT(IN) :: dtD(NDIM_STRESS_Voigt)
            REAL(DP),          INTENT(IN) :: Cel(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)

            TYPE(Stress) :: stress_k
            REAL(DP)     :: eps_k
            REAL(DP)     :: Ep_k(NDIM_STRESS_Voigt)
            REAL(DP)     :: strain_k(NDIM_STRESS_Voigt)
            
            stress_k = stress_trial
            eps_k    = igp%eps
            strain_k = igp%strain
            Ep_k     = igp%Ep
        
            IF (plasticty_type .EQ. 'J2') THEN
                
                IF (ALLOCATED(material_logs)) THEN
                    CALL add_log_message(material_logs, &
                    &   info='J2 plasticity - radial return algo used for integration')
                END IF

                CALL radial_return_algo(Cep, stress_k, Ep_k, eps_k, Cel)
                igp%sigma  = stress_k
                igp%strain = strain_k + dtD
                igp%Ep     = Ep_k
                igp%eps    = eps_k 
                
                RETURN 

            ELSE IF (plasticty_type .EQ. 'BW') THEN
                
                ! Need to be developed in radial return algo ...

            END IF

        END SUBROUTINE implicit_integration

        SUBROUTINE radial_return_algo(Cep, stress_k, Ep, eps, Cel)

            USE matlib_module,      ONLY: gets_yield_dd, get_yield_function, &
                                    &     gets_yield_derivative, hardeningfunc
            USE parameter_module,   ONLY: TOLERANCE
            USE stress_module,      ONLY: OPERATOR(-), extract

            REAL(DP),       INTENT(OUT) :: Cep(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            TYPE(Stress), INTENT(INOUT) :: stress_k
            REAL(DP),     INTENT(INOUT) :: Ep(NDIM_STRESS_Voigt)
            REAL(DP),     INTENT(INOUT) :: eps
            REAL(DP),        INTENT(IN) :: Cel(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)

            REAL(DP) :: b, a
            REAL(DP) :: mu
            REAL(DP) :: dellambda
            REAL(DP) :: r_stress(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: Ihat(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: Ctilda(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: stress_VM
            ! REAL(DP) :: stress_VM0

            REAL(DP) :: f_k
            REAL(DP) :: lambda_k
            REAL(DP) :: f_q
            REAL(DP) :: f_stress(NDIM_STRESS_Voigt)
            ! REAL(DP) :: sigY
            ! REAL(DP) :: H
            
            REAL(DP) :: f_stressT(1,NDIM_STRESS_Voigt)
            REAL(DP) :: Ctildar(NDIM_STRESS_Voigt,1)
            REAL(DP) :: rCtilda(1,NDIM_STRESS_Voigt)
            REAL(DP) :: denom

            INTEGER  :: Niter, i
            
            CHARACTER(len=MSGLENGTH) :: logmsg

            logmsg = ''

            dellambda = ZERO
            mu = Cel(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            
            Niter = 0

            ! CALL extract(stress_k, stress_VM=stress_VM0)

            DO 
                Niter = Niter + 1
                
                CALL get_yield_function(f_k, stress_k, eps)
                CALL gets_yield_derivative(f_q, f_stress, stress_k, eps)
                ! f_k = stress_VM0 - THREE * mu * dellambda -  sigY

                IF (ABS(f_k) .LT. TOLERANCE) EXIT
                    
                lambda_k = f_k / (THREE * mu - f_q)
                
                WRITE(logmsg,'(A,I4,2(A,ES15.7))') &
                    &   'For iter no : ', Niter, ', lamda_k : ', lambda_k, &
                    &   ' f_k : ', f_k
                
                IF (ALLOCATED(material_logs)) THEN
                    CALL add_log_message(material_logs, info=logmsg)
                END IF

                dellambda = dellambda + lambda_k
                eps       = eps + lambda_k
                Ep        = Ep + lambda_k * f_stress
                stress_k  = stress_k - lambda_k * MATMUL(Cel, f_stress)
                
            END DO
            
            CALL gets_yield_derivative(f_q, f_stress, stress_k, eps)

            IF (ALLOCATED(material_logs)) THEN
                CALL add_log_message(material_logs, &
                    &   info='Integration successful')
            END IF
            
            CALL gets_yield_dd(r_stress=r_stress, thisstress=stress_k, &
                                &   eps=eps)

            CALL extract(stress_k, stress_VM=stress_VM)
            a = THREE * dellambda / (TWO * stress_VM) 
            b = TWO * mu * a / (ONE + TWO * mu * a)
            
            Ihat = (TWO * stress_VM / THREE) * r_stress
            Ctilda = Cel - TWO * mu * b * Ihat

            f_stressT(1,:) = f_stress
            Ctildar(:,1) = MATMUL(Ctilda, f_stress)
            rCtilda      = MATMUL(f_stressT, Ctilda)

            denom = - f_q + SUM(f_stress * MATMUL(Ctilda, f_stress))

            Cep = Ctilda - (ONE/denom) * MATMUL(Ctildar, rCtilda)
             
        END SUBROUTINE radial_return_algo

        SUBROUTINE backward_return_algo(Amat_k, stress_k, Cel, eps_k, dellambda)

            USE matlib_module,      ONLY: gets_yield_dd
            USE utility_module,     ONLY: inv, Inverse3

            REAL(DP),    INTENT(OUT) :: Amat_k(NDIM_STRESS_Voigt+1, NDIM_STRESS_Voigt+1)
            TYPE(Stress), INTENT(IN) :: stress_k
            REAL(DP),     INTENT(IN) :: Cel(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP),     INTENT(IN) :: eps_k
            REAL(DP),     INTENT(IN) :: dellambda

            REAL(DP) :: invAmat(NDIM_STRESS_Voigt+1,NDIM_STRESS_Voigt+1)
            REAL(DP) :: Cinv(NDIM_STRESS_Voigt,NDIM_STRESS_Voigt)
            REAL(DP) :: r_stress(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: r_q(NDIM_STRESS_Voigt)
            REAL(DP) :: h_stress(NDIM_STRESS_Voigt)
            REAL(DP) :: h_q

            INTEGER :: info

            Amat_k = ZERO

            CALL inv(Cel, Cinv, info)
            
            CALL gets_yield_dd(r_stress=r_stress, r_q=r_q, h_stress=h_stress, &
                            &   h_q=h_q, thisStress=stress_k, eps=eps_k)

            invAmat(1:NDIM_STRESS_Voigt,1:NDIM_STRESS_Voigt) = Cinv + &
                                                &   dellambda * r_stress
                                                
            invAmat(1:NDIM_STRESS_Voigt,NDIM_STRESS_Voigt+1) = r_q
            invAmat(NDIM_STRESS_Voigt+1,1:NDIM_STRESS_Voigt) = h_stress
            invAmat(NDIM_STRESS_Voigt+1,NDIM_STRESS_Voigt+1) = - ONE + h_q

            CALL inv(invAmat, Amat_k, info)

        END SUBROUTINE backward_return_algo

        ! SUBROUTINE explicit_integration_check(tau, tau_np1, hardening, eps, newdt, stat, logmsg)
        !     !Purpose:
        !     ! The subroutine checks if the current effect stress value is within the bounds of explicit
        !     ! integration scheme

        !     REAL(DP),                  INTENT(IN) :: tau(NDIM_STRESS_Voigt)
        !     REAL(DP),                  INTENT(IN) :: tau_np1(NDIM_STRESS_Voigt)
        !     TYPE(SH_const),            INTENT(IN) :: hardening
        !     REAL(DP),                  INTENT(IN) :: eps
        !     REAL(DP),                 INTENT(OUT) :: newdt
        !     INTEGER,                  INTENT(OUT) :: stat
        !     CHARACTER(len=MSGLENGTH), INTENT(OUT) :: logmsg
            
        !     REAL(DP) :: tauY
        !     REAL(DP) :: tau_eff
        !     REAL(DP) :: tau_trial_eff
        !     REAL(DP) :: tauT

        !     CHARACTER(len=MSGLENGTH) :: msg

        !     stat = 0

        !     logmsg = ''

        !     ! CALL hardeningfunc(eps, sigY=tauY)

        !     tauY = hardening%A * (hardening%e0 + eps)**hardening%n

        !     CALL get_Stress_State(tau, tau_eff=tau_eff, stat=stat)

        !     IF (stat .NE. 0) GOTO 10
                
        !     IF ((tau_eff .GT. EXPLICIT_MARGIN*tauY) .AND. (tau_eff .LT. (TWO - EXPLICIT_MARGIN)*tauY)) THEN
        !         newdt = ONE
        !     ELSE IF (tau_eff .GT. tauY) THEN
        !         stat = 1 
        !         msg = 'Effecive stress value outside the yield envelop: '
        !         WRITE(logmsg, '(A, 2ES15.7)') TRIM(msg), tau_eff, tauY
        !         RETURN
        !     ELSE
        !         CALL get_Stress_State(tau_np1, tau_eff=tau_trial_eff, stat=stat)
                
        !         IF (stat .NE. 0) GOTO 10

        !         tauT = (ONE + EXPLICIT_MARGIN) * tauY/TWO
        !         newdt = (tauT - tau_eff)/(tau_trial_eff - tau_eff)

        !         RETURN
        !     END IF 

        ! 10  msg = 'Eigenvalue calculation Failed!! (LAPACK ERROR), loc:explicit integration check\n'
        !     WRITE(logmsg, '(A)') TRIM(msg)
        !     RETURN
                
        ! END SUBROUTINE explicit_integration_check

        PURE FUNCTION Jaumann_to_Tressdell(CepJ, S, Stype, detF) RESULT(CepT)
            ! Purpose:
            !   Converst the Jaumann Stress tangent moduli to the to the Tressdell
            !   tangent moduli
            ! Inputs:
            !   CepJ : Jaumann tangent moduli (in the Voigt notation)
            !      S : Stress tensor(in the Voigt Notation)
            !  Stype : Integer specifying type of Stress used in plastcity integration
            !          1 - Cauchy stress
            !          2 - Kirchoff stress
            !   detF : Jacobian determinant
            !
            
            USE utility_module, ONLY : KDelta

            REAL(DP),           INTENT(IN) :: CepJ(:,:)
            REAL(DP),           INTENT(IN) :: S(:)
            INTEGER,            INTENT(IN) :: Stype
            REAL(DP), OPTIONAL, INTENT(IN) :: detF
             

            REAL(DP), ALLOCATABLE :: CepT(:,:)
            REAL(DP), ALLOCATABLE :: C_dash(:,:)
            REAL(DP), ALLOCATABLE :: sigotI(:,:)
            
            INTEGER :: i, n_dim, n_stress

            n_dim = SIZE(CepJ,1)
            n_stress = SIZE(S)

            ALLOCATE(CepT(n_dim,n_dim))
            ALLOCATE(C_dash(n_dim,n_dim))
            ALLOCATE(sigotI(n_dim,n_dim))

            C_dash = ZERO
            sigotI = ZERO

            IF (n_dim .EQ. 3) THEN
                IF (n_stress .EQ. 6) THEN
                    C_dash(1,1) = TWO*S(1)
                    C_dash(1,3) = S(6)
                    C_dash(2,2) = TWO*S(2)
                    C_dash(2,3) = S(6)
                    C_dash(3,3) = HALF * (S(1) + S(2))
                    C_dash(3,2) = C_dash(2,3)
                    C_dash(3,1) = C_dash(1,3)

                    IF (Stype .EQ. 1) THEN
                        DO i = 1,n_dim-1
                            sigotI(i,1) = S(i)
                            sigotI(i,2) = S(i)
                        END DO
                        sigotI(3,1) = S(6)
                        sigotI(3,2) = S(6)
                    END IF
                ELSE IF (n_stress .EQ. 3) THEN
                    C_dash(1,1) = TWO*S(1)
                    C_dash(1,3) = S(3)
                    C_dash(2,2) = TWO*S(2)
                    C_dash(2,3) = S(3)
                    C_dash(3,3) = HALF * (S(1) + S(2))
                    C_dash(3,2) = C_dash(2,3)
                    C_dash(3,1) = C_dash(1,3)

                    IF (Stype .EQ. 1) THEN
                        DO i = 1,n_dim
                            sigotI(i,1) = S(i)
                            sigotI(i,2) = S(i)
                        END DO
                    END IF
                END IF
            END IF

            IF (n_dim .EQ. 6) THEN
                C_dash(1,1) = TWO*S(1)
                C_dash(2,2) = TWO*S(2)
                C_dash(3,3) = TWO*S(3)
                C_dash(4,4) = HALF * (S(2) + S(3))
                C_dash(5,5) = HALF * (S(1) + S(3))
                C_dash(6,6) = HALF * (S(1) + S(2))
                C_dash(1,5) = S(5)
                C_dash(5,1) = S(5)
                C_dash(1,6) = S(6)
                C_dash(6,1) = S(6)
                C_dash(2,4) = S(4)
                C_dash(4,2) = S(4)
                C_dash(2,6) = S(6)
                C_dash(6,2) = S(6)
                C_dash(3,4) = S(4)
                C_dash(4,3) = S(4)
                C_dash(3,5) = S(5)
                C_dash(5,3) = S(5)
                C_dash(4,5) = HALF*S(6)
                C_dash(5,4) = HALF*S(6)
                C_dash(4,6) = HALF*S(5)
                C_dash(6,4) = HALF*S(5)
                C_dash(5,6) = HALF*S(4)
                C_dash(6,5) = HALF*S(4)

                IF (Stype .EQ. 1) THEN
                    DO i = 1,n_dim
                        sigotI(i,1) = S(i)
                        sigotI(i,2) = S(i)
                        sigotI(i,3) = S(i)
                    END DO
                END IF
            END IF

            IF (Stype .EQ. 1) THEN
                CepT = CepJ - C_dash + sigotI
            ELSE IF (Stype .EQ. 2) THEN
                IF (PRESENT(detF)) CepT = (CepJ - C_dash)/detF
            END IF
            
        END FUNCTION Jaumann_to_Tressdell

END MODULE material_module
