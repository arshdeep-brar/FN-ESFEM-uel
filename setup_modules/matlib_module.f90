MODULE matlib_module

    USE parameter_module,   ONLY: DP, ZERO, ONE, SIX, ONE_THIRD, PLASTICITY_FLAG 

    IMPLICIT NONE

    PRIVATE

    !Emperical conststants for BW plasricity model 
    TYPE, PRIVATE :: BWM_yield_const
        PRIVATE
        REAL(DP) :: Cs
        REAL(DP) :: Ct
        REAL(DP) :: Cc
        REAL(DP) :: Ceta
        REAL(DP) :: m
        REAL(DP) :: eta0
    END TYPE

    TYPE, PRIVATE :: Elasticity_const
        PRIVATE
        REAL(DP) :: E
        REAL(DP) :: nu
    END TYPE

    !Strain hardening constants for power rule
    ! sigma = A(e0 + e)^n
    TYPE, PRIVATE :: SH_const_exp
        PRIVATE
        REAL(DP) :: A  
        REAL(DP) :: e0 
        REAL(DP) :: n
    END TYPE SH_const_exp

    TYPE, PRIVATE :: SH_const_soft
        PRIVATE
        REAL(DP) :: sig0 
        REAL(DP) :: sigINF 
        REAL(DP) :: delta 
        REAL(DP) :: hm
    END TYPE SH_const_soft

    CHARACTER(len=10),      ALLOCATABLE, SAVE :: mat_name
    TYPE(Elasticity_const), ALLOCATABLE, SAVE :: elastic_const
    TYPE(BWM_yield_const),  ALLOCATABLE, SAVE :: BW_const
    TYPE(SH_const_exp),     ALLOCATABLE, SAVE :: hardening_const_exp
    TYPE(SH_const_soft),    ALLOCATABLE, SAVE :: hardening_const_soft
    
    CHARACTER(len=2), PUBLIC, ALLOCATABLE, SAVE :: plasticty_type

    PUBLIC :: set_material, get_elasticity_moduli, get_yield_function, &
            & gets_yield_derivative, gets_yield_dd, hardeningfunc,  &
            & display_material_status, deallocate_material

    CONTAINS

        SUBROUTINE set_material(name)
            ! Purpose
            ! Sets the material from the library of following materials
            !   name : Test-EL (Elastic material with no plasticity)
            !       Elastic properties : 
            !       ====================
            !       |    E    |   nu   |
            !       |  70 GPa |   0.3  |
            !       ====================
            !
            !   name : Test-EPP (Elastic prefectly plastic material)
            !       Elastic properties : 
            !       ====================
            !       |    E    |   nu   |
            !       | 210 GPa |   0.3  |
            !       ====================
            !       Hardening function sigY = A*(e0 + eps)^n
            !       =============================
            !       |     A     |  e0   |   n   |
            !       |  240 MPa  |  1.0  |  0.0  |
            !       =============================
            !
            !   name : Test-1 (Elastoplastic material with J2 plasticity)
            !       Elastic properties :
            !       ========================
            !       |      E      |   nu   |
            !       |  71.119 GPa |   0.3  |
            !       ========================
            !       Hardening function sigY = A*(e0 + eps)^n
            !       ===================================
            !       |     A     |    e0    |     n    |
            !       |  908 MPa  |  0.0058  |  0.1742  |
            !       ===================================
            !
            !   name : AL-BW (Aluminium elastoplastic model with BW plasticity)
            !       Elastic properties :
            !       ========================
            !       |      E      |   nu   |
            !       |  71.119 GPa |   0.3  |
            !       ========================
            !       BW constants :
            !       ===================================================
            !       |    Cs   |  Ct   |  Cc   |  Ceta  |  eta0  |  n  |
            !       |  0.885  |  1.0  |  0.9  |  0.09  |  0.33  |  6  |
            !       ===================================================  
            !       Hardening function sigY = A*(e0 + eps)^n
            !       ===================================
            !       |     A     |    e0    |     n    |
            !       |  908 MPa  |  0.0058  |  0.1742  |
            !       ===================================
            !
            !   name : Test-2 (Elastoplastic material with J2 plasticity)
            !       Elastic properties :
            !       =========================
            !       |      E      |    nu   |
            !       |  206.9 GPa  |   0.29  |
            !       =========================
            !       Hardening function sigY = A*(e0 + eps)^n
            !       ===================================
            !       |      A       |    e0    |   n   |
            !       |  239.67 MPa  |  0.8345  |  1.0  |
            !       ===================================
            !
            !   name : Test-3 (Elastoplastic material with J2 plasticity)
            !       Elastic properties :
            !       ======================
            !       |     E     |   nu   |
            !       |  210 GPa  |   0.3  |
            !       ======================
            !       Hardening function sigY = A*(e0 + eps)^n
            !       ===================================
            !       |      A       |    e0    |   n   |
            !       |  239.67 MPa  |  0.8345  |  1.0  |
            !       ===================================
            !
            !   name : Test-4 (Elastoplastic material with J2 plasticity)
            !       Elastic properties :
            !       ========================
            !       |      E      |   nu   |
            !       |  206.9 GPa  |   0.3  |
            !       ========================
            !       Exponential with linear hardening : 
            !           sigY = sig0 + (sigINF - sig0)(1-exp(-d*eps))+K*eps
            !       ============================================
            !       |  sig0   |  sigINF   |   d   |     K      | 
            !       | 450 MPa |  715 MPa  | 16.93 | 129.24 MPa |
            !       ============================================
            !
            USE parameter_module,   ONLY: EXIT_FUNCTION

            CHARACTER(len=*), INTENT(IN) :: name
            
            ALLOCATE(mat_name)
            mat_name = name

            IF (name .EQ. 'Test-EL') THEN
                CALL set_elasticity_consts(E=70000._DP, nu=0.3_DP)
            ELSE IF (name .EQ. 'Test-EPP') THEN
                CALL set_elasticity_consts(E=210000._DP, nu=0.3_DP)
                CALL set_hardening_const_exp(A=240._DP, e0=1._DP, &
                                        &    n=0._DP)

                ALLOCATE(plasticty_type)
                plasticty_type = 'J2'
            ELSE IF (name .EQ. 'Test-1') THEN
                CALL set_elasticity_consts(E=71119._DP, nu=0.3_DP)
                CALL set_hardening_const_exp(A=908._DP, e0=0.0058_DP, &
                                            &   n=0.1742_DP)
                ALLOCATE(plasticty_type)
                plasticty_type = 'J2'
            
            ELSE IF(name .EQ. 'AL-BW') THEN
                CALL set_elasticity_consts(E=71119._DP, nu=0.3_DP)
                CALL set_const_BW(Cs=0.885_DP, Ct=1._DP, Cc=0.9_DP, &
                            &     Ceta=0.09_DP, m=SIX, eta0=ONE_THIRD)
                CALL set_hardening_const_exp(A=908._DP, e0=0.0058_DP, &
                                            &   n=0.1742_DP)
                ALLOCATE(plasticty_type)
                plasticty_type = 'BW'
            
            ELSE IF(name .EQ. 'Test-2') THEN
                CALL set_elasticity_consts(E=206900._DP, nu=0.29_DP)
                CALL set_hardening_const_exp(A=239.673227_DP, e0=0.8344695_DP, &
                                            &    n=1._DP)
                ALLOCATE(plasticty_type)
                plasticty_type = 'J2'

            ELSE IF(name .EQ. 'Test-3') THEN
                CALL set_elasticity_consts(E=210000._DP, nu=0.3_DP)
                CALL set_hardening_const_exp(A=239.673227_DP, e0=0.8344695_DP, &
                                            &    n=1._DP)
                ALLOCATE(plasticty_type)
                plasticty_type = 'J2'
            ELSE IF(name .EQ. 'Test-4') THEN
                ! sig0 = 450._DP
                ! sigINF = 715._DP
                ! delta = 16.93_DP
                ! hm = 129.24_DP
                CALL set_elasticity_consts(E=206900._DP, nu=0.29_DP)
                CALL set_hardening_const_soft(sig0=450._DP, sigINF=715._DP,    &
                                            &   delta=16.93_DP, hm=129.24_DP)
                ALLOCATE(plasticty_type)
                plasticty_type = 'J2'
            ELSE
                WRITE(*,*) 'Material ', name, ' not present in library'
                CALL EXIT_FUNCTION
            END IF

        END SUBROUTINE set_material

        SUBROUTINE set_elasticity_consts(E, nu)
            REAL(DP), INTENT(IN) :: E
            REAL(DP), INTENT(IN) :: nu
            
            ALLOCATE(elastic_const)

            elastic_const%E  = E
            elastic_const%nu = nu

        END SUBROUTINE set_elasticity_consts
        
        SUBROUTINE set_const_BW(Cs, Ct, Cc, Ceta, m, eta0)
            REAL(DP), INTENT(IN) :: Cs
            REAL(DP), INTENT(IN) :: Ct
            REAL(DP), INTENT(IN) :: Cc
            REAL(DP), INTENT(IN) :: Ceta
            REAL(DP), INTENT(IN) :: m
            REAL(DP), INTENT(IN) :: eta0

            ALLOCATE(BW_const)

            BW_const%Cs   = Cs
            BW_const%Ct   = Ct
            BW_const%Cc   = Cc
            BW_const%Ceta = Ceta
            BW_const%m    = m
            BW_const%eta0 = eta0

        END SUBROUTINE set_const_BW

        SUBROUTINE set_hardening_const_exp(A, e0, n)
            REAL(DP), INTENT(IN) :: A  
            REAL(DP), INTENT(IN) :: e0 
            REAL(DP), INTENT(IN) :: n

            ALLOCATE(hardening_const_exp)

            hardening_const_exp%A  = A
            hardening_const_exp%e0 = e0
            hardening_const_exp%n  = n

        END SUBROUTINE set_hardening_const_exp

        SUBROUTINE set_hardening_const_soft(sig0, sigINF, delta, hm)
            REAL(DP), INTENT(IN) :: sig0 
            REAL(DP), INTENT(IN) :: sigINF 
            REAL(DP), INTENT(IN) :: delta 
            REAL(DP), INTENT(IN) :: hm
            
            ALLOCATE(hardening_const_soft)

            hardening_const_soft%sig0   = sig0
            hardening_const_soft%sigINF = sigINF
            hardening_const_soft%delta  = delta
            hardening_const_soft%hm     = hm
            
        END SUBROUTINE set_hardening_const_soft

        SUBROUTINE display_material_status()

            USE parameter_module, ONLY:EXIT_FUNCTION
            
            WRITE(*,*) ''
            WRITE(*,*) '*************MATERIAL CONSTANTS INFORMATION*************'
            WRITE(*,*) '========================================================'
            
            IF (ALLOCATED(mat_name)) WRITE(*,'(2X,A)') 'Material name : ', TRIM(mat_name)

            IF (ALLOCATED(elastic_const)) THEN
                WRITE(*,*) ''
                WRITE(*,'(5X,A,ES16.6)') "Young's Modulus : ", elastic_const%E
                WRITE(*,'(5X,A,ES16.6)') "Poisson's ratio : ", elastic_const%nu
                WRITE(*,*) ''
            ELSE 
                WRITE(*,*) 'ERROR: No elastic constants provided'
                CALL EXIT_FUNCTION
            END IF

            IF (ALLOCATED(BW_const)) THEN
                WRITE(*,'(2X,3A)') 'PLASTICITY : ', plasticty_type, &
                                &   ' yield suface used for analysis'
                WRITE(*, '(5X,A,ES16.6)') "C_s   : ", BW_const%Cs
                WRITE(*, '(5X,A,ES16.6)') "C_t   : ", BW_const%Ct
                WRITE(*, '(5X,A,ES16.6)') "C_c   : ", BW_const%Cc
                WRITE(*, '(5X,A,ES16.6)') "C_eta : ", BW_const%Ct
                WRITE(*, '(5X,A,ES16.6)') "eta_0 : ", BW_const%eta0
                WRITE(*, '(5X,A,ES16.6)') "m     : ", BW_const%m
                WRITE(*,*) ''
            ELSE
                IF (ALLOCATED(plasticty_type)) THEN 
                    WRITE(*,'(2X,3A)') 'PLASTICITY : ', plasticty_type, &
                                    &   ' yield surface used for analysis'
                    WRITE(*,*) ''
                ELSE
                    WRITE(*,'(2X,A)') 'ELASTIC MATERIAL'
                END IF
            END IF

            IF (ALLOCATED(hardening_const_exp)) THEN
                WRITE(*,'(2X,A)') 'Exponential hardening function:'
                WRITE(*,'(5X,A)') 'sigY = A(e0+eps)^n'
                WRITE(*,'(5X,A,ES16.6)') 'A  : ', hardening_const_exp%A
                WRITE(*,'(5X,A,ES16.6)') 'e0 : ', hardening_const_exp%e0
                WRITE(*,'(5X,A,ES16.6)') 'n  : ', hardening_const_exp%n

            ELSE IF (ALLOCATED(hardening_const_soft)) THEN
                WRITE(*,'(2X,A)') 'Exponential + linear hardening function:'
                WRITE(*,'(5X,A)') 'sigY = sig0 + (sigINF - sig0)(1-exp(-d*eps))+K*eps'
                WRITE(*,'(5X,A,ES16.6)') 'sig0   : ', hardening_const_soft%sig0 
                WRITE(*,'(5X,A,ES16.6)') 'sigINF : ', hardening_const_soft%sigINF    
                WRITE(*,'(5X,A,ES16.6)') 'd      : ', hardening_const_soft%delta    
                WRITE(*,'(5X,A,ES16.6)') 'K      : ', hardening_const_soft%hm    
            ELSE
                IF (ALLOCATED(plasticty_type)) THEN
                    WRITE(*,*) 'NOTE: A new hardening function should be defined in'
                    WRITE(*,'(2X,A)') 'SUBROITINE : hardeningfunc in MODULE : matlib_module'
                END IF
            END IF

            WRITE(*,*) '' 
            WRITE(*,*) '========================================================'
            WRITE(*,*) ''

        END SUBROUTINE display_material_status

        SUBROUTINE deallocate_material()
            
            IF (ALLOCATED(mat_name))             DEALLOCATE(mat_name)
            IF (ALLOCATED(elastic_const))        DEALLOCATE(elastic_const)
            IF (ALLOCATED(BW_const))             DEALLOCATE(BW_const)
            IF (ALLOCATED(hardening_const_exp))  DEALLOCATE(hardening_const_exp)
            IF (ALLOCATED(hardening_const_soft)) DEALLOCATE(hardening_const_soft)
            IF (ALLOCATED(plasticty_type))       DEALLOCATE(plasticty_type)

        END SUBROUTINE deallocate_material

        PURE SUBROUTINE hardeningfuc_exp(hardening, eps, sigY, H)
            TYPE(SH_const_exp),  INTENT(IN) :: hardening
            REAL(DP),            INTENT(IN) :: eps
            REAL(DP), OPTIONAL, INTENT(OUT) :: sigY
            REAL(DP), OPTIONAL, INTENT(OUT) :: H

            IF (PRESENT(sigY)) THEN
                sigY = hardening%A * (hardening%e0 + eps)**hardening%n
            END IF

            IF (PRESENT(H)) THEN
                H = hardening%A * hardening%n * (hardening%e0 + eps) ** &
                                                &   (hardening%n - ONE)
            END IF

        END SUBROUTINE hardeningfuc_exp
    
        PURE SUBROUTINE hardeningfunc_soft(hardening, eps, sigY, H)
            TYPE(SH_const_soft), INTENT(IN) :: hardening
            REAL(DP),            INTENT(IN) :: eps
            REAL(DP), OPTIONAL, INTENT(OUT) :: sigY
            REAL(DP), OPTIONAL, INTENT(OUT) :: H
            
            IF (PRESENT(sigY)) THEN
                sigY = hardening%sig0 + (hardening%sigINF - hardening%sig0) * &
                        &   (1 - EXP(-hardening%delta*eps)) + hardening%hm*eps
            END IF

            IF (PRESENT(H)) THEN
                H = (hardening%sigINF - hardening%sig0) * &
                    & hardening%delta * EXP(-hardening%delta*eps) + hardening%hm
            END IF
            
        END SUBROUTINE hardeningfunc_soft

        PURE SUBROUTINE hardeningfunc(eps, sigY, H)
            REAL(DP),            INTENT(IN) :: eps
            REAL(DP), OPTIONAL, INTENT(OUT) :: sigY
            REAL(DP), OPTIONAL, INTENT(OUT) :: H

            REAL(DP) :: stressY
            REAL(DP) :: dsigYde

            IF (ALLOCATED(hardening_const_exp)) THEN
                CALL hardeningfuc_exp(hardening_const_exp, eps=eps, &
                                    &    sigY=stressY, H=dsigYde)

            ELSE IF (ALLOCATED(hardening_const_soft)) THEN
                CALL hardeningfunc_soft(hardening_const_soft, eps=eps, &
                                    &    sigY=stressY, H=dsigYde)
            ELSE
                ! New hardening function define 
            END IF

            IF (PRESENT(sigY)) sigY = stressY
            IF (PRESENT(H))    H    = dsigYde
            
        END SUBROUTINE

        FUNCTION get_elasticity_moduli(dim) RESULT(Cel)
            
            USE parameter_module,   ONLY: HALF, TWO

            INTEGER, INTENT(IN) :: dim

            REAL(DP) :: Cel(dim*(dim+1)/2, dim*(dim+1)/2)

            REAL(DP) :: E, nu
            REAL(DP) :: lambda, mu

            Cel(:,:) = ZERO

            E = elastic_const%E
            nu = elastic_const%nu
            
            IF (dim .EQ. 3) THEN
                lambda = E*nu/((1+nu)*(1-2*nu))
                mu = E/(TWO*(1+nu))
                
                Cel(1,1) = lambda + TWO * mu
                Cel(2,2) = lambda + TWO * mu
                Cel(3,3) = lambda + TWO * mu
                Cel(4,4) = mu
                Cel(5,5) = mu
                Cel(6,6) = mu
                Cel(1,2) = lambda
                Cel(1,3) = lambda
                Cel(2,3) = lambda
                Cel(2,1) = Cel(1,2)
                Cel(3,1) = Cel(1,3)
                Cel(3,2) = Cel(2,3)
            
            ELSE IF (dim .EQ. 2) THEN
                Cel(1,1) = E/(1-nu**2)
                Cel(2,2) = E/(1-nu**2)
                Cel(1,2) = nu*E/(1-nu**2)
                Cel(2,1) = nu*E/(1-nu**2)
                Cel(3,3) = HALF*(1-nu)*E/(1-nu**2)
            END IF

        END FUNCTION get_elasticity_moduli 

        SUBROUTINE get_yield_function(yield_func, thisStress, eps)
            
            USE stress_module,      ONLY: Stress, extract
            USE parameter_module,   ONLY: PI, ONE_SIXTH

            REAL(DP),    INTENT(OUT) :: yield_func
            TYPE(Stress), INTENT(IN) :: thisStress
            REAL(DP),     INTENT(IN) :: eps

            REAL(DP) :: stress_VM
            REAL(DP) :: eta
            REAL(DP) :: sigY
            REAL(DP) :: f_eta
            REAL(DP) :: f_theta
            REAL(DP) :: gamma
            REAL(DP) :: GAMMACONST
            REAL(DP) :: lode_angle
            REAL(DP) :: lode_angle_parameter
            
            stress_VM            = ZERO
            sigY                 = ZERO
            f_eta                = ZERO
            f_theta              = ZERO
            gamma                = ZERO
            GAMMACONST           = ZERO
            lode_angle           = ZERO
            lode_angle_parameter = ZERO

            CALL hardeningfunc(eps=eps, sigY=sigY)            
            
            IF (ALLOCATED(BW_const)) THEN

                CALL extract(thisStress, stress_VM=stress_VM,                    &
                            &    stress_triaxiality=eta, lode_angle=lode_angle,  &
                            &    lode_angle_parameter=lode_angle_parameter)

                GAMMACONST = COS(PI/SIX)/(ONE - COS(PI/SIX))

                gamma = GAMMACONST * ((ONE/COS(lode_angle - PI * ONE_SIXTH)) - ONE)

                IF (lode_angle_parameter .GE. ZERO) THEN
                    f_theta = BW_const%Cs + (BW_const%Ct - BW_const%Cs) *         &
                        &     (gamma - (gamma**(BW_const%m + ONE))/(BW_const%m + ONE))
                    
                ELSE
                    f_theta = BW_const%Cs + (BW_const%Cc - BW_const%Cs) *         &
                        &     (gamma - (gamma**(BW_const%m + ONE))/(BW_const%m + ONE))
                
                END IF

                f_eta = ONE - BW_const%Ceta * (eta - BW_const%eta0)

                yield_func = stress_VM - sigY * f_theta * f_eta
            
            ELSE

                CALL extract(thisStress, stress_VM=stress_VM)
                
                yield_func = stress_VM - sigY

            END IF

        END SUBROUTINE get_yield_function

        SUBROUTINE gets_yield_derivative(yield_q, yield_stress, thisStress, eps)

            USE stress_module,      ONLY: Stress, extract
            USE parameter_module,   ONLY: PI, ONE_SIXTH, TWO, THREE, ERROR_SMALL, &
                                    &     NDIM_STRESS, NDIM_STRESS_Voigt
            USE utility_module,     ONLY: Identity, Voigt2

            REAL(DP),    INTENT(OUT) :: yield_q
            REAL(DP),    INTENT(OUT) :: yield_stress(NDIM_STRESS_Voigt)
            TYPE(Stress), INTENT(IN) :: thisStress
            REAL(DP),     INTENT(IN) :: eps

            ! Internal variables 
            ! f_eta : Correction w.r.t. hydrostatic pressure 
            ! f_theta : correction w.r.t. load angle parameter 
            ! gamma : internal variable for f_theta
            ! Dstress_VMDstress : Derivative of effective stress w.r.t. stress tensor
            ! DetaDstress : Derivative of stress triaxiality w.r.t. stress tensor
            ! DgammaDstress : Derivative of gamma w.r.t. stress tensor
            ! angle : supporting variable 
            ! Id : indentity matrix
            ! H : plastic modulus 
            ! C_ax : BW model yield function constant
            REAL(DP) :: f_eta
            REAL(DP) :: f_theta
            REAL(DP) :: gamma
            REAL(DP) :: Dstress_VMDstress(NDIM_STRESS,NDIM_STRESS)
            REAL(DP) :: DetaDstress(NDIM_STRESS,NDIM_STRESS)
            REAL(DP) :: DgammaDstress(NDIM_STRESS,NDIM_STRESS)
            REAL(DP) :: angle
            REAL(DP) :: Id(NDIM_STRESS,NDIM_STRESS)
            REAL(DP) :: sigY, H
            REAL(DP) :: C_ax

            REAL(DP) :: GAMMACONST
            REAL(DP) :: lode_angle
            REAL(DP) :: lode_angle_parameter
            REAL(DP) :: stress_VM
            REAL(DP) :: eta
            REAL(DP) :: stress_dev(NDIM_STRESS, NDIM_STRESS)
            REAL(DP) :: f_stress(NDIM_STRESS, NDIM_STRESS)

            f_eta                = ZERO
            f_theta              = ZERO
            gamma                = ZERO
            Dstress_VMDstress    = ZERO
            DetaDstress          = ZERO
            DgammaDstress        = ZERO
            angle                = ZERO
            Id                   = ZERO
            sigY                 = ZERO
            H                    = ZERO
            C_ax                 = ZERO
            GAMMACONST           = ZERO
            lode_angle           = ZERO
            lode_angle_parameter = ZERO
            stress_VM            = ZERO
            eta                  = ZERO
            stress_dev           = ZERO

            CALL hardeningfunc(eps=eps, sigY=sigY, H=H)
            
            IF (ALLOCATED(BW_const)) THEN
                CALL extract(thisStress, stress_VM=stress_VM,                    &
                            &    stress_triaxiality=eta, lode_angle=lode_angle,  &
                            &    lode_angle_parameter=lode_angle_parameter,      &
                            &    stress_deviatoric=stress_dev)
                
                GAMMACONST = COS(PI/SIX)/(ONE - COS(PI/SIX))

                gamma = GAMMACONST * ((ONE/COS(lode_angle - PI * ONE_SIXTH)) - ONE)

                IF (lode_angle_parameter .GE. ZERO) THEN
                    f_theta = BW_const%Cs + (BW_const%Ct - BW_const%Cs) *         &
                        &     (gamma - (gamma**(BW_const%m + ONE))/(BW_const%m + ONE))
                    C_ax = BW_const%Ct

                ELSE
                    f_theta = BW_const%Cs + (BW_const%Cc - BW_const%Cs) *         &
                        &     (gamma - (gamma**(BW_const%m + ONE))/(BW_const%m + ONE))
                    C_ax = BW_const%Cc
                
                END IF

                f_eta = ONE - BW_const%Ceta * (eta - BW_const%eta0)

                Dstress_VMDstress = (THREE/(TWO*stress_VM)) * stress_dev

                DetaDstress = - (THREE*eta/(TWO*stress_VM**2)) * stress_dev

                angle = lode_angle - ONE_SIXTH * PI

                Id = Identity(NDIM_STRESS)

                DgammaDstress = THREE * GAMMACONST * (TAN(angle)/COS(angle)) *    &
                        &   (ONE/(stress_VM*SIN(THREE*lode_angle+ERROR_SMALL))) * &
                        &   (ONE_THIRD * Id +                                     &
                        &   (COS(THREE*lode_angle)/(2*stress_VM))*stress_dev -    &
                        &   (THREE/(TWO*stress_VM**2)) * MATMUL(stress_dev, stress_dev))
                
                f_stress = Dstress_VMDstress + sigY * f_theta * BW_const%Ceta * &
                        &   DetaDstress - sigY * f_eta * (C_ax - BW_const%Cs) *     &
                        &   (1 - gamma**BW_const%m) * DgammaDstress

                yield_stress = Voigt2(f_stress, 2)
                yield_q = -H * f_eta * f_theta
            
            ELSE 
                
                CALL extract(thisStress, stress_VM=stress_VM,    &
                            &    stress_deviatoric=stress_dev)

                Dstress_VMDstress = (THREE/(TWO*stress_VM)) * stress_dev

                f_stress = Dstress_VMDstress

                yield_stress = Voigt2(f_stress, 2)
                yield_q = -H 

            END IF
        END SUBROUTINE gets_yield_derivative

        SUBROUTINE gets_yield_dd(r_stress, r_q, h_stress, h_q, thisStress, eps)

            USE stress_module,      ONLY: Stress, extract
            USE parameter_module,   ONLY: PI, ONE_SIXTH, TWO, THREE, ERROR_SMALL, &
                                    &     NDIM_STRESS, NDIM_STRESS_Voigt
            USE utility_module,     ONLY: KDelta, Identity4s, IdoxId

            REAL(DP), OPTIONAL, INTENT(OUT) :: r_stress(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP), OPTIONAL, INTENT(OUT) :: r_q(NDIM_STRESS_Voigt)
            REAL(DP), OPTIONAL, INTENT(OUT) :: h_q
            REAL(DP), OPTIONAL, INTENT(OUT) :: h_stress(NDIM_STRESS_Voigt)
            TYPE(Stress),        INTENT(IN) :: thisStress
            REAL(DP),            INTENT(IN) :: eps

            REAL(DP) :: sigY
            REAL(DP) :: H
            REAL(DP) :: ddsVMdds(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: stress_dev(NDIM_STRESS, NDIM_STRESS)
            REAL(DP) :: stress_VM
            
            CALL hardeningfunc(eps=eps, sigY=sigY, H=H)

            IF (ALLOCATED(BW_const)) THEN
                ! Yield function double derivative for BW plasticity

            ELSE
                ! Yield function double derivative for J2 plasticity
                CALL extract(thisStress, stress_VM=stress_VM,       &
                            &      stress_deviatoric=stress_dev)
                CALL ddstress_VMddstress(ddsVMdds, stress_dev, stress_VM)

                IF (PRESENT(r_stress)) r_stress = ddsVMdds
                IF (PRESENT(r_q))      r_q      = ZERO
                IF (PRESENT(h_q))      h_q      = ZERO
                IF (PRESENT(h_stress)) h_stress = ZERO

            END IF

        END SUBROUTINE gets_yield_dd

        PURE SUBROUTINE ddstress_VMddstress(ddsVMdds, stress_dev, stress_VM)

            USE parameter_module,   ONLY: NDIM_STRESS_Voigt, NDIM_STRESS, TWO, THREE
            USE utility_module,     ONLY: KDelta, Identity4s, IdoxId, TensorCrossProduct

            REAL(DP), INTENT(OUT) :: ddsVMdds(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP),  INTENT(IN) :: stress_dev(NDIM_STRESS, NDIM_STRESS)
            REAL(DP),  INTENT(IN) :: stress_VM

            REAL(DP) :: Idev(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)
            REAL(DP) :: SoxS(NDIM_STRESS_Voigt, NDIM_STRESS_Voigt)

            REAL(DP) :: Const1 
            REAL(DP) :: Const2

            Idev = Identity4s(NDIM_STRESS) - ONE_THIRD * IdoxId(NDIM_STRESS)

            SoxS = TensorCrossProduct(stress_dev, stress_dev)

            Const1 = THREE/(TWO * stress_VM)
            Const2 = THREE/(TWO * stress_VM**2)

            ddsVMdds = Const1 * (Idev - Const2 * SoxS) 
            
        END SUBROUTINE ddstress_VMddstress

END MODULE matlib_module
