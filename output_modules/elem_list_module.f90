MODULE elem_list_module 

    USE parameter_module, ONLY: DP, ZERO, NDIM

    IMPLICIT NONE
    
    PRIVATE

    TYPE, PUBLIC :: elem_vars
        PRIVATE
        INTEGER              :: Nnodes
        INTEGER, ALLOCATABLE :: connec(:)
        REAL(DP)             :: Stress(6)      = ZERO
        REAL(DP)             :: sig_eff        = ZERO
        REAL(DP)             :: Sig12(3)       = ZERO
        REAL(DP)             :: eta            = ZERO
        REAL(DP)             :: la             = ZERO
        REAL(DP)             :: lap            = ZERO
        REAL(DP)             :: trialStress(6) = ZERO
        REAL(DP)             :: PlStress(6)    = ZERO
        REAL(DP)             :: totStrain(6)   = ZERO
        REAL(DP)             :: plStrain(6)    = ZERO
        REAL(DP)             :: eps            = ZERO
        INTEGER              :: YCond          = 0
    END TYPE elem_vars

    INTERFACE update
        MODULE PROCEDURE update_elem_output
    END INTERFACE update

    INTERFACE extract
        MODULE PROCEDURE extract_elem_output
    END INTERFACE extract

    INTERFACE empty
        MODULE PROCEDURE empty_elem_list
    END INTERFACE empty

    INTERFACE set
        MODULE PROCEDURE set_elem_list
    END INTERFACE set

    INTERFACE display
        MODULE PROCEDURE display_elem_output
    END INTERFACE display

    TYPE(elem_vars), ALLOCATABLE, PUBLIC, SAVE :: elem_list(:)

    PUBLIC :: update, extract, empty, set, display, igp_to_elemList
    
    CONTAINS 

        SUBROUTINE update_elem_output(thiselem, Nnodes, connec, Stress, Sig12, sig_eff, &
                                &     eta, la, lap, trialStress, PlStress, totStrain,   &
                                &     plStrain, eps, YCond)

            TYPE(elem_vars), INTENT(INOUT) :: thiselem
            INTEGER,  OPTIONAL, INTENT(IN) :: Nnodes
            INTEGER,  OPTIONAL, INTENT(IN) :: connec(:)
            REAL(DP), OPTIONAL, INTENT(IN) :: Stress(6)
            REAL(DP), OPTIONAL, INTENT(IN) :: Sig12(3)
            REAL(DP), OPTIONAL, INTENT(IN) :: sig_eff
            REAL(DP), OPTIONAL, INTENT(IN) :: eta
            REAL(DP), OPTIONAL, INTENT(IN) :: la
            REAL(DP), OPTIONAL, INTENT(IN) :: lap  
            REAL(DP), OPTIONAL, INTENT(IN) :: trialStress(6)
            REAL(DP), OPTIONAL, INTENT(IN) :: PlStress(6)
            REAL(DP), OPTIONAL, INTENT(IN) :: totStrain(6)
            REAL(DP), OPTIONAL, INTENT(IN) :: plStrain(6)
            REAL(DP), OPTIONAL, INTENT(IN) :: eps
            INTEGER,  OPTIONAL, INTENT(IN) :: YCond          

            IF (PRESENT(Nnodes)) thiselem%Nnodes = Nnodes
            
            IF (.NOT. ALLOCATED(thiselem%connec)) ALLOCATE(thiselem%connec(Nnodes))
            IF (PRESENT(connec) .AND. ALLOCATED(thiselem%connec)) thiselem%connec = connec
 
            IF (PRESENT(Stress))      thiselem%Stress      = Stress
            IF (PRESENT(Sig12))       thiselem%Sig12       = Sig12
            IF (PRESENT(sig_eff))     thiselem%sig_eff     = sig_eff
            IF (PRESENT(eta))         thiselem%eta         = eta
            IF (PRESENT(la))          thiselem%la          = la
            IF (PRESENT(lap))         thiselem%lap         = lap
            IF (PRESENT(trialStress)) thiselem%trialStress = trialStress
            IF (PRESENT(PlStress))    thiselem%PlStress    = PlStress
            IF (PRESENT(totStrain))   thiselem%totStrain   = totStrain
            IF (PRESENT(plStrain))    thiselem%plStrain    = plStrain
            IF (PRESENT(eps))         thiselem%eps         = eps
            IF (PRESENT(YCond))       thiselem%YCond       = YCond

        END SUBROUTINE update_elem_output
        
        PURE SUBROUTINE extract_elem_output(thiselem, Nnodes, connec, Stress, Sig12, sig_eff, &
                                      &     eta, la, lap, trialStress, PlStress, totStrain,   &
                                      &     plStrain, eps, YCond)

            TYPE(elem_vars),     INTENT(IN) :: thiselem
            INTEGER,  OPTIONAL, INTENT(OUT) :: Nnodes
            INTEGER,  OPTIONAL, INTENT(OUT) :: connec(:)
            REAL(DP), OPTIONAL, INTENT(OUT) :: Stress(6)
            REAL(DP), OPTIONAL, INTENT(OUT) :: Sig12(3)
            REAL(DP), OPTIONAL, INTENT(OUT) :: sig_eff
            REAL(DP), OPTIONAL, INTENT(OUT) :: eta
            REAL(DP), OPTIONAL, INTENT(OUT) :: la
            REAL(DP), OPTIONAL, INTENT(OUT) :: lap
            REAL(DP), OPTIONAL, INTENT(OUT) :: trialStress(6)
            REAL(DP), OPTIONAL, INTENT(OUT) :: PlStress(6)
            REAL(DP), OPTIONAL, INTENT(OUT) :: totStrain(6)
            REAL(DP), OPTIONAL, INTENT(OUT) :: plStrain(6)
            REAL(DP), OPTIONAL, INTENT(OUT) :: eps
            INTEGER,  OPTIONAL, INTENT(OUT) :: YCond

            IF (PRESENT(Nnodes))      Nnodes      = thiselem%Nnodes 
            IF (PRESENT(connec))      connec      = thiselem%connec
            IF (PRESENT(Stress))      Stress      = thiselem%Stress 
            IF (PRESENT(Sig12))       Sig12       = thiselem%Sig12
            IF (PRESENT(sig_eff))     sig_eff     = thiselem%sig_eff
            IF (PRESENT(eta))         eta         = thiselem%eta 
            IF (PRESENT(la))          la          = thiselem%la
            IF (PRESENT(lap))         lap         = thiselem%lap
            IF (PRESENT(trialStress)) trialStress = thiselem%trialStress
            IF (PRESENT(PlStress))    PlStress    = thiselem%PlStress
            IF (PRESENT(totStrain))   totStrain   = thiselem%totStrain
            IF (PRESENT(plStrain))    plStrain    = thiselem%plStrain
            IF (PRESENT(eps))         eps         = thiselem%eps 
            IF (PRESENT(YCond))       YCond       = thiselem%YCond

        END SUBROUTINE extract_elem_output

        SUBROUTINE igp_to_elemList(this, igp)
            
            USE material_module,    ONLY: IG_point, extract
            USE stress_module,      ONLY: Stress, extract
            USE parameter_module,   ONLY: NDIM_STRESS_Voigt, EXIT_FUNCTION, &
                                    &     ERROR

            TYPE(elem_vars), INTENT(INOUT) :: this
            TYPE(IG_point),     INTENT(IN) :: igp

            TYPE(Stress) :: sigma
            REAL(DP)     :: stress_v(NDIM_STRESS_Voigt)
            REAL(DP)     :: stress_principal(3)
            REAL(DP)     :: stress_VM
            REAL(DP)     :: stress_triaxiality
            REAL(DP)     :: lode_angle
            REAL(DP)     :: lode_angle_parameter
            REAL(DP)     :: trialStress(NDIM_STRESS_Voigt)
            REAL(DP)     :: PlStress(NDIM_STRESS_Voigt)
            REAL(DP)     :: totStrain(NDIM_STRESS_Voigt)
            REAL(DP)     :: plStrain(NDIM_STRESS_Voigt)
            REAL(DP)     :: eps

            INTEGER :: info 

            CALL extract(igp, sigma=sigma, Ep=plStrain, strain=totStrain, &
                    &    trialStress=trialStress, eps=eps)

            CALL extract(sigma, stress_v=stress_v, stress_VM=stress_VM,  &
                    &    stress_principal=stress_principal,              &
                    &    stress_triaxiality=stress_triaxiality,          &
                    &    lode_angle=lode_angle, lode_angle_parameter=    &
                    &    lode_angle_parameter, info=info)
            
            IF (info .NE. 0) THEN
                WRITE(*,*) "Unable to calculate Principal Stress"
                CALL EXIT_FUNCTION
            END IF
            
            PlStress = trialStress - stress_v

            IF (ABS(SUM(PlStress*PlStress)-ERROR) .EQ. ZERO) THEN 
                this%YCond = 0
            ELSE
                this%YCond = 1 
            END IF

            this%Sig12   = stress_principal
            this%sig_eff = stress_VM
            this%eta     = stress_triaxiality
            this%la      = lode_angle
            this%lap     = lode_angle_parameter
            this%eps     = eps

            IF (NDIM_STRESS_Voigt .EQ. 3) THEN
                this%Stress      = ZERO
                this%Stress(1:2) = stress_v(1:2)
                this%Stress(6)   = stress_v(NDIM_STRESS_Voigt)

                this%trialStress      = ZERO
                this%trialStress(1:2) = trialStress(1:2)
                this%trialStress(6)   = trialStress(NDIM_STRESS_Voigt)
                
                this%PlStress      = ZERO
                this%PlStress(1:2) = PlStress(1:2)
                this%PlStress(6)   = PlStress(NDIM_STRESS_Voigt)
                
                this%totStrain      = ZERO
                this%totStrain(1:2) = totStrain(1:2)
                this%totStrain(6)   = totStrain(NDIM_STRESS_Voigt)
                
                this%plStrain      = ZERO
                this%plStrain(1:2) = plStrain(1:2)
                this%plStrain(6)   = plStrain(NDIM_STRESS_Voigt)
                
            ELSE IF (NDIM_STRESS_Voigt .EQ. 6) THEN
                this%Stress(1:NDIM_STRESS_Voigt)      = stress_v (1:NDIM_STRESS_Voigt)
                this%trialStress(1:NDIM_STRESS_Voigt) = trialStress(1:NDIM_STRESS_Voigt)
                this%PlStress(1:NDIM_STRESS_Voigt)    = PlStress(1:NDIM_STRESS_Voigt)
                this%totStrain(1:NDIM_STRESS_Voigt)   = totStrain(1:NDIM_STRESS_Voigt)
                this%plStrain(1:NDIM_STRESS_Voigt)    = plStrain(1:NDIM_STRESS_Voigt)
            END IF

        END SUBROUTINE igp_to_elemList 
         
        SUBROUTINE display_elem_output(thiselem)

            TYPE(elem_vars), INTENT(IN) :: thiselem

            CHARACTER(len=10) :: fmtconn
            INTEGER :: i

            fmtconn = ''

            WRITE(fmtconn, '(A, I1, A)') '(', thiselem%Nnodes, 'I5)'

            WRITE(*,'(1X, A, I5)') 'Number of Nodes : ', thiselem%Nnodes
            
            WRITE(*,'(1X, A)') 'Element Topology Connectivity'

            WRITE(*,fmtconn) (thiselem%connec(i), i=1,SIZE(thiselem%connec,1)) 

        END SUBROUTINE display_elem_output
        
        PURE SUBROUTINE set_elem_list(thislist, Nelem)

            TYPE(elem_vars), ALLOCATABLE, INTENT(INOUT) :: thislist(:)
            INTEGER,                           INTENT(IN) :: Nelem

            IF (.NOT. ALLOCATED(thislist)) ALLOCATE(thislist(Nelem))

        END SUBROUTINE set_elem_list
        
        SUBROUTINE empty_elem_list(thislist)

            TYPE(elem_vars), INTENT(INOUT) :: thislist(:)

            INTEGER ::  i

            DO i = 1, SIZE(thislist, 1)
                thislist(i)%Stress      = ZERO
                thislist(i)%sig_eff     = ZERO
                thislist(i)%Sig12       = ZERO
                thislist(i)%eta         = ZERO
                thislist(i)%la          = ZERO
                thislist(i)%lap         = ZERO
                thislist(i)%trialStress = ZERO
                thislist(i)%PlStress    = ZERO
                thislist(i)%totStrain   = ZERO
                thislist(i)%plStrain    = ZERO
                thislist(i)%eps         = ZERO
                thislist(i)%YCond       = 0
            END DO

        END SUBROUTINE empty_elem_list

END MODULE elem_list_module
