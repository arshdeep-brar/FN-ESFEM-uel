INCLUDE 'global_modules/parameter_module.f90'
INCLUDE 'global_modules/utility_module.f90'
INCLUDE 'support_modules/log_module.f90'
INCLUDE 'support_modules/stress_module.f90' 
INCLUDE 'setup_modules/input_module.f90'
INCLUDE 'setup_modules/matlib_module.f90'
INCLUDE 'analysis_modules/EdgeSFEM.f90'
INCLUDE 'analysis_modules/material_module.f90'
INCLUDE 'analysis_modules/element_module.f90'
INCLUDE 'output_modules/node_list_module.f90'
INCLUDE 'output_modules/elem_list_module.f90'
INCLUDE 'output_modules/output_module.f90'

SUBROUTINE UEXTERNALDB(LOP, LRESTART, TIME, DTIME, KSTEP, KINC)
    USE parameter_module,   ONLY: DP, DIRLENGTH, EXIT_FUNCTION, MSG_FILE,    &
                            &     STRESS_TYPE
    USE input_module,       ONLY: MeshInfo, indir, ReadMeshInfo, extract
    USE matlib_module,      ONLY: set_material, display_material_status,     &
                            &     deallocate_material
    USE elem_list_module,   ONLY: elem_list, update, empty, set
    USE node_list_module,   ONLY: node_list, update, empty, set
    USE output_module,      ONLY: write_output, outdir
    USE log_module,         ONLY: logsdir
    USE utility_module,     ONLY: newunit     

    INCLUDE 'ABA_PARAM.INC'

    REAL(DP), INTENT(IN) :: TIME(2)
    REAL(DP), INTENT(IN) :: DTIME
    INTEGER,  INTENT(IN) :: LOP, LRESTART, KSTEP, KINC

    ! Internal variable 
    CHARACTER(len=DIRLENGTH) :: cwd
    INTEGER                  :: lenworkdir 

    TYPE(MeshInfo) :: thismesh
    
    ! Variables for mesh information
    INTEGER               :: Npoints
    INTEGER               :: Nelems
    INTEGER,  ALLOCATABLE :: ElemConn(:,:)
    REAL(DP), ALLOCATABLE :: PointCoord(:,:)
    
    ! Number of Nodes for each element
    INTEGER :: Nnodes

    ! Loop Variable 
    INTEGER :: i

    ! Output request
    CHARACTER(len=10) :: outrequest(8)

    ! Initaiting outrequest variable 
    DO i = 1, SIZE(outrequest)
        outrequest(i) = ''
    END DO

    !Initiating variables 
    cwd        = ''
    lenworkdir = 0
    Nnodes     = 0
    Npoints    = 0
    Nelems     = 0

    !--------------------------------------------------------------------
    !                   Element Output Request
    !--------------------------------------------------------------------
    !Add the string to outrequest ARRAY:
    !   STRESS   : For the stress output
    !   PLSTRAIN : For Plastic Strain 
    !   STRAIN   : For total Strain
    !   PSTRESS  : Principal Stress
    !   VMSTRESS : VonMises Stress
    !   ETA      : Stres Triaxiality 
    !   LA       : Lode Angle
    !   LAP      : Lode Angle Parameter
    ! Start here 

    outrequest(1) = 'STRESS'
    outrequest(2) = 'PLSTRAIN'
    outrequest(3) = 'VMSTRESS'
    outrequest(4) = 'TRSTRESS'
    outrequest(5) = 'PSTRESS'
    outrequest(6) = 'EPS'
    outrequest(7) = 'YCOND'

    SELECT CASE (LOP)

        ! At the start of analysis
        CASE (0) 

            indir = ''
            outdir = ''

            CALL GETOUTDIR(cwd, lenworkdir)

            ! Check if the DIRLENGHT is sufficient 
            IF (lenworkdir .GT. DIRLENGTH-30) THEN
                WRITE(MSG_FILE, '(A, I3)') 'DIRError : Increare the DIRLENGHT to ', lenworkdir + 30
                CALL EXIT_FUNCTION
            END IF


            outdir = TRIM(cwd) // '/Outputs/'
            indir  = TRIM(cwd) // '/Input/'
            logsdir = TRIM(cwd) // '/'

            ! ---------------------------------------------------------------------------- !
            !                         READING MESH INFO FILE
            ! ---------------------------------------------------------------------------- !

            CALL ReadMeshInfo(indir, thismesh)

            ! Set Material 
            CALL set_material('Test-4')
            CALL display_material_status

            IF (STRESS_TYPE .EQ. 'K') THEN
                WRITE(*,*) 'NOTE : Kirchoff stress formulation used for analysis'
                WRITE(*,*) ''
            ELSE IF (STRESS_TYPE .EQ. 'C') THEN
                WRITE(*,*) 'NOTE : Cauchy stress formulation used for analysis'
                WRITE(*,*) ''
            ELSE
                WRITE(*,*) 'Unrecognized value of the parameter STRESS_TYPE '
                CALL EXIT_FUNCTION
            END IF

            ! Extract Information from Mesh Info
            CALL extract(thismesh, Npoints=Npoints, Nelem=Nelems)

            ! Allocate ElemConn and PointCoord
            ALLOCATE(ElemConn(Nelems, 4))
            ALLOCATE(PointCoord(Npoints, 2))

            ! Set node and elem data list
            CALL set(node_list, Npoints)
            CALL set(elem_list, Nelems)

            ! Extract ElemConn and PointCoord data 
            CALL extract(thismesh, PointCoord=PointCoord, ElemConn=ElemConn)

            ! Update the Static Data Point data 
            DO i = 1, Npoints
                CALL update(node_list(i), X=PointCoord(i,:))
            END DO

            DO i = 1, Nelems
                Nnodes = COUNT(ElemConn(i,:) .GT. 0) 
                CALL update(elem_list(i), Nnodes=Nnodes, &
                        &   connec=PACK(ElemConn(i,:), ElemConn(i,:).GT.0))
            END DO

        ! Start of each increment
        CASE(1)
            CALL empty(elem_list)
            CALL empty(node_list)

        ! End of each increment
        CASE(2)
            CALL write_output(KINC, KSTEP, outdir, outrequest)
        ! End of analysis
        CASE(3)
            CALL deallocate_material
    END SELECT
    RETURN
END

SUBROUTINE UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS,           &
    &   PROPS, NPROPS, COORD, MCRD, NNODE, U, DU, V, A, JTYPE, TIME,       &
    &   DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF, &
    &   NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROPS,   &
    &   PERIOD)

    USE parameter_module,   ONLY: DP, ZERO, NDIM, NDIM_STRESS_Voigt, ONE_THIRD, &
                            &     NDIM_Voigt, NDIM_STRESS, MSG_FILE, MSGLENGTH, &
                            &     EXIT_FUNCTION, ONE, STRESS_TYPE
    USE node_list_module,   ONLY: node_list, update, extract
    USE elem_list_module,   ONLY: elem_list, update, extract, igp_to_elemList
    USE material_module,    ONLY: IG_point, ASSIGNMENT(=), start_material_logging, &
                            &     end_material_logging
    USE element_module,     ONLY: ESFEM_elem, update, integrate, display, extract
    USE log_module,         ONLY: Logs, set_log_level, set_log_elem, add_log_message, &
                            &     write_logs

    INCLUDE 'ABA_PARAM.INC'

    DIMENSION RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), PROPS(*), SVARS(*),       &
        &   ENERGY(8), COORD(MCRD,NNODE), U(NDOFEL), DU(MLVARX,*), V(NDOFEL), &
        &   TIME(2), PARAMS(*), JDLTYP(MDLOAD,*), ADLMAG(MDLOAD,*),           &
        &   DDLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE), LFLAGS(*), JPROPS(*)
    
    TYPE(ESFEM_elem) :: thisElem
    TYPE(Ig_point)   :: igp

    ! TYPE(Logs)       :: elemlogs

    ! Local copy of the element variables
    REAL(DP) :: NodesCoord(NNODE,2)
    REAL(DP) :: K(NDOFEL, NDOFEL)
    REAL(DP) :: fvec(NDOFEL)
    REAL(DP) :: defU(NDOFEL)
    REAL(DP) :: incDefU(NDOFEL)

    REAL(DP) :: Disp(2)
    INTEGER  :: connec(NNODE) 

    INTEGER :: stat

    !Initiating varibales
    NodesCoord(:,:) = ZERO
    K(:,:)          = ZERO
    fvec(:)         = ZERO
    defU(:)         = ZERO 
    incDefU         = ZERO
    
    Disp(:) = ZERO

    connec(:) = 0
    istat     = 0

    IF (NSVARS .NE. 3*NDIM_STRESS_Voigt+1) THEN
        WRITE(*,*) 'Incorrect dimension of SVARS array'
        WRITE(*,*) ' NSVARS = 3*NDIM_STRESS_Voigt+1'
        CALL EXIT_FUNCTION
    END IF

    ! CALL set_log_level(elemlogs, level=2)
    ! CALL set_log_elem(elemlogs, elem=JELEM)

    ! CALL add_log_message(elemlogs, debug="UEL initialized")

    !Copying Node coordinates to the local copy
    DO i = 1,NNODE
        NodesCoord(i,1) = COORD(1,i)
        NodesCoord(i,2) = COORD(2,i)
    END DO

    DO i = 1,NDOFEL
        incDefU(i) = DU(i,1)
        defU(i)    = U(i) - DU(i,1) 
    END DO

    CALL extract(elem_list(JELEM), connec=connec)
    
    ! Updating the integration point varibales with the SVAR array
    igp = SVARS(1:NSVARS)
    
    ! Update element deformation object 
    CALL update(thisElem, NodesCoord, defU, incDefU, igp, stat)

    IF (stat .NE. 0) THEN
        WRITE(*,*) 'Unable to update the element'
    !     CALL add_log_message(elemlogs, &
    ! &   err='Eigenvalue calculation Failed in updating material variables(LAPACK ERROR!!)')
        CALL EXIT_FUNCTION
    END IF

    ! CALL add_log_message(elemlogs, debug="Object update completed")

    IF (JELEM .EQ. 1) THEN
        CALL start_material_logging(level=2, elem=JELEM, inc=KINC)
    END IF
    CALL integrate(thisElem, K, fvec, STRESS_TYPE) 
    IF (JELEM .EQ. 1) THEN
        CALL end_material_logging
    END IF

    CALL extract(thisElem, igp=igp)

    SVARS(1:NSVARS) = igp
    
    
    ! ------------------------------------------------------------------ !
    !                     Updating the output
    ! ------------------------------------------------------------------ !

    CALL igp_to_elemList(elem_list(JELEM), igp)

    IF (JTYPE .EQ. 2003) THEN
        CALL update(node_list(connec(1)), u=U(1:2))
        CALL update(node_list(connec(3)), u=U(3:4))
        Disp = ZERO
        DO i = 1, 3
            Disp(1) = Disp(1) + ONE_THIRD * U(2*i-1)
            Disp(2) = Disp(2) + ONE_THIRD * U(2*i)
        END DO
        CALL update(node_list(connec(2)), u=Disp)
    
    ELSE IF (JTYPE .EQ. 2004) THEN
        CALL update(node_list(connec(1)), u = U(1:2))
        Disp=ZERO
        Disp(1) = ONE_THIRD * (U(1) + U(3) + U(5))
        Disp(2) = ONE_THIRD * (U(2) + U(4) + U(6))
        CALL update(node_list(connec(2)), u = Disp)
        CALL update(node_list(connec(3)), u = U(5:6))
        Disp=ZERO
        Disp(1) = ONE_THIRD * (U(1) + U(5) + U(7))
        Disp(2) = ONE_THIRD * (U(2) + U(6) + U(8))
        CALL update(node_list(connec(4)), u = Disp)
    END IF
    ! ----------------------------------------------------------------------------!
    !                        End of the Output update
    ! ----------------------------------------------------------------------------!
    
    ! CALL add_log_message(elemlogs, debug='Element output transfered')
    
    RHS(:,1) = -fvec(:)
    AMATRX = K
    
    ! CALL write_logs(elemlogs)

    RETURN
END
