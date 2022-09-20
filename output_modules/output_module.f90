MODULE output_module
    USE parameter_module, ONLY: DP, DIRLENGTH, NDIM, ZERO, MSG_FILE, &  
                            &   NUM_SMALL

    IMPLICIT NONE
    PRIVATE

    CHARACTER(len=DIRLENGTH), PUBLIC, SAVE :: outdir

    PUBLIC :: write_output
    
    
    CONTAINS

        SUBROUTINE write_output(Kinc, KStep, dir, outrequest)
            
            USE elem_list_module, ONLY: elem_vars, elem_list, extract
            USE node_list_module, ONLY: node_list, extract
            USE utility_module,   ONLY: newunit

            INTEGER,                     INTENT(IN) :: Kinc
            INTEGER,                     INTENT(IN) :: KStep
            CHARACTER(len=DIRLENGTH),    INTENT(IN) :: dir
            CHARACTER(len=10), OPTIONAL, INTENT(IN) :: outrequest(:)

            ! Local Variables 
            INTEGER :: thisunit
            INTEGER :: i

            ! File vartiables 
            CHARACTER(len=DIRLENGTH) :: outfile
            CHARACTER(len=DIRLENGTH) :: outnum
            
            ! Format variables
            CHARACTER(len=10)        :: FMTNNODE, FMTNELEM, FMTKINC, FMTFLOAT

            ! Integer variables for output data
            INTEGER :: Npoints, Nelems
            INTEGER :: Nsize

            ! Nodal Coordinates
            REAL(DP) :: X(NDIM)
            ! Nodal displacements
            REAL(DP) :: disp(NDIM)

            ! File check
            LOGICAL :: file_exist

            thisunit = 0
            Nsize = 0
            Nelems = 0
            Npoints = 0
            outnum = ''
            outfile = ''
            FMTNNODE = ''
            FMTNELEM = ''
            FMTKINC = ''
            FMTFLOAT = ''

            FMTNNODE = 'I5' 
            FMTNELEM = 'I5'
            FMTKINC = 'I1.1,I5.5'
            FMTFLOAT = 'ES20.10'
            
            thisunit = newunit()

            ! outnum of the file based on the increment number with 
            ! 10 placeholder for the increment number
            WRITE(outnum, '('//TRIM(FMTKINC)//')') KStep, Kinc

            ! Full address of the output file 
            outfile = TRIM(dir) // 'esfem_'//TRIM(outnum) // '.vtk'

            INQUIRE(FILE=outfile, EXIST=file_exist)

            IF (file_exist) THEN
                OPEN(UNIT=thisunit, FILE=outfile, STATUS='REPLACE', ACTION='WRITE')
            ELSE
                OPEN(UNIT=thisunit, FILE=outfile, STATUS='NEW', ACTION='WRITE')
            END IF

            ! VTK file header
            WRITE(thisunit, '(A)') '# vtk DataFile Version 3.1'
            WRITE(thisunit, '(A)') 'ESFEM Model output'

            ! File format
            WRITE(thisunit, '(A)') 'ASCII'

            ! vtk data type
            WRITE(thisunit, '(A)') 'DATASET UNSTRUCTURED_GRID' 

            !-------------------------------------------------------------------
            !               Writing points data (Nodes+Centroids)
            !-------------------------------------------------------------------
            ! Total Number of points
            Npoints = SIZE(node_list)
            WRITE(thisunit, '(A, '//TRIM(FMTNNODE) //', A)')'POINTS', Npoints, ' double'
            DO i = 1, Npoints
                CALL extract(node_list(i), x=X)
                IF (NDIM .EQ. 2) THEN
                    ! Check the dimension of the problem for 2d problems use
                    WRITE(thisunit, '(3' // TRIM(FMTFLOAT) // ')') X(1), X(2), ZERO
                ! ELSE IF(NDIM .EQ. 3) THEN
                !     WRITE(thisunit, '(3' // TRIM(FMTFLOAT) // ')') X(1), X(2), X(3)
                ELSE
                    WRITE(MSG_FILE, '(A)') 'OUTPUT ERROR : Incorrent dimension'
                END IF
            END DO
            WRITE(thisunit, '(A)') ''


            !----------------------------------------------------------------------
            !               Writing element topology data                           
            !----------------------------------------------------------------------
            Nelems = SIZE(elem_list)

            ! Calculating the size of nodes used for Element output
            CALL get_size(elem_list, Nsize)
            ! Element connections with the points, vtk format
            WRITE(thisunit, '(A, ' //TRIM(FMTNELEM) // ', I10)') 'CELLS ', Nelems, Nsize
            CALL writeElementOutput('connec') ! flag : connec, writes the point connectivity data 
            
            ! Type of cell for each element
            WRITE(thisunit, '(A, ' // TRIM(FMTNELEM) // ')') 'CELL_TYPES', Nelems
            CALL writeElementOutput('eltype') ! flag : eltype, writes the type of element based on the 
                                              !  number of nodes
            
            !-----------------------------------------------------------------------
            !                       Point Data begins
            !-----------------------------------------------------------------------

            WRITE(thisunit, '(A, ' //TRIM(FMTNNODE) // ')') 'POINT_DATA', Npoints

            !-----------------------------------------------------------------------
            !                   Writing Point Displacement
            !-----------------------------------------------------------------------
            WRITE(thisunit, '(A)') 'VECTORS displacement double'
            DO i = 1, Npoints
                CALL extract(node_list(i), u=disp)
                IF (NDIM .EQ. 2) THEN
                    ! Check the dimension of the problem for 2d problems use
                    IF (ABS(disp(1)) .LT. NUM_SMALL) disp(1) = ZERO
                    IF (ABS(disp(2)) .LT. NUM_SMALL) disp(2) = ZERO
                    WRITE(thisunit, '(3' // TRIM(FMTFLOAT) // ')') disp(1), disp(2), ZERO
                ! ELSE IF(NDIM .EQ. 3) THEN
                !     WRITE(thisunit, '(3' // TRIM(FMTFLOAT) // ')') disp(1), disp(2), disp(3)
                ELSE
                    WRITE(MSG_FILE, '(A)') 'OUTPUT ERROR : Incorrent dimension'
                END IF
            END DO
            WRITE(thisunit, '(A)') ''

            !------------------------------------------------------------------------
            !                      Element Data Begins
            !------------------------------------------------------------------------

            WRITE(thisunit, '(A, ' // TRIM(FMTNELEM) // ')') 'CELL_DATA', Nelems
    
            IF (PRESENT(outrequest)) THEN
                DO i = 1, SIZE(outrequest)
                    IF (outrequest(i) .NE. '') THEN
                        CALL writeElementOutput(outrequest(i))
                    END IF
                END DO
            END IF

            CLOSE(thisunit)

            CONTAINS
                PURE SUBROUTINE get_size(elemlist, N)
                
                    TYPE(elem_vars), INTENT(IN) :: elemlist(:)
                    INTEGER,        INTENT(OUT) :: N

                    INTEGER :: Nnodes
                    INTEGER :: iter

                    N = 0

                    DO iter = 1, SIZE(elemlist)
                        CALL extract(elemlist(iter), Nnodes=Nnodes)
                        N = N + Nnodes + 1
                    END DO

                END SUBROUTINE get_size

                SUBROUTINE writeElementOutput(outflag)

                    CHARACTER(len=*), INTENT(IN) :: outflag

                    INTEGER              :: Nnodes
                    INTEGER, ALLOCATABLE :: connec(:)
                    REAL(DP)             :: Stress(6)
                    REAL(DP)             :: trialStress(6)
                    REAL(DP)             :: plStrain(6)
                    REAL(DP)             :: PlStress(6)
                    REAL(DP)             :: totStrain(6)
                    REAL(DP)             :: Sig12(3)
                    REAL(DP)             :: sig_eff
                    REAL(DP)             :: eps
                    REAL(DP)             :: eta
                    REAL(DP)             :: la
                    REAL(DP)             :: lap
                    INTEGER              :: YCond

                    INTEGER :: thiselem
                    INTEGER :: j
                    REAL(DP) :: PStress(6)
                    PStress = ZERO

                    SELECT CASE (TRIM(ADJUSTL(outflag)))

                    CASE ('STRESS')
                        WRITE(thisunit, '(A)') 'TENSORS stress double'
                    CASE ('PLSTRAIN')
                        WRITE(thisunit, '(A)') 'TENSORS PlasticStrain double'
                    CASE ('TRSTRESS')
                        WRITE(thisunit, '(A)') 'TENSORS TrialStress double'
                    CASE('STRAIN')
                        WRITE(thisunit, '(A)') 'TENSORS TotalStrain double'
                    CASE ('PLSTRESS')
                        WRITE(thisunit, '(A)') 'TENSORS PlStress double'
                    CASE('PSTRESS')
                        WRITE(thisunit, '(A)') 'TENSORS PrincipalStress double'
                    CASE('VMSTRESS')
                        WRITE(thisunit, '(A)') 'SCALARS VMStress double'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    CASE('EPS')
                        WRITE(thisunit, '(A)') 'SCALARS EffectivePlasticStrain double'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    CASE('ETA')
                        WRITE(thisunit, '(A)') 'SCALARS StressTriaxaility double'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    CASE('LA')
                        WRITE(thisunit, '(A)') 'SCALARS LodeAngle double'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    CASE('LAP')
                        WRITE(thisunit, '(A)') 'SCALARS LodeAngleParameter double'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    CASE('YCOND')
                        WRITE(thisunit, '(A)') 'SCALARS YieldCondition int'
                        WRITE(thisunit, '(A)') 'LOOKUP_TABLE default'
                    END SELECT

                    DO thiselem = 1, SIZE(elem_list)
                        
                        CALL extract(elem_list(thiselem), Nnodes = Nnodes)

                        ALLOCATE(connec(Nnodes))

                        SELECT CASE (TRIM(ADJUSTL(outflag)))

                        CASE('connec')
                            CALL extract(elem_list(thiselem), connec=connec)

                            WRITE(thisunit, '(' // TRIM(FMTNNODE) // ')', ADVANCE='no') Nnodes
                            DO j = 1, Nnodes
                                WRITE(thisunit, '(' // TRIM(FMTNNODE) // ')', ADVANCE='no') connec(j)-1
                            END DO
                            WRITE(thisunit, '(A)') ''

                        CASE('eltype')
                            IF (Nnodes .EQ. 3) THEN
                                WRITE(thisunit, '(I2)') 5 ! For triangle elements
                            ELSE IF (Nnodes .EQ. 4) THEN
                                WRITE(thisunit, '(I2)') 9 ! For quad element
                            ELSE
                                WRITE(MSG_FILE, '(A)') 'OUTPUT ERROR : Incorrect number of nodes' 
                            END IF
                        
                        CASE('STRESS')
                            CALL extract(elem_list(thiselem), Stress=Stress)
                            CALL writeTensor(Stress)

                        CASE('PLSTRAIN')
                            CALL extract(elem_list(thiselem), plStrain= plStrain)
                            CALL writeTensor(plStrain)
                        
                        CASE('PLSTRESS')
                            CALL extract(elem_list(thiselem), PlStress= PlStress)
                            CALL writeTensor(PlStress)
                        
                        CASE('TRSTRESS')
                            CALL extract(elem_list(thiselem), trialStress=trialStress)
                            CALL writeTensor(trialStress)    

                        CASE('STRAIN')
                            CALL extract(elem_list(thiselem), totStrain = totStrain)
                            CALL writeTensor(totStrain)

                        CASE('PSTRESS')
                            CALL extract(elem_list(thiselem), Sig12=Sig12)
                            PStress(1:3) = Sig12
                            CALL writeTensor(PStress)
                            
                        CASE('VMSTRESS')
                            CALL extract(elem_list(thiselem), sig_eff=sig_eff)
                            WRITE(thisunit, '(' // TRIM(FMTFLOAT) // ')') sig_eff

                        CASE('EPS')
                            CALL extract(elem_list(thiselem), eps=eps)
                            WRITE(thisunit, '(' // TRIM(FMTFLOAT) // ')') eps
                        
                        CASE('ETA')
                            CALL extract(elem_list(thiselem), eta=eta)
                            WRITE(thisunit, '(' // TRIM(FMTFLOAT) // ')') eta
                        
                        CASE('LA')
                            CALL extract(elem_list(thiselem), la=la)
                            WRITE(thisunit, '(' // TRIM(FMTFLOAT) // ')') la
                        
                        CASE('LAP')
                            CALL extract(elem_list(thiselem), lap=lap)
                            WRITE(thisunit, '(' // TRIM(FMTFLOAT) // ')') lap
                        
                        CASE('YCOND')
                            CALL extract(elem_list(thiselem), YCond=YCond)
                            WRITE(thisunit, '(I2)') YCond
                        
                        END SELECT

                        DEALLOCATE(connec)

                    END DO

                    WRITE(thisunit, '(A)') ''

                END SUBROUTINE writeElementOutput

                
                SUBROUTINE writeTensor(thisvect)
                    
                    REAL(DP), INTENT(IN) :: thisvect(6)
                    
                    REAL(DP) :: tensor(3,3)

                    INTEGER :: l

                    tensor = ZERO

                    tensor(1,1) = thisvect(1)
                    tensor(2,2) = thisvect(2)
                    tensor(3,3) = thisvect(3)
                    tensor(2,3) = thisvect(4)
                    tensor(1,3) = thisvect(5)
                    tensor(1,2) = thisvect(6)
                    tensor(3,2) = tensor(2,3)
                    tensor(3,1) = tensor(1,3)
                    tensor(2,1) = tensor(1,2)
                    
                    DO l = 1,3
                        WRITE(thisunit, '(3' // TRIM(FMTFLOAT) // ')') tensor(l,1), tensor(l,2), tensor(l,3)
                    END DO
                    WRITE(thisunit, '(A)') ' '

                END SUBROUTINE writeTensor

        END SUBROUTINE write_output

END MODULE output_module
