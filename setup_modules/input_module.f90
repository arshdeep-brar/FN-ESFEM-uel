MODULE input_module
    USE parameter_module, ONLY: DP, ZERO, ONE, ONE_THIRD, SIX, DIRLENGTH, &
                            &   MSG_FILE, EXIT_FUNCTION

    IMPLICIT NONE
    PRIVATE

    TYPE, PUBLIC :: MeshInfo
        INTEGER               :: Npoints
        INTEGER               :: Nelem
        REAL(DP), ALLOCATABLE :: PointCoord(:,:)
        INTEGER,  ALLOCATABLE :: ElemConn(:,:)
    END TYPE MeshInfo

    INTERFACE extract
        MODULE PROCEDURE extract_MeshInfo
    END INTERFACE extract
    
    CHARACTER(len=DIRLENGTH), PUBLIC, SAVE :: indir

    PUBLIC :: extract, ReadMeshInfo

    CONTAINS

        PURE SUBROUTINE extract_MeshInfo(thismesh, Npoints, Nelem, PointCoord, ElemConn)

            TYPE(MeshInfo),     INTENT(IN) :: thismesh
            INTEGER,  OPTIONAL, INTENT(OUT) :: Npoints
            INTEGER,  OPTIONAL, INTENT(OUT) :: Nelem
            REAL(DP), OPTIONAL, INTENT(OUT) :: PointCoord(:,:)
            INTEGER,  OPTIONAL, INTENT(OUT) :: ElemConn(:,:)

            IF (PRESENT(Npoints)) Npoints = thismesh%Npoints

            IF (PRESENT(Nelem)) Nelem = thismesh%Nelem

            IF ((PRESENT(PointCoord)) .AND. (ALLOCATED(thismesh%PointCoord))) &
                    &   PointCoord = thismesh%PointCoord

            IF ((PRESENT(ElemConn)) .AND. (ALLOCATED(thismesh%ElemConn))) & 
                    &   ElemConn = thismesh%ElemConn
   
        END SUBROUTINE extract_MeshInfo

        SUBROUTINE ReadMeshInfo(dir, Topo)
            USE utility_module, ONLY: newunit

            CHARACTER(len=DIRLENGTH), INTENT(IN) :: dir
            TYPE(MeshInfo),          INTENT(OUT) :: Topo

            ! Local varaibles
            CHARACTER(len=DIRLENGTH) :: filename
            CHARACTER(len=13) :: shortchar
            CHARACTER(len=60) :: longchar

            INTEGER  :: Nnodes
            INTEGER  :: Ncentroid
            INTEGER  :: Nelem
            INTEGER  :: Npoints
            INTEGER  :: thisunit
            INTEGER  :: iostat
            INTEGER  :: i
            INTEGER  :: elemid, pointid
            INTEGER  :: ival1, ival2, ival3, ival4
            REAL(DP) :: rval1, rval2    

            thisunit = newunit()            

            filename = TRIM(dir) // 'MeshInfo.txt'

            OPEN(UNIT=thisunit, FILE=filename, STATUS='OLD', IOSTAT=iostat)

            IF (iostat .NE. 0) THEN
                WRITE(MSG_FILE, *) 'FILE ERROR : Error in opening MeshInfo file'
                CALL EXIT_FUNCTION
                RETURN
            END IF
            
            READ(thisunit, *)shortchar, Nelem
            READ(thisunit, *)shortchar, Nnodes
            READ(thisunit, *)shortchar, Ncentroid

            Npoints = Nnodes + Ncentroid

            ALLOCATE(Topo%PointCoord(Npoints,2))
            ALLOCATE(Topo%ElemConn(Nelem,4))

            Topo%Nelem = Nelem
            Topo%Npoints = Npoints

            READ(thisunit, *) longchar

            DO i =  1, Npoints
                READ(thisunit, *) pointid, rval1, rval2
                Topo%PointCoord(pointid, 1) = rval1
                Topo%PointCoord(pointid, 2) = rval2
            END DO

            READ(thisunit, *) longchar

            DO i = 1, Nelem
                READ(thisunit,*) elemid, ival1, ival2, ival3, ival4
                Topo%ElemConn(elemid, 1) = ival1
                Topo%ElemConn(elemid, 2) = ival2
                Topo%ElemConn(elemid, 3) = ival3
                Topo%ElemConn(elemid, 4) = ival4 
            END DO

            CLOSE(thisunit)

        END SUBROUTINE ReadMeshInfo

END MODULE input_module