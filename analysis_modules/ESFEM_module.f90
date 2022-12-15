MODULE ESFEM_module
    
    USE parameter_module,   ONLY: DP
    
    IMPLICIT NONE

    !Defining a point in the domain
    TYPE, PRIVATE :: Point2D
        PRIVATE
        REAL(DP) :: x
        REAl(DP) :: y
    END TYPE Point2D

    !Data Structure to store a length and normal of a line 
    ! segment in 2D space
    TYPE, PRIVATE :: lineseg
        PRIVATE
        REAL(DP) :: l
        REAL(DP) :: nx
        REAL(DP) :: ny
    END TYPE lineseg

    TYPE, PUBLIC :: SmoothedDomain
        PRIVATE
        INTEGER                    :: Npoints 
        TYPE(Point2D), ALLOCATABLE :: Points(:)
        TYPE(Point2D), ALLOCATABLE :: IGPoints(:)
        TYPE(lineseg), ALLOCATABLE :: Lines(:)
    END TYPE SmoothedDomain
    
    CONTAINS
    
END MODULE ESFEM_module