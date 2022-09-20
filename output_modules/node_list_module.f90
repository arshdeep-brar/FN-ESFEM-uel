MODULE node_list_module

    USE parameter_module, ONLY: DP, ZERO, NDIM

    IMPLICIT NONE

    PRIVATE

    TYPE, PUBLIC ::node_output
        PRIVATE
        REAL(DP) :: x(NDIM)
        REAL(DP) :: u(NDIM) 
    END TYPE node_output

    INTERFACE update
        MODULE PROCEDURE update_node_output
    END INTERFACE update

    INTERFACE extract
        MODULE PROCEDURE extract_node_output
    END INTERFACE extract

    INTERFACE set
        MODULE PROCEDURE set_node_list
    END INTERFACE set

    INTERFACE empty
        MODULE PROCEDURE empty_node_list
    END INTERFACE empty

    TYPE(node_output), ALLOCATABLE, PUBLIC, SAVE :: node_list(:)

    PUBLIC :: update, set, empty, extract

    CONTAINS

        PURE SUBROUTINE set_node_list(thislist, Npoints)
            !Purpose:
            ! The subroutine sets up the node_list for the model
            ! and allocates the size for the node_list array

            TYPE(node_output), ALLOCATABLE, INTENT(INOUT) :: thislist(:)
            INTEGER,                           INTENT(IN) :: Npoints

            IF (.NOT. ALLOCATED(thislist)) ALLOCATE(thislist(Npoints))

        END SUBROUTINE set_node_list

        PURE SUBROUTINE update_node_output(node, x, u)
            
            TYPE(node_output), INTENT(INOUT) :: node
            REAL(DP),   OPTIONAL, INTENT(IN) :: x(NDIM)
            REAL(DP),   OPTIONAL, INTENT(IN) :: u(NDIM)

            IF (PRESENT(x)) node%x = x

            IF (PRESENT(u)) node%u = u
 
        END SUBROUTINE update_node_output 

        PURE SUBROUTINE extract_node_output(node, x, u)
            
            TYPE(node_output),     INTENT(IN) :: node
            REAL(DP),   OPTIONAL, INTENT(OUT) :: x(NDIM)
            REAL(DP),   OPTIONAL, INTENT(OUT) :: u(NDIM)

            IF (PRESENT(x)) x = node%x

            IF (PRESENT(u)) u = node%u
 
        END SUBROUTINE extract_node_output 

        PURE SUBROUTINE empty_node_list(thislist)
            
            TYPE(node_output), INTENT(INOUT) :: thislist(:)

            INTEGER :: i
            REAL(DP) :: ZEROARRAYN(NDIM)
            
            ZEROARRAYN = ZERO

            DO i = 1, SIZE(thislist, 1)
                CALL update(thislist(i), u = ZEROARRAYN)
            END DO

        END SUBROUTINE empty_node_list

END MODULE node_list_module