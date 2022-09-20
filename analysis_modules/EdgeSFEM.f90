MODULE EdgeSFEM
    !Purpose:
    ! Define a Edge based finite element Method that can be used to create the
    ! 
    ! 
    ! This verson of the code deals with the non-linear analysis for a triangular
    ! element and account for Geometric Nonlinear problems.
    ! 
    ! Programmer : A. S. Brar
    ! Date       : 07/01/2021
    
    USE parameter_module,   ONLY: DP, THREE, HALF, TWO_THIRD, ONE_SIXTH, ZERO, &
                                & ONE, ONE_THIRD, MSG_FILE 

    IMPLICIT NONE

    !Defining a point in the domain
    TYPE, PUBLIC :: Point
        REAL(DP) :: x
        REAl(DP) :: y
    END TYPE Point

    !Data structure to store Smoothed Domainomain in the element  
    TYPE, PUBLIC :: FreeEdgeDomain
        PRIVATE
        TYPE(Point) :: P1
        TYPE(Point) :: P2
        TYPE(Point) :: P3
    END TYPE FreeEdgeDomain

    TYPE, EXTENDS(FreeEdgeDomain) :: InsideEdgeDomain
        PRIVATE
        TYPE(Point) :: P4
    END TYPE InsideEdgeDomain


    !Data Structure to store a length and normal of a line segment in 2D space
    TYPE, PUBLIC :: lineseg
        PRIVATE
        REAL(DP) :: l
        REAL(DP) :: nx
        REAL(DP) :: ny
    END TYPE lineseg

    !Data Structure to store the value of shape function in natural coordinate frame
    TYPE, PUBLIC :: freeedge_IP
        PRIVATE
        REAL(DP) :: N1
        REAL(DP) :: N2
        REAL(DP) :: N3
    END TYPE

    TYPE, EXTENDS(freeedge_IP) :: InsideEdge_IP
        PRIVATE
        REAL(DP) :: N4
    END TYPE

    PUBLIC :: calc_SF_Derivative_FE, calc_SF_Derivative_IE

    CONTAINS

        PURE SUBROUTINE findCentroid(Node1, Node2, Node3, Centroid)
            !Purpose:
            !Find the centroid of the triangular Element

            TYPE(Point), INTENT(IN) :: Node1
            TYPE(Point), INTENT(IN) :: Node2
            TYPE(Point), INTENT(IN) :: Node3
            TYPE(Point), INTENT(OUT) :: Centroid

            Centroid%x = (Node1%x + Node2%x + Node3%x)/THREE

            Centroid%y = (Node1%y + Node2%y + Node3%y)/THREE

        END SUBROUTINE findCentroid

        PURE SUBROUTINE calc_Triangle_area(Node1, Node2, Node3, Area)
            !Purpose:
            ! Calculates the area of the triangle 

            TYPE(Point), INTENT(IN) :: Node1
            TYPE(Point), INTENT(IN) :: Node2
            TYPE(Point), INTENT(IN) :: Node3
            REAL(DP), INTENT(OUT) :: Area

            REAL(DP) :: term1
            REAL(DP) :: term2
            REAL(DP) :: term3

            term1 = Node2%x * Node3%y - Node3%x * Node2%y
            term2 = Node1%x * Node3%y - Node3%x * Node1%y
            term3 = Node1%x * Node2%y - Node2%x * Node1%y

            Area =  HALF * ABS(term1 - term2 + term3)  
 
        END SUBROUTINE


        PURE SUBROUTINE free_edge_domain(Node1, Node2, Centroid, Domain)
            !Purpose:
            !Create smoothed domain  for a free Edge
            !defined by Node1 and Node2
            TYPE(Point), INTENT(IN) :: Node1
            TYPE(Point), INTENT(IN) :: Node2
            TYPE(Point), INTENT(IN) :: Centroid
            TYPE(FreeEdgeDomain), INTENT(OUT) :: Domain
            
            Domain%P1 = Node1
            Domain%P2 = Centroid
            Domain%P3 = Node2

        END SUBROUTINE free_edge_domain
        
        PURE SUBROUTINE Inside_edge_domain(Node1, Node2, Centroid1, Centroid2, Domain)
            !Purpose:
            !Create smoothed domain  for an inside Edge
            !defined by Node1 and Node2
            TYPE(Point), INTENT(IN) :: Node1
            TYPE(Point), INTENT(IN) :: Node2
            TYPE(Point), INTENT(IN) :: Centroid1
            TYPE(Point), INTENT(IN) :: Centroid2
            TYPE(InsideEdgeDomain), INTENT(OUT) :: Domain
            
            Domain%P1 = Node1
            Domain%P2 = Centroid1
            Domain%P3 = Node2
            Domain%P4 = Centroid2

        END SUBROUTINE Inside_edge_domain

        PURE SUBROUTINE def_lineseg(linesegID, FE_domain, IE_domain, line)
            !Purpose:
            !Defines the boundary of the smoothed domain and calculate the length
            !of the line segment and normal to the line segment
            !Inputs:
            !   FE_domain: Smoothed domain for a free edge defined by nodes and centroid
            !
            !   IE_domain: Smoothed domain for an inside edge defined by nodes and centroid
            !
            !   linesegID: ID given to line segment
            !           (for a free Edge Domain)
            !           {1 if Node 1 to centroid
            !            2 if Centroid to Node 2
            !            3 if Node 2 to Node 1}
            !           (for an Inside Edge)
            !           {1 if Node 1 to Centroid1
            !            2 if Centroid1 to Node2
            !            3 if Node2 to Centroid2
            !            4 if Centroid2 to Node1}
            ! 
            !Note: Only one of the FE_domain or IE_domain can be present
            ! 
            !Output:
            !   line: lineseg data structure containing information about the normal
            !        and length of the line 

            INTEGER, INTENT(IN) :: linesegID
            TYPE(FreeEdgeDomain), OPTIONAL, INTENT(IN) :: FE_domain
            TYPE(InsideEdgeDomain), OPTIONAL, INTENT(IN) :: IE_domain
            TYPE(lineseg), INTENT(OUT) :: line

            ! Temorary variable to store length of the line segment 
            REAL(DP) :: length

            ! For the free edge domain calculating line segments lenghts and normals
            IF (PRESENT(FE_domain)) THEN
                IF (linesegID == 1) THEN
                    
                    length = SQRT((FE_domain%P2%x - FE_domain%P1%x)**2 + (FE_domain%P2%y - FE_domain%P1%y)**2)
                    line%l = length
                    line%nx = -(1/length) * (FE_domain%P2%y -  FE_domain%P1%y)
                    line%ny = (1/length) * (FE_domain%P2%x -  FE_domain%P1%x)

                ELSE IF (linesegID == 2) THEN

                    length = SQRT((FE_domain%P3%x - FE_domain%P2%x)**2 + (FE_domain%P3%y - FE_domain%P2%y)**2)
                    line%l = length
                    line%nx = -(1/length) * (FE_domain%P3%y -  FE_domain%P2%y)
                    line%ny = (1/length) * (FE_domain%P3%x -  FE_domain%P2%x)

                ELSE IF (linesegID == 3) THEN

                    length = SQRT((FE_domain%P1%x - FE_domain%P3%x)**2 + (FE_domain%P1%y - FE_domain%P3%y)**2)
                    line%l = length
                    line%nx = -(1/length) * (FE_domain%P1%y -  FE_domain%P3%y)
                    line%ny = (1/length) * (FE_domain%P1%x -  FE_domain%P3%x)

                END IF
            
            !For an inside edge domain, line segments normals and length calculation
            ELSE IF (PRESENT(IE_domain)) THEN
                IF (linesegID == 1) THEN
                    
                    length = SQRT((IE_domain%P2%x - IE_domain%P1%x)**2 + (IE_domain%P2%y - IE_domain%P1%y)**2)
                    line%l = length
                    line%nx = (1/length) * (IE_domain%P2%y -  IE_domain%P1%y)
                    line%ny = -(1/length) * (IE_domain%P2%x -  IE_domain%P1%x)

                ELSE IF (linesegID == 2) THEN

                    length = SQRT((IE_domain%P3%x - IE_domain%P2%x)**2 + (IE_domain%P3%y - IE_domain%P2%y)**2)
                    line%l = length
                    line%nx = (1/length) * (IE_domain%P3%y -  IE_domain%P2%y)
                    line%ny = -(1/length) * (IE_domain%P3%x -  IE_domain%P2%x)
                
                ELSE IF (linesegID == 3) THEN

                    length = SQRT((IE_domain%P4%x - IE_domain%P3%x)**2 + (IE_domain%P4%y - IE_domain%P3%y)**2)
                    line%l = length
                    line%nx = (1/length) * (IE_domain%P4%y -  IE_domain%P3%y)
                    line%ny = -(1/length) * (IE_domain%P4%x -  IE_domain%P3%x)

                ELSE IF (linesegID == 4) THEN

                    length = SQRT((IE_domain%P1%x - IE_domain%P4%x)**2 + (IE_domain%P1%y - IE_domain%P4%y)**2)
                    line%l = length
                    line%nx = (1/length) * (IE_domain%P1%y -  IE_domain%P4%y)
                    line%ny = -(1/length) * (IE_domain%P1%x -  IE_domain%P4%x)

                END IF

            END IF 

        END SUBROUTINE def_lineseg

    
        PURE SUBROUTINE def_IP_freeEdge(linesegID, NxG)
            !Purpose:
            ! Calculates the value of shape funtions at integration points
            ! Note: Local triangle should be defined such that Free Edge is 
            !       given by Node1 and Node2 
            !
            !Inputs:
            !   lindesegID : As defined in def_lineseg
        
            !Outputs:
            !   NxG : Values of shape functions N1, N2, N3 at integration points
            !

            INTEGER, INTENT(IN) :: linesegID
            TYPE(freeedge_IP), INTENT(OUT) :: NxG

            ! Value of shape function for Smoothed domain 1
            IF (linesegID == 1) THEN
                ! Shape function value for line segment 1
                NxG%N1 = TWO_THIRD
                NxG%N2 = ONE_SIXTH
                NxG%N3 = ONE_SIXTH

            ELSEIF (linesegID == 2) THEN
                ! Shape function value for line segment 2
                NxG%N1 = ONE_SIXTH
                NxG%N2 = TWO_THIRD
                NxG%N3 = ONE_SIXTH

            ELSEIF (linesegID == 3) THEN
                ! Shape function value for line segment 3
                NxG%N1 = HALF 
                NxG%N2 = HALF
                NxG%N3 = ZERO
            END IF

        END SUBROUTINE def_IP_freeEdge
        
        PURE SUBROUTINE def_IP_insideEdge(linesegID, NxG)
            !Purpose
            ! Calculates the value of shape funtions at integration points
            ! Note: Local triangles should be defined such that Edge is 
            !       given by Node1 and Node3  
            !
            !Inputs:
            !   lindesegID : As defined in def_lineseg
        
            !Outputs:
            !   NxG : Values of shape functions N1, N2, N3 at integration points
            !

            INTEGER, INTENT(IN) :: linesegID
            TYPE(InsideEdge_IP), INTENT(OUT) :: NxG

            IF (linesegID == 1) THEN

                NxG%N1 = TWO_THIRD
                NxG%N2 = ONE_SIXTH
                NxG%N3 = ONE_SIXTH
                NxG%N4 = ZERO

            ELSE IF (linesegID == 2) THEN

                NxG%N1 = ONE_SIXTH
                NxG%N2 = ONE_SIXTH
                NxG%N3 = TWO_THIRD
                NxG%N4 = ZERO
            
            ELSE IF (linesegID == 3) THEN

                NxG%N1 = ONE_SIXTH
                NxG%N2 = ZERO
                NxG%N3 = TWO_THIRD
                NxG%N4 = ONE_SIXTH
                
            ELSE IF (linesegID == 4) THEN

                NxG%N1 = TWO_THIRD
                NxG%N2 = ZERO
                NxG%N3 = ONE_SIXTH
                NxG%N4 = ONE_SIXTH
                
            END IF

        END SUBROUTINE def_IP_insideEdge
        
        PURE SUBROUTINE calc_SF_Derivative_FE(Nodes, dNdX, Vsk)
            !Purpose:
            ! Calculates the shape function deritive for free edge 
            ! 
            ! Inputs: 
            !   Nodes: An array (3,2) containing the information about the coordinates
            !          of nodes of the element
            !  
            ! Outputs: 
            !   dNdX : An array (3,2) containing the value of derivative of the shape function 
            !          for a free edge smoothed domain
            !          


            REAL(DP), DIMENSION(3,2), INTENT(IN) :: Nodes
            REAL(DP), DIMENSION(3,2), INTENT(OUT) :: dNdX
            REAL(DP), INTENT(OUT) :: Vsk
            
            ! Internal variables for the subroutine
            TYPE(Point) :: Centroid
            TYPE(Point) :: Node1
            TYPE(Point) :: Node2
            TYPE(Point) :: Node3
            REAL(DP) :: Area
            INTEGER :: I
            TYPE(lineseg) :: line
            TYPE(freeedge_IP) :: NxG
            TYPE(FreeEdgeDomain) :: Domain

            ! Initializing B matrix
            dNdX(:,:) = ZERO

            ! Aloocating the nodal coordinates to local copy of node points
            Node1%x = Nodes(1,1)
            Node1%y = Nodes(1,2)
            Node2%x = Nodes(2,1)
            Node2%y = Nodes(2,2)
            Node3%x = Nodes(3,1)
            Node3%y = Nodes(3,2)

            !Finding the centroid of the triangle
            CALL findCentroid(Node1, Node2, Node3, Centroid)
            !Finding the area of the triangle
            CALL calc_Triangle_area(Node1, Node2, Node3, Area)
            
            Vsk = Area/THREE

            ! Creating a domain for the current edge
            CALL free_edge_domain(Node1, Node2, Centroid, Domain)

            DO I = 1,3
                CALL def_lineseg(linesegID=I, FE_domain=Domain, line=line)
                CALL def_IP_freeEdge(linesegID=I, NxG=NxG)
                
                dNdX(1,1) = dNdX(1,1) + (ONE/Vsk) * NxG%N1 * line%nx * line%l 
                dNdX(1,2) = dNdX(1,2) + (ONE/Vsk) * NxG%N1 * line%ny * line%l
                
                dNdX(2,1) = dNdX(2,1) + (ONE/Vsk) * NxG%N2 * line%nx * line%l
                dNdX(2,2) = dNdX(2,2) + (ONE/Vsk) * NxG%N2 * line%ny * line%l
 
                dNdX(3,1) = dNdX(3,1) + (ONE/Vsk) * NxG%N3 * line%nx * line%l
                dNdX(3,2) = dNdX(3,2) + (ONE/Vsk) * NxG%N3 * line%ny * line%l

            END DO

        END SUBROUTINE calc_SF_Derivative_FE

        PURE SUBROUTINE calc_SF_Derivative_IE(Nodes, dNdX, Vsk)
            !Purpose:
            ! Calculate shape function derivative for an internal Edge 
            ! 
            ! Inputs: 
            !   Nodes: An array (4,2) containing the information about the coordinates
            !          of nodes of the element
            !  
            ! Outputs: 
            !   dNdX : An array (4,2) containing derivative of shape function for an inner edge
            !          smoothed domain

            REAL(DP), DIMENSION(4,2), INTENT(IN) :: Nodes
            REAL(DP), DIMENSION(4,2), INTENT(OUT) :: dNdX
            REAL(DP), INTENT(OUT) :: Vsk
            
            ! Internal variables for the subroutine
            TYPE(Point) :: Centroid1
            TYPE(Point) :: Centroid2
            TYPE(Point) :: Node1
            TYPE(Point) :: Node2
            TYPE(Point) :: Node3
            TYPE(Point) :: Node4
            REAL(DP) :: Area1
            REAL(DP) :: Area2
            INTEGER :: I
            TYPE(lineseg) :: line
            TYPE(InsideEdge_IP) :: NxG
            TYPE(InsideEdgeDomain) :: Domain

            ! Initializing B matrix
            dNdX(:,:) = ZERO

            Node1%x = Nodes(1,1)
            Node1%y = Nodes(1,2)
            Node2%x = Nodes(2,1)
            Node2%y = Nodes(2,2)
            Node3%x = Nodes(3,1)
            Node3%y = Nodes(3,2)
            Node4%x = Nodes(4,1)
            Node4%y = Nodes(4,2)

            CALL findCentroid(Node1, Node2, Node3, Centroid1)
            CALL calc_Triangle_area(Node1, Node2, Node3, Area1)
            CALL findCentroid(Node1, Node3, Node4, Centroid2)
            CALL calc_Triangle_area(Node1, Node3, Node4, Area2)

            Vsk = ONE_THIRD * (Area1 + Area2)

            CALL Inside_edge_domain(Node1, Node3, Centroid1, Centroid2, Domain)

            DO I = 1,4

                CALL def_lineseg(linesegID=I, IE_domain=Domain, line=line)
                CALL def_IP_insideEdge(linesegID=I, NxG=NxG)

                dNdX(1,1) = dNdX(1,1) + (ONE/Vsk) * NxG%N1 * line%nx * line%l 
                dNdX(1,2) = dNdX(1,2) + (ONE/Vsk) * NxG%N1 * line%ny * line%l
                
                dNdX(2,1) = dNdX(2,1) + (ONE/Vsk) * NxG%N2 * line%nx * line%l
                dNdX(2,2) = dNdX(2,2) + (ONE/Vsk) * NxG%N2 * line%ny * line%l
                
                dNdX(3,1) = dNdX(3,1) + (ONE/Vsk) * NxG%N3 * line%nx * line%l
                dNdX(3,2) = dNdX(3,2) + (ONE/Vsk) * NxG%N3 * line%ny * line%l
                
                dNdX(4,1) = dNdX(4,1) + (ONE/Vsk) * NxG%N4 * line%nx * line%l
                dNdX(4,2) = dNdX(4,2) + (ONE/Vsk) * NxG%N4 * line%ny * line%l
                
            END DO

        END SUBROUTINE calc_SF_Derivative_IE
        
END MODULE EdgeSFEM
