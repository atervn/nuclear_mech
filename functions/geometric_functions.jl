function get_volume!(shellStruct)
    #=
    based on https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    and http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    =#

    # array to store the volumes of each triangle
    volumes = zeros(length(shellStruct.tri));

    for i = eachindex(shellStruct.tri)

        # calculate the volume of each tetrahedron using the dot product and cross product
        volumes[i] = 1/6*dot(shellStruct.vert[shellStruct.tri[i][1]],cross(shellStruct.vert[shellStruct.tri[i][2]],shellStruct.vert[shellStruct.tri[i][3]]))
    end

    return sum(volumes)

end

function get_area!(shellStruct)

    # array to store the areas of each triangle
    areas = zeros(length(shellStruct.tri));

    # iterate through the triangles
    for i = eachindex(shellStruct.tri)

        # dot product of the two edge vectors
        dotProduct = dot(shellStruct.edgeVectors[shellStruct.triEdge1[i]],shellStruct.edgeVectors[shellStruct.triEdge2[i]])

        # angle between the two edge vectors
        temp = dotProduct/(shellStruct.edgeVectorNorms[shellStruct.triEdge1[i]]*shellStruct.edgeVectorNorms[shellStruct.triEdge2[i]])
        if temp > 1
            temp = 1
        elseif temp < -1
            temp = -1
        end

        θ = acos(temp);
        
        # area of the triangle
        areas[i] = 1/2*shellStruct.edgeVectorNorms[shellStruct.triEdge1[i]]*shellStruct.edgeVectorNorms[shellStruct.triEdge2[i]]*sin(θ);

    end

    return areas;

end

function get_local_curvatures!(shellStruct)

    # arrays to store angles
    angles1 = zeros(Float64,length(shellStruct.edges));
    angles2 = zeros(Float64,length(shellStruct.edges));

    # iterate through the edges
    for i = eachindex(shellStruct.edges)
        
        # define vectors between the vectors of one of the neihboring triangles edges
        AB = shellStruct.vert[shellStruct.edges[i][1]] - shellStruct.vert[shellStruct.edges3Vertex[i][1]]   
        AC = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges3Vertex[i][1]]

        # calculate the vector norms and calculate dot product
        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        # calculate the angle
        angles1[i] = acosd(dotProduct/(ABnorm*ACnorm));
        
        # define vectors between the vectors of one of the neihboring triangles edges
        AB = shellStruct.vert[shellStruct.edges[i][1]] - shellStruct.vert[shellStruct.edges3Vertex[i][2]]        
        AC = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges3Vertex[i][2]]

        # calculate the vector norms and calculate dot product
        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        # calculate the angle
        angles2[i] = acosd(dotProduct/(ABnorm*ACnorm));

    end

    # initialize vector for the curvatures
    shellStruct.curvatures = zeros(length(shellStruct.vert));

    # https://computergraphics.stackexchange.com/a/1721

    # iterate through the vertices
    for i = 1:length(shellStruct.vert)

        # init a vector for the cotangent sums
        cotangentSum = zeros(Float64,3);

        # iterate through the neibhoring edges
        for j = shellStruct.vertexEdges[i]

            # Compute the sum of cotangent weights multiplied by edge vectors to approximate the Laplace-Beltrami operator.
            cotangentSum .+= (cotd(angles1[j]) + cotd(angles2[j])).*shellStruct.edgeVectors[j]
            
        end

        # Compute the curvatures using the Laplace-Beltrami operator
        LaplaceBeltrami = 1/(2*sum(shellStruct.voronoiAreas[i]))*cotangentSum;
        LaplaceBeltramiNorm = norm(LaplaceBeltrami)

        # Local curvature for the vertex
        shellStruct.curvatures[i] = LaplaceBeltramiNorm/2;

    end

end

function get_voronoi_areas!(shellStruct)

    # create an array to store the Voronoi areas for each vertex
    shellStruct.voronoiAreas = Vector{Float64}(undef,length(shellStruct.vert))

    for i = 1:length(shellStruct.vert)

        # initialize a temporary variable to store the sum of triangle areas
        voronoiAreasTemp = 0;

        for j = 1:length(shellStruct.vertexTri[i])

            # add one-third of the triangle area to the temporary variable
            voronoiAreasTemp += shellStruct.triangleAreas[shellStruct.vertexTri[i][j]]/3

        end

        # assign the calculated Voronoi area to the corresponding vertex
        shellStruct.voronoiAreas[i] = voronoiAreasTemp;
        
    end
end

function get_area_unit_vectors!(shellStruct)

    # calculate barycenters for each triangle
    baryocenters = Vector{Vec{3,Float64}}(undef, length(shellStruct.tri));
    for i = eachindex(shellStruct.tri)
        
        # calculate the barycenter as the mean of the vertex positions of the triangle
        baryocenters[i] = mean(shellStruct.vert[shellStruct.tri[i]]);
    
    end

    # create an array to store area unit vectors for each vertex
    shellStruct.areaUnitVectors = Vector{Vector{Vec{3,Float64}}}(undef, length(shellStruct.vert))


    for i = 1:length(shellStruct.vert)

        # create an array to store area unit vectors for the current vertex
        areaUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.vertexTri[i]))
        j2 = 1;
        for j = shellStruct.vertexTri[i]

            # calculate the vector from the current vertex to the barycenter of the adjacent triangle
            vectorTemp = shellStruct.vert[i] - baryocenters[j]

            # calculate the norm of the vector
            vectorNorm = norm(vectorTemp);

            # calculate the area unit vector by dividing the vector by its norm
            areaUnitVectors[j2] = vectorTemp/vectorNorm;
            
            # increase the index for the vertices triagnles
            j2 += 1
        end

        # assign the calculated area unit vectors to the corresponding vertex
        shellStruct.areaUnitVectors[i] = areaUnitVectors;

    end
end

function get_shell_normals!(shellStruct)

    # calculate triangle normal unit vectors
    shellStruct.triangleNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.tri));
    for i = eachindex(shellStruct.tri)
        
        # calculate the normal vector for the triangle using the cross product of its edge vectors
        normalVector = cross(shellStruct.edgeVectors[shellStruct.triEdge1[i]],shellStruct.edgeVectors[shellStruct.triEdge2[i]])
        
        # normalize the normal vector to obtain the triangle's normal unit vector
        shellStruct.triangleNormalUnitVectors[i] = normalVector/norm(normalVector);

    end

    # calculate edge normal unit vectors
    shellStruct.edgeNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)

        # calculate the sum of the triangle normal unit vectors adjacent to the edge
        vector = shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]] + shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]]
        
        # normalize the sum vector to obtain the edge's normal unit vector
        shellStruct.edgeNormalUnitVectors[i] = vector./norm(vector);

    end

    # calculate vertex normal unit vectors
    shellStruct.vertexNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.vert));
    
    for i = 1:length(shellStruct.vert)

        # sum the triangle normal unit vectors associated with the vertex
        vector = sum(shellStruct.triangleNormalUnitVectors[shellStruct.vertexTri[i]]);

        # normalize the sum vector to obtain the vertex's normal unit vector
        shellStruct.vertexNormalUnitVectors[i] = vector./norm(vector);

    end

end

function get_triangle_angles(shellStruct)

    # init a angles vector
    angles = zeros(Float64,length(shellStruct.edges));
    
    # iterate through the angles
    for i = eachindex(angles)

        # calculate the angle between the triangle edges using the dot product of their unit normal vectors
        angles[i] = acos(max(-1.0,min(1.0,dot(shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]],shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]]))))
        
        # check the sign of the angle based on the cross product of the edge vector and the unit normal vectors
        if dot(-shellStruct.edgeVectors[i],cross(shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]],shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]])) > 0
            angles[i] = -angles[i]
        end
        
    end

    return angles

end

function get_edge_vectors!(shellStruct)

    # iterate through the shell edges
    for i = eachindex(shellStruct.edges)

        # check if the edge is the first edge of a vertex pair
        if shellStruct.firstEdges[i] == 1
            
            # calculate the edge vector as the difference between the second and first vertex of the edge
            shellStruct.edgeVectors[i] = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges[i][1]];
            
            # set the mirror edge's vector as the negation of the original edge vector
            shellStruct.edgeVectors[shellStruct.mirrorEdges[i]] = -shellStruct.edgeVectors[i]
            
            # calculate the norm of the original edge vector
            shellStruct.edgeVectorNorms[i] = norm(shellStruct.edgeVectors[i])
            
            # set the norm of the mirror edge as the same as the original edge's norm
            shellStruct.edgeVectorNorms[shellStruct.mirrorEdges[i]] = shellStruct.edgeVectorNorms[i]
            
            # calculate the unit vector of the original edge by dividing the edge vector by its norm
            shellStruct.edgeUnitVectors[i] = shellStruct.edgeVectors[i]/shellStruct.edgeVectorNorms[i]
            
            # set the unit vector of the mirror edge as the negation of the original edge's unit vector
            shellStruct.edgeUnitVectors[shellStruct.mirrorEdges[i]] = -shellStruct.edgeUnitVectors[i]
        end
    end
end

function line_point_distance(AB::Vec{3,Float64},AC::Vec{3,Float64})

    # calculate the cross product of AB and AC
    crossProduct = cross(AB, AC)
    
    # calculate the magnitude of AB
    ABMagnitude = norm(AB)
    
    # calculate the normalized cross product divided by the magnitude of AB
    normalizedCrossProduct = crossProduct / ABMagnitude
    
    # calculate the distance between the line segment and the point
    distance = norm(normalizedCrossProduct)
    
    # return the distance
    return distance

end

function vertex_triangle_distance(shell, vertex::Vec3, tri::Int)

    # get the triangle vertices from the pipetteType object
    tri = shell.tri[tri];

    # extract the base vertex and edge vectors of the triangle
    B  = shell.vert[tri[1]];
    E0 = shell.vert[tri[2]] - B;
    E1 = shell.vert[tri[3]] - B;   

    # calculate the difference vector between B and the vertex of the triangle
    D = B - vertex

    # call the vertex_triangle_distance_sub function to perform the distance calculation
    dist, PP0, vertices = vertex_triangle_distance_sub(B,E0,E1,D,tri)
    
    # return the distance, closest point on the triangle, and triangle vertices
    return dist, PP0, vertices
    
end

function vertex_triangle_distance_sub(B,E0,E1,D,tri)

    # calculate dot products of edge vectors and differences
    a = dot(E0,E0);
    b = dot(E0,E1);
    c = dot(E1,E1);
    d = dot(E0,D);
    e = dot(E1,D);;
    f = dot(D,D); 
    
    # calculate the determinant and barycentric coordinates
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e
    
    if (s + t) <= det
        if s < 0
            if t < 0
                # region4: minimum on the vertex region of the triangle
                if d < 0
                    t = 0
                    if -d >= a
                        s = 1;
                        sqrdistance = a + 2*d + f;
                    else
                        s = -d / a;
                        sqrdistance = d*s + f;
                    end
                else
                    s = 0;
                    if e >= 0;
                        t = 0;
                        sqrdistance = f;
                    else
                        if -e >= c
                            t = 1;
                            sqrdistance = c + 2*e + f;
                        else
                            t = -e/c;
                            sqrdistance = e*t + f;
                        end
                    end
                end

                vertices = [tri[1]]

            else
                # region 3: minimum on the edge region E1
                s = 0;
                if e >= 0
                    t = 0;
                    sqrdistance = f;
                else
                    if -e >= c
                        t = 1;
                        sqrdistance = c + 2*e + f;
                    else
                        t = -e/c;
                        sqrdistance = e*t + f;
                    end
                end

                vertices = tri[[1,3]]

            end
        else
            if t < 0
                # region 5: minimum on the edge region E0
                t = 0
                if d >= 0
                    s = 0;
                    sqrdistance = f;
                else
                    if -d >= a
                        s = 1;
                        sqrdistance = a + 2*d + f;
                    else
                        s = -d/a;
                        sqrdistance = d*s + f;
                    end
                end

                vertices = tri[[1,2]]

            else
                # region 0: minimum inside the triangle
                invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;

                vertices = tri[[1,2,3]]

            end
        end
    else
        if s < 0
            # region 2: minimum on the edge region E1+E2
            tmp0 = b + d;
            tmp1 = c + e;
            if tmp1 > tmp0  # minimum on edge s+t=1
                numer = tmp1 - tmp0;
                denom = a - 2*b + c;
                if numer >= denom
                    s = 1;
                    t = 0;
                    sqrdistance = a + 2*d + f;
                else
                    s = numer/denom;
                    t = 1 - s;
                    sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                end
    
            else  # minimum on edge s=0
                s = 0.0
                if tmp1 <= 0
                    t = 1;
                    sqrdistance = c + 2*e + f;
                else
                    if e >= 0
                        t = 0;
                        sqrdistance = f;
                    else
                        t = -e/c;
                        sqrdistance = e*t + f;
                    end
                end
            end
            vertices = [tri[3]]
        else
            if t < 0
                # region 6: minimum on the edge region E0+E2
                tmp0 = b + e;
                tmp1 = a + d;
                if tmp1 > tmp0
                    numer = tmp1 - tmp0;
                    denom = a - 2*b + c;
                    if numer >= denom
                        t = 1;
                        s = 0;
                        sqrdistance = c + 2*e + f;
                    else
                        t = numer / denom;
                        s = 1 - t;
                        sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    end
                else
                    t = 0;
                    if tmp1 <= 0
                        s = 1;
                        sqrdistance = a + 2*d + f;
                    else
                        if d >= 0
                            s = 0;
                            sqrdistance = f;
                        else
                            s = -d/a;
                            sqrdistance = d*s + f;
                        end
                    end
                end

                vertices = [tri[2]]
            
            else
                # region 1: minimum on the edge region E0+E1
                numer = c + e - b - d;
                if numer <= 0
                    s = 0
                    t = 1
                    sqrdistance = c + 2*e + f;
                else
                    denom = a - 2*b + c;
                    if numer >= denom
                        s = 1;
                        t = 0;
                        sqrdistance = a + 2*d + f;
                    else
                        s = numer/denom;
                        t = 1 - s;
                        sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    end
                end

                vertices = tri[[2,3]]

            end
        end
    end
    
    # account for numerical round-off error
    if sqrdistance < 0
        sqrdistance = 0
    end
    
    # calculate the final distance
    dist = sqrt(sqrdistance);

    # calculate the closest point on the triangle to the given point B
    PP0 = B + s*E0 + t*E1

    return dist, PP0, vertices

end