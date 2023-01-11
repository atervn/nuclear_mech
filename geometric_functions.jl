function get_volume!(shellStruct)
    #=
    based on https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    and http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    =#
    volumes = zeros(length(shellStruct.tri));

    for i = eachindex(shellStruct.tri)
        volumes[i] = 1/6*dot(shellStruct.vert[shellStruct.tri[i][1]],cross(shellStruct.vert[shellStruct.tri[i][2]],shellStruct.vert[shellStruct.tri[i][3]]))
    end

    return sum(volumes)

end

function get_area!(shellStruct)

    areas = zeros(length(shellStruct.tri));

    for i = eachindex(shellStruct.tri)

        dotProduct = dot(shellStruct.edgeVectors[shellStruct.triEdge1[i]],shellStruct.edgeVectors[shellStruct.triEdge2[i]])

        θ = acos(dotProduct/(shellStruct.edgeVectorNorms[shellStruct.triEdge1[i]]*shellStruct.edgeVectorNorms[shellStruct.triEdge2[i]]));
        
        areas[i] = 1/2*shellStruct.edgeVectorNorms[shellStruct.triEdge1[i]]*shellStruct.edgeVectorNorms[shellStruct.triEdge2[i]]*sin(θ);
        
    end

    return areas;

end

function get_local_curvatures!(shellStruct)

    angles1 = zeros(Float64,length(shellStruct.edges));
    angles2 = zeros(Float64,length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)
        

        AB = shellStruct.vert[shellStruct.edges[i][1]] - shellStruct.vert[shellStruct.edges3Vertex[i][1]]   
        AC = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges3Vertex[i][1]]

        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        angles1[i] = acosd(dotProduct/(ABnorm*ACnorm));
        
        AB = shellStruct.vert[shellStruct.edges[i][1]] - shellStruct.vert[shellStruct.edges3Vertex[i][2]]        
        AC = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges3Vertex[i][2]]

        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        angles2[i] = acosd(dotProduct/(ABnorm*ACnorm));
    end

    curvatures = zeros(length(shellStruct.vert));

    # https://computergraphics.stackexchange.com/a/1721

    for i = 1:length(shellStruct.vert)

        cotangentSum = zeros(Float64,3);

        for j = shellStruct.vertexEdges[i]

            cotangentSum .+= (cotd(angles1[j]) + cotd(angles2[j])).*shellStruct.edgeVectors[j]

        end

        LaplaceBeltrami = 1/(2*sum(shellStruct.voronoiAreas[i]))*cotangentSum;

        LaplaceBeltramiNorm = norm(LaplaceBeltrami)

        curvatures[i] = LaplaceBeltramiNorm/2;
    end

    shellStruct.curvatures = curvatures;

end

function get_voronoi_areas!(shellStruct)

    voronoiAreas = Vector{Float64}(undef,length(shellStruct.vert))

    for i = 1:length(shellStruct.vert)

        voronoiAreasTemp = 0;

        for j = 1:length(shellStruct.vertexTri[i])

            voronoiAreasTemp += shellStruct.triangleAreas[shellStruct.vertexTri[i][j]]/3

        end

        voronoiAreas[i] = voronoiAreasTemp;

    end

    shellStruct.voronoiAreas = voronoiAreas;

end

function get_area_unit_vectors!(shellStruct)


    baryocenters = Vector{Vec{3,Float64}}(undef, length(shellStruct.tri));

    for i = eachindex(shellStruct.tri)
        baryocenters[i] = mean(shellStruct.vert[shellStruct.tri[i]]);
    end

    shellStruct.areaUnitVectors = Vector{Vector{Vec{3,Float64}}}(undef, length(shellStruct.vert))

    for i = 1:length(shellStruct.vert)

        areaUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.vertexTri[i]))
        j2 = 1;
        for j = shellStruct.vertexTri[i]

            vectorTemp = shellStruct.vert[i] - baryocenters[j]

            vectorNorm = norm(vectorTemp);

            areaUnitVectors[j2] = vectorTemp/vectorNorm;
            
            j2 += 1
        end

        shellStruct.areaUnitVectors[i] = areaUnitVectors;

    end

end

function get_triangle_normals!(shellStruct)

    #triangleNormalUnitVectors = zeros(Float64,size(shellStruct.tri,1),3);
    triangleNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.tri));

    for i = eachindex(shellStruct.tri)

        # tri = @view shellStruct.tri[i,:];

        normalVector = cross(shellStruct.edgeVectors[shellStruct.triEdge1[i]],shellStruct.edgeVectors[shellStruct.triEdge2[i]])

        #normalVector = cross_product(p1,p2,p3);

        triangleNormalUnitVectors[i] = normalVector/norm(normalVector);

    end

    shellStruct.triangleNormalUnitVectors = triangleNormalUnitVectors;

    edgeNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)

        vector = shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]] + shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]]
        edgeNormalUnitVectors[i] = vector./norm(vector);

    end
    
    shellStruct.edgeNormalUnitVectors = edgeNormalUnitVectors;


    vertexNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.vert));
    
    for i = 1:length(shellStruct.vert)

        vector = sum(shellStruct.triangleNormalUnitVectors[shellStruct.vertexTri[i]]);
        vertexNormalUnitVectors[i] = vector./norm(vector);

    end

    shellStruct.vertexNormalUnitVectors = vertexNormalUnitVectors;

end

function get_triangle_angles(shellStruct)

    angles = zeros(Float64,length(shellStruct.edges));
    
    for i = eachindex(angles)
        angles[i] = acos(min(1.0,dot(shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]],shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]])))
        # println(shellStruct.edgeVectors[i])
        if dot(shellStruct.edgeVectors[i],cross(shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][1]],shellStruct.triangleNormalUnitVectors[shellStruct.edgesTri[i][2]])) > 0
            angles[i] = -angles[i]
        end
    end

    return angles

end

function get_edge_vectors!(shellStruct)

    for i = eachindex(shellStruct.edges)
        if shellStruct.firstEdges[i] == 1
            shellStruct.edgeVectors[i] = shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges[i][1]];
            shellStruct.edgeVectors[shellStruct.mirrorEdges[i]] = -shellStruct.edgeVectors[i]
            shellStruct.edgeVectorNorms[i] = norm(shellStruct.edgeVectors[i])
            shellStruct.edgeVectorNorms[shellStruct.mirrorEdges[i]] = shellStruct.edgeVectorNorms[i]
            shellStruct.edgeUnitVectors[i] = shellStruct.edgeVectors[i]/shellStruct.edgeVectorNorms[i]
            shellStruct.edgeUnitVectors[shellStruct.mirrorEdges[i]] = -shellStruct.edgeUnitVectors[i]
        end
    end

end

function line_point_distance(AB::Vec{3,Float64},AC::Vec{3,Float64})

    return norm(cross(AB,AC)/norm(AB))
 
end

function vertex_triangle_distance(shellStruct, vertex::Vec3, tri::Int)

    # based on https://gist.github.com/joshuashaffer/99d58e4ccbd37ca5d96e

    # if pip === nothing
        tri = shellStruct.tri[tri];
        B  = shellStruct.vert[tri[1]];
        E0 = shellStruct.vert[tri[2]] - B;
        E1 = shellStruct.vert[tri[3]] - B;   
    # else
    #     tri = pip.tri[tri,:];
    #     B  = pip.vert[tri[1]];
    #     E0 = pip.vert[tri[2]] - B;
    #     E1 = pip.vert[tri[3]] - B;
    # end

    D = B - vertex
    
    a = dot(E0,E0);
    b = dot(E0,E1);
    c = dot(E1,E1);
    d = dot(E0,D);
    e = dot(E1,D);;
    f = dot(D,D); 
    
    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e
    
    if (s + t) <= det
        if s < 0
            if t < 0
                # region4
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
                # region 3
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
                # region 5
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
                # region 0
                invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;

                vertices = tri[[1,2,3]]

            end
        end
    else
        if s < 0
            # region 2
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
                # region6
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
                # region 1
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
    
    dist = sqrt(sqrdistance);


    PP0 = B + s*E0 + t*E1
    return dist, PP0, vertices
    
end

function vertex_triangle_distance(vertex::Vec3, tri::Int, pip::pipetteType)

    # based on https://gist.github.com/joshuashaffer/99d58e4ccbd37ca5d96e

    tri = pip.tri[tri];
    B  = pip.vert[tri[1]];
    E0 = pip.vert[tri[2]] - B;
    E1 = pip.vert[tri[3]] - B;

    D = B - vertex
    
    a = dot(E0,E0);
    b = dot(E0,E1);
    c = dot(E1,E1);
    d = dot(E0,D);
    e = dot(E1,D);;
    f = dot(D,D); 
    
    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e
    
    if (s + t) <= det
        if s < 0
            if t < 0
                # region4
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
                # region 3
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
                # region 5
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
                # region 0
                invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;

                vertices = tri[[1,2,3]]

            end
        end
    else
        if s < 0
            # region 2
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
                # region6
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
                # region 1
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
    
    dist = sqrt(sqrdistance);


    PP0 = B + s*E0 + t*E1
    return dist, PP0, vertices
    
end