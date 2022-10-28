function get_volume!(nuc)
    #=
    based on https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    and http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    =#
    volumes = zeros(length(nuc.tri));

    for i = eachindex(nuc.tri)
        volumes[i] = 1/6*dot(nuc.vert[nuc.tri[i][1]],cross(nuc.vert[nuc.tri[i][2]],nuc.vert[nuc.tri[i][3]]))
    end

    return sum(volumes)

end

function get_area!(nuc)

    areas = zeros(length(nuc.tri));

    for i = eachindex(nuc.tri)

        dotProduct = dot(nuc.edgeVectors[nuc.triEdge1[i]],nuc.edgeVectors[nuc.triEdge2[i]])

        θ = acos(dotProduct/(nuc.edgeVectorNorms[nuc.triEdge1[i]]*nuc.edgeVectorNorms[nuc.triEdge2[i]]));
        
        areas[i] = 1/2*nuc.edgeVectorNorms[nuc.triEdge1[i]]*nuc.edgeVectorNorms[nuc.triEdge2[i]]*sin(θ);
        
    end

    return areas;

end

function get_local_curvatures!(nuc)

    angles1 = zeros(Float64,length(nuc.edges));
    angles2 = zeros(Float64,length(nuc.edges));

    for i = eachindex(nuc.edges)
        

        AB = nuc.vert[nuc.edges[i][1]] - nuc.vert[nuc.edges3Vertex[i][1]]   
        AC = nuc.vert[nuc.edges[i][2]] - nuc.vert[nuc.edges3Vertex[i][1]]

        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        angles1[i] = acosd(dotProduct/(ABnorm*ACnorm));
        
        AB = nuc.vert[nuc.edges[i][1]] - nuc.vert[nuc.edges3Vertex[i][2]]        
        AC = nuc.vert[nuc.edges[i][2]] - nuc.vert[nuc.edges3Vertex[i][2]]

        ABnorm = norm(AB)
        ACnorm = norm(AC)
        dotProduct = dot(AB,AC)

        angles2[i] = acosd(dotProduct/(ABnorm*ACnorm));
    end

    curvatures = zeros(length(nuc.vert));

    # https://computergraphics.stackexchange.com/a/1721

    for i = 1:length(nuc.vert)

        cotangentSum = zeros(Float64,3);

        for j = nuc.vertexEdges[i]

            cotangentSum .+= (cotd(angles1[j]) + cotd(angles2[j])).*nuc.edgeVectors[j]

        end

        LaplaceBeltrami = 1/(2*sum(nuc.voronoiAreas[i]))*cotangentSum;

        LaplaceBeltramiNorm = norm(LaplaceBeltrami)

        curvatures[i] = LaplaceBeltramiNorm/2;
    end

    nuc.curvatures = curvatures;

end

function get_voronoi_areas!(nuc)

    voronoiAreas = Vector{Float64}(undef,length(nuc.vert))

    for i = 1:length(nuc.vert)

        voronoiAreasTemp = 0;

        for j = 1:length(nuc.vertexTri[i])

            voronoiAreasTemp += nuc.triangleAreas[nuc.vertexTri[i][j]]/3

        end

        voronoiAreas[i] = voronoiAreasTemp;

    end

    nuc.voronoiAreas = voronoiAreas;

end

function get_area_unit_vectors!(nuc)


    baryocenters = Vector{Vec{3,Float64}}(undef, length(nuc.tri));

    for i = eachindex(nuc.tri)
        baryocenters[i] = mean(nuc.vert[nuc.tri[i]]);
    end

    nuc.areaUnitVectors = Vector{Vector{Vec{3,Float64}}}(undef, length(nuc.vert))

    for i = 1:length(nuc.vert)

        areaUnitVectors = Vector{Vec{3,Float64}}(undef, length(nuc.vertexTri[i]))
        j2 = 1;
        for j = nuc.vertexTri[i]

            vectorTemp = nuc.vert[i] - baryocenters[j]

            vectorNorm = norm(vectorTemp);

            areaUnitVectors[j2] = vectorTemp/vectorNorm;
            
            j2 += 1
        end

        nuc.areaUnitVectors[i] = areaUnitVectors;

    end

end

function get_triangle_normals!(nuc)

    #triangleNormalUnitVectors = zeros(Float64,size(nuc.tri,1),3);
    triangleNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(nuc.tri));

    for i = eachindex(nuc.tri)

        # tri = @view nuc.tri[i,:];

        normalVector = cross(nuc.edgeVectors[nuc.triEdge1[i]],nuc.edgeVectors[nuc.triEdge2[i]])

        #normalVector = cross_product(p1,p2,p3);

        triangleNormalUnitVectors[i] = normalVector/norm(normalVector);

    end

    nuc.triangleNormalUnitVectors = triangleNormalUnitVectors;

    edgeNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(nuc.edges));

    for i = eachindex(nuc.edges)

        vector = nuc.triangleNormalUnitVectors[nuc.edgesTri[i][1]] + nuc.triangleNormalUnitVectors[nuc.edgesTri[i][2]]
        edgeNormalUnitVectors[i] = vector./norm(vector);

    end
    
    nuc.edgeNormalUnitVectors = edgeNormalUnitVectors;


    vertexNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    
    for i = 1:length(nuc.vert)

        vector = sum(nuc.triangleNormalUnitVectors[nuc.vertexTri[i]]);
        vertexNormalUnitVectors[i] = vector./norm(vector);

    end

    nuc.vertexNormalUnitVectors = vertexNormalUnitVectors;

end

function get_triangle_angles(nuc)

    angles = zeros(Float64,length(nuc.edges));
    
    for i = eachindex(angles)
        angles[i] = acos(min(1.0,dot(nuc.triangleNormalUnitVectors[nuc.edgesTri[i][1]],nuc.triangleNormalUnitVectors[nuc.edgesTri[i][2]])))
    end

    return angles

end

function get_edge_vectors!(nuc)

    for i = eachindex(nuc.edges)
        if nuc.firstEdges[i] == 1
            nuc.edgeVectors[i] = nuc.vert[nuc.edges[i][2]] - nuc.vert[nuc.edges[i][1]];
            nuc.edgeVectors[nuc.mirrorEdges[i]] = -nuc.edgeVectors[i]
            nuc.edgeVectorNorms[i] = norm(nuc.edgeVectors[i])
            nuc.edgeVectorNorms[nuc.mirrorEdges[i]] = nuc.edgeVectorNorms[i]
            nuc.edgeUnitVectors[i] = nuc.edgeVectors[i]/nuc.edgeVectorNorms[i]
            nuc.edgeUnitVectors[nuc.mirrorEdges[i]] = -nuc.edgeUnitVectors[i]
        end
    end

end

function line_point_distance(AB::Vec{3,Float64},AC::Vec{3,Float64})

    return norm(cross(AB,AC)/norm(AB))
 
end

function vertex_triangle_distance(nuc, vertex, tri, pip = nothing)

    # based on https://gist.github.com/joshuashaffer/99d58e4ccbd37ca5d96e

    # if pip === nothing
        tri = nuc.tri[tri];
        B  = nuc.vert[tri[1]];
        E0 = nuc.vert[tri[2]] - B;
        E1 = nuc.vert[tri[3]] - B;   
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
