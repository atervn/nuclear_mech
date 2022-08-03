function get_volume!(nuc)
    #=
    based on https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    and http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    =#
    volumes = zeros(size(nuc.tri,1),1);

    for i = 1:size(nuc.tri,1)
        t = nuc.tri[i,:];

        v321 = nuc.x[t[3]]*nuc.y[t[2]]*nuc.z[t[1]];
        v231 = nuc.x[t[2]]*nuc.y[t[3]]*nuc.z[t[1]];
        v312 = nuc.x[t[3]]*nuc.y[t[1]]*nuc.z[t[2]];
        v132 = nuc.x[t[1]]*nuc.y[t[3]]*nuc.z[t[2]];
        v213 = nuc.x[t[2]]*nuc.y[t[1]]*nuc.z[t[3]];
        v123 = nuc.x[t[1]]*nuc.y[t[2]]*nuc.z[t[3]];
        
        volumes[i] = 1/6*(-v321 + v231 + v312 - v132 - v213 + v123);

    end

    return sum(volumes)

end

function get_area!(nuc)

    areas = zeros(size(nuc.tri,1),1);

    for i = 1:1:size(nuc.tri,1)
        ABx = nuc.x[nuc.tri[i,1]] - nuc.x[nuc.tri[i,2]];
        ABy = nuc.y[nuc.tri[i,1]] - nuc.y[nuc.tri[i,2]];
        ABz = nuc.z[nuc.tri[i,1]] - nuc.z[nuc.tri[i,2]];
        
        ACx = nuc.x[nuc.tri[i,1]] - nuc.x[nuc.tri[i,3]];
        ACy = nuc.y[nuc.tri[i,1]] - nuc.y[nuc.tri[i,3]];
        ACz = nuc.z[nuc.tri[i,1]] - nuc.z[nuc.tri[i,3]];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        θ = acosd(dotProduct/(ABnorm*ACnorm));
        
        areas[i] = 1/2*ABnorm*ACnorm*sind(θ);
        
    end

    return sum(areas);

end

function get_local_curvatures!(nuc)

    angles1 = zeros(Float64,size(nuc.edges,1));
    angles2 = zeros(Float64,size(nuc.edges,1));

    nuc.edges3vertex = zeros(Int64,size(nuc.edges,1),2);

    for i = 1:size(nuc.edges,1)
        
        neighboringTriangles = findall(sum(nuc.tri .== nuc.edges[i,1],dims=2) .> 0 .&& sum(nuc.tri .== nuc.edges[i,2],dims=2) .> 0);
        neighboringTriangles = [j[1] for j in neighboringTriangles];

        thirdVertex1 = nuc.tri[neighboringTriangles[1],.!(nuc.tri[neighboringTriangles[1],:] .== nuc.edges[i,1] .|| nuc.tri[neighboringTriangles[1],:] .== nuc.edges[i,2])][1];
        thirdVertex2 = nuc.tri[neighboringTriangles[2],.!(nuc.tri[neighboringTriangles[2],:] .== nuc.edges[i,1] .|| nuc.tri[neighboringTriangles[2],:] .== nuc.edges[i,2])][1];

        nuc.edges3vertex[i,1] = thirdVertex1;
        nuc.edges3vertex[i,2] = thirdVertex2; 

        ABx = nuc.x[nuc.edges[i,1]] - nuc.x[thirdVertex1];
        ABy = nuc.y[nuc.edges[i,1]] - nuc.y[thirdVertex1];
        ABz = nuc.z[nuc.edges[i,1]] - nuc.z[thirdVertex1];
        
        ACx = nuc.x[nuc.edges[i,2]] - nuc.x[thirdVertex1];
        ACy = nuc.y[nuc.edges[i,2]] - nuc.y[thirdVertex1];
        ACz = nuc.z[nuc.edges[i,2]] - nuc.z[thirdVertex1];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        angles1[i] = acosd(dotProduct/(ABnorm*ACnorm));
        
        ABx = nuc.x[nuc.edges[i,1]] - nuc.x[thirdVertex2];
        ABy = nuc.y[nuc.edges[i,1]] - nuc.y[thirdVertex2];
        ABz = nuc.z[nuc.edges[i,1]] - nuc.z[thirdVertex2];
        
        ACx = nuc.x[nuc.edges[i,2]] - nuc.x[thirdVertex2];
        ACy = nuc.y[nuc.edges[i,2]] - nuc.y[thirdVertex2];
        ACz = nuc.z[nuc.edges[i,2]] - nuc.z[thirdVertex2];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        angles2[i] = acosd(dotProduct/(ABnorm*ACnorm));
    end

    curvatures = zeros(length(nuc.x));
    vertexNormalUnitVectors = zeros(length(nuc.x),3);

    for i = 1:length(nuc.x)

        cotangentSum = [0;0;0];

        for j = nuc.vertexEdges[i]

            cotangentSum = cotangentSum .+ (cotd(angles1[j]) + cotd(angles2[j])).*[nuc.x[nuc.edges[j,1]] - nuc.x[nuc.edges[j,2]] ; nuc.y[nuc.edges[j,1]] - nuc.y[nuc.edges[j,2]] ; nuc.z[nuc.edges[j,1]] - nuc.z[nuc.edges[j,2]]]

        end

        LaplaceBeltrami = 1/(2*sum(nuc.voronoiAreas[i])).*cotangentSum;

        LaplaceBeltramiNorm = sqrt(sum(LaplaceBeltrami.^2));

        curvatures[i] = LaplaceBeltramiNorm./2;
        vertexNormalUnitVectors[i,:] = LaplaceBeltrami./LaplaceBeltramiNorm;

    end

    nuc.curvatures = curvatures;
    nuc.vertexNormalUnitVectors = vertexNormalUnitVectors;

end

function get_voronoi_areas!(nuc)

    voronoiAreas = fill(Float64[], length(nuc.x));

    for i = 1:length(nuc.x)

        voronoiAreasTemp = zeros(Float64,length(nuc.vertexTri[i,1]));

        for j = 1:length(nuc.vertexTri[i,1])

            v = circshift(nuc.tri[nuc.vertexTri[i,1][j],:],(-nuc.vertexTri[i,2][j]));
            centerX3 = sum(nuc.x[v])/3;
            centerY3 = sum(nuc.y[v])/3;
            centerZ3 = sum(nuc.z[v])/3;

            tempArea = 0;

            for k = 2:3

                centerX2 = (nuc.x[i] + nuc.x[v[k]])/2;
                centerY2 = (nuc.y[i] + nuc.y[v[k]])/2;
                centerZ2 = (nuc.z[i] + nuc.z[v[k]])/2;

                dist1 = sqrt((centerX3 - centerX2).^2 + (centerY3 - centerY2).^2 + (centerZ3 - centerZ2).^2);
                dist2 = sqrt((nuc.x[i] - centerX2).^2 + (nuc.y[i] - centerY2).^2 + (nuc.z[i] - centerZ2).^2);

                tempArea += dist1*dist2/2;

            end

            voronoiAreasTemp[j] = tempArea;

        end

        voronoiAreas[i] = voronoiAreasTemp;

    end

    nuc.voronoiAreas = voronoiAreas;

end

function get_area_unit_vectors!(nuc)


    baryocenters = zeros(Float64,size(nuc.tri,1,),3);

    for i = 1:size(nuc.tri,1)
        baryocenters[i,:] = [sum(nuc.x[nuc.tri[i,:]]) sum(nuc.y[nuc.tri[i,:]]) sum(nuc.z[nuc.tri[i,:]])]./3;
    end

    nuc.areaUnitVectors = fill(Float64[], length(nuc.x));

    for i = 1:length(nuc.x)

        areaUnitVectors = zeros(Float64,length(nuc.vertexTri[i,1]),3);

        j2 = 1;
        for j = nuc.vertexTri[i,1]

            vectorTemp = [nuc.x[i] - baryocenters[j,1] nuc.y[i] - baryocenters[j,2] nuc.z[i] - baryocenters[j,3]];

            vectorNorm = sqrt(sum(vectorTemp.^2));

            areaUnitVectors[j2,:] = vectorTemp./vectorNorm;
            
            j2 += 1
        end

        nuc.areaUnitVectors[i] = areaUnitVectors;

    end

end

function get_triangle_normals!(nuc)

    triangleNormalUnitVectors = zeros(Float64,size(nuc.tri,1),3);

    for i = 1:size(nuc.tri,1)

        tri = nuc.tri[i,:];
        p1 = [nuc.x[tri[1]], nuc.y[tri[1]], nuc.z[tri[1]]];
        p2 = [nuc.x[tri[2]], nuc.y[tri[2]], nuc.z[tri[2]]];
        p3 = [nuc.x[tri[3]], nuc.y[tri[3]], nuc.z[tri[3]]];

        normalVector = cross_product(p1,p2,p3);

        triangleNormalUnitVectors[i,:] = normalVector./sqrt(sum(normalVector.^2));

    end

    nuc.triangleNormalUnitVectors = triangleNormalUnitVectors;

end

function cross_product(p1,p2,p3)

    return [((p2[2] - p1[2]) * (p3[3] - p1[3])) - ((p2[3] - p1[3]) * (p3[2] - p1[2]))
            ((p2[3] - p1[3]) * (p3[1] - p1[1])) - ((p2[1] - p1[1]) * (p3[3] - p1[3]))
            ((p2[1] - p1[1]) * (p3[2] - p1[2])) - ((p2[2] - p1[2]) * (p3[1] - p1[1]))];

end

function get_triangle_angles(nuc)

    angles = zeros(Float64,size(nuc.edges,1),1);
    

    for i = 1:size(nuc.edges,1)
        normalVector1 = nuc.triangleNormalUnitVectors[nuc.edgesTri[i,1],:];
        normalVector2 = nuc.triangleNormalUnitVectors[nuc.edgesTri[i,2],:];

        dotProduct = normalVector1[1]*normalVector2[1] + normalVector1[2]*normalVector2[2] + normalVector1[3]*normalVector2[3];
        angles[i] = acosd(dotProduct);
    end

    return angles

end

function line_point_distance(p1,p2,p3)

    return sqrt(sum(cross_product(p1,p2,p3).^2))/sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2 + (p2[3] - p1[3])^2)

end

function vertex_triangle_distance(nuc, p, tri)

    tri = nuc.tri[tri,:];

    B = [nuc.x[tri[1]], nuc.y[tri[1]], nuc.z[tri[1]]];
    E0 = [nuc.x[tri[2]], nuc.y[tri[2]], nuc.z[tri[2]]] .- B;
    E1 = [nuc.x[tri[3]], nuc.y[tri[3]], nuc.z[tri[3]]] .- B;
    D = B .- [nuc.x[p], nuc.y[p], nuc.z[p]];
    
    
    a = E0[1]^2 + E0[2]^2 + E0[3]^2;
    b = E0[1]*E1[1] + E0[2]*E1[2] + E0[3]*E1[3];
    c = E1[1]^2 + E1[2]^2 + E1[3]^2;
    d = E0[1]*D[1] + E0[2]*D[2] + E0[3]*D[3];
    e = E1[1]*D[1] + E1[2]*D[2] + E1[3]*D[3];
    f = D[1]^2 + D[2]^2 + D[3]^2;
    
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
            else
                # region 0
                invDet = 1/det;
                s = s*invDet;
                t = t*invDet;
                sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
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
            end
        end
    end
    
    # account for numerical round-off error
    if sqrdistance < 0
        sqrdistance = 0
    end
    
    dist = sqrt(sqrdistance);
    
    PP0 = Float64.(B + s*E0 + t*E1)
    return dist, PP0
    
end