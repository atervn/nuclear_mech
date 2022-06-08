function get_volume!(nucleus)
    #=
    based on https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
    and http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    =#
    volumes = zeros(size(nucleus.tri,1),1);

    for i = 1:size(nucleus.tri,1)
        t = nucleus.tri[i,:];

        v321 = nucleus.x[t[3]]*nucleus.y[t[2]]*nucleus.z[t[1]];
        v231 = nucleus.x[t[2]]*nucleus.y[t[3]]*nucleus.z[t[1]];
        v312 = nucleus.x[t[3]]*nucleus.y[t[1]]*nucleus.z[t[2]];
        v132 = nucleus.x[t[1]]*nucleus.y[t[3]]*nucleus.z[t[2]];
        v213 = nucleus.x[t[2]]*nucleus.y[t[1]]*nucleus.z[t[3]];
        v123 = nucleus.x[t[1]]*nucleus.y[t[2]]*nucleus.z[t[3]];
        
        volumes[i] = 1/6*(-v321 + v231 + v312 - v132 - v213 + v123);

    end

    return sum(volumes)

end

function get_area!(nucleus)

    areas = zeros(size(nucleus.tri,1),1);

    for i = 1:1:size(nucleus.tri,1)
        ABx = nucleus.x[nucleus.tri[i,1]] - nucleus.x[nucleus.tri[i,2]];
        ABy = nucleus.y[nucleus.tri[i,1]] - nucleus.y[nucleus.tri[i,2]];
        ABz = nucleus.z[nucleus.tri[i,1]] - nucleus.z[nucleus.tri[i,2]];
        
        ACx = nucleus.x[nucleus.tri[i,1]] - nucleus.x[nucleus.tri[i,3]];
        ACy = nucleus.y[nucleus.tri[i,1]] - nucleus.y[nucleus.tri[i,3]];
        ACz = nucleus.z[nucleus.tri[i,1]] - nucleus.z[nucleus.tri[i,3]];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        θ = acosd(dotProduct/(ABnorm*ACnorm));
        
        areas[i] = 1/2*ABnorm*ACnorm*sind(θ);
        
    end

    return sum(areas);

end

function get_local_curvatures!(nucleus)

    angles1 = zeros(Float64,size(nucleus.edges,1));
    angles2 = zeros(Float64,size(nucleus.edges,1));

    for i = 1:size(nucleus.edges,1)
        
        neighboringTriangles = findall(sum(nucleus.tri .== nucleus.edges[i,1],dims=2) .> 0 .&& sum(nucleus.tri .== nucleus.edges[i,2],dims=2) .> 0);
        neighboringTriangles = [j[1] for j in neighboringTriangles];

        thirdVertex1 = nucleus.tri[neighboringTriangles[1],.!(nucleus.tri[neighboringTriangles[1],:] .== nucleus.edges[i,1] .|| nucleus.tri[neighboringTriangles[1],:] .== nucleus.edges[i,2])][1];
        thirdVertex2 = nucleus.tri[neighboringTriangles[2],.!(nucleus.tri[neighboringTriangles[2],:] .== nucleus.edges[i,1] .|| nucleus.tri[neighboringTriangles[2],:] .== nucleus.edges[i,2])][1];

        ABx = nucleus.x[nucleus.edges[i,1]] - nucleus.x[thirdVertex1];
        ABy = nucleus.y[nucleus.edges[i,1]] - nucleus.y[thirdVertex1];
        ABz = nucleus.z[nucleus.edges[i,1]] - nucleus.z[thirdVertex1];
        
        ACx = nucleus.x[nucleus.edges[i,2]] - nucleus.x[thirdVertex1];
        ACy = nucleus.y[nucleus.edges[i,2]] - nucleus.y[thirdVertex1];
        ACz = nucleus.z[nucleus.edges[i,2]] - nucleus.z[thirdVertex1];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        angles1[i] = acosd(dotProduct/(ABnorm*ACnorm));
        
        ABx = nucleus.x[nucleus.edges[i,1]] - nucleus.x[thirdVertex2];
        ABy = nucleus.y[nucleus.edges[i,1]] - nucleus.y[thirdVertex2];
        ABz = nucleus.z[nucleus.edges[i,1]] - nucleus.z[thirdVertex2];
        
        ACx = nucleus.x[nucleus.edges[i,2]] - nucleus.x[thirdVertex2];
        ACy = nucleus.y[nucleus.edges[i,2]] - nucleus.y[thirdVertex2];
        ACz = nucleus.z[nucleus.edges[i,2]] - nucleus.z[thirdVertex2];

        ABnorm = sqrt(ABx^2 + ABy^2 + ABz^2);
        ACnorm = sqrt(ACx^2 + ACy^2 + ACz^2);
        dotProduct = ABx*ACx + ABy*ACy + ABz*ACz;

        angles2[i] = acosd(dotProduct/(ABnorm*ACnorm));
    end

    curvatures = zeros(length(nucleus.x));
    vertexNormalUnitVectors = zeros(length(nucleus.x),3);

    for i = 1:length(nucleus.x)

        cotangentSum = [0;0;0];

        for j = nucleus.vertexEdges[i]

            cotangentSum = cotangentSum .+ (cotd(angles1[j]) + cotd(angles2[j])).*[nucleus.x[nucleus.edges[j,1]] - nucleus.x[nucleus.edges[j,2]] ; nucleus.y[nucleus.edges[j,1]] - nucleus.y[nucleus.edges[j,2]] ; nucleus.z[nucleus.edges[j,1]] - nucleus.z[nucleus.edges[j,2]]]

        end

        LaplaceBeltrami = 1/(2*sum(nucleus.voronoiAreas[i])).*cotangentSum;

        LaplaceBeltramiNorm = sqrt(sum(LaplaceBeltrami.^2));

        curvatures[i] = LaplaceBeltramiNorm./2;
        vertexNormalUnitVectors[i,:] = LaplaceBeltrami./LaplaceBeltramiNorm;

    end

    nucleus.curvatures = curvatures;
    nucleus.vertexNormalUnitVectors = vertexNormalUnitVectors;

end

function get_voronoi_areas!(nucleus)

    voronoiAreas = fill(Float64[], length(nucleus.x));

    for i = 1:length(nucleus.x)

        voronoiAreasTemp = zeros(Float64,length(nucleus.vertexTri[i,1]));

        for j = 1:length(nucleus.vertexTri[i,1])

            v = circshift(nucleus.tri[nucleus.vertexTri[i,1][j],:],(-nucleus.vertexTri[i,2][j]));
            centerX3 = sum(nucleus.x[v])/3;
            centerY3 = sum(nucleus.y[v])/3;
            centerZ3 = sum(nucleus.z[v])/3;

            tempArea = 0;

            for k = 2:3

                centerX2 = (nucleus.x[i] + nucleus.x[v[k]])/2;
                centerY2 = (nucleus.y[i] + nucleus.y[v[k]])/2;
                centerZ2 = (nucleus.z[i] + nucleus.z[v[k]])/2;

                dist1 = sqrt((centerX3 - centerX2).^2 + (centerY3 - centerY2).^2 + (centerZ3 - centerZ2).^2);
                dist2 = sqrt((nucleus.x[i] - centerX2).^2 + (nucleus.y[i] - centerY2).^2 + (nucleus.z[i] - centerZ2).^2);

                tempArea += dist1*dist2/2;

            end

            voronoiAreasTemp[j] = tempArea;

        end

        voronoiAreas[i] = voronoiAreasTemp;

    end

    nucleus.voronoiAreas = voronoiAreas;

end

function get_area_unit_vectors!(nucleus)


    baryocenters = zeros(Float64,size(nucleus.tri,1,),3);

    for i = 1:size(nucleus.tri,1)
        baryocenters[i,:] = [sum(nucleus.x[nucleus.tri[i,:]]) sum(nucleus.y[nucleus.tri[i,:]]) sum(nucleus.z[nucleus.tri[i,:]])]./3;
    end

    nucleus.areaUnitVectors = fill(Float64[], length(nucleus.x));

    for i = 1:length(nucleus.x)

        areaUnitVectors = zeros(Float64,length(nucleus.vertexTri[i,1]),3);

        j2 = 1;
        for j = nucleus.vertexTri[i,1]

            vectorTemp = [nucleus.x[i] - baryocenters[j,1] nucleus.y[i] - baryocenters[j,2] nucleus.z[i] - baryocenters[j,3]];

            vectorNorm = sqrt(sum(vectorTemp.^2));

            areaUnitVectors[j2,:] = vectorTemp./vectorNorm;
            
        end

        nucleus.areaUnitVectors[i] = areaUnitVectors;

    end

end

function get_triangle_normals!(nucleus)

    triangleNormalUnitVectors = zeros(Float64,size(nucleus.tri,1),3);

    for i = 1:size(nucleus.tri,1)

        tri = nucleus.tri[i,:];
        p1 = [nucleus.x[tri[1]], nucleus.y[tri[1]], nucleus.z[tri[1]]];
        p2 = [nucleus.x[tri[2]], nucleus.y[tri[2]], nucleus.z[tri[2]]];
        p3 = [nucleus.x[tri[3]], nucleus.y[tri[3]], nucleus.z[tri[3]]];

        normalVector = cross_product(p1,p2,p3);

        triangleNormalUnitVectors[i,:] = normalVector./sqrt(sum(normalVector.^2));

    end

    nucleus.triangleNormalUnitVectors = triangleNormalUnitVectors;

end

function cross_product(p1,p2,p3)

    return [((p2[2] - p1[2]) * (p3[3] - p1[3])) - ((p2[3] - p1[3]) * (p3[2] - p1[2]))
            ((p2[3] - p1[3]) * (p3[1] - p1[1])) - ((p2[1] - p1[1]) * (p3[3] - p1[3]))
            ((p2[1] - p1[1]) * (p3[2] - p1[2])) - ((p2[2] - p1[2]) * (p3[1] - p1[1]))];

end

function get_triangle_angles(nucleus)

    angles = zeros(Float64,size(nucleus.edges,1),1);
    
    for i = 1:size(nucleus.edges,1)

        normalVector1 = nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,1],:];
        normalVector2 = nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,2],:];

        dotProduct = normalVector1[1]*normalVector2[1] + normalVector1[2]*normalVector2[2] + normalVector1[3]*normalVector2[3];
        angles[i] = acosd(dotProduct);
    end

    return angles

end