function get_volume_forces!(nuc,spar)

    nucleusVolume = get_volume!(nuc);

    pressure = -spar.bulkModulus*log10(nucleusVolume/nuc.normalVolume);

    if length(nuc.forces.volume) == 0
        nuc.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end

    for i = eachindex(nuc.vert)
        if nuc.vert[i][1]>0 && sqrt(nuc.vert[i][2]^2 + nuc.vert[i][3]^2) < 0.3
            forceMagnitude = 20*sum(nuc.voronoiAreas[i]);
        else
            forceMagnitude = pressure*sum(nuc.voronoiAreas[i]);
        end
        nuc.forces.volume[i] = forceMagnitude*nuc.vertexNormalUnitVectors[i]
    end
end

function get_area_forces!(nuc,spar)

    triangleAreas = get_area!(nuc);
    nucleusArea = sum(triangleAreas);

    forceMagnitude = spar.areaCompressionStiffness*(nucleusArea - nuc.normalArea)/nuc.normalArea;

    if length(nuc.forces.area) == 0
        nuc.forces.area = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end

    for i = 1:length(nuc.vert)
        nuc.forces.area[i] = -forceMagnitude*mean(nuc.areaUnitVectors[i])
    end

    for i = 1:size(nuc.tri,1)

        baryocenter = mean(nuc.trii[i]);
        magnitude = 0.1*spar.areaCompressionStiffness*(triangleAreas[i] - nuc.normalTriangleAreas[i])/nuc.normalTriangleAreas[i];
        
        for j = 1:3

            vector = nuc.trii[i][j] - baryocenter;
            unitVector = vector/norm(vector)

            nuc.forces.area[nuc.tri[i,j]] += -magnitude*unitVector
        end

    end

end

function get_bending_forces!(nuc,spar)

    if length(nuc.forces.bending) == 0
        nuc.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end

    angles = get_triangle_angles(nuc);

    for i = 1:size(nuc.edges,1)
    
        if nuc.firstEdges[i] == 1
            moment = spar.bendingStiffness*sind(angles[i] - nuc.normalAngle);

            
            distance1 = line_point_distance(nuc.ep1[i][1] - nuc.ep2[i][1], nuc.ep1[i][1] - nuc.ep31[i][1])
            
            force = moment*distance1*nuc.triangleNormalUnitVectors[nuc.edgesTri[i,1]]

            nuc.forces.bending[nuc.edges3vertex[i,1]] += force;
            
            nuc.forces.bending[nuc.edges[i,1]] += -0.5.*force;
            
            nuc.forces.bending[nuc.edges[i,2]] += -0.5.*force;

            distance2 = line_point_distance(nuc.ep1[i][1] - nuc.ep2[i][1], nuc.ep1[i][1] - nuc.ep32[i][1])
            
            force = moment*distance2.*nuc.triangleNormalUnitVectors[nuc.edgesTri[i,2]]

            nuc.forces.bending[nuc.edges3vertex[i,2]] += force
            
            nuc.forces.bending[nuc.edges[i,1]] += -0.5.*force
            
            nuc.forces.bending[nuc.edges[i,2]] += -0.5.*force

        end
    end
end

function get_elastic_forces!(nuc,spar)

    nuc.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    for i = eachindex(nuc.vert)
        nuc.forces.elastic[i] = Vec(0.,0.,0.);
    end

    for i = 1:size(nuc.edges,1)
        if nuc.firstEdges[i] == 1

            vector = nuc.vert[nuc.edges[i,2]] - nuc.vert[nuc.edges[i,1]]
            vectorNorm = norm(vector);

            unitVector = vector/vectorNorm;

            force = spar.laminaStiffness*(vectorNorm - nuc.normalLengths[i])*unitVector;

            

            nuc.forces.elastic[nuc.edges[i,1]] += force
            nuc.forces.elastic[nuc.edges[i,2]] -= force
        end
    end
end

function get_repulsion_forces!(nuc,spar)

    if length(nuc.forces.repulsion) == 0
        nuc.forces.repulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end

    #idx = collect(1:length(nuc.x));

    tree = KDTree(nuc.vert);

    for i = 1:length(nuc.vert)

        #idxTemp = copy(idx);

#        sqrtDistances = (nuc.x[i] .- nuc.x).^2 .+ (nuc.y[i] .- nuc.y).^2 .+ (nuc.z[i] .- nuc.z).^2;

 #       selfAndNeighbors = sort([i; nuc.neighbors[i]])

  #      deleteat!(sqrtDistances,selfAndNeighbors);
   #     deleteat!(idxTemp,selfAndNeighbors);

    #    closest = findmin(sqrtDistances);
     #   closest = idxTemp[closest[2]];

        closest = knn(tree, nuc.vert[i],1,true,j -> any(j .== [i; nuc.neighbors[i]]))

        nTri = length(nuc.vertexTri[closest[1],1]);

        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);

        for j = 1:nTri

            tri = nuc.vertexTri[closest[1],1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc, i, tri)
        end
        
        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < spar.repulsionDistance

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = (nuc.vert[i] - closeCoords)/closestDistance;
            
            forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)

            nuc.forces.repulsion[i] = forceMagnitude*unitVector;

        end
    end
end

function flat_repulsion_forces(nuc,spar)

    repulsionForces = zeros(Float64,length(nuc.x),3);

    for i = 1:length(nuc.x)

        if nuc.z[i] < spar.repulsionDistance

            distance = nuc.z[i];
            repulsionForces[i,3] = spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2);

        end
    end
end

function flat_repulsion_forces2(nuc,spar,t)

    repulsionForces = zeros(Float64,length(nuc.x),3);

    planePos = -0.002*(t) + 1.5;
    if planePos < 0.4
        planePos = 0.4;
    end

    for i = 1:length(nuc.x)

        if nuc.z[i] > planePos - spar.repulsionDistance

            distance = abs(planePos - nuc.z[i]);
            repulsionForces[i,3] = -spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2);
        elseif nuc.z[i] < -planePos + spar.repulsionDistance
            distance = abs(-planePos - nuc.z[i]);
            repulsionForces[i,3] = spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2);
        end

    end

    return repulsionForces

end

function get_aspiration_repulsion_forces(nuc,pip,spar)


    repulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    tree = KDTree(pip.vert);

    for i = 1:length(nuc.vert)

        closest = knn(tree, nuc.vert[i],1,true)

        nTri = length(pip.vertexTri[closest[1]][1]);
        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Any}(undef,nTri);

        for j = 1:nTri

            tri = pip.vertexTri[closest[1]][1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc, i, tri, pip);
        end

        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < spar.repulsionDistance

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = (nuc.vert[i] - closeCoords)./closestDistance;
            
            forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2);

            repulsion[i] = forceMagnitude*unitVector;
        else
            repulsion[i] = Vec(0.,0.,0.)
        end
    end

    return repulsion

end