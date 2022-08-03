function get_volume_forces!(nuc,spar)

    nucleusVolume = get_volume!(nuc);

    pressure = -spar.bulkModulus*log10(nucleusVolume/nuc.normalVolume);

    if length(nuc.forces.volumeX) == 0
        nuc.forces.volumeX = zeros(Float64,length(nuc.x));
        nuc.forces.volumeY = zeros(Float64,length(nuc.x));
        nuc.forces.volumeZ = zeros(Float64,length(nuc.x));
    end

    for i = 1:length(nuc.x)
        forceMagnitude = pressure*sum(nuc.voronoiAreas[i]);
        nuc.forces.volumeX[i] = forceMagnitude.*nuc.vertexNormalUnitVectors[i,1];
        nuc.forces.volumeY[i] = forceMagnitude.*nuc.vertexNormalUnitVectors[i,2];
        nuc.forces.volumeZ[i] = forceMagnitude.*nuc.vertexNormalUnitVectors[i,3];
    end

end

function get_area_forces!(nuc,spar)

    nucleusArea = get_area!(nuc);

    forceMagnitude = spar.areaCompressionStiffness*(nucleusArea - nuc.normalArea)/nuc.normalArea;

    if length(nuc.forces.areaX) == 0
        nuc.forces.areaX = zeros(Float64,length(nuc.x));
        nuc.forces.areaY = zeros(Float64,length(nuc.x));
        nuc.forces.areaZ = zeros(Float64,length(nuc.x));
    end

    for i = 1:length(nuc.x)

        nuc.forces.areaX[i] = -forceMagnitude.*sum(nuc.areaUnitVectors[i][:,1]);
        nuc.forces.areaY[i] = -forceMagnitude.*sum(nuc.areaUnitVectors[i][:,2]);
        nuc.forces.areaZ[i] = -forceMagnitude.*sum(nuc.areaUnitVectors[i][:,3]);

    end

end

function get_bending_forces!(nuc,spar)

    if length(nuc.forces.bendingX) == 0
        nuc.forces.bendingX = zeros(Float64,length(nuc.x));
        nuc.forces.bendingY = zeros(Float64,length(nuc.x));
        nuc.forces.bendingZ = zeros(Float64,length(nuc.x));
    end

    angles = get_triangle_angles(nuc);

    for i = 1:size(nuc.edges,1)
    
        if nuc.firstEdges[i] == 1
            moment = spar.bendingStiffness*sind(angles[i] - nuc.normalAngle);

            p1 = [nuc.x[nuc.edges[i,1]], nuc.y[nuc.edges[i,1]], nuc.z[nuc.edges[i,1]]];
            p2 = [nuc.x[nuc.edges[i,2]], nuc.y[nuc.edges[i,2]], nuc.z[nuc.edges[i,2]]];
            p3 = [nuc.x[nuc.edges3vertex[i,1]], nuc.y[nuc.edges3vertex[i,1]], nuc.z[nuc.edges3vertex[i,1]]];

            distance1 = line_point_distance(p1,p2,p3)
            
            force = moment*distance1.*nuc.triangleNormalUnitVectors[nuc.edgesTri[i,1],:]

            nuc.forces.bendingX[nuc.edges3vertex[i,1]] += force[1];
            nuc.forces.bendingY[nuc.edges3vertex[i,1]] += force[2];
            nuc.forces.bendingZ[nuc.edges3vertex[i,1]] += force[3];
            
            nuc.forces.bendingX[nuc.edges[i,1]] += -0.5.*force[1];
            nuc.forces.bendingY[nuc.edges[i,1]] += -0.5.*force[2];
            nuc.forces.bendingZ[nuc.edges[i,1]] += -0.5.*force[3];
            
            nuc.forces.bendingX[nuc.edges[i,2]] += -0.5.*force[1];
            nuc.forces.bendingY[nuc.edges[i,2]] += -0.5.*force[2];
            nuc.forces.bendingZ[nuc.edges[i,2]] += -0.5.*force[3];

            p3 = [nuc.x[nuc.edges3vertex[i,2]], nuc.y[nuc.edges3vertex[i,2]], nuc.z[nuc.edges3vertex[i,2]]];

            distance2 = line_point_distance(p1,p2,p3)

            force = moment*distance2.*nuc.triangleNormalUnitVectors[nuc.edgesTri[i,2],:]

            nuc.forces.bendingX[nuc.edges3vertex[i,2]] += force[1];
            nuc.forces.bendingY[nuc.edges3vertex[i,2]] += force[2];
            nuc.forces.bendingZ[nuc.edges3vertex[i,2]] += force[3];
            
            nuc.forces.bendingX[nuc.edges[i,1]] += -0.5.*force[1];
            nuc.forces.bendingY[nuc.edges[i,1]] += -0.5.*force[2];
            nuc.forces.bendingZ[nuc.edges[i,1]] += -0.5.*force[3];
            
            nuc.forces.bendingX[nuc.edges[i,2]] += -0.5.*force[1];
            nuc.forces.bendingY[nuc.edges[i,2]] += -0.5.*force[2];
            nuc.forces.bendingZ[nuc.edges[i,2]] += -0.5.*force[3];

        end
    end
end

function get_elastic_forces!(nuc,spar)

    nuc.forces.elasticX = zeros(Float64,length(nuc.x));
    nuc.forces.elasticY = zeros(Float64,length(nuc.x));
    nuc.forces.elasticZ = zeros(Float64,length(nuc.x));

    for i = 1:size(nuc.edges,1)
        if nuc.firstEdges[i] == 1

            vector = [nuc.x[nuc.edges[i,2]] - nuc.x[nuc.edges[i,1]],
                      nuc.y[nuc.edges[i,2]] - nuc.y[nuc.edges[i,1]],
                      nuc.z[nuc.edges[i,2]] - nuc.z[nuc.edges[i,1]]];
                      
            vectorNorm = norm(vector);

            unitVector = vector./vectorNorm;

            force = spar.laminaStiffness*(vectorNorm - nuc.normalLengths[i]).*unitVector;

            nuc.forces.elasticX[nuc.edges[i,1]] += force[1];
            nuc.forces.elasticY[nuc.edges[i,1]] += force[2];
            nuc.forces.elasticZ[nuc.edges[i,1]] += force[3];
            
            nuc.forces.elasticX[nuc.edges[i,2]] -= force[1];
            nuc.forces.elasticY[nuc.edges[i,2]] -= force[2];
            nuc.forces.elasticZ[nuc.edges[i,2]] -= force[3];
        end
    end
end

function get_repulsion_forces!(nuc,spar)

    if length(nuc.forces.repulsionX) == 0
        nuc.forces.repulsionX = zeros(Float64,length(nuc.x));
        nuc.forces.repulsionY = zeros(Float64,length(nuc.x));
        nuc.forces.repulsionZ = zeros(Float64,length(nuc.x));
    end

    idx = collect(1:length(nuc.x));

    for i = 1:length(nuc.x)

        idxTemp = copy(idx);

        sqrtDistances = (nuc.x[i] .- nuc.x).^2 .+ (nuc.y[i] .- nuc.y).^2 .+ (nuc.z[i] .- nuc.z).^2;

        selfAndNeighbors = sort([i; nuc.neighbors[i]])

        deleteat!(sqrtDistances,selfAndNeighbors);
        deleteat!(idxTemp,selfAndNeighbors);

        closest = findmin(sqrtDistances);
        closest = idxTemp[closest[2]];

        nTri = length(nuc.vertexTri[closest,1]);

        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Any}(undef,nTri);

        for j = 1:nTri

            tri = nuc.vertexTri[closest,1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc, i, tri)
        end

        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < spar.repulsionDistance

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = [nuc.x[i] - closeCoords[1], nuc.y[i] - closeCoords[2], nuc.z[i] - closeCoords[3]]./closestDistance;
            
            forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)

            nuc.forces.repulsionX[i] = forceMagnitude*unitVector[1];
            nuc.forces.repulsionY[i] = forceMagnitude*unitVector[2];
            nuc.forces.repulsionZ[i] = forceMagnitude*unitVector[3];

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
    if planePos < 0.2
        planePos = 0.2;
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