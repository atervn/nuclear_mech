function get_volume_forces(nucleus)

    nucleusVolume = get_volume!(nucleus);

    pressure = -log10(nucleusVolume/nucleus.normalVolume);

    volumeForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:length(nucleus.x)
        forceMagnitude = pressure*sum(nucleus.voronoiAreas[i]);
        volumeForces[i,:] = forceMagnitude.*nucleus.vertexNormalUnitVectors[i,:];
    end

    return volumeForces

end

function get_area_forces(nucleus)

    nucleusArea = get_area!(nucleus);

    forceMagnitude = (nucleusArea - nucleus.normalArea)/nucleus.normalArea;

    areaForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:length(nucleus.x)

        areaForces[i,:] = -forceMagnitude.*sum(nucleus.areaUnitVectors[i],dims=1);

    end

    return areaForces

end

function get_bending_forces(nucleus)

    bendingForces = zeros(Float64,length(nucleus.x),3);
    
    angles = get_triangle_angles(nucleus);

    for i = 1:size(nucleus.edges,1)
    
        if nucleus.firstEdges[i] == 1
            moment = sind(angles[i] - nucleus.normalAngle);

            p1 = [nucleus.x[nucleus.edges[i,1]] nucleus.y[nucleus.edges[i,1]] nucleus.z[nucleus.edges[i,1]]];
            p2 = [nucleus.x[nucleus.edges[i,2]] nucleus.y[nucleus.edges[i,2]] nucleus.z[nucleus.edges[i,2]]];
            p3 = [nucleus.x[nucleus.edges3vertex[i,1]] nucleus.y[nucleus.edges3vertex[i,1]] nucleus.z[nucleus.edges3vertex[i,1]]];

            distance1 = line_point_distance(p1,p2,p3)
            
            force = moment*distance1.*nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,1],:]

            bendingForces[nucleus.edges3vertex[i,1],:] += force;
            bendingForces[nucleus.edges[i,1],:] += -0.5.*force;
            bendingForces[nucleus.edges[i,2],:] += -0.5.*force;

            p3 = [nucleus.x[nucleus.edges3vertex[i,2]] nucleus.y[nucleus.edges3vertex[i,2]] nucleus.z[nucleus.edges3vertex[i,2]]];

            distance2 = line_point_distance(p1,p2,p3)

            force = moment*distance2.*nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,2],:]

            bendingForces[nucleus.edges3vertex[i,2],:] += force;
            bendingForces[nucleus.edges[i,1],:] += -0.5.*force;
            bendingForces[nucleus.edges[i,2],:] += -0.5.*force;
        end
    end

    return bendingForces

end

function get_elastic_forces(nucleus)

    elasticForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:size(nucleus.edges,1)
        if nucleus.firstEdges[i] == 1

            vector = [nucleus.x[nucleus.edges[i,2]] - nucleus.x[nucleus.edges[i,1]],
                      nucleus.y[nucleus.edges[i,2]] - nucleus.y[nucleus.edges[i,1]],
                      nucleus.z[nucleus.edges[i,2]] - nucleus.z[nucleus.edges[i,1]]];
                      
            vectorNorm = norm(vector);

            unitVector = vector./vectorNorm;

            force = -(vectorNorm - nucleus.normalLength).*unitVector;

            elasticForces[nucleus.edges[i,1],:] += force;
            elasticForces[nucleus.edges[i,2],:] += -force;

        end
    end
    return elasticForces
end

function get_repulsion_forces(nucleus)

    repulsionForces = zeros(Float64,length(nucleus.x),3);

    idx = collect(1:length(nucleus.x));

    for i = 1:length(nucleus.x)

        idxTemp = copy(idx);

        sqrtDistances = (nucleus.x[i] .- nucleus.x).^2 .+ (nucleus.y[i] .- nucleus.y).^2 .+ (nucleus.z[i] .- nucleus.z).^2;

        selfAndNeighbors = sort([i; nucleus.neighbors[i]])

        deleteat!(sqrtDistances,selfAndNeighbors);
        deleteat!(idxTemp,selfAndNeighbors);

        closest = findmin(sqrtDistances);
        closest = idxTemp[closest[2]];

        nTri = length(nucleus.vertexTri[closest,1]);

        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Any}(undef,nTri);

        for j = 1:nTri

            tri = nucleus.vertexTri[closest,1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nucleus, i, tri)
        end

        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < 0.5

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = [nucleus.x[i] - closeCoords[1], nucleus.y[i] - closeCoords[2], nucleus.z[i] - closeCoords[3]]./closestDistance;
            
            repulsionForces[i,:] = (0.5 - closestDistance)^(3/2).*unitVector;

        end
    end

    return repulsionForces

end