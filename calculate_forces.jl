function get_volume_forces!(nuc,spar)

    nucleusVolume = get_volume!(nuc);

    pressure = -spar.bulkModulus*log10(nucleusVolume/nuc.normalVolume);

    if length(nuc.forces.volume) == 0
        nuc.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end

    for i = eachindex(nuc.vert)

        forceMagnitude = pressure*sum(nuc.voronoiAreas[i]);

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

            vector = nuc.edgeVectors[i];
            vectorNorm = norm(vector);

            unitVector = vector/vectorNorm;

            force = spar.laminaStiffness*(vectorNorm - nuc.normalLengths[i])*unitVector;

            nuc.forces.elastic[nuc.edges[i,1]] += force
            nuc.forces.elastic[nuc.edges[i,2]] -= force
        end
    end
end

function get_repulsion_forces!(nuc,spar,envelopeTree)

    if length(nuc.forces.envelopeRepulsion) == 0
        nuc.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end
    
    for i = 1:length(nuc.vert)

        closest,distance = knn(envelopeTree, nuc.vert[i],1,true,j -> any(j .== [i; nuc.neighbors[i]]))

        if distance[1] < mean(nuc.normalLengths)*2

            nTri = length(nuc.vertexTri[closest[1],1]);

            closePointDistances = Vector{Float64}(undef,nTri);
            closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);

            for j = 1:nTri

                tri = nuc.vertexTri[closest[1],1][j];

                closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc, nuc.vert[i], tri)
            end
            
            closestPoint = findmin(closePointDistances);
            closestDistance = closestPoint[1];

            if closestDistance < spar.repulsionDistance

                closeCoords = closePointCoordinates[closestPoint[2]];

                unitVector = (nuc.vert[i] - closeCoords)/closestDistance;
                
                forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)
                nuc.forces.envelopeRepulsion[i] = forceMagnitude*unitVector;
            else
                nuc.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
            end
        else
            nuc.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
        end
    end
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

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc, nuc.vert[i], tri, pip);
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

function get_aspiration_forces(nuc,pip,spar)

    force = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    for i = eachindex(nuc.vert)
        if nuc.vert[i][1]>0 && sqrt(nuc.vert[i][2]^2 + nuc.vert[i][3]^2) < 3
            forceMagnitude = 1*sum(nuc.voronoiAreas[i]);
            force[i] = forceMagnitude*nuc.vertexNormalUnitVectors[i]
        else
            force[i] = Vec(0.,0.,0.);
        end
    end

    return force

end

function get_linear_chromatin_forces!(chro,spar)

    for i = 1:spar.chromatinNumber
        chro.forces.strandLinear[i][1:end-1] += 10*(chro.vectorNorms[i] .- spar.chroVertexDistance).*chro.vectors[i]./chro.vectorNorms[i];
        chro.forces.strandLinear[i][2:end] -= chro.forces.strandLinear[i][1:end-1];
    end
end

function get_bending_chromatin_forces!(chro,spar,normAngle)
    
    for i = 1:spar.chromatinNumber

        vectors1 = -chro.vectors[i][1:end-1]
        vectors2 = chro.vectors[i][2:end]

        angles = dot.(vectors1, vectors2)./(chro.vectorNorms[i][1:end-1].*chro.vectorNorms[i][2:end]);

        angles[angles .>= 1] .= 1.; 

        angles = acosd.(angles);
        
        unitVectors1 = cross.(vectors1,cross.(vectors1,vectors2));
        unitVectors2 = cross.(vectors2,cross.(vectors2,vectors1));

        unitVectors1 = unitVectors1./norm.(unitVectors1);
        unitVectors2 = unitVectors2./norm.(unitVectors2);
        
        chro.forces.strandBending[i][1:end-2] += -1e-2.*(angles .- normAngle).*unitVectors1;
        chro.forces.strandBending[i][3:end] += -1e-2.*(angles .- normAngle).*unitVectors2;
    end
end

function get_chromation_chromation_repulsion_forces!(chro,spar,chromatinTree)

    for i = 1:spar.chromatinNumber

        closeVertices = inrange(chromatinTree, chro.vert[chro.strandIdx[i]], spar.repulsionDistance, false)

        for j = 1:spar.chromatinLength
            for k = eachindex(closeVertices[j])
                if chro.strandIdx[i][j] != closeVertices[j][k]
                    vector = chro.strandVert[i][j] - chro.vert[closeVertices[j][k]]
                    vectorNorm = norm(vector)
                    chro.forces.strandChroRepulsion[i][j] += -spar.repulsionConstant*(vectorNorm - spar.repulsionDistance)*vector/vectorNorm ;#24*0.01/vectorNorm^2*(2*0.12^12/vectorNorm^12 - 0.12^6/vectorNorm^6)*vector;
                end
            end
        end
    end
end

function get_random_fluctuations(spar,dt)

    strength = sqrt(2*spar.viscosity*1.38e-23*300/(dt))/spar.viscosity/spar.scalingLength*spar.scalingTime;

    fluctuationForces = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = 1:spar.chromatinLength*spar.chromatinNumber
        fluctuationForces[i] = strength.*Vec(randn(),randn(),randn())
    end

    return fluctuationForces

end

function initialize_chromatin_forces!(chro)

    for i = 1:length(chro.vert)

        chro.forces.linear[i] = Vec(0.,0.,0.)
        chro.forces.bending[i] = Vec(0.,0.,0.)
        chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
        chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)

    end
end

function get_envelope_chromatin_repulsion_forces!(nuc,chro,spar,envelopeTree)

    nuc.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = eachindex(nuc.vert)
        nuc.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef, length(chro.vert));

    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(envelopeTree, chro.vert[i],1,true)

             if distance[1] < mean(nuc.normalLengths)*2

            nTri = length(nuc.vertexTri[closest[1],1]);

            closePointDistances = Vector{Float64}(undef,nTri);
            closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);
            closePointVertices = Vector{Vector{Int64}}(undef, nTri);

            for j = 1:nTri

                tri = nuc.vertexTri[closest[1],1][j];

                closePointDistances[j],closePointCoordinates[j],closePointVertices[j] = vertex_triangle_distance(nuc, chro.vert[i], tri)
            end
            
            closestPoint = findmin(closePointDistances);
            closestDistance = closestPoint[1];

            if closestDistance < spar.repulsionDistance

                closeCoords = closePointCoordinates[closestPoint[2]];

                unitVector = (chro.vert[i] - closeCoords)/closestDistance;
                
                forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)
                chro.forces.enveRepulsion[i] = forceMagnitude*unitVector;

                closeVertices = closePointVertices[closestPoint[2]]; 

                if length(closeVertices) == 1

                    nuc.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - nuc.vert[closeVertices[2]][1])/(nuc.vert[closeVertices[2]][1] - nuc.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    nuc.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    nuc.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(nuc.vert[closeVertices[1]] - nuc.vert[closeVertices[2]])*norm(nuc.vert[closeVertices[1]] - nuc.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - nuc.vert[closeVertices[2]])*norm(closeCoords - nuc.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - nuc.vert[closeVertices[1]])*norm(closeCoords - nuc.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - nuc.vert[closeVertices[1]])*norm(closeCoords - nuc.vert[closeVertices[2]])

                    nuc.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    nuc.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    nuc.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;

                    closeCoords
                end

            else
                chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
            end
        else
            chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
        end
    end
end

function get_micromanipulation_forces(nuc,mm,spar)

    micromanipulation = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = eachindex(nuc.vert)
        micromanipulation[i] = Vec(0.,0.,0.);
    end

    micromanipulation[mm.leftmostVertex] = 100*(mm.leftmostVertexPosition .- nuc.vert[mm.leftmostVertex])
    micromanipulation[mm.leftNeighbors] = 50*(mm.leftNeigborPositions .- nuc.vert[mm.leftNeighbors])

    forceVector = 50*Vec(1.,0.,0.);
    micromanipulation[mm.rightmostVertex] = forceVector
    for i = eachindex(mm.rightNeighbors)
        micromanipulation[mm.rightNeighbors[i]] = 0.5*forceVector
    end

    return micromanipulation
end

function get_lad_forces(nuc,chro,spar)

    ladEnveForces = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    ladChroForces = Vector{Vec{3,Float64}}(undef, length(chro.vert));

    for i = 1:length(nuc.vert)
        ladEnveForces[i] = Vec(0.,0.,0.) 
    end
    for i = 1:length(chro.vert)
        ladChroForces[i] = Vec(0.,0.,0.) 
    end

    for i = 1:spar.chromatinNumber

        for j = 1:length(chro.lads[i])
            vector = nuc.vert[nuc.lads[i][j]] - chro.vert[chro.strandIdx[i][chro.lads[i][j]]]
            distance = norm(vector)
            magnitude = -1*(distance - 0.5)
            ladEnveForces[nuc.lads[i][j]] = magnitude*vector/distance
            ladChroForces[chro.strandIdx[i][chro.lads[i][j]]] = -magnitude*vector/distance

        end
    end

    return ladEnveForces, ladChroForces

end