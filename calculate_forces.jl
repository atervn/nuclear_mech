function get_volume_forces!(nuc,spar)

    nucleusVolume = get_volume!(nuc.enve);

    pressure = -spar.bulkModulus*log10(nucleusVolume/nuc.enve.normalVolume);

    nuc.enve.forces.volume = pressure*nuc.enve.voronoiAreas.*nuc.enve.vertexNormalUnitVectors

end

function get_area_forces!(nuc,spar)

    nucleusArea = sum(nuc.enve.triangleAreas);

    forceMagnitude = spar.areaCompressionStiffness*(nucleusArea - nuc.enve.normalArea)/nuc.enve.normalArea;

    for i = 1:length(nuc.enve.vert)
        nuc.enve.forces.area[i] = -forceMagnitude*mean(nuc.enve.areaUnitVectors[i])
    end

    for i = eachindex(nuc.enve.tri)

        baryocenter = mean(nuc.enve.vert[nuc.enve.tri[i]]);
        magnitude = 0.1*spar.areaCompressionStiffness*(nuc.enve.triangleAreas[i] - nuc.enve.normalTriangleAreas[i])/nuc.enve.normalTriangleAreas[i];
        
        for j = 1:3

            vector = nuc.enve.vert[nuc.enve.tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            nuc.enve.forces.area[nuc.enve.tri[i][j]] += -magnitude*unitVector
        end

    end

end

function get_bending_forces!(nuc,spar)

    nuc.enve.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));
    for i = 1:length(nuc.enve.vert)
        nuc.enve.forces.bending[i] = Vec(0.,0.,0.)
    end

    angles = get_triangle_angles(nuc.enve);

    moment = spar.bendingStiffness*(angles .- nuc.enve.normalAngle);

    for i = eachindex(nuc.enve.edges)
        
        if nuc.enve.firstEdges[i] == 1
            
            distance1 = line_point_distance(nuc.enve.edgeVectors[i], nuc.enve.edgeVectors[nuc.enve.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(nuc.enve.edgeVectors[i], nuc.enve.edgeVectors[nuc.enve.edgeThirdVertices[i][2]])
            
            force1 = -moment[i]/distance1*nuc.enve.triangleNormalUnitVectors[nuc.enve.edgesTri[i][1]]
            force2 = -moment[i]/distance2*nuc.enve.triangleNormalUnitVectors[nuc.enve.edgesTri[i][2]]

            counterForces = 0.5*force1 + 0.5*force2

            nuc.enve.forces.bending[nuc.enve.edges3Vertex[i][1]] += force1;
            nuc.enve.forces.bending[nuc.enve.edges3Vertex[i][2]] += force2;
            nuc.enve.forces.bending[nuc.enve.edges[i][1]] -= counterForces;
            nuc.enve.forces.bending[nuc.enve.edges[i][2]] -= counterForces;

        end
    end
end

function get_elastic_forces!(nuc,spar)

    nuc.enve.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));

    for i = eachindex(nuc.enve.vert)
        nuc.enve.forces.elastic[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(nuc.enve.edges)
        if nuc.enve.firstEdges[i] == 1

            force = spar.laminaStiffness*(nuc.enve.edgeVectorNorms[i] - nuc.enve.normalLengths[i])*nuc.enve.edgeUnitVectors[i];

            nuc.enve.forces.elastic[nuc.enve.edges[i][1]] += force
            nuc.enve.forces.elastic[nuc.enve.edges[i][2]] -= force
        end
    end
end

function get_repulsion_forces!(nuc,spar,envelopeTree)

    if length(nuc.enve.forces.envelopeRepulsion) == 0
        nuc.enve.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));
    end
    
    for i = 1:length(nuc.enve.vert)

        closest,distance = knn(envelopeTree, nuc.enve.vert[i],1,true,j -> any(j .== [i; nuc.enve.neighbors[i]]))

        if distance[1] < mean(nuc.enve.normalLengths)*1.5

            nTri = length(nuc.enve.vertexTri[closest[1]]);

            closePointDistances = Vector{Float64}(undef,nTri);
            closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);

            for j = 1:nTri

                tri = nuc.enve.vertexTri[closest[1]][j];

                closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc.enve, nuc.enve.vert[i], tri)
            end
            
            closestPoint = findmin(closePointDistances);
            closestDistance = closestPoint[1];

            if closestDistance < spar.repulsionDistance

                closeCoords = closePointCoordinates[closestPoint[2]];

                unitVector = (nuc.enve.vert[i] - closeCoords)/closestDistance;
                
                forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)
                nuc.enve.forces.envelopeRepulsion[i] = forceMagnitude*unitVector;
            else
                nuc.enve.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
            end
        else
            nuc.enve.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
        end
    end
end

function get_aspiration_repulsion_forces(nuc,pip,spar)


    repulsion = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));

    tree = KDTree(pip.vert);

    for i = 1:length(nuc.enve.vert)

        closest = knn(tree, nuc.enve.vert[i],1,true)

        nTri = length(pip.vertexTri[closest[1]][1]);
        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Any}(undef,nTri);

        for j = 1:nTri

            tri = pip.vertexTri[closest[1]][1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(nuc.enve, nuc.enve.vert[i], tri, pip);
        end

        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < spar.repulsionDistance

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = (nuc.enve.vert[i] - closeCoords)./closestDistance;
            
            forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2);

            repulsion[i] = forceMagnitude*unitVector;
        else
            repulsion[i] = Vec(0.,0.,0.)
        end
    end

    return repulsion

end

function get_aspiration_forces(nuc,pip,spar)

    force = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));

    for i = eachindex(nuc.enve.vert)
        if nuc.enve.vert[i][1]>0 && sqrt(nuc.enve.vert[i][2]^2 + nuc.enve.vert[i][3]^2) < 3
            forceMagnitude = 1*sum(nuc.enve.voronoiAreas[i]);
            force[i] = forceMagnitude*nuc.enve.vertexNormalUnitVectors[i]
        else
            force[i] = Vec(0.,0.,0.);
        end
    end

    return force

end

function get_linear_chromatin_forces!(nuc,spar)

    for i = 1:length(nuc.chro.vert)
        nuc.chro.forces.linear[i] = Vec(0.,0.,0.)
    end

    for i = 1:spar.chromatinNumber
        nuc.chro.forces.strandLinear[i][1:end-1] += spar.chromatinStiffness*(nuc.chro.vectorNorms[i] .- spar.chroVertexDistance).*nuc.chro.vectors[i]./nuc.chro.vectorNorms[i];
        nuc.chro.forces.strandLinear[i][2:end] -= nuc.chro.forces.strandLinear[i][1:end-1];
    end
end

function get_bending_chromatin_forces!(nuc,spar)
    
    for i = 1:spar.chromatinNumber

        vectors1 = -nuc.chro.vectors[i][1:end-1]
        vectors2 = nuc.chro.vectors[i][2:end]

        angles = dot.(vectors1, vectors2)./(nuc.chro.vectorNorms[i][1:end-1].*nuc.chro.vectorNorms[i][2:end]);

        angles[angles .>= 1] .= 1.; 

        angles = acos.(angles);
        
        unitVectors1 = cross.(vectors1,cross.(vectors1,vectors2));
        unitVectors2 = cross.(vectors2,cross.(vectors2,vectors1));

        unitVectors1 = unitVectors1./norm.(unitVectors1);
        unitVectors2 = unitVectors2./norm.(unitVectors2);
        
        nuc.chro.forces.strandBending[i][1:end-2] = -spar.chromatinBendingModulus./nuc.chro.vectorNorms[i][1:end-1].*(angles .- spar.chromatinNormalAngle).*unitVectors1;
        nuc.chro.forces.strandBending[i][3:end] = -spar.chromatinBendingModulus./nuc.chro.vectorNorms[i][2:end].*(angles .- spar.chromatinNormalAngle).*unitVectors2;
    end
end

function get_chromation_chromation_repulsion_forces!(nuc,spar,chromatinTree)

    for i = 1:length(nuc.chro.vert)
        nuc.chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinNumber

        closeVertices = inrange(chromatinTree, nuc.chro.vert[nuc.chro.strandIdx[i]], spar.repulsionDistance, false)

        for j = 1:spar.chromatinLength
            for k = eachindex(closeVertices[j])
                if nuc.chro.strandIdx[i][j] != closeVertices[j][k]
                    vector = nuc.chro.strandVert[i][j] - nuc.chro.vert[closeVertices[j][k]]
                    vectorNorm = norm(vector)
                    nuc.chro.forces.strandChroRepulsion[i][j] += -spar.repulsionConstant*(vectorNorm - spar.repulsionDistance)*vector/vectorNorm ;#24*0.01/vectorNorm^2*(2*0.12^12/vectorNorm^12 - 0.12^6/vectorNorm^6)*vector;
                end
            end
        end
    end
end

function get_random_fluctuations(spar,dt)

    strength = sqrt(2*spar.boltzmannConst*spar.temperature/(dt));

    fluctuationForces = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = 1:spar.chromatinLength*spar.chromatinNumber
        fluctuationForces[i] = strength.*Vec(randn(),randn(),randn())
    end

    return fluctuationForces

end

function get_random_enve_fluctuations(spar,nuc,dt)

    strength = sqrt(0*spar.boltzmannConst*spar.temperature/(dt));

    fluctuationForces = Vector{Vec{3,Float64}}(undef,length(nuc.enve.vert))
    for i = 1:length(nuc.enve.vert)
        fluctuationForces[i] = strength.*Vec(randn(),randn(),randn())
    end

    return fluctuationForces

end

function initialize_chromatin_forces!(nuc)

    for i = 1:length(nuc.chro.vert)

        nuc.chro.forces.linear[i] = Vec(0.,0.,0.)
        nuc.chro.forces.bending[i] = Vec(0.,0.,0.)
        nuc.chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
        nuc.chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
        nuc.chro.forces.ladChroForces[i] = Vec(0.,0.,0.)
    end
end

function get_envelope_chromatin_repulsion_forces!(nuc,spar,envelopeTree)

    nuc.enve.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));
    for i = eachindex(nuc.enve.vert)
        nuc.enve.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(nuc.chro.vert)
        nuc.chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(envelopeTree, nuc.chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(nuc.enve.vertexTri[closest[1]]);

            neighbors = nuc.enve.neighbors[closest[1]]
            neigborsCoords = nuc.enve.vert[neighbors]

            triangles = nuc.enve.tri[nuc.enve.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(nuc.enve.vert[triangles[j][k]] - nuc.chro.vert[i])
                end
            end
            
            tri = nuc.enve.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(nuc.enve, nuc.chro.vert[i], tri)

            getForce = false
            
            unitVector = (nuc.chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = nuc.enve.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) >= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (nuc.chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = nuc.enve.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(nuc.enve.edges,1) .== closeVertices[1] .&& getindex.(nuc.enve.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = nuc.enve.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = nuc.enve.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) >= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                nuc.chro.forces.enveRepulsion[i] = forceMagnitude*unitVector;

                if length(closeVertices) == 1

                    nuc.enve.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - nuc.enve.vert[closeVertices[2]][1])/(nuc.enve.vert[closeVertices[2]][1] - nuc.enve.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    nuc.enve.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    nuc.enve.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(nuc.enve.vert[closeVertices[1]] - nuc.enve.vert[closeVertices[2]])*norm(nuc.enve.vert[closeVertices[1]] - nuc.enve.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - nuc.enve.vert[closeVertices[2]])*norm(closeCoords - nuc.enve.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - nuc.enve.vert[closeVertices[1]])*norm(closeCoords - nuc.enve.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - nuc.enve.vert[closeVertices[1]])*norm(closeCoords - nuc.enve.vert[closeVertices[2]])

                    nuc.enve.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    nuc.enve.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    nuc.enve.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;

                    closeCoords
                end
            end
        end
    end
end

function get_micromanipulation_forces(nuc,mm,spar)

    micromanipulation = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));
    for i = eachindex(nuc.enve.vert)
        micromanipulation[i] = Vec(0.,0.,0.);
    end


    micromanipulation[mm.leftmostVertex] = 2*spar.pullingForce*(mm.leftmostVertexPosition .- nuc.enve.vert[mm.leftmostVertex])
    # micromanipulation[mm.leftNeighbors] = spar.pullingForce*(mm.leftNeigborPositions .- nuc.enve.vert[mm.leftNeighbors])

    forceVector = spar.pullingForce*Vec(1.,0.,0.);
    micromanipulation[mm.rightmostVertex] = forceVector
    # for i = eachindex(mm.rightNeighbors)
    #     micromanipulation[mm.rightNeighbors[i]] = 0.5*1/length(mm.rightNeighbors)*forceVector
    # end
    

    return micromanipulation
end

function get_lad_forces!(nuc,spar)

    nuc.enve.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));

    for i = 1:length(nuc.enve.vert)
        nuc.enve.forces.ladEnveForces[i] = Vec(0.,0.,0.) 
    end

    for i = 1:spar.chromatinNumber

        for j = 1:length(nuc.chro.lads[i])
            vector = nuc.enve.vert[nuc.enve.lads[i][j]] - nuc.chro.vert[nuc.chro.strandIdx[i][nuc.chro.lads[i][j]]]
            distance = norm(vector)
            magnitude = -spar.ladStrength*(distance - 0.2)
            nuc.enve.forces.ladEnveForces[nuc.enve.lads[i][j]] = magnitude*vector/distance
            nuc.chro.forces.ladChroForces[nuc.chro.strandIdx[i][nuc.chro.lads[i][j]]] = -magnitude*vector/distance

        end
    end

end

function get_plane_repulsion(nuc,ext,spar)

    planeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert));
    for i = eachindex(nuc.enve.vert)
        planeRepulsion[i] = Vec(0.,0.,0.);
    end

    touchTop = zeros(Bool,length(nuc.enve.vert))

    for i = eachindex(nuc.enve.vert)
        if ext[1] - nuc.enve.vert[i][3] < spar.repulsionDistance
            
            touchTop[i] = true

            if ext[1] - nuc.enve.vert[i][3] < 0
                planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*0.1)
            else
                distance = ext[1] - nuc.enve.vert[i][3];
                planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
                ext[3][i] = true
            end
        end
    end

    for i = eachindex(nuc.enve.vert)
        if nuc.enve.vert[i][3] - ext[2] < spar.repulsionDistance
            if nuc.enve.vert[i][3] - ext[2] < 0
                planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*0.1)
            else
                distance = nuc.enve.vert[i][3] - ext[2];
                planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
            end
        end
    end

    return planeRepulsion, touchTop

end

function get_forces!(nuc,spar,ext,simset)

    get_volume_forces!(nuc,spar);
    get_area_forces!(nuc,spar);
    get_bending_forces!(nuc,spar);
    get_elastic_forces!(nuc,spar);

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);

    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(nuc,ext[1],spar)

    elseif cmp(simset.simType,"PC") == 0
        planeRepulsion = get_plane_repulsion(nuc,ext,spar)
    end

    get_linear_chromatin_forces!(nuc,spar);
    get_bending_chromatin_forces!(nuc,spar)
    get_chromation_chromation_repulsion_forces!(nuc,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(nuc,spar,simset.envelopeTree)
    get_crosslink_forces!(nuc,spar)
    get_lad_forces!(nuc,spar)

    nuc.enve.forces.total = nuc.enve.forces.volume .+ nuc.enve.forces.area .+ nuc.enve.forces.elastic .+ nuc.enve.forces.bending .+ nuc.enve.forces.chromationRepulsion .+ nuc.enve.forces.ladEnveForces;  # .+ nuc.enve.forces.envelopeRepulsion
    
    if cmp(simset.simType,"MA") == 0 
        nuc.enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        nuc.enve.forces.total .+= micromanipulation
    elseif cmp(simset.simType,"PC") == 0

        ext[4] = 0;
        for i = findall(ext[3])
            if nuc.enve.forces.total[i][3] > 0
                ext[4] += nuc.enve.forces.total[i][3]
            end
        end

        nuc.enve.forces.total .+= planeRepulsion
    end

    if simset.simType == "VRC"
        get_repl_comp_elastic_forces!(ext,spar)
        get_repl_comp_volume_forces!(ext,spar)
        get_repl_comp_area_forces!(ext,spar)
        get_repl_comp_bending_forces!(ext,spar)        
        replCompRepulsion = get_repl_comp_chromatin_repulsion_forces!(nuc,spar,nuc.repl.tree)
        nuc.repl.forces.total = nuc.repl.forces.elastic .+ nuc.repl.forces.volume .+ nuc.repl.forces.area .+ nuc.repl.forces.bending .+ nuc.repl.forces.chromationRepulsion;
    end

    nuc.chro.forces.total = nuc.chro.forces.linear .+ nuc.chro.forces.bending .+ nuc.chro.forces.chroRepulsion .+ nuc.chro.forces.enveRepulsion .+ nuc.chro.forces.crosslink .+ nuc.chro.forces.ladChroForces

    if simset.simType == "VRC"
        nuc.chro.forces.total = nuc.chro.forces.total .+ replCompRepulsion;
    end

end

function get_crosslink_forces!(nuc,spar)

    nuc.chro.forces.crosslink = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = eachindex(nuc.chro.forces.crosslink)
        nuc.chro.forces.crosslink[i] = Vec(0.,0.,0.)
    end

    for i = 1:length(nuc.chro.crosslinks)

        vector = nuc.chro.vert[nuc.chro.crosslinks[i][1]] - nuc.chro.vert[nuc.chro.crosslinks[i][2]]
        vectorNorm = norm(vector);
        nuc.chro.forces.crosslink[nuc.chro.crosslinks[i][1]] = -spar.ladStrength*(vectorNorm - 0.4).*vector./vectorNorm;
        nuc.chro.forces.crosslink[nuc.chro.crosslinks[i][2]] = -nuc.chro.forces.crosslink[nuc.chro.crosslinks[i][1]]

    end
end

function get_repl_comp_elastic_forces!(replComp,spar)

    replComp.forces.elastic = Vector{Vec{3,Float64}}(undef, length(replComp.vert));

    for i = eachindex(replComp.vert)
        replComp.forces.elastic[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(replComp.edges)
        if replComp.firstEdges[i] == 1

            force = spar.laminaStiffness*(replComp.edgeVectorNorms[i] - replComp.normalLengths[i])*replComp.edgeUnitVectors[i];

            replComp.forces.elastic[replComp.edges[i][1]] += force
            replComp.forces.elastic[replComp.edges[i][2]] -= force
        end
    end

end

function get_repl_comp_volume_forces!(replComp,spar)

    replCompVolume = get_volume!(replComp);

    pressure = -100*spar.bulkModulus*log10(replCompVolume/replComp.normalVolume);

    replComp.forces.volume = pressure*replComp.voronoiAreas.*replComp.vertexNormalUnitVectors

end

function get_repl_comp_area_forces!(replComp,spar)

    for i = 1:length(replComp.vert)
        replComp.forces.area[i] = Vec(0.,0.,0.)
    end

    for i = eachindex(replComp.tri)

        baryocenter = mean(replComp.vert[replComp.tri[i]]);
        magnitude = 0.1*spar.areaCompressionStiffness*(replComp.triangleAreas[i] - replComp.normalTriangleAreas[i])/replComp.normalTriangleAreas[i];
                
        for j = 1:3
        
            vector = replComp.vert[replComp.tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            replComp.forces.area[replComp.tri[i][j]] += -magnitude*unitVector
        end

    end

end

function get_repl_comp_bending_forces!(replComp,spar)

    replComp.forces.bending = Vector{Vec{3,Float64}}(undef, length(replComp.vert));
    for i = 1:length(replComp.vert)
        replComp.forces.bending[i] = Vec(0.,0.,0.)
    end

    angles = get_triangle_angles(replComp);

    moment = 0.01*spar.bendingStiffness*(angles .- replComp.normalAngle);

    for i = eachindex(replComp.edges)
        
        if replComp.firstEdges[i] == 1
            
            distance1 = line_point_distance(replComp.edgeVectors[i], replComp.edgeVectors[replComp.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(replComp.edgeVectors[i], replComp.edgeVectors[replComp.edgeThirdVertices[i][2]])
            
            force1 = -moment[i]/distance1*replComp.triangleNormalUnitVectors[replComp.edgesTri[i][1]]
            force2 = -moment[i]/distance2*replComp.triangleNormalUnitVectors[replComp.edgesTri[i][2]]

            counterForces = 0.5*force1 + 0.5*force2

            replComp.forces.bending[replComp.edges3Vertex[i][1]] += force1;
            replComp.forces.bending[replComp.edges3Vertex[i][2]] += force2;
            replComp.forces.bending[replComp.edges[i][1]] -= counterForces;
            replComp.forces.bending[replComp.edges[i][2]] -= counterForces;

        end
    end
end

function get_repl_comp_chromatin_repulsion_forces!(nuc,spar,replCompTree)

    nuc.repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert));
    for i = eachindex(nuc.repl.vert)
        nuc.repl.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    replCompRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.chro.vert));
    for i = eachindex(nuc.chro.vert)
        replCompRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(replCompTree, nuc.chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(nuc.repl.vertexTri[closest[1]]);

            neighbors = nuc.repl.neighbors[closest[1]]
            neigborsCoords = nuc.repl.vert[neighbors]

            triangles = nuc.repl.tri[nuc.repl.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(nuc.repl.vert[triangles[j][k]] - nuc.chro.vert[i])
                end
            end
            
            tri = nuc.repl.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(nuc.repl, nuc.chro.vert[i], tri)

            getForce = false
            
            unitVector = (nuc.chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = nuc.repl.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) <= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (nuc.chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = nuc.repl.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(nuc.repl.edges,1) .== closeVertices[1] .&& getindex.(nuc.repl.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = nuc.repl.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = nuc.repl.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) <= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                replCompRepulsion[i] = forceMagnitude*unitVector;

                if length(closeVertices) == 1

                    nuc.repl.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - nuc.repl.vert[closeVertices[2]][1])/(nuc.repl.vert[closeVertices[2]][1] - nuc.repl.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    nuc.repl.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    nuc.repl.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(nuc.repl.vert[closeVertices[1]] - nuc.repl.vert[closeVertices[2]])*norm(nuc.repl.vert[closeVertices[1]] - nuc.repl.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - nuc.repl.vert[closeVertices[2]])*norm(closeCoords - nuc.repl.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - nuc.repl.vert[closeVertices[1]])*norm(closeCoords - nuc.repl.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - nuc.repl.vert[closeVertices[1]])*norm(closeCoords - nuc.repl.vert[closeVertices[2]])

                    nuc.repl.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    nuc.repl.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    nuc.repl.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;
                end
            end
        end
    end

    return replCompRepulsion

end

function get_forces_adh_init!(nuc,spar,ext)

    get_volume_forces!(nuc,spar);
    get_area_forces!(nuc,spar);
    get_bending_forces!(nuc,spar);
    get_elastic_forces!(nuc,spar);
    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);

    planeRepulsion, touchTop = get_plane_repulsion(nuc,ext,spar)

    nuc.enve.forces.total = nuc.enve.forces.volume .+ nuc.enve.forces.area .+ nuc.enve.forces.elastic .+ nuc.enve.forces.bending .+ nuc.enve.forces.chromationRepulsion .+ nuc.enve.forces.ladEnveForces;  # .+ nuc.enve.forces.envelopeRepulsion

    ext[4] = 0;
    for i = findall(ext[3])
        if nuc.enve.forces.total[i][3] > 0
            ext[4] += nuc.enve.forces.total[i][3]
        end
    end

    nuc.enve.forces.total .+= planeRepulsion

    return planeRepulsion, touchTop

end