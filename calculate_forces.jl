function get_volume_forces!(nuc,spar)

    nucleusVolume = get_volume!(nuc);

    pressure = -spar.bulkModulus*log10(nucleusVolume/nuc.normalVolume);

    nuc.forces.volume = pressure*nuc.voronoiAreas.*nuc.vertexNormalUnitVectors

end

function get_area_forces!(nuc,spar)

    nucleusArea = sum(nuc.triangleAreas);

    forceMagnitude = spar.areaCompressionStiffness*(nucleusArea - nuc.normalArea)/nuc.normalArea;

    for i = 1:length(nuc.vert)
        nuc.forces.area[i] = -forceMagnitude*mean(nuc.areaUnitVectors[i])
    end

    for i = eachindex(nuc.tri)

        baryocenter = mean(nuc.vert[nuc.tri[i]]);
        magnitude = 0.1*spar.areaCompressionStiffness*(nuc.triangleAreas[i] - nuc.normalTriangleAreas[i])/nuc.normalTriangleAreas[i];
        
        for j = 1:3

            vector = nuc.vert[nuc.tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            nuc.forces.area[nuc.tri[i][j]] += -magnitude*unitVector
        end

    end

end

function get_bending_forces!(nuc,spar)

    nuc.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = 1:length(nuc.vert)
        nuc.forces.bending[i] = Vec(0.,0.,0.)
    end

    angles = get_triangle_angles(nuc);

    moment = spar.bendingStiffness*(angles .- nuc.normalAngle);

    for i = eachindex(nuc.edges)
        
        if nuc.firstEdges[i] == 1
            
            distance1 = line_point_distance(nuc.edgeVectors[i], nuc.edgeVectors[nuc.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(nuc.edgeVectors[i], nuc.edgeVectors[nuc.edgeThirdVertices[i][2]])
            
            force1 = -moment[i]/distance1*nuc.triangleNormalUnitVectors[nuc.edgesTri[i][1]]
            force2 = -moment[i]/distance2*nuc.triangleNormalUnitVectors[nuc.edgesTri[i][2]]

            counterForces = 0.5*force1 + 0.5*force2

            nuc.forces.bending[nuc.edges3Vertex[i][1]] += force1;
            nuc.forces.bending[nuc.edges3Vertex[i][2]] += force2;
            nuc.forces.bending[nuc.edges[i][1]] -= counterForces;
            nuc.forces.bending[nuc.edges[i][2]] -= counterForces;

        end
    end
end

function get_elastic_forces!(nuc,spar)

    nuc.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    for i = eachindex(nuc.vert)
        nuc.forces.elastic[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(nuc.edges)
        if nuc.firstEdges[i] == 1

            force = spar.laminaStiffness*(nuc.edgeVectorNorms[i] - nuc.normalLengths[i])*nuc.edgeUnitVectors[i];

            nuc.forces.elastic[nuc.edges[i][1]] += force
            nuc.forces.elastic[nuc.edges[i][2]] -= force
        end
    end
end

function get_repulsion_forces!(nuc,spar,envelopeTree)

    if length(nuc.forces.envelopeRepulsion) == 0
        nuc.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    end
    
    for i = 1:length(nuc.vert)

        closest,distance = knn(envelopeTree, nuc.vert[i],1,true,j -> any(j .== [i; nuc.neighbors[i]]))

        if distance[1] < mean(nuc.normalLengths)*1.5

            nTri = length(nuc.vertexTri[closest[1]]);

            closePointDistances = Vector{Float64}(undef,nTri);
            closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);

            for j = 1:nTri

                tri = nuc.vertexTri[closest[1]][j];

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

    for i = 1:length(chro.vert)
        chro.forces.linear[i] = Vec(0.,0.,0.)
    end

    for i = 1:spar.chromatinNumber
        chro.forces.strandLinear[i][1:end-1] += spar.chromatinStiffness*(chro.vectorNorms[i] .- spar.chroVertexDistance).*chro.vectors[i]./chro.vectorNorms[i];
        chro.forces.strandLinear[i][2:end] -= chro.forces.strandLinear[i][1:end-1];
    end
end

function get_bending_chromatin_forces!(chro,spar)
    
    for i = 1:spar.chromatinNumber

        vectors1 = -chro.vectors[i][1:end-1]
        vectors2 = chro.vectors[i][2:end]

        angles = dot.(vectors1, vectors2)./(chro.vectorNorms[i][1:end-1].*chro.vectorNorms[i][2:end]);

        angles[angles .>= 1] .= 1.; 

        angles = acos.(angles);
        
        unitVectors1 = cross.(vectors1,cross.(vectors1,vectors2));
        unitVectors2 = cross.(vectors2,cross.(vectors2,vectors1));

        unitVectors1 = unitVectors1./norm.(unitVectors1);
        unitVectors2 = unitVectors2./norm.(unitVectors2);
        
        chro.forces.strandBending[i][1:end-2] = -spar.chromatinBendingModulus./chro.vectorNorms[i][1:end-1].*(angles .- spar.chromatinNormalAngle).*unitVectors1;
        chro.forces.strandBending[i][3:end] = -spar.chromatinBendingModulus./chro.vectorNorms[i][2:end].*(angles .- spar.chromatinNormalAngle).*unitVectors2;
    end
end

function get_chromation_chromation_repulsion_forces!(chro,spar,chromatinTree)

    for i = 1:length(chro.vert)
        chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
    end


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

    strength = sqrt(2*spar.boltzmannConst*spar.temperature/(dt));

    fluctuationForces = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = 1:spar.chromatinLength*spar.chromatinNumber
        fluctuationForces[i] = strength.*Vec(randn(),randn(),randn())
    end

    return fluctuationForces

end

function get_random_enve_fluctuations(spar,nuc,dt)

    strength = sqrt(0*spar.boltzmannConst*spar.temperature/(dt));

    fluctuationForces = Vector{Vec{3,Float64}}(undef,length(nuc.vert))
    for i = 1:length(nuc.vert)
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
        chro.forces.ladChroForces[i] = Vec(0.,0.,0.)

    end
end

function get_envelope_chromatin_repulsion_forces!(nuc,chro,spar,envelopeTree)

    nuc.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = eachindex(nuc.vert)
        nuc.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(chro.vert)
        chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(envelopeTree, chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(nuc.vertexTri[closest[1]]);

            neighbors = nuc.neighbors[closest[1]]
            neigborsCoords = nuc.vert[neighbors]

            triangles = nuc.tri[nuc.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(nuc.vert[triangles[j][k]] - chro.vert[i])
                end
            end
            
            tri = nuc.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(nuc, chro.vert[i], tri)

            getForce = false
            
            unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = nuc.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) >= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = nuc.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(nuc.edges,1) .== closeVertices[1] .&& getindex.(nuc.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = nuc.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = nuc.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) >= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                chro.forces.enveRepulsion[i] = forceMagnitude*unitVector;

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
            end
        end
    end
end

function get_micromanipulation_forces(nuc,mm,spar)

    micromanipulation = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = eachindex(nuc.vert)
        micromanipulation[i] = Vec(0.,0.,0.);
    end


    micromanipulation[mm.leftmostVertex] = 2*spar.pullingForce*(mm.leftmostVertexPosition .- nuc.vert[mm.leftmostVertex])
    # micromanipulation[mm.leftNeighbors] = spar.pullingForce*(mm.leftNeigborPositions .- nuc.vert[mm.leftNeighbors])

    forceVector = spar.pullingForce*Vec(1.,0.,0.);
    micromanipulation[mm.rightmostVertex] = forceVector
    # for i = eachindex(mm.rightNeighbors)
    #     micromanipulation[mm.rightNeighbors[i]] = 0.5*1/length(mm.rightNeighbors)*forceVector
    # end
    

    return micromanipulation
end

function get_lad_forces!(nuc,chro,spar)

    nuc.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    for i = 1:length(nuc.vert)
        nuc.forces.ladEnveForces[i] = Vec(0.,0.,0.) 
    end

    for i = 1:spar.chromatinNumber

        for j = 1:length(chro.lads[i])
            vector = nuc.vert[nuc.lads[i][j]] - chro.vert[chro.strandIdx[i][chro.lads[i][j]]]
            distance = norm(vector)
            magnitude = -spar.ladStrength*(distance - 0.2)
            nuc.forces.ladEnveForces[nuc.lads[i][j]] = magnitude*vector/distance
            chro.forces.ladChroForces[chro.strandIdx[i][chro.lads[i][j]]] = -magnitude*vector/distance

        end
    end

end

function get_plane_repulsion(nuc,ext,spar)

    planeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    for i = eachindex(nuc.vert)
        planeRepulsion[i] = Vec(0.,0.,0.);
    end

    touchTop = zeros(Bool,length(nuc.vert))

    for i = eachindex(nuc.vert)
        if ext[1] - nuc.vert[i][3] < spar.repulsionDistance
            
            touchTop[i] = true

            if ext[1] - nuc.vert[i][3] < 0
                planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*0.1)
            else
                distance = ext[1] - nuc.vert[i][3];
                planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
                ext[3][i] = true
            end
        end
    end

    for i = eachindex(nuc.vert)
        if nuc.vert[i][3] - ext[2] < spar.repulsionDistance
            if nuc.vert[i][3] - ext[2] < 0
                planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*0.1)
            else
                distance = nuc.vert[i][3] - ext[2];
                planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
            end
        end
    end

    return planeRepulsion, touchTop

end

function get_forces!(nuc,chro,spar,ext,simset)

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
    get_linear_chromatin_forces!(chro,spar);
    get_bending_chromatin_forces!(chro,spar)
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(nuc,chro,spar,simset.envelopeTree)
    crosslinkForces = get_crosslink_forces!(chro,spar,simset)
    get_lad_forces!(nuc,chro,spar)


    nuc.forces.total = nuc.forces.volume .+ nuc.forces.area .+ nuc.forces.elastic .+ nuc.forces.bending .+ nuc.forces.chromationRepulsion .+ nuc.forces.ladEnveForces;  # .+ nuc.forces.envelopeRepulsion
    
    if cmp(simset.simType,"MA") == 0 
        nuc.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        nuc.forces.total .+= micromanipulation
    elseif cmp(simset.simType,"PC") == 0

        ext[4] = 0;
        for i = findall(ext[3])
            if nuc.forces.total[i][3] > 0
                ext[4] += nuc.forces.total[i][3]
            end
        end

        nuc.forces.total .+= planeRepulsion
    end

    if simset.simType == "VRC"
        get_repl_comp_elastic_forces!(ext,spar)
        get_repl_comp_volume_forces!(ext,spar)
        get_repl_comp_area_forces!(ext,spar)
        get_repl_comp_bending_forces!(ext,spar)        
        replCompRepulsion = get_repl_comp_chromatin_repulsion_forces!(ext,chro,spar,ext.tree)
        ext.forces.total = ext.forces.elastic .+ ext.forces.volume .+ ext.forces.area .+ ext.forces.bending .+ ext.forces.chromationRepulsion;
    end

    chro.forces.total = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ chro.forces.enveRepulsion .+ crosslinkForces .+ chro.forces.ladChroForces

    if simset.simType == "VRC"
        chro.forces.total = chro.forces.total .+ replCompRepulsion;
    end


end

function get_crosslink_forces!(chro,spar,simset)

    crosslinkForces = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = eachindex(crosslinkForces)
        crosslinkForces[i] = Vec(0.,0.,0.)
    end

    for i = 1:length(chro.crosslinks)

        vector = chro.vert[chro.crosslinks[i][1]] - chro.vert[chro.crosslinks[i][2]]
        vectorNorm = norm(vector);
        crosslinkForces[chro.crosslinks[i][1]] = -spar.ladStrength*(vectorNorm - 0.4).*vector./vectorNorm;
        crosslinkForces[chro.crosslinks[i][2]] = -crosslinkForces[chro.crosslinks[i][1]]

    end

    return crosslinkForces
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

function get_repl_comp_chromatin_repulsion_forces!(ext,chro,spar,replCompTree)

    ext.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(ext.vert));
    for i = eachindex(ext.vert)
        ext.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    replCompRepulsion = Vector{Vec{3,Float64}}(undef, length(chro.vert));
    for i = eachindex(chro.vert)
        replCompRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(replCompTree, chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(ext.vertexTri[closest[1]]);

            neighbors = ext.neighbors[closest[1]]
            neigborsCoords = ext.vert[neighbors]

            triangles = ext.tri[ext.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(ext.vert[triangles[j][k]] - chro.vert[i])
                end
            end
            
            tri = ext.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(ext, chro.vert[i], tri)

            getForce = false
            
            unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = ext.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) <= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = ext.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(ext.edges,1) .== closeVertices[1] .&& getindex.(ext.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = ext.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = ext.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) <= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                replCompRepulsion[i] = forceMagnitude*unitVector;

                if length(closeVertices) == 1

                    ext.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - ext.vert[closeVertices[2]][1])/(ext.vert[closeVertices[2]][1] - ext.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    ext.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    ext.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(ext.vert[closeVertices[1]] - ext.vert[closeVertices[2]])*norm(ext.vert[closeVertices[1]] - ext.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - ext.vert[closeVertices[2]])*norm(closeCoords - ext.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - ext.vert[closeVertices[1]])*norm(closeCoords - ext.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - ext.vert[closeVertices[1]])*norm(closeCoords - ext.vert[closeVertices[2]])

                    ext.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    ext.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    ext.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;
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

    nuc.forces.total = nuc.forces.volume .+ nuc.forces.area .+ nuc.forces.elastic .+ nuc.forces.bending .+ nuc.forces.chromationRepulsion .+ nuc.forces.ladEnveForces;  # .+ nuc.forces.envelopeRepulsion

    ext[4] = 0;
    for i = findall(ext[3])
        if nuc.forces.total[i][3] > 0
            ext[4] += nuc.forces.total[i][3]
        end
    end

    nuc.forces.total .+= planeRepulsion

    return planeRepulsion, touchTop

end