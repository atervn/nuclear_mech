function get_volume_forces!(enve,spar)

    nucleusVolume = get_volume!(enve);

    pressure = -spar.bulkModulus*log10(nucleusVolume/enve.normalVolume);

    enve.forces.volume = pressure*enve.voronoiAreas.*enve.vertexNormalUnitVectors

end

function get_area_forces!(enve, spar)

    nucleusArea = sum(enve.triangleAreas);

    forceMagnitude = spar.areaCompressionStiffness*(nucleusArea - enve.normalArea)/enve.normalArea;

    for i = 1:length(enve.vert)
        enve.forces.area[i] = -forceMagnitude*mean(enve.areaUnitVectors[i])
    end

    for i = eachindex(enve.tri)

        baryocenter = mean(enve.vert[enve.tri[i]]);
        magnitude = 0.1*spar.areaCompressionStiffness*(enve.triangleAreas[i] - enve.normalTriangleAreas[i])/enve.normalTriangleAreas[i];
        
        for j = 1:3

            vector = enve.vert[enve.tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            enve.forces.area[enve.tri[i][j]] += -magnitude*unitVector
        end

    end

end

function get_bending_forces!(enve,spar)

    enve.forces.bending = Vector{Vec{3,Float64}}(undef, length(enve.vert));
    for i = 1:length(enve.vert)
        enve.forces.bending[i] = Vec(0.,0.,0.)
    end

    angles = get_triangle_angles(enve);

    moment = spar.bendingStiffness*(angles .- enve.normalAngle);

    for i = eachindex(enve.edges)
        
        if enve.firstEdges[i] == 1
            
            distance1 = line_point_distance(enve.edgeVectors[i], enve.edgeVectors[enve.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(enve.edgeVectors[i], enve.edgeVectors[enve.edgeThirdVertices[i][2]])
            
            force1 = -moment[i]/distance1*enve.triangleNormalUnitVectors[enve.edgesTri[i][1]]
            force2 = -moment[i]/distance2*enve.triangleNormalUnitVectors[enve.edgesTri[i][2]]

            counterForces = 0.5*force1 + 0.5*force2

            enve.forces.bending[enve.edges3Vertex[i][1]] += force1;
            enve.forces.bending[enve.edges3Vertex[i][2]] += force2;
            enve.forces.bending[enve.edges[i][1]] -= counterForces;
            enve.forces.bending[enve.edges[i][2]] -= counterForces;

        end
    end
end

function get_elastic_forces!(enve, spar)

    enve.forces.elastic = Vector{Vec{3,Float64}}(undef, length(enve.vert));

    for i = eachindex(enve.vert)
        enve.forces.elastic[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(enve.edges)
        if enve.firstEdges[i] == 1

            force = spar.laminaStiffness*(enve.edgeVectorNorms[i] - enve.normalLengths[i])*enve.edgeUnitVectors[i];

            enve.forces.elastic[enve.edges[i][1]] += force
            enve.forces.elastic[enve.edges[i][2]] -= force
        end
    end
end

function get_repulsion_forces!(enve,spar,envelopeTree)

    if length(enve.forces.envelopeRepulsion) == 0
        enve.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert));
    end
    
    for i = 1:length(enve.vert)

        closest,distance = knn(envelopeTree, enve.vert[i],1,true,j -> any(j .== [i; enve.neighbors[i]]))

        if distance[1] < mean(enve.normalLengths)*1.5

            nTri = length(enve.vertexTri[closest[1]]);

            closePointDistances = Vector{Float64}(undef,nTri);
            closePointCoordinates = Vector{Vec{3,Float64}}(undef, nTri);

            for j = 1:nTri

                tri = enve.vertexTri[closest[1]][j];

                closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(enve, enve.vert[i], tri)
            end
            
            closestPoint = findmin(closePointDistances);
            closestDistance = closestPoint[1];

            if closestDistance < spar.repulsionDistance

                closeCoords = closePointCoordinates[closestPoint[2]];

                unitVector = (enve.vert[i] - closeCoords)/closestDistance;
                
                forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2)
                enve.forces.envelopeRepulsion[i] = forceMagnitude*unitVector;
            else
                enve.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
            end
        else
            enve.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)
        end
    end
end

function get_aspiration_repulsion_forces(enve,pip,spar)


    repulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert));

    tree = KDTree(pip.vert);

    for i = 1:length(enve.vert)

        closest = knn(tree, enve.vert[i],1,true)

        nTri = length(pip.vertexTri[closest[1]][1]);
        closePointDistances = Vector{Float64}(undef,nTri);
        closePointCoordinates = Vector{Any}(undef,nTri);

        for j = 1:nTri

            tri = pip.vertexTri[closest[1]][1][j];

            closePointDistances[j],closePointCoordinates[j] = vertex_triangle_distance(enve, enve.vert[i], tri, pip);
        end

        closestPoint = findmin(closePointDistances);
        closestDistance = closestPoint[1];

        if closestDistance < spar.repulsionDistance

            closeCoords = closePointCoordinates[closestPoint[2]];

            unitVector = (enve.vert[i] - closeCoords)./closestDistance;
            
            forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closestDistance)^(3/2);

            repulsion[i] = forceMagnitude*unitVector;
        else
            repulsion[i] = Vec(0.,0.,0.)
        end
    end

    return repulsion

end

function get_aspiration_forces(enve,pip,spar)

    force = Vector{Vec{3,Float64}}(undef, length(enve.vert));

    for i = eachindex(enve.vert)
        if enve.vert[i][1]>0 && sqrt(enve.vert[i][2]^2 + enve.vert[i][3]^2) < 3
            forceMagnitude = 1*sum(enve.voronoiAreas[i]);
            force[i] = forceMagnitude*enve.vertexNormalUnitVectors[i]
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

function get_random_enve_fluctuations(spar,enve,dt)

    strength = sqrt(0*spar.boltzmannConst*spar.temperature/(dt));

    fluctuationForces = Vector{Vec{3,Float64}}(undef,length(enve.vert))
    for i = 1:length(enve.vert)
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

function get_envelope_chromatin_repulsion_forces!(enve,chro,spar,envelopeTree)

    enve.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert));
    for i = eachindex(enve.vert)
        enve.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    for i = eachindex(chro.vert)
        chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(envelopeTree, chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(enve.vertexTri[closest[1]]);

            neighbors = enve.neighbors[closest[1]]
            neigborsCoords = enve.vert[neighbors]

            triangles = enve.tri[enve.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(enve.vert[triangles[j][k]] - chro.vert[i])
                end
            end
            
            tri = enve.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(enve, chro.vert[i], tri)

            getForce = false
            
            unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = enve.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) >= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = enve.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(enve.edges,1) .== closeVertices[1] .&& getindex.(enve.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = enve.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = enve.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) >= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                chro.forces.enveRepulsion[i] = forceMagnitude*unitVector;

                if length(closeVertices) == 1

                    enve.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - enve.vert[closeVertices[2]][1])/(enve.vert[closeVertices[2]][1] - enve.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    enve.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    enve.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(enve.vert[closeVertices[1]] - enve.vert[closeVertices[2]])*norm(enve.vert[closeVertices[1]] - enve.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - enve.vert[closeVertices[2]])*norm(closeCoords - enve.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - enve.vert[closeVertices[1]])*norm(closeCoords - enve.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - enve.vert[closeVertices[1]])*norm(closeCoords - enve.vert[closeVertices[2]])

                    enve.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    enve.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    enve.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;

                    closeCoords
                end
            end
        end
    end
end

function get_micromanipulation_forces(enve,mm,spar)

    micromanipulation = Vector{Vec{3,Float64}}(undef, length(enve.vert));
    for i = eachindex(enve.vert)
        micromanipulation[i] = Vec(0.,0.,0.);
    end


    micromanipulation[mm.leftmostVertex] = 2*spar.pullingForce*(mm.leftmostVertexPosition .- enve.vert[mm.leftmostVertex])
    # micromanipulation[mm.leftNeighbors] = spar.pullingForce*(mm.leftNeigborPositions .- enve.vert[mm.leftNeighbors])

    forceVector = spar.pullingForce*Vec(1.,0.,0.);
    micromanipulation[mm.rightmostVertex] = forceVector
    # for i = eachindex(mm.rightNeighbors)
    #     micromanipulation[mm.rightNeighbors[i]] = 0.5*1/length(mm.rightNeighbors)*forceVector
    # end
    

    return micromanipulation
end

function get_lad_forces!(enve,chro,spar)

    enve.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(enve.vert));

    for i = 1:length(enve.vert)
        enve.forces.ladEnveForces[i] = Vec(0.,0.,0.) 
    end

    for i = 1:spar.chromatinNumber

        for j = 1:length(chro.lads[i])
            vector = enve.vert[enve.lads[i][j]] - chro.vert[chro.strandIdx[i][chro.lads[i][j]]]
            distance = norm(vector)
            magnitude = -spar.ladStrength*(distance - 0.2)
            enve.forces.ladEnveForces[enve.lads[i][j]] = magnitude*vector/distance
            chro.forces.ladChroForces[chro.strandIdx[i][chro.lads[i][j]]] = -magnitude*vector/distance

        end
    end

end

function get_plane_repulsion(enve,simset,spar)

    enve.forces.planeRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert));
    for i = eachindex(enve.vert)
        enve.forces.planeRepulsion[i] = Vec(0.,0.,0.);
    end

    simset.adh.touchingTop = zeros(Bool,length(enve.vert))

    for i = eachindex(enve.vert)
        if simset.adh.topPlane - enve.vert[i][3] < spar.repulsionDistance
            
            simset.adh.touchingTop[i] = true

            if simset.adh.topPlane - enve.vert[i][3] < 0
                enve.forces.planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*0.1)
            else
                distance = simset.adh.topPlane - enve.vert[i][3];
                enve.forces.planeRepulsion[i] = Vec(0.,0.,-spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
            end
        end
    end

    for i = eachindex(enve.vert)
        if enve.vert[i][3] - simset.adh.bottomPlane < spar.repulsionDistance
            if enve.vert[i][3] - simset.adh.bottomPlane < 0
                enve.forces.planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*0.1)
            else
                distance = enve.vert[i][3] - simset.adh.bottomPlane;
                enve.forces.planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
            end
        end
    end
end

function get_forces!(enve, chro , spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);
    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)
    end

    get_linear_chromatin_forces!(chro,spar);
    get_bending_chromatin_forces!(chro,spar)
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(enve,chro,spar,simset.envelopeTree)
    get_crosslink_forces!(chro,spar)
    get_lad_forces!(enve,chro,spar)

    if simset.adh.adherent
        get_plane_repulsion(enve,simset,spar)
    end

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion
    
    if simset.adh.adherent
        enve.forces.total .+= enve.forces.planeRepulsion;
    end

    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    end

    chro.forces.total = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ chro.forces.enveRepulsion .+ chro.forces.crosslink .+ chro.forces.ladChroForces
end

function get_forces!(enve, chro, repl, spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);

    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)

    elseif cmp(simset.simType,"PC") == 0
        planeRepulsion = get_plane_repulsion(nuc,ext,spar)
    end

    get_linear_chromatin_forces!(chro,spar);
    get_bending_chromatin_forces!(chro,spar)
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(enve,chro,spar,simset.envelopeTree)
    get_crosslink_forces!(chro,spar)
    get_lad_forces!(enve,chro,spar)

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion
    
    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    elseif cmp(simset.simType,"PC") == 0

        ext[4] = 0;
        for i = findall(ext[3])
            if enve.forces.total[i][3] > 0
                ext[4] += enve.forces.total[i][3]
            end
        end

        enve.forces.total .+= planeRepulsion
    end

    get_repl_comp_elastic_forces!(ext,spar)
    get_repl_comp_volume_forces!(ext,spar)
    get_repl_comp_area_forces!(ext,spar)
    get_repl_comp_bending_forces!(ext,spar)
    replCompRepulsion = get_repl_comp_chromatin_repulsion_forces!(nuc,spar,repl.tree)
    repl.forces.total = repl.forces.elastic .+ repl.forces.volume .+ repl.forces.area .+ repl.forces.bending .+ repl.forces.chromationRepulsion;

    chro.forces.total = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ chro.forces.enveRepulsion .+ chro.forces.crosslink .+ chro.forces.ladChroForces .+ replCompRepulsion
end

function get_forces!(enve, spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);
    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)
    end

    if simset.adh.adherent
        get_plane_repulsion(enve,simset,spar)
    end

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion
    
    if simset.adh.adherent
        enve.forces.total .+= enve.forces.planeRepulsion;
    end

    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    end

end

function get_crosslink_forces!(chro,spar)

    chro.forces.crosslink = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
    for i = eachindex(chro.forces.crosslink)
        chro.forces.crosslink[i] = Vec(0.,0.,0.)
    end

    for i = 1:length(chro.crosslinks)

        vector = chro.vert[chro.crosslinks[i][1]] - chro.vert[chro.crosslinks[i][2]]
        vectorNorm = norm(vector);
        chro.forces.crosslink[chro.crosslinks[i][1]] = -spar.ladStrength*(vectorNorm - 0.4).*vector./vectorNorm;
        chro.forces.crosslink[chro.crosslinks[i][2]] = -chro.forces.crosslink[chro.crosslinks[i][1]]

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

function get_repl_comp_chromatin_repulsion_forces!(chro,spar,replCompTree)

    repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert));
    for i = eachindex(repl.vert)
        repl.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    replCompRepulsion = Vector{Vec{3,Float64}}(undef, length(chro.vert));
    for i = eachindex(chro.vert)
        replCompRepulsion[i] = Vec(0.,0.,0.)
    end


    for i = 1:spar.chromatinLength*spar.chromatinNumber
        closest,distance = knn(replCompTree, chro.vert[i],1,true)

        if distance[1] < spar.meanLaminaLength*2

            nTri = length(repl.vertexTri[closest[1]]);

            neighbors = repl.neighbors[closest[1]]
            neigborsCoords = repl.vert[neighbors]

            triangles = repl.tri[repl.vertexTri[closest[1]]]
            distanceVector = zeros(size(triangles,1))
            for j = 1:nTri
                for k = 1:3
                    distanceVector[j] += norm(repl.vert[triangles[j][k]] - chro.vert[i])
                end
            end
            
            tri = repl.vertexTri[closest[1]][argmin(distanceVector)]

            closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(repl, chro.vert[i], tri)

            getForce = false
            
            unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
            if closePointDistance < spar.repulsionDistance
                getForce = true
                
            else
                enveUnitVector = repl.vertexNormalUnitVectors[closest[1]];
                if dot(unitVector,enveUnitVector) <= 0
                    getForce = true
                end
            end

            if getForce

                # unitVector = (chro.vert[i] - closeCoords)/closePointDistance;
                
                if length(closeVertices) == 1
                    enveUnitVector = repl.vertexNormalUnitVectors[closeVertices[1]];
                elseif length(closeVertices) == 2
                    edge = findall((getindex.(repl.edges,1) .== closeVertices[1] .&& getindex.(repl.edges,2) .== closeVertices[2]))[1]
                    enveUnitVector = repl.edgeNormalUnitVectors[edge]
                else
                    enveUnitVector = repl.triangleNormalUnitVectors[tri]
                end

                if dot(unitVector,enveUnitVector) <= 0
                    unitVector = -unitVector
                    forceMagnitude = 0.5*spar.repulsionConstant
                else
                    forceMagnitude = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)
                end

                replCompRepulsion[i] = forceMagnitude*unitVector;

                if length(closeVertices) == 1

                    repl.forces.chromationRepulsion[closeVertices[1]] += -forceMagnitude*unitVector;

                elseif length(closeVertices) == 2

                    w1 = (closeCoords[1] - repl.vert[closeVertices[2]][1])/(repl.vert[closeVertices[2]][1] - repl.vert[closeVertices[2]][2])
                    w2 = 1 - w1;

                    repl.forces.chromationRepulsion[closeVertices[1]] += -w1*forceMagnitude*unitVector;
                    repl.forces.chromationRepulsion[closeVertices[2]] += -w2*forceMagnitude*unitVector;

                else

                    # get baryocentric weights
                    # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
                    # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

                    fullArea = 0.5*norm(repl.vert[closeVertices[1]] - repl.vert[closeVertices[2]])*norm(repl.vert[closeVertices[1]] - repl.vert[closeVertices[3]])
                    area1 = 0.5*norm(closeCoords - repl.vert[closeVertices[2]])*norm(closeCoords - repl.vert[closeVertices[3]])
                    area2 = 0.5*norm(closeCoords - repl.vert[closeVertices[1]])*norm(closeCoords - repl.vert[closeVertices[3]])
                    area3 = 0.5*norm(closeCoords - repl.vert[closeVertices[1]])*norm(closeCoords - repl.vert[closeVertices[2]])

                    repl.forces.chromationRepulsion[closeVertices[1]] += -area1/fullArea*forceMagnitude*unitVector;
                    repl.forces.chromationRepulsion[closeVertices[2]] += -area2/fullArea*forceMagnitude*unitVector;
                    repl.forces.chromationRepulsion[closeVertices[3]] += -area3/fullArea*forceMagnitude*unitVector;
                end
            end
        end
    end

    return replCompRepulsion

end

function get_forces_adh_init!(enve,spar,ext)

    get_volume_forces!(nuc,spar);
    get_area_forces!(nuc,spar);
    get_bending_forces!(nuc,spar);
    get_elastic_forces!(nuc,spar);
    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);

    planeRepulsion, touchTop = get_plane_repulsion(nuc,ext,spar)

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion

    ext[4] = 0;
    for i = findall(ext[3])
        if enve.forces.total[i][3] > 0
            ext[4] += enve.forces.total[i][3]
        end
    end

    enve.forces.total .+= planeRepulsion

    return planeRepulsion, touchTop

end