function get_volume_forces!(enve::envelopeType,spar::scaledParametersType,simset::simulationSettingsType)

    # compute the volume of the nucleus
    nucleusVolume = get_volume!(enve);

    # if the nucleus is smaller than the normal volume
    if nucleusVolume < enve.normalVolume || simset.newVolumeSimulation

        # compute the volume pressure based on the difference between nucleus volume and normal volume
        volPressure = -spar.bulkModulus*log10(nucleusVolume/(enve.normalVolume));

    # otherwise, set to 0
    else
        volPressure = 0
    end

    # compute the volume forces on envelope vertices
    enve.forces.volume = (spar.osmoticPressure + volPressure)*enve.voronoiAreas.*enve.vertexNormalUnitVectors

end

function get_volume_forces!(enve,repl::replicationCompartmentType,spar,simset::simulationSettingsType)

    # compute the volume as the different between the nuclear volume the replication compartment volume
    nucleusVolume = get_volume!(enve);# - get_volume!(repl);

    # if the nucleus is smaller than the normal volume
    if nucleusVolume < enve.normalVolume  #|| simset.newVolumeSimulation # - repl.normalVolume

        # compute the volume pressure based on the difference between nucleus volume and normal volume
        volPressure = 0*-spar.bulkModulus*log10(nucleusVolume/(enve.normalVolume));# - repl.normalVolume));

    # otherwise, set to 0
    else

        volPressure = 0

    end

    # compute the volume forces on envelope vertices
    enve.forces.volume = (spar.osmoticPressure + volPressure)*enve.voronoiAreas.*enve.vertexNormalUnitVectors

end

function get_area_forces!(enve, spar)

    for i = eachindex(enve.vert)

        enve.forces.area[i] = Vec(0.,0.,0.)

    end

    nucleusArea = sum(enve.triangleAreas);

    globalMagnitude = 10e-4/spar.viscosity*spar.scalingTime*(nucleusArea - enve.normalArea)/enve.normalArea;

    # compute area forces for each triangle
    for i = eachindex(enve.tri)

        # get the baryocenter
        centroid = mean(enve.vert[enve.tri[i]]);

        vectors = Vector{Any}(undef, 3)

        vectors[1] = centroid - enve.vert[enve.tri[i][1]]
        vectors[2] = centroid - enve.vert[enve.tri[i][2]]
        vectors[3] = centroid - enve.vert[enve.tri[i][3]]

        vectorSum = norm(vectors[1])^2 + norm(vectors[2])^2 + norm(vectors[3])^2

        # get the force magnitude
        localMagnitude = spar.areaCompressionStiffness*(enve.triangleAreas[i] - enve.normalTriangleAreas[i])/vectorSum

        globalMagnitudeTemp = globalMagnitude*enve.triangleAreas[i]/vectorSum

        # for each triangle vertex
        for j = 1:3

            # calculate the force on the vertex
            enve.forces.area[enve.tri[i][j]] += localMagnitude*vectors[j] + globalMagnitudeTemp*vectors[j]

        end
    end
end

function get_bending_forces!(enve,spar)

    # set the forces to zero
    for i = 1:length(enve.vert)
        enve.forces.bending[i] = Vec(0.,0.,0.)
    end

    # compute angles of triangles in the envelope
    angles = get_triangle_angles(enve);

    # compute bending moments based on difference between angles and normal angles
    moments = spar.laminaBendingStiffness*enve.laminaDisintegrationMultipliers.*(angles .- enve.normalAngle);

    # for each edge
    for i = eachindex(enve.edges)
        
        # if this is a first edge
        if enve.firstEdges[i] == 1
            
            # compute distances between current edge and third vertex of the adjacent triangles
            distance1 = line_point_distance(enve.edgeVectors[i], enve.edgeVectors[enve.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(enve.edgeVectors[i], enve.edgeVectors[enve.edgeThirdVertices[i][2]])
            
            # compute bending forces for adjacent vertices
            force1 = moments[i]/distance1*enve.triangleNormalUnitVectors[enve.edgesTri[i][1]]
            force2 = moments[i]/distance2*enve.triangleNormalUnitVectors[enve.edgesTri[i][2]]

            # compute counter forces for edge vertices
            counterForces = 0.5*force1 + 0.5*force2

            # add to the bending forces for vertices
            enve.forces.bending[enve.edges3Vertex[i][1]] += force1;
            enve.forces.bending[enve.edges3Vertex[i][2]] += force2;
            enve.forces.bending[enve.edges[i][1]] -= counterForces;
            enve.forces.bending[enve.edges[i][2]] -= counterForces;

        end
    end
end

function get_elastic_forces!(enve, spar)

    # set the forces to zero
    for i = eachindex(enve.vert)
        enve.forces.elastic[i] = Vec(0.,0.,0.);
    end

    # compute elastic forces for each edge
    for i = eachindex(enve.edges)

        # if first edge
        if enve.firstEdges[i] == 1

            # compute elastic force based on stiffness, edge vector norm difference, and edge unit vector
            force = enve.laminaDisintegrationMultipliers[i]*enve.envelopeMultipliers[i]*spar.laminaStiffness*(enve.edgeVectorNorms[i] - enve.normalLengths[i])*enve.edgeUnitVectors[i];

            # update elastic forces for vertices
            enve.forces.elastic[enve.edges[i][1]] += force
            enve.forces.elastic[enve.edges[i][2]] -= force

        end
    end
end

function get_repulsion_forces!(enve,spar,simset)

    # compute repulsion forces for each vertex
    for i = 1:length(enve.vert)

        # get the intra envelope interactions
        getForce, unitVector, closePointDistance, _ = vertex_shell_interaction(enve.vert[i],enve,simset.envelopeTree,spar,[i; enve.neighbors[i]],false,false,spar.repulsionDistance)

        # check if the closest distance is within the repulsion distance
        if getForce
            
            # calculate the force
            enve.forces.envelopeRepulsion[i] = spar.repulsionConstant*(spar.repulsionDistance - closePointDistance)*unitVector
        
        # if farther then the repulsion distance, set force to 0
        else

            # assign zero force
            enve.forces.envelopeRepulsion[i] = Vec(0.,0.,0.)

        end
    end
end

function get_aspiration_repulsion_forces!(enve,pip,spar)

    # compute aspiration repulsion forces for each envelope vertex
    for i = 1:length(enve.vert)

        # get the interaction with the pipette
        getForce, unitVector, closePointDistance, closeCoords, closeVertices,tri = vertex_shell_interaction(enve.vert[i],pip,pip.pipetteTree,spar,[],true,false,spar.repulsionDistance)

        # if the repulsion is calculated
        if getForce

            # get the interaction force with the pipette
            unitVector, forceMagnitude,acrossBoundary = get_vertex_shell_interaction_vertex_force(pip,unitVector,closePointDistance, closeVertices, tri, spar,false,spar.repulsionDistance)

            # calculate the force
            enve.forces.pipetteRepulsion[i] = forceMagnitude*unitVector;

        else

            # assign zero force
            enve.forces.pipetteRepulsion[i] = Vec(0.,0.,0.)

        end
    end
end

function get_aspiration_forces!(enve,spar)

    # init variable to sum all the surface normal x-components for the vertices within the pipette
    xComponents = 0

    # init a vector to collect the enve vertices within the pipette
    inPipette = []

    # compute total x-components and identify vertices in the pipette
    for i = eachindex(enve.vert)

        # check if within the pipette
        if enve.vert[i][1]>0 && sqrt(enve.vert[i][2]^2 + enve.vert[i][3]^2) < spar.pipetteRadius

            # add the vertex to the suction region
            push!(inPipette,i)

            # add the xComponent
            xComponents += enve.vertexNormalUnitVectors[i][1]*enve.voronoiAreas[i]

        end
    end

    # get the total force
    totalForce = spar.aspirationPressure*pi*spar.pipetteRadius^2

    # get the multiplier based on the xComponents of the vertices within the pipette
    multiplier = totalForce/xComponents;

    # iterate through the vertices
    for i = eachindex(enve.vert)

        # if within the pipette
        if any(i .== inPipette)

            # calculate the force
            enve.forces.aspiration[i] = multiplier*enve.vertexNormalUnitVectors[i]*enve.voronoiAreas[i]

        # outside the pipette
        else

            # set force to 0
            enve.forces.aspiration[i] = Vec(0.,0.,0.)

        end
    end
end

function get_linear_chromatin_forces!(chro,spar)

    # initialize linear chromatin forces array
    for i = 1:length(chro.vert)
        chro.forces.linear[i] = Vec(0.,0.,0.)
    end

    # iterate through the chromosomes
    for i = 1:spar.chromatinNumber

        # compute linear forces for each strand
        chro.forces.strandLinear[i][1:end-1] += spar.chromatinStiffness*(chro.vectorNorms[i] .- spar.chroVertexDistance).*chro.vectors[i]./chro.vectorNorms[i];
        chro.forces.strandLinear[i][2:end] -= chro.forces.strandLinear[i][1:end-1];
    end
end

function get_bending_chromatin_forces!(chro,spar)
    
    # iterate through the chromosomes
    for i = 1:spar.chromatinNumber

        # compute vectors for bending forces calculation
        vectors1 = -chro.vectors[i][1:end-1]
        vectors2 = chro.vectors[i][2:end]

        # calculate angles between adjacent vectors
        angles = dot.(vectors1, vectors2)./(chro.vectorNorms[i][1:end-1].*chro.vectorNorms[i][2:end]);

         # limit the angles to the range [-1, 1] to avoid NaN errors in acos
        angles[angles .>= 1] .= 1.; 

        # calculate actual angles using arccosine
        angles = acos.(angles);
        
        # compute unit vectors for bending forces
        unitVectors1 = cross.(vectors1,cross.(vectors1,vectors2));
        unitVectors2 = cross.(vectors2,cross.(vectors2,vectors1));

        # normalize unit vectors
        unitVectors1 = unitVectors1./norm.(unitVectors1);
        unitVectors2 = unitVectors2./norm.(unitVectors2);
        
        # calculate bending forces
        chro.forces.strandBending[i][1:end-2] = -spar.chromatinBendingModulus./chro.vectorNorms[i][1:end-1].*(angles .- spar.chromatinNormalAngle).*unitVectors1;
        chro.forces.strandBending[i][3:end] = -spar.chromatinBendingModulus./chro.vectorNorms[i][2:end].*(angles .- spar.chromatinNormalAngle).*unitVectors2;

    end
end

function get_chromation_chromation_repulsion_forces!(chro,spar,chromatinTree)

    # initialize repulsion forces to zero
    for i = 1:length(chro.vert)
        chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
    end

    # calculate repulsion forces between chromatin vertices
    for i = 1:spar.chromatinNumber

        # find close vertices within repulsion distance for the vertices in chromosome i
        closeVertices = inrange(chromatinTree, chro.vert[chro.strandIdx[i]], spar.repulsionDistance, false)

        # compute repulsion forces for each vertex in the strand
        for j = 1:spar.chromatinLength

            # iterate through the close vertices for vertex j
            for k = eachindex(closeVertices[j])

                # if this is not the vertex itself
                if chro.strandIdx[i][j] != closeVertices[j][k]

                    # get the vector adn vector norm
                    vector = chro.strandVert[i][j] - chro.vert[closeVertices[j][k]]
                    vectorNorm = norm(vector)

                    # add to the repulsion force
                    chro.forces.strandChroRepulsion[i][j] += -spar.repulsionConstant*(vectorNorm - spar.repulsionDistance)*vector/vectorNorm

                end
            end
        end
    end
end

function get_random_fluctuations(spar,simset,nVertices)

    # calculate the strength of random fluctuations
    strength = sqrt(2*spar.boltzmannConst*spar.temperature/(spar.dt*simset.timeStepMultiplier));

    # init the forces
    fluctuations = Vector{Vec{3,Float64}}(undef,nVertices)

    # generate random fluctuation forces
    for i = 1:nVertices
        fluctuations[i] = strength.*Vec(randn(),randn(),randn())
    end

    return fluctuations

end

function get_envelope_chromatin_repulsion_forces!(enve,chro,spar,envelopeTree)

    # initialize envelope repulsion forces
    for i = eachindex(enve.vert)
        enve.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    # initialize chromatin repulsion force
    for i = eachindex(chro.vert)
        chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
    end

    # iterate through the chromatin vertices
    for i = 1:spar.chromatinLength*spar.chromatinNumber

        # get the interactions between vertex and shell
        getForce, unitVector, closePointDistance, closeCoords, closeVertices, tri = vertex_shell_interaction(chro.vert[i],enve,envelopeTree,spar,[],true,false,spar.repulsionDistance)

        # if the repulsion is calculated
        if getForce

            # get the envelope interaction
            unitVector, forceMagnitude, acrossBoundary = get_vertex_shell_interaction_vertex_force(enve,unitVector,closePointDistance, closeVertices, tri, spar,false,spar.repulsionDistance)

            # calculate the force
            chro.forces.enveRepulsion[i] = forceMagnitude*unitVector;

            if !acrossBoundary

                # get the forces on the envelope
                forces = get_vertex_shell_interaction_shell_force(enve,forceMagnitude,unitVector,closeVertices,closeCoords,1)

                # if the closest point in the enve is a vertex
                if length(closeVertices) == 1

                    # assign force
                    enve.forces.chromationRepulsion[closeVertices[1]] += forces[1]

                # if the closest point is on an edge 
                elseif length(closeVertices) == 2

                    # assign forces
                    enve.forces.chromationRepulsion[closeVertices[1]] += forces[1]
                    enve.forces.chromationRepulsion[closeVertices[2]] += forces[2]

                # if the closest is on the triangle
                else

                    # assign forces
                    enve.forces.chromationRepulsion[closeVertices[1]] += forces[1]
                    enve.forces.chromationRepulsion[closeVertices[2]] += forces[2]
                    enve.forces.chromationRepulsion[closeVertices[3]] += forces[3]

                end
            end
        end
    end
end

function get_micromanipulation_forces!(enve,mm,spar)

    # initialize the forces
    for i = eachindex(enve.vert)
        enve.forces.micromanipulation[i] = Vec(0.,0.,0.);
    end

    # calculate the force from the static pipette
    enve.forces.micromanipulation[mm.leftmostVertex] = spar.staticPullingForceMultiplier*spar.pullingForce*(mm.leftmostVertexPosition .- enve.vert[mm.leftmostVertex])
 
    # calculate the force from the pipette being moved
    enve.forces.micromanipulation[mm.rightmostVertex] = spar.pullingForce*(mm.pipettePosition - enve.vert[mm.rightmostVertex])
    
end

function get_lad_forces!(enve,chro,spar)

    # iterate through the chromsomes
    for i = 1:spar.chromatinNumber

        # iterate through the LADs in the chromosome
        for j = 1:length(chro.lads[i])

            # get the vector between the envelope and chromosome vertices
            vector = enve.vert[enve.lads[i][j]] - chro.vert[chro.strandIdx[i][chro.lads[i][j]]]

            # get the vector norm
            vectorNorm = norm(vector)

            # get the magnitude
            force = -spar.ladStrength*(vectorNorm - spar.ladLength)*vector/vectorNorm

            # get the forces
            enve.forces.ladEnveForces[enve.lads[i][j]] = force
            chro.forces.ladChroForces[chro.strandIdx[i][chro.lads[i][j]]] = -force

        end
    end
end

function get_plane_repulsion!(enve,simset,spar)

    # initiate the plane repulsion forces
    for i = eachindex(enve.vert)
        enve.forces.planeRepulsion[i] = Vec(0.,0.,0.);
    end

    # initiate the vector collecting the vertices that are in contact with the top plane
    simset.adh.touchingTop = zeros(Bool,length(enve.vert))

    # define the top plane center (of the sphere)
    topPlaneCenter = Vec(0.,0.,simset.adh.topPlane)

    # go through the vertices
    for i = eachindex(enve.vert)

        # check if the vertex is close to the top plane
        if norm(topPlaneCenter - enve.vert[i]) > spar.cytoskeletonPlaneRadius - spar.repulsionDistance 

            # set the touchingTop to true
            simset.adh.touchingTop[i] = true

            # get the unit vector for the interaction
            unitVector = (topPlaneCenter - enve.vert[i])/norm(topPlaneCenter - enve.vert[i])

            # check if the enve vertex is above the top plane
            if norm(topPlaneCenter - enve.vert[i]) > spar.cytoskeletonPlaneRadius

                # calculate the force
                enve.forces.planeRepulsion[i] = spar.planeRepulsionConstant*unitVector

            # otherwise
            else

                # calculate the distance from the top plane
                distance = spar.cytoskeletonPlaneRadius - norm(topPlaneCenter - enve.vert[i]);

                # calculate the force
                enve.forces.planeRepulsion[i] = spar.repulsionConstant*(spar.repulsionDistance - distance)*unitVector
            end
        end

        # check if the vertex is close to the bottom plane
        if enve.vert[i][3] - simset.adh.bottomPlane < spar.repulsionDistance

            # if the vertex is below the bottom plane
            if enve.vert[i][3] - simset.adh.bottomPlane < 0

                # calculate the force
                enve.forces.planeRepulsion[i] = Vec(0.,0.,spar.planeRepulsionConstant)

            # otherwise
            else

                # calculate the distance from the top plane
                distance = enve.vert[i][3] - simset.adh.bottomPlane;

                # calculate the force
                enve.forces.planeRepulsion[i] = Vec(0.,0.,spar.repulsionConstant*(spar.repulsionDistance - distance))

            end
        end

        if simset.adh.stickyBottom && enve.vert[i][3] - simset.adh.bottomPlane < 2*spar.repulsionDistance

            forceX = 100000*(simset.adh.originalCoordinates[i][1] - enve.vert[i][1])
            forceY = 100000*(simset.adh.originalCoordinates[i][2] - enve.vert[i][2])

            enve.forces.planeRepulsion[i] += Vec(forceX,forceY,0)

        end
    end
end

function get_crosslink_forces!(chro,spar)

    # initiate the forces
    for i = eachindex(chro.forces.crosslink)
        chro.forces.crosslink[i] = Vec(0.,0.,0.)
    end

    # iterate through the crosslinks
    for i = 1:length(chro.crosslinks)

        # get the vector and the vector norm
        vector = chro.vert[chro.crosslinks[i][1]] - chro.vert[chro.crosslinks[i][2]]
        vectorNorm = norm(vector);

        # get the force
        force = -spar.crosslinkStrength*(vectorNorm - spar.crosslinkLength).*vector./vectorNorm

        # assign the forces
        chro.forces.crosslink[chro.crosslinks[i][1]] = force
        chro.forces.crosslink[chro.crosslinks[i][2]] = -force

    end
end

function get_repl_elastic_forces!(repl,spar)

    # initialize the forces
    for i = eachindex(repl.vert)
        repl.forces.elastic[i] = Vec(0.,0.,0.);
    end

    # calculate the normal length as the average of the current distances
    normalLength = mean(repl.edgeVectorNorms)

    # go through the edges
    for i = eachindex(repl.edges)

        # check if first edge
        if repl.firstEdges[i] == 1

            # calculate the force
            force = spar.replStiffness*(repl.edgeVectorNorms[i] - normalLength)*repl.edgeUnitVectors[i];

            # assign the forces
            repl.forces.elastic[repl.edges[i][1]] += force
            repl.forces.elastic[repl.edges[i][2]] -= force
        end
    end

end

function get_repl_volume_forces!(repl,spar)
    
    replVolume = get_volume!(repl)

    if spar.replPressure != 0 && !repl.growth

        if replVolume < spar.replTargetVolume
            try
                volPressure = -spar.replBulkModulus*log10(replVolume/(repl.normalVolume));
            catch
                volPressure = -spar.replBulkModulus*1e-10;
            end
        else    
            volPressure = 0
        end
    else
        volPressure = 0
    end

    if replVolume < spar.replTargetVolume && repl.growth
        # calculate the volume force based on the growth pressure
        repl.forces.volume = (spar.replPressure + volPressure)*repl.voronoiAreas.*repl.vertexNormalUnitVectors
    else
        repl.forces.volume = volPressure*repl.voronoiAreas.*repl.vertexNormalUnitVectors
    end

end

function get_repl_area_forces!(repl,spar)

    # initialize the forces
    for i = 1:length(repl.vert)
        repl.forces.area[i] = Vec(0.,0.,0.)
    end

    # iterate through the triangles
    for i = eachindex(repl.tri)

        # get the baryocenter
        baryocenter = mean(repl.vert[repl.tri[i]]);

        # get the force magnitude
        magnitude = 0*spar.areaCompressionStiffness*(repl.triangleAreas[i] - repl.normalTriangleAreas[i])/repl.normalTriangleAreas[i];
                
        # iterate through the vertices
        for j = 1:3
        
            # get the vector and unitVector
            vector = repl.vert[repl.tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            # add the force
            repl.forces.area[repl.tri[i][j]] += -magnitude*unitVector

        end
    end
end

function get_repl_bending_forces!(repl,spar)

    # set the forces to zero
    for i = 1:length(repl.vert)
        repl.forces.bending[i] = Vec(0.,0.,0.)
    end

    # compute angles of triangles in the envelope
    angles = get_triangle_angles(repl);

    # compute bending moments based on difference between angles and normal angles
    moments = spar.replBendingStiffness*(angles .- repl.normalAngle);

    # for each edge
    for i = eachindex(repl.edges)
         
        # if this is a first edge
        if repl.firstEdges[i] == 1
            
            # compute distances between current edge and third vertex of the adjacent triangles
            distance1 = line_point_distance(repl.edgeVectors[i], repl.edgeVectors[repl.edgeThirdVertices[i][1]])
            distance2 = line_point_distance(repl.edgeVectors[i], repl.edgeVectors[repl.edgeThirdVertices[i][2]])
            
            # compute bending forces for adjacent vertices
            force1 = -moments[i]/distance1*repl.triangleNormalUnitVectors[repl.edgesTri[i][1]]
            force2 = -moments[i]/distance2*repl.triangleNormalUnitVectors[repl.edgesTri[i][2]]

            # compute counter forces for edge vertices
            counterForces = 0.5*force1 + 0.5*force2

            # add to the bending forces for vertices
            repl.forces.bending[repl.edges3Vertex[i][1]] -= force1;
            repl.forces.bending[repl.edges3Vertex[i][2]] -= force2;
            repl.forces.bending[repl.edges[i][1]] += counterForces;
            repl.forces.bending[repl.edges[i][2]] += counterForces;

        end
    end
end

function get_repl_chromatin_repulsion_forces!(chro, repl, spar)

    # initialize repl forces
    for i = eachindex(repl.vert)
        repl.forces.chromationRepulsion[i] = Vec(0.,0.,0.);
    end

    # initialize chro force
    for i = eachindex(chro.vert)
        chro.forces.replRepulsion[i] = Vec(0.,0.,0.)
    end

    # iterate through the chromatin vertices
    for i = 1:spar.chromatinLength*spar.chromatinNumber
        
        getForce, unitVector, closePointDistance, closeCoords, closeVertices, tri = vertex_shell_interaction(chro.vert[i],repl,repl.tree,spar,[],true,true,spar.repulsionDistance)

        if getForce
            
            # get the repl interaction forces
            unitVector, forceMagnitude,acrossBoundary = get_vertex_shell_interaction_vertex_force(repl,unitVector,closePointDistance, closeVertices, tri, spar,true,spar.repulsionDistance)

            # calculate the force
            chro.forces.replRepulsion[i] = forceMagnitude*unitVector

            if !acrossBoundary

                # get the forces on the envelope
                forces = get_vertex_shell_interaction_shell_force(repl,forceMagnitude,unitVector,closeVertices,closeCoords,5)

                # if the closest point in the enve is a vertex
                if length(closeVertices) == 1

                    # assign force
                    repl.forces.chromationRepulsion[closeVertices[1]] += forces[1]

                # if the closest point is on an edge    
                elseif length(closeVertices) == 2

                    # assign forces
                    repl.forces.chromationRepulsion[closeVertices[1]] += forces[1]
                    repl.forces.chromationRepulsion[closeVertices[2]] += forces[2]

                # if the closest is on the triangle
                else

                    # assign forces
                    repl.forces.chromationRepulsion[closeVertices[1]] += forces[1]
                    repl.forces.chromationRepulsion[closeVertices[2]] += forces[2]
                    repl.forces.chromationRepulsion[closeVertices[3]] += forces[3]

                end
            end
        end
    end
end

function get_repl_comp_enve_repulsion_forces!(enve, repl, spar)

    # set the force to 0
    for i = eachindex(repl.vert)
        repl.forces.envelopeRepulsion[i] = Vec(0.,0.,0.);
    end

    # set the force to 0
    for i = eachindex(enve.vert)
        enve.forces.replRepulsion[i] = Vec(0.,0.,0.)
    end

    # iterate through the vertices
    for i = eachindex(enve.vert)

        # get the interactions between vertex and shell
        getForce, unitVector, closePointDistance, closeCoords, closeVertices, tri = vertex_shell_interaction(enve.vert[i],repl,repl.tree,spar,[],true,true,spar.repulsionDistance*3)
        
        # if the repulsion is calculated
        if getForce

            # get the repl interaction
            unitVector, forceMagnitude,acrossBoundary = get_vertex_shell_interaction_vertex_force(repl,unitVector,closePointDistance, closeVertices,tri,spar,true,spar.repulsionDistance*3)
            # calculate the force   
            # enve.forces.replRepulsion[i] = forceMagnitude*unitVector;

            if !acrossBoundary

                # get the forces on the repl
                forces = get_vertex_shell_interaction_shell_force(repl,forceMagnitude,unitVector,closeVertices,closeCoords,1)

                # if the closest point in the enve is a vertex
                if length(closeVertices) == 1

                    # assign force
                    repl.forces.envelopeRepulsion[closeVertices[1]] += forces[1]

                # if the closest point is on an edge 
                elseif length(closeVertices) == 2

                    # assign forces
                    repl.forces.envelopeRepulsion[closeVertices[1]] += forces[1]
                    repl.forces.envelopeRepulsion[closeVertices[2]] += forces[2]

                # if the closest is on the triangle
                else
                    # assign forces
                    repl.forces.envelopeRepulsion[closeVertices[1]] += forces[1]
                    repl.forces.envelopeRepulsion[closeVertices[2]] += forces[2]
                    repl.forces.envelopeRepulsion[closeVertices[3]] += forces[3]

                end
            end
        end
    end

end

function vertex_shell_interaction(vertex,shell,tree,spar,neighbors,checkOutside,isRepl,repulsionDistance)

    # get the closest shell vertex
    if length(neighbors) == 0
        closest,distance = knn(tree, vertex, 1, true)
    else
        closest,distance = knn(tree, vertex, 1, true, j -> any(j .== neighbors))
    end

    # flag for calculating the force
    getForce = false

    # check if the closest point is within a certain distance
    if distance[1] < spar.meanLaminaLength*2

        # get the number of triangles of the closest vertex
        nTri = length(shell.vertexTri[closest[1]]);

        # get the neigboring triangles
        triangles = shell.tri[shell.vertexTri[closest[1]]]

        # init a vector for distances
        distanceVector = zeros(size(triangles,1))

        # iterate through the neighboring triangles
        for j = 1:nTri
            
            #  go through the triangle vertices
            for k = 1:3
                
                # calculate the sum of distances for each triangle
                distanceVector[j] += norm(shell.vert[triangles[j][k]] - vertex)

            end
        end
        
        # select the closest of the neighboring triangles based on the smallest total distance to its vertices
        tri = shell.vertexTri[closest[1]][argmin(distanceVector)]

        # compute distances and coordinates for the neighboring triangle
        closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(shell, vertex, tri)

        # get the unit vector between the chromatin vertex and the closest point in the envelope
        unitVector = (vertex - closeCoords)/closePointDistance;

        # if the closest point is within the repulsion distance
        if closePointDistance < repulsionDistance

            # calculate the force
            getForce = true
        
        # otherwise, check if the pipette surface unit vector and unitVector are somewhat in the same direction
        elseif checkOutside

            # get the normal surface vector for the closest envelope vertex
            if isRepl
                shellUnitVector = -shell.vertexNormalUnitVectors[closest[1]];
            else
                shellUnitVector = shell.vertexNormalUnitVectors[closest[1]];
            end

            # if dot product is larger than zero (more than 90 degrees difference in direction)
            if dot(unitVector,shellUnitVector) >= 0

                # calculate the force
                getForce = true

            end

        end

        return getForce, unitVector, closePointDistance, closeCoords, closeVertices, tri

    end

    return getForce, [], [], [], [], []

end

function get_vertex_shell_interaction_vertex_force(shell,unitVector,closePointDistance, closeVertices, tri, spar, isRepl, repulsionDistance)

    # if the closest point in the enve is a vertex
    if length(closeVertices) == 1
        
        # get the vertex normal unit vector
        shellUnitVector = shell.vertexNormalUnitVectors[closeVertices[1]];
    
    # if the closest point is on an edge
    elseif length(closeVertices) == 2

        # get the edge index and get the normal unit vector
        edge = findall((getindex.(shell.edges,1) .== closeVertices[1] .&& getindex.(shell.edges,2) .== closeVertices[2]))[1]
        shellUnitVector = shell.edgeNormalUnitVectors[edge]

    # if the closest point is on a triangle
    else

        # get the triangle normal unit vector
        shellUnitVector = shell.triangleNormalUnitVectors[tri]

    end

    if isRepl
        shellUnitVector = -shellUnitVector
    end

    # check the dot product again with the more accurate unit vector to check if the angle between the vectors is more than 90 degrees (has passed through the envelope surface)
    if dot(unitVector,shellUnitVector) >= 0

        # set the unit vector and force magnitude
        unitVector = -unitVector
        forceMagnitude = spar.outsideRepulsionMultiplier*spar.repulsionConstant
        acrossBoundary = true

    # outside the pipette surface
    else

        # get the magnitude
        forceMagnitude = spar.repulsionConstant*(repulsionDistance - closePointDistance)
        acrossBoundary = false

    end

    return unitVector, forceMagnitude, acrossBoundary

end

function get_vertex_shell_interaction_shell_force(shell,forceMagnitude,unitVector,closeVertices,closeCoords,multiplier)

    # get full force magnitude
    fullForceVector = -multiplier*forceMagnitude*unitVector

    # if the closest point in the enve is a vertex
    if length(closeVertices) == 1

        # assign force
        forces = tuple(fullForceVector)

    # if the closest point is on an edge 
    elseif length(closeVertices) == 2

        vert1 = shell.vert[closeVertices[1]]
        vert2 = shell.vert[closeVertices[2]]

        # calculate the force weight on each edge vertex based on the distance
        if vert1[1] != vert2[1]
            w1 = (closeCoords[1] - vert1[1])/(vert1[1] - vert2[1])
        elseif  vert1[2] != vert2[2]
            w1 = (closeCoords[2] - vert1[2])/(vert1[2] - vert2[2])
        else
            w1 = (closeCoords[3] - vert1[3])/(vert1[3] - vert2[3])
        end
        w2 = 1 - w1;

        # assign the counter forces to the edge vertices based on the weights
        forces = (w1*fullForceVector, w2*fullForceVector)

    # if the closest is on the triangle
    else

        # get baryocentric weights
        # https://answers.unity.com/questions/383804/calculate-uv-coordinates-of-3d-point-on-plane-of-m.html
        # https://math.stackexchange.com/questions/1727200/compute-weight-of-a-point-on-a-3d-triangle

        # calculate the triangle area
        fullArea = 0.5*norm(shell.vert[closeVertices[1]] - shell.vert[closeVertices[2]])*norm(shell.vert[closeVertices[1]] - shell.vert[closeVertices[3]])
        
        # calculate the part areas for each vertex (weights)
        w1 = 0.5*norm(closeCoords - shell.vert[closeVertices[2]])*norm(closeCoords - shell.vert[closeVertices[3]])/fullArea
        w2 = 0.5*norm(closeCoords - shell.vert[closeVertices[1]])*norm(closeCoords - shell.vert[closeVertices[3]])/fullArea
        w3 = 0.5*norm(closeCoords - shell.vert[closeVertices[1]])*norm(closeCoords - shell.vert[closeVertices[2]])/fullArea

        # assign the counter forces to the edge vertices based on the weights calculated from the relative areas
        forces = (w1*fullForceVector,w2*fullForceVector,w3*fullForceVector)

    end

    return forces

end

function get_afm_forces!(enve,ext,spar)

    # initiate the plane repulsion forces
    for i = eachindex(enve.vert)
        enve.forces.afmRepulsion[i] = Vec(0.,0.,0.);
    end

    # initiate the vector collecting the vertices that are in contact with the top plane
    ext.touchingIdx = zeros(Bool,length(enve.vert))

    # go through the vertices
    for i = eachindex(enve.vert)

        # check if the vertex is close to the top plane
        if norm(ext.beadPosition - enve.vert[i]) < 3.31 + spar.repulsionDistance 

            # set the touchingTop to true
            ext.touchingIdx[i] = true

            # get the unit vector for the interaction
            unitVector = (enve.vert[i]-ext.beadPosition)/norm(enve.vert[i]-ext.beadPosition)

            # check if the enve vertex is above the top plane
            if norm(enve.vert[i]-ext.beadPosition) < 3.31

                # calculate the force
                enve.forces.afmRepulsion[i] = 10*spar.repulsionConstant*unitVector

            # otherwise
            else

                # calculate the distance from the top plane
                distance = norm(ext.beadPosition - enve.vert[i]) - 3.31;

                # calculate the force
                enve.forces.afmRepulsion[i] = 10*spar.repulsionConstant*(spar.repulsionDistance - distance)*unitVector
            end
        end
    end

end