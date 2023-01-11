function create_all_chromsomes(enve,chro,spar,ladCenterIdx)

    chro.vert = Vector{Vec{3,Float64}}(undef, spar.chromatinNumber*spar.chromatinLength);
    chro.strandIdx = Vector{Vector{Int64}}(undef,spar.chromatinNumber)
    chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    envelopeTree = KDTree(enve.vert);

    limits = [minimum(getindex.(enve.vert,1)) maximum(getindex.(enve.vert,1));
              minimum(getindex.(enve.vert,2)) maximum(getindex.(enve.vert,2));
              minimum(getindex.(enve.vert,3)) maximum(getindex.(enve.vert,3))]

    reDoAllCounter = 0;    
    while true

        reDoAll = false

        for k = 1:spar.chromatinNumber

            startInd = 1 + (k-1)*spar.chromatinLength
            endInd = spar.chromatinLength + (k-1)*spar.chromatinLength

            reDoCounter = 0;
            while true

                vert = create_chromatin_polymer(enve,chro,spar,startInd,ladCenterIdx,k,limits,envelopeTree)

                if vert isa Bool
                    reDoCounter += 1
                else
                    chro.vert[startInd:endInd] = vert;
                    chro.strandIdx[k] = collect(startInd:endInd);
                    break
                end
                if reDoCounter > 1000
                    reDoAll = true
                    break
                end
            end

            if reDoAll
                reDoAllCounter += 1
                break
            end
        end

        if !reDoAll

            chro.forces.linear = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            chro.forces.bending = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,spar.chromatinNumber)
            chro.vectorNorms = Vector{Vector{Float64}}(undef,spar.chromatinNumber)
            chro.forces.strandLinear = Vector{Any}(undef,spar.chromatinNumber)
            chro.forces.strandBending = Vector{Any}(undef,spar.chromatinNumber)
            chro.forces.strandChroRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            chro.forces.strandEnveRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            chro.forces.strandLadChroForces = Vector{Any}(undef,spar.chromatinNumber)
            initialize_chromatin_forces!(chro);
            
            for i = 1:spar.chromatinNumber
                chro.strandVert[i] = @view chro.vert[chro.strandIdx[i]];
                chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];
                chro.vectorNorms[i] = norm.(chro.vectors[i])
                chro.forces.strandLinear[i] = @view chro.forces.linear[chro.strandIdx[i]];
                chro.forces.strandBending[i] = @view chro.forces.bending[chro.strandIdx[i]];
                chro.forces.strandChroRepulsion[i] = @view chro.forces.chroRepulsion[chro.strandIdx[i]];
                chro.forces.strandEnveRepulsion[i] = @view chro.forces.enveRepulsion[chro.strandIdx[i]];
                chro.forces.strandLadChroForces[i] = @view chro.forces.ladChroForces[chro.strandIdx[i]];
            end


            return chro
        elseif reDoAllCounter > 100
            return false
        end

    end

    return chro

end

function create_chromatin_polymer(enve,chro, spar, startInd, ladCenterIdx, k, limits, envelopeTree)
    
    vert = Vector{Vec{3,Float64}}(undef,spar.chromatinLength)

    newVertex = get_first_vertex(enve, chro, spar, startInd, ladCenterIdx[k],  limits, envelopeTree)

    if newVertex isa Bool
        return false
    else
        vert[1] = newVertex;
    end

    # get the other vertices
    for i = 2:spar.chromatinLength
        
        newVertex = get_other_vertices(enve, chro, spar, vert, startInd, i, ladCenterIdx, k, envelopeTree)

        if newVertex isa Bool
            return false
        else
            vert[i] = newVertex;
        end

    end

    return vert
end

function get_first_vertex(enve, chro, spar, startInd, ladCenterIdx, limits, envelopeTree)
    
    reDoCounter = 0;

    # get the starting vertex
    while true
            
        reDoStart = false
        
        # convert to cartesian coordinates
        vertex = Vec((rand(1)*(limits[1,2]-limits[1,1]))[1] + limits[1,1],
                    (rand(1)*(limits[2,2]-limits[2,1]))[1] + limits[2,1],
                    (rand(1)*(limits[3,2]-limits[3,1]))[1] + limits[3,1]);

        # check if this is too close to any existing vertex
        for j = 1:startInd-1
            if norm(vertex - chro.vert[j]) < spar.repulsionDistance
                reDoStart = true
                break
            end
        end

        if check_if_outside_envelope(enve,envelopeTree,vertex,spar)
            reDoStart = true
        end

        if norm(vertex - ladCenterIdx) > spar.freeNucleusRadius/3
            reDoStart = true
        end
        if !reDoStart
            return vertex
        else
            reDoCounter += 1
        end

        if reDoCounter > 100
            return false
        end
    end
end

function get_other_vertices(enve, chro, spar, vert, startInd, currentInd, ladCenterIdx, k, envelopeTree)

    reDoCounter = 0;
    
    while true

        reDoVert = false;
        
        # get new direction
        polarAngle = rand()*180;
        azimuthAngle = rand()*360;

        # get new vertex
        newPoint = vert[currentInd-1] + Vec(spar.chroVertexDistance*sind(polarAngle)*cosd(azimuthAngle), spar.chroVertexDistance*sind(polarAngle)*sind(azimuthAngle),spar.chroVertexDistance*cosd(polarAngle))

        # check if other vertices are too close (in other chromatins)
        for j = 1:startInd-1
            if norm(newPoint - chro.vert[j]) < spar.repulsionDistance
                reDoVert = true
                break
            end
        end

        if !reDoVert
            # check if other vertices are too close (in the same chromatin)
            for j = 1:currentInd-1
                if norm(newPoint - vert[j]) < spar.repulsionDistance
                    reDoVert = true
                    break
                end
            end
            
            if !reDoVert

                # check if too close to nuclear membrane
                if check_if_outside_envelope(enve,envelopeTree,newPoint,spar)
                    reDoVert = true
                end
                
                if !reDoVert

                    distanceFromOwn = norm(newPoint - ladCenterIdx[k])

                    for j = 1:spar.chromatinNumber

                        if j != k && (norm(newPoint - ladCenterIdx[j]) - distanceFromOwn < -1)
                            reDoVert = true
                            break
                        end

                    end

                    # if norm(newPoint - ladCenterIdx) > spar.freeNucleusRadius/3
                    #     reDoVert = true
                    # end

                    if !reDoVert
                        return newPoint
                    end
                end

            end
        end

        reDoCounter += 1

        if reDoCounter > 1000
            return false
        end
    end

end

function check_if_outside_envelope(enve,envelopeTree,vertex,spar)

    closest,distance = knn(envelopeTree, vertex,1,true)

    nTri = length(enve.vertexTri[closest[1]]);

    neighbors = enve.neighbors[closest[1]]
    neigborsCoords = enve.vert[neighbors]

    triangles = enve.tri[enve.vertexTri[closest[1]]]
    distanceVector = zeros(size(triangles,1))
    for j = 1:nTri
        for k = 1:3
            distanceVector[j] += norm(enve.vert[triangles[j][k]] - vertex)
        end
    end

    tri = enve.vertexTri[closest[1]][argmin(distanceVector)]

    closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(enve, vertex, tri)

    unitVector = (vertex - closeCoords)/closePointDistance;

    if length(closeVertices) == 1
        enveUnitVector = enve.vertexNormalUnitVectors[closeVertices[1]];
    elseif length(closeVertices) == 2
        edge = findall((getindex.(enve.edges,1) .== closeVertices[1] .&& getindex.(enve.edges,2) .== closeVertices[2]))[1]
        enveUnitVector = enve.edgeNormalUnitVectors[edge]
    else
        enveUnitVector = enve.triangleNormalUnitVectors[tri]
    end

    if dot(unitVector,enveUnitVector) <= 0
        return false
    elseif closePointDistance < spar.repulsionDistance
        return false
    else
        return true
    end

end