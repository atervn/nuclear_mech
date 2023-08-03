function create_all_chromsomes(enve,chro,spar,ladCenterIdx)

    # Initialize the chromosome's vertex coordinates and vertex and strand index arrays
    chro.vert = Vector{Vec{3,Float64}}(undef, spar.chromatinNumber*spar.chromatinLength);
    chro.strandIdx = Vector{Vector{Int64}}(undef,spar.chromatinNumber)
    chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    # create a KDTree of the envelope vertices
    envelopeTree = KDTree(enve.vert);

    # create a vector of [min, max] values for each dimension of the envelope
    limits = [minimum(getindex.(enve.vert,1)) maximum(getindex.(enve.vert,1));
              minimum(getindex.(enve.vert,2)) maximum(getindex.(enve.vert,2));
              minimum(getindex.(enve.vert,3)) maximum(getindex.(enve.vert,3))]


    # counter to see if chromatins have to be redone
    reDoAllCounter = 0;

    # iterate until the chromosomes are created successfully
    while true

        # flag to indicate if all chromosomes need to be recreated
        reDoAll = false
        # iterate over each chromosome
        for k = 1:spar.chromatinNumber

            # start and end indices of the current chromosome
            startInd = 1 + (k-1)*spar.chromatinLength
            endInd = spar.chromatinLength + (k-1)*spar.chromatinLength

            # Iterate until a successful chromosome is created
            reDoCounter = 0;
            while true

                # create a polymer of chromatin vertices
                vert = create_chromatin_polymer(enve,chro,spar,startInd,ladCenterIdx,k,limits,envelopeTree)

                # if the polymer is not valid, recreate the chromosome
                if vert isa Bool
                    reDoCounter += 1
                
                # otherwise, store the vertices and strand indices for the current chromosome
                else
                    chro.vert[startInd:endInd] = vert;
                    chro.strandIdx[k] = collect(startInd:endInd);
                    break
                end

                # if the maximum number of retries is exceeded, recreate all chromosomes
                if reDoCounter > 1000
                    reDoAll = true
                    break
                end
            end

            # if all chromosomes need to be recreated, break out of the loop
            if reDoAll
                reDoAllCounter += 1
                break
            end
        end
        
        # if chromatin creating has failed, return false
        if reDoAllCounter > 100 && reDoAll
            return false
        elseif !reDoAll
            break
        end
    end

    chro = initialize_chromatin_properties(chro,spar)

    return chro

end

function create_chromatin_polymer(enve,chro, spar, startInd, ladCenterIdx, k, limits, envelopeTree)
    
    # initialize the polymer of chromatin vertices
    vert = Vector{Vec{3,Float64}}(undef,spar.chromatinLength)

    # get the first vertex
    newVertex = get_first_vertex(enve, chro, spar, startInd, ladCenterIdx[k],  limits, envelopeTree)

    # if the vertex is not valid, return false
    if newVertex isa Bool
        return false

    # otherwise, add the vertex to the polymer
    else
        vert[1] = newVertex;
    end

    # get the other vertices
    for i = 2:spar.chromatinLength
        
        # get the next vertex
        newVertex = get_other_vertices(enve, chro, spar, vert, startInd, i, ladCenterIdx, k, envelopeTree)

        # if the vertex is not valid, return false
        if newVertex isa Bool
            return false

        # otherwise, add the vertex to the polymer
        else
            vert[i] = newVertex;
        end

    end

    return vert
end

function get_first_vertex(enve, chro, spar, startInd, ladCenterIdx, limits, envelopeTree)
    
    # iterate until a valid starting vertex is found
    reDoCounter = 0;
    while true
        
        # flag to indicate if the vertex needs to be retried
        reDoStart = false
        
        # generate a random vertex in the limits of the envelope
        vertex = Vec((rand(1)*(limits[1,2]-limits[1,1]))[1] + limits[1,1],
                    (rand(1)*(limits[2,2]-limits[2,1]))[1] + limits[2,1],
                    (rand(1)*(limits[3,2]-limits[3,1]))[1] + limits[3,1]);

        # check if the vertex is too close to any existing vertices
        for j = 1:startInd-1
            if norm(vertex - chro.vert[j]) < spar.repulsionDistance
                reDoStart = true
                break
            end
        end

        # check if the vertex is outside the envelope
        if check_if_outside_envelope(enve,envelopeTree,vertex,spar)
            reDoStart = true
        end

        # check if the vertex is too far from the lad center
        if norm(vertex - ladCenterIdx) > spar.freeNucleusRadius*spar.maxLadDistanceMultiplier
            reDoStart = true
        end

        # if the vertex is valid, return it
        if !reDoStart
            return vertex

        # otherwise, try again
        else
            reDoCounter += 1
        end

        # if the maximum number of retries is exceeded, return false
        if reDoCounter > 100
            return false
        end
    end
end

function get_other_vertices(enve, chro, spar, vert, startInd, currentInd, ladCenterIdx, k, envelopeTree)

    # iterate until a valid vertex is found.
    reDoCounter = 0;
    while true

        # flag to indicate if the vertex needs to be retried
        reDoVert = false;
        
        # generate a new vertex direction
        polarAngle = rand()*180;
        azimuthAngle = rand()*360;

        # generate a new vertex position
        newPoint = vert[currentInd-1] + Vec(spar.chroVertexDistance*sind(polarAngle)*cosd(azimuthAngle), spar.chroVertexDistance*sind(polarAngle)*sind(azimuthAngle),spar.chroVertexDistance*cosd(polarAngle))

        # check if the new vertex is too close to any other vertices (in other chromatins)
        for j = 1:startInd-1
            if norm(newPoint - chro.vert[j]) < spar.repulsionDistance
                reDoVert = true
                break
            end
        end

        # if everything is fine
        if !reDoVert
            # check if the new vertex is too close to any other vertices (in the same chromatin)
            for j = 1:currentInd-1
                if norm(newPoint - vert[j]) < spar.repulsionDistance
                    reDoVert = true
                    break
                end
            end
            
            # if everything is fine
            if !reDoVert

                # check if too close or ouside of the envelope
                if check_if_outside_envelope(enve,envelopeTree,newPoint,spar)
                    reDoVert = true
                end
                
                # if everything is fine
                if !reDoVert

                    # Check if the new vertex is too close to the other lad centers
                    distanceFromOwn = norm(newPoint - ladCenterIdx[k])
                    for j = 1:spar.chromatinNumber
                        if j != k && (norm(newPoint - ladCenterIdx[j]) - distanceFromOwn < -spar.maxChromosomeOverlap)
                            reDoVert = true
                            break
                        end
                    end

                    # if the new vertex is valid, return it
                    if !reDoVert
                        return newPoint
                    end
                end

            end
        end

        reDoCounter += 1

        # if the maximum number of retries is exceeded, return false
        if reDoCounter > 1000
            return false
        end
    end
end

function check_if_outside_envelope(enve,envelopeTree,vertex,spar)

    # find the closest vertex in the envelope to the input vertex
    closest,distance = knn(envelopeTree, vertex,1,true)

    # get the number of triangles that the closest vertex belongs to
    nTri = length(enve.vertexTri[closest[1]]);

    # get the neighbors of the closest vertex
    neighbors = enve.neighbors[closest[1]]

    # the triangle vertices for the closest triangles
    triangles = enve.tri[enve.vertexTri[closest[1]]]

    # Create a vector that stores the distances between the input vertex and
    # the vertices of the triangles that the closest vertex belongs to
    distanceVector = zeros(size(triangles,1))
    for j = 1:nTri
        for k = 1:3
            distanceVector[j] += norm(enve.vert[triangles[j][k]] - vertex)
        end
    end

    # find the triangle that is closest to the input vertex
    tri = enve.vertexTri[closest[1]][argmin(distanceVector)]

    # calculate the distance between the input vertex and the closest point on the triangle
    closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(enve, vertex, tri)

    # Calculate the unit vector pointing from the closest point on the triangle to the input vertex
    unitVector = (vertex - closeCoords)/closePointDistance;

    # if the closest point in the triangle is a vertex, get its unit normal vector
    if length(closeVertices) == 1
        enveUnitVector = enve.vertexNormalUnitVectors[closeVertices[1]];

    # if the closest point in the triangle is an edge, get its unit normal vector
    elseif length(closeVertices) == 2
        edge = findall((getindex.(enve.edges,1) .== closeVertices[1] .&& getindex.(enve.edges,2) .== closeVertices[2]))[1]
        enveUnitVector = enve.edgeNormalUnitVectors[edge]
    
    # if the closest point in the triangle is on the triangle surface, get the triangle's unit normal vector
    else
        enveUnitVector = enve.triangleNormalUnitVectors[tri]
    end

    # if the dot product is below zero (the vector from surface to the vertex is pointing oppposite to the 
    # surface normal, indicating that the vertex is most likely inside)
    if dot(unitVector,enveUnitVector) >= 0
        return true

    # otherwise, if the vertex is farther than repulsion distance, return false (not outside)
    elseif closePointDistance < spar.repulsionDistance
        return true

    # otherwise return true (outside or too close)
    else
        return false
    end

end

function initialize_chromatin_forces!(chro,spar)

    # init the force vectors for all chromosomes
    chro.forces.linear = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.bending = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.crosslink = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.replRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)

    # init the force vectors for each chromosome strand
    chro.forces.strandLinear = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandBending = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandCrosslink = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandChroRepulsion = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandEnveRepulsion = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandReplRepulsion = Vector{Any}(undef,spar.chromatinNumber)
    chro.forces.strandLadChroForces = Vector{Any}(undef,spar.chromatinNumber)

    # set each force to zeros
    for i = 1:length(chro.vert)
        chro.forces.linear[i] = Vec(0.,0.,0.)
        chro.forces.bending[i] = Vec(0.,0.,0.)
        chro.forces.crosslink[i] = Vec(0.,0.,0.)
        chro.forces.chroRepulsion[i] = Vec(0.,0.,0.)
        chro.forces.enveRepulsion[i] = Vec(0.,0.,0.)
        chro.forces.replRepulsion[i] = Vec(0.,0.,0.)
        chro.forces.ladChroForces[i] = Vec(0.,0.,0.)
    end

    # set each forces for each strand as the view on the all forces
    for i = 1:spar.chromatinNumber
        chro.forces.strandLinear[i] = @view chro.forces.linear[chro.strandIdx[i]];
        chro.forces.strandBending[i] = @view chro.forces.bending[chro.strandIdx[i]];
        chro.forces.strandCrosslink[i] = @view chro.forces.crosslink[chro.strandIdx[i]];
        chro.forces.strandChroRepulsion[i] = @view chro.forces.chroRepulsion[chro.strandIdx[i]];
        chro.forces.strandEnveRepulsion[i] = @view chro.forces.enveRepulsion[chro.strandIdx[i]];
        chro.forces.strandReplRepulsion[i] = @view chro.forces.replRepulsion[chro.strandIdx[i]];
        chro.forces.strandLadChroForces[i] = @view chro.forces.ladChroForces[chro.strandIdx[i]];
    end
end

function initialize_chromatin_properties(chro,spar)
    
    # initialize vectors for chromatin vectors, vector norms and strand vectors
    chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,spar.chromatinNumber)
    chro.vectorNorms = Vector{Vector{Float64}}(undef,spar.chromatinNumber)
    chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    # for each chromosome, set the strand vertex coordinate views and calculate the vectors and norms
    for i = 1:spar.chromatinNumber
        chro.strandVert[i] = @view chro.vert[chro.strandIdx[i]];
        chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];
        chro.vectorNorms[i] = norm.(chro.vectors[i])
    end

    # initialize chromatin force vectors
    initialize_chromatin_forces!(chro,spar)

    return chro

end

function get_chromatin_neighbors!(chro,spar)

    chro.neighbors = Vector{Vector{Int64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    ind = 1
    for _ = 1:spar.chromatinNumber
        for j = 1:spar.chromatinLength

            # set neighbors for current vertex
            if j == 1
                chro.neighbors[ind] = [ind, ind+1, ind+2]
            elseif j == 2
                chro.neighbors[ind] = [ind-1, ind, ind+1, ind+2]
            elseif j == spar.chromatinLength-1
                chro.neighbors[ind] = [ind-2, ind-1, ind, ind+1]
            elseif j == spar.chromatinLength
                chro.neighbors[ind] = [ind-2, ind-1, ind]
            else
                chro.neighbors[ind] = [ind-2, ind-1, ind, ind+1, ind+2]
            end
            
            ind += 1
        end
    end

end