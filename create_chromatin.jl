function create_all_chromsomes(chro,spar,ladCenterIdx)

    chro.vert = Vector{Vec{3,Float64}}(undef, spar.chromatinNumber*spar.chromatinLength);
    chro.strandIdx = Vector{Vector{Int64}}(undef,spar.chromatinNumber)
    chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    reDoAllCounter = 0;    
    while true

        reDoAll = false

        for k = 1:spar.chromatinNumber

            startInd = 1 + (k-1)*spar.chromatinLength
            endInd = spar.chromatinLength + (k-1)*spar.chromatinLength

            reDoCounter = 0;
            while true

                vert = create_chromatin_polymer(chro,spar,startInd,ladCenterIdx,k)

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
end

function create_chromatin_polymer(chro,spar,startInd,ladCenterIdx,k)
    
    vert = Vector{Vec{3,Float64}}(undef,spar.chromatinLength)

    newVertex = get_first_vertex(chro,spar,startInd,ladCenterIdx[k])

    if newVertex isa Bool
        return false
    else
        vert[1] = newVertex;
    end

    # get the other vertices
    for i = 2:spar.chromatinLength
        
        newVertex = get_other_vertices(chro,spar,vert,startInd,i,ladCenterIdx,k)

        if newVertex isa Bool
            return false
        else
            vert[i] = newVertex;
        end

    end

    return vert
end

function get_first_vertex(chro,spar,startInd,ladCenterIdx)
    
    reDoCounter = 0;

    # get the starting vertex
    while true
            
        reDoStart = false
        
        # get spherical coodinates
        polarAngle = rand()*180;
        azimuthAngle = rand()*360;
        orgDist = rand()*(spar.freeNucleusRadius-spar.repulsionDistance);

        # convert to cartesian coordinates
        vertex = Vec(orgDist*sind(polarAngle)*cosd(azimuthAngle), orgDist*sind(polarAngle)*sind(azimuthAngle),orgDist*cosd(polarAngle))

        # check if this is too close to any existing vertex
        for j = 1:startInd-1
            if norm(vertex - chro.vert[j]) < spar.repulsionDistance
                reDoStart = true
                break
            end
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

function get_other_vertices(chro,spar,vert,startInd,currentInd,ladCenterIdx,k)

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
                if norm(newPoint) >= (spar.freeNucleusRadius - spar.repulsionDistance)
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