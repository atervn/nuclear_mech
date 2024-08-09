function get_lad_centers(enve,spar)

    # calculate the nucleus area
    nucleusArea = 4*pi*spar.freeNucleusRadius^2;

    # calculate min distance between lad centers
    minDistance = 1.4*sqrt(nucleusArea/spar.chromatinNumber/pi)

    # init list to store lad center indices
    ladCenterIdx = zeros(Int64,spar.chromatinNumber)

    # init counter to track num of loop iterations
    reDoAllCounter = 0;

    # repeat until solution found or max num of iterations reached
    while true

        # set flag to indicate if all lad centers are redone
        reDoAll = false

        # iterate over lad centers
        for i = 1:spar.chromatinNumber

            # init counter to track num of loop iterations for single chromatin
            reDoChromosomeCounter = 0

            # repeat until solution found or max num of iterations reached
            while true
                
                # generate random vertex index
                tempIdx = rand(1:length(enve.vert))

                # set flag to indicate if all no other lad centers are not too close
                allGood = true

                # check if vertex far enough from all other lad centers
                for j = 1:i-1

                    # calculate distance between current vertex and jth lad center, if 
                    # distance < min distance, then vertex not far enough away
                    if norm(enve.vert[tempIdx] - enve.vert[ladCenterIdx[j]]) < minDistance
                        allGood = false
                        break
                    end
                end

                # if vertex far enough away, then add it to list of lad centers
                if allGood
                    ladCenterIdx[i] = tempIdx
                    break
                end
                
                # add to the counter to redo chromatin
                reDoChromosomeCounter += 1

                # stop if maximum iteration number is reached and redo all
                if reDoChromosomeCounter == 100
                    reDoAll = true
                    break
                end
            end

            # if all chromatins are redone, add to the counter
            if reDoAll
                reDoAllCounter += 1
                break
            end
        end

        # if flag not set, then break out of loop
        if !reDoAll
            break
        end

        # if max num of iterations reached, then break out of loop
        if reDoAllCounter == 100
            break
        end

        # otherwise, decrease min distance by 5% and try again
        minDistance = 0.95*minDistance;
    end

    return ladCenterIdx

end

function get_lad_enve_vertices(ladCenterIdx,enve,spar)

    # initialize a list to store the lad vertices
    ladVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    # get the number of lad vertices for each lad region
    ladNumbers = rand(spar.minLadNumber:spar.maxLadNumber,spar.chromatinNumber)

    # initialize a list to store the close vertices for each lad region
    closeVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    for i = 1:spar.chromatinNumber
        closeVertices[i] = []
    end

    # For each vertex in the envelope, find the closest lad center
    for i = 1:length(enve.vert)

        # initialize the closest distance to a large value
        closestDistance = 1000;
        
        # initialize the closest vertex
        closestVertex = 0;

        # for each lad center, check if the current vertex is closer than the closest distance
        for j = 1:spar.chromatinNumber

            # Calculate the distance between the current vertex and the jth lad center and
            # if the distance is less than the closest distance, then update the closest distance and closest vertex
            if norm(enve.vert[i] - enve.vert[ladCenterIdx[j]]) < closestDistance
                closestDistance = norm(enve.vert[i] - enve.vert[ladCenterIdx[j]])
                closestVertex = j
            end
        end

        # add the current vertex to the list of close vertices for the closest lad center
        push!(closeVertices[closestVertex],i)
    end

    # for each lad region, sample the close vertices to get the lad vertices
    for i = 1:spar.chromatinNumber

        # if the number of close vertices is greater than or equal to the number of lad vertices for the current lad region,
        # get the specified number of lad vertices
        if length(closeVertices[i]) >= ladNumbers[i]
            ladVertices[i] = sample(closeVertices[i],ladNumbers[i],replace = false)
        
        # otherwise, get all close vertices in a random order
        else
            ladVertices[i] = sample(closeVertices[i],length(closeVertices[i]),replace = false)
        end
    end

    return ladVertices

end

function get_lad_chro_vertices(enve,chro,spar)

    # int a vector to store chromatin lad vertices
    ladVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    # for each lad region
    for i = 1:spar.chromatinNumber

        # get number of lad vertices in current lad region
        nLads = length(enve.lads[i]);
        
        # get number of vertices per lad
        vertPerLad = floor(spar.chromatinLength/nLads)

        # init a vector for the possible lad vertices for current lad region
        possibleLadVertices = Vector{Vector{Int64}}(undef, nLads);

        # initialize start index to 1
        startIdx = 1;

        # for each lad in current lad region
        for j = 1:nLads

            # if start index + number of vertices per lad is greater than chromatin length, then set possible lad
            # vertices to be from start index to chromatin length
            if startIdx + vertPerLad > spar.chromatinLength
                possibleLadVertices[j] = startIdx:spar.chromatinLength

            # otherwise, set possible lad vertices to be from start index to start index + number of vertices per lad
            else
                possibleLadVertices[j] = startIdx:startIdx+vertPerLad

                # increment start index by number of vertices per lad + 1
                startIdx += vertPerLad + 1;
            end
        end

        allGood = true
        
        while true

            # initialize lad vertices for current lad region to be empty
            ladVertices[i] = [];

            # for each lad in current lad region
            for j = 1:nLads
                
                tryCounter = 0

                if allGood

                    # if current lad is the first lad, then randomly select a vertex from possible lad vertices
                    if j == 1
                        tempVertex = rand(possibleLadVertices[j])
                    
                    # otherwise, randomly select a vertex from possible lad vertices until the distance between
                    # the current vertex and the previous vertex is greater than 8
                    else
                        while true

                            tempVertex = rand(possibleLadVertices[j])
                            
                            if tempVertex - ladVertices[i][end] >= spar.minLadVertexDistance

                                # calculate distance between current vertex and previous vertex in envelope
                                distanceEnve = norm(enve.vert[enve.lads[i][j]] - enve.vert[enve.lads[i][j-1]])

                                # calculate distance between current vertex and previous vertex in chromatin
                                distanceChro = (tempVertex - ladVertices[i][end])*spar.chroVertexDistance

                                # if distance in chromatin is greater than distance in envelope, then break out of loop
                                if distanceChro > distanceEnve
                                    break
                                end
                            end

                            tryCounter += 1

                            if tryCounter > 20
                                allGood = false
                                break
                            end
                        end
                    end

                end

                # add current vertex to lad vertices for current lad region
                push!(ladVertices[i],tempVertex)
            end

            if allGood
                break
            else
                allGood = true
            end
        end
    end

    # set lads attribute of chro to ladVertices
    chro.lads = ladVertices;

    return chro
    
end