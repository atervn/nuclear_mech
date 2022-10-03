function create_lads(chro,nuc,spar)

    ladCenterIdx = get_lad_centers(nuc,spar)

    nuc.lads = get_lad_enve_vertices(ladCenterIdx,nuc,spar)

    chro.lads = get_lad_chro_vertices(nuc,spar)

    return nuc,chro

end

function get_lad_centers(nuc,spar)

    nucleusArea = 4*pi*spar.freeNucleusRadius^2;
    minDistance = 1.4*sqrt(nucleusArea/spar.chromatinNumber/pi)

    ladCenterIdx = zeros(Int64,46)

    reDoAllCounter = 0;
    while true
        reDoAll = false
        for i = 1:spar.chromatinNumber
            reDoChromosomeCounter = 0
            while true
                
                tempIdx = rand(1:length(nuc.vert))
                allGood = true

                for j = 1:i-1

                    if norm(nuc.vert[tempIdx] - nuc.vert[ladCenterIdx[j]]) < minDistance
                        allGood = false
                        break
                    end

                end
                if allGood
                    ladCenterIdx[i] = tempIdx
                    break
                end
                
                reDoChromosomeCounter += 1

                if reDoChromosomeCounter == 100
                    reDoAll = true
                    break
                end
            end

            if reDoAll
                reDoAllCounter += 1
                break
            end
        end

        if reDoAllCounter == 100
            break
        elseif !reDoAll
            break
        else
            minDistance = 0.95*minDistance;
        end
    end

    return ladCenterIdx

end

function get_lad_enve_vertices(ladCenterIdx,nuc,spar)

    ladVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    usedVertices = [];
    ladNumbers = rand(8:13,spar.chromatinNumber)

    closeVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    for i = 1:spar.chromatinNumber
        closeVertices[i] = []
    end

    for i = 1:length(nuc.vert)
        closestDistance = 1000;
        closestVertex = 0;
        for j = 1:spar.chromatinNumber
            if norm(nuc.vert[i] - nuc.vert[ladCenterIdx[j]]) < closestDistance
                closestDistance = norm(nuc.vert[i] - nuc.vert[ladCenterIdx[j]])
                closestVertex = j
            end
        end

        push!(closeVertices[closestVertex],i)
    end

    for i = 1:spar.chromatinNumber
        ladVertices[i] = sample(closeVertices[i],ladNumbers[i])
    end


    return ladVertices

end

function get_lad_chro_vertices(nuc,spar)

    ladVertices = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    for i = 1:spar.chromatinNumber

        nLads = length(nuc.lads[i]);
        
        vertPerLad = round(spar.chromatinLength/nLads)

        possibleLadVertices = Vector{Vector{Int64}}(undef, nLads);

        startIdx = 1;
        for j = 1:nLads

            if startIdx + vertPerLad > spar.chromatinLength
                possibleLadVertices[j] = startIdx:spar.chromatinLength
            else
                possibleLadVertices[j] = startIdx:startIdx+vertPerLad
                startIdx += vertPerLad + 1;
            end
        end

        ladVertices[i] = [];

        for j = 1:nLads
            
            push!(ladVertices[i],rand(possibleLadVertices[j]))
            
        end


    end

    return ladVertices

end
