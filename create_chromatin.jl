function create_all_chromsomes(nuc,spar,ladCenterIdx)

    nuc.chro.vert = Vector{Vec{3,Float64}}(undef, spar.chromatinNumber*spar.chromatinLength);
    nuc.chro.strandIdx = Vector{Vector{Int64}}(undef,spar.chromatinNumber)
    nuc.chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    reDoAllCounter = 0;    
    while true

        reDoAll = false

        for k = 1:spar.chromatinNumber

            startInd = 1 + (k-1)*spar.chromatinLength
            endInd = spar.chromatinLength + (k-1)*spar.chromatinLength

            reDoCounter = 0;
            while true

                vert = create_chromatin_polymer(nuc,spar,startInd,ladCenterIdx,k)

                if vert isa Bool
                    reDoCounter += 1
                else
                    nuc.chro.vert[startInd:endInd] = vert;
                    nuc.chro.strandIdx[k] = collect(startInd:endInd);
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

            nuc.chro.forces.linear = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.bending = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,spar.chromatinNumber)
            nuc.chro.vectorNorms = Vector{Vector{Float64}}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandLinear = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandBending = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandChroRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandEnveRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandLadChroForces = Vector{Any}(undef,spar.chromatinNumber)
            initialize_chromatin_forces!(nuc);
            
            for i = 1:spar.chromatinNumber
                nuc.chro.strandVert[i] = @view nuc.chro.vert[nuc.chro.strandIdx[i]];
                nuc.chro.vectors[i] = nuc.chro.strandVert[i][2:end] .- nuc.chro.strandVert[i][1:end-1];
                nuc.chro.vectorNorms[i] = norm.(nuc.chro.vectors[i])
                nuc.chro.forces.strandLinear[i] = @view nuc.chro.forces.linear[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandBending[i] = @view nuc.chro.forces.bending[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandChroRepulsion[i] = @view nuc.chro.forces.chroRepulsion[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandEnveRepulsion[i] = @view nuc.chro.forces.enveRepulsion[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandLadChroForces[i] = @view nuc.chro.forces.ladChroForces[nuc.chro.strandIdx[i]];
            end


            return nuc
        elseif reDoAllCounter > 100
            return false
        end

    end

    return nuc

end

function create_chromatin_polymer(nuc,spar,startInd,ladCenterIdx,k)
    
    vert = Vector{Vec{3,Float64}}(undef,spar.chromatinLength)

    newVertex = get_first_vertex(nuc,spar,startInd,ladCenterIdx[k])

    if newVertex isa Bool
        return false
    else
        vert[1] = newVertex;
    end

    # get the other vertices
    for i = 2:spar.chromatinLength
        
        newVertex = get_other_vertices(nuc,spar,vert,startInd,i,ladCenterIdx,k)

        if newVertex isa Bool
            return false
        else
            vert[i] = newVertex;
        end

    end

    return vert
end

function get_first_vertex(nuc,spar,startInd,ladCenterIdx)
    
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
            if norm(vertex - nuc.chro.vert[j]) < spar.repulsionDistance
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

function get_other_vertices(nuc,spar,vert,startInd,currentInd,ladCenterIdx,k)

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
            if norm(newPoint - nuc.chro.vert[j]) < spar.repulsionDistance
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

function add_chromatin_adh_init(importFolder,exportFolder; parameterFile::String="./parameters_adh_init.txt")

    ipar = inputParametersType();
    ipar = read_parameters(ipar,parameterFile);

    nuc.enve = setup_nucleus(ipar,"load",importFolder)

    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);

    nuc.chro = chromatinType();

    ladCenterIdx = get_lad_centers(nuc,spar);
    nuc.enve.lads = get_lad_enve_vertices(ladCenterIdx,nuc,spar);

    nuc.chro = create_all_chromsomes_adh_init(nuc.chro,spar,nuc.enve.vert[ladCenterIdx],nuc);

    nuc.chro.lads = get_lad_chro_vertices(nuc,spar)

    nuc.chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

    ex = setup_export(exportFolder,nuc,chro,spar,true,true)

    simset = simulationSettingsType()
    simset.simType = "INIT";
    
    ext = check_adhesion_file!(ex,"load",importFolder,simset)

    export_data(nuc,chro,spar,ex,ext,0,simset)

    return ex.folderName

end

function create_all_chromsomes_adh_init(chro,spar,ladCenterIdx,nuc)

    nuc.chro.vert = Vector{Vec{3,Float64}}(undef, spar.chromatinNumber*spar.chromatinLength);
    nuc.chro.strandIdx = Vector{Vector{Int64}}(undef,spar.chromatinNumber)
    nuc.chro.strandVert = Vector{Any}(undef,spar.chromatinNumber)

    envelopeTree = KDTree(nuc.enve.vert);

    reDoAllCounter = 0;    
    while true

        reDoAll = false

        for k = 1:spar.chromatinNumber

            startInd = 1 + (k-1)*spar.chromatinLength
            endInd = spar.chromatinLength + (k-1)*spar.chromatinLength

            reDoCounter = 0;
            while true

                vert = create_chromatin_polymer_adh_init(chro,spar,startInd,ladCenterIdx,k,envelopeTree,nuc)

                if vert isa Bool
                    reDoCounter += 1
                else
                    nuc.chro.vert[startInd:endInd] = vert;
                    nuc.chro.strandIdx[k] = collect(startInd:endInd);
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
            println("chromatin $k done")
        end

        if !reDoAll

            nuc.chro.forces.linear = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.bending = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,spar.chromatinNumber*spar.chromatinLength)
            nuc.chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,spar.chromatinNumber)
            nuc.chro.vectorNorms = Vector{Vector{Float64}}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandLinear = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandBending = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandChroRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandEnveRepulsion = Vector{Any}(undef,spar.chromatinNumber)
            nuc.chro.forces.strandLadChroForces = Vector{Any}(undef,spar.chromatinNumber)
            initialize_chromatin_forces!(chro);
            
            for i = 1:spar.chromatinNumber
                nuc.chro.strandVert[i] = @view nuc.chro.vert[nuc.chro.strandIdx[i]];
                nuc.chro.vectors[i] = nuc.chro.strandVert[i][2:end] .- nuc.chro.strandVert[i][1:end-1];
                nuc.chro.vectorNorms[i] = norm.(nuc.chro.vectors[i])
                nuc.chro.forces.strandLinear[i] = @view nuc.chro.forces.linear[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandBending[i] = @view nuc.chro.forces.bending[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandChroRepulsion[i] = @view nuc.chro.forces.chroRepulsion[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandEnveRepulsion[i] = @view nuc.chro.forces.enveRepulsion[nuc.chro.strandIdx[i]];
                nuc.chro.forces.strandLadChroForces[i] = @view nuc.chro.forces.ladChroForces[nuc.chro.strandIdx[i]];
            end


            return chro
        elseif reDoAllCounter > 100
            return false
        end

    end
end

function create_chromatin_polymer_adh_init(chro,spar,startInd,ladCenterIdx,k,envelopeTree,nuc)
    
    vert = Vector{Vec{3,Float64}}(undef,spar.chromatinLength)

    newVertex = get_first_vertex_adh_init(chro,spar,startInd,ladCenterIdx[k],envelopeTree,nuc)

    if newVertex isa Bool
        return false
    else
        vert[1] = newVertex;
    end

    # get the other vertices
    for i = 2:spar.chromatinLength
        
        newVertex = get_other_vertices_adh_init(chro,spar,vert,startInd,i,ladCenterIdx,k,envelopeTree,nuc)

        if newVertex isa Bool
            return false
        else
            vert[i] = newVertex;
        end
    end

    return vert
end

function get_first_vertex_adh_init(chro,spar,startInd,ladCenterIdx,envelopeTree,nuc)
    
    xLimits = [minimum(getindex.(nuc.enve.vert,1)) maximum(getindex.(nuc.enve.vert,1))]
    yLimits = [minimum(getindex.(nuc.enve.vert,2)) maximum(getindex.(nuc.enve.vert,2))]
    zLimits = [minimum(getindex.(nuc.enve.vert,3)) maximum(getindex.(nuc.enve.vert,3))]
    
    reDoCounter = 0;

    # get the starting vertex
    while true
            
        reDoStart = false
        

        # get spherical coodinates
        # polarAngle = rand()*180;
        # azimuthAngle = rand()*360;
        # orgDist = rand()*(maxDistance-spar.repulsionDistance);

        # convert to cartesian coordinates
        # vertex = Vec(orgDist*sind(polarAngle)*cosd(azimuthAngle), orgDist*sind(polarAngle)*sind(azimuthAngle),orgDist*cosd(polarAngle))

        vertex = Vec((rand(1)*(xLimits[2]-xLimits[1]))[1] + xLimits[1],(rand(1)*(yLimits[2]-yLimits[1]))[1] + yLimits[1],(rand(1)*(zLimits[2]-zLimits[1]))[1] + zLimits[1])

        # check if this is too close to any existing vertex
        for j = 1:startInd-1
            if norm(vertex - nuc.chro.vert[j]) < spar.repulsionDistance
                reDoStart = true
                break
            end
        end

        if check_if_outside_envelope(nuc,envelopeTree,vertex,spar)
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

function check_if_outside_envelope(nuc,envelopeTree,vertex,spar)

    closest,distance = knn(envelopeTree, vertex,1,true)

    nTri = length(nuc.enve.vertexTri[closest[1]]);

    neighbors = nuc.enve.neighbors[closest[1]]
    neigborsCoords = nuc.enve.vert[neighbors]

    triangles = nuc.enve.tri[nuc.enve.vertexTri[closest[1]]]
    distanceVector = zeros(size(triangles,1))
    for j = 1:nTri
        for k = 1:3
            distanceVector[j] += norm(nuc.enve.vert[triangles[j][k]] - vertex)
        end
    end

    tri = nuc.enve.vertexTri[closest[1]][argmin(distanceVector)]

    closePointDistance,closeCoords,closeVertices = vertex_triangle_distance(nuc.enve, vertex, tri)

    unitVector = (vertex - closeCoords)/closePointDistance;

    if length(closeVertices) == 1
        enveUnitVector = nuc.enve.vertexNormalUnitVectors[closeVertices[1]];
    elseif length(closeVertices) == 2
        edge = findall((getindex.(nuc.enve.edges,1) .== closeVertices[1] .&& getindex.(nuc.enve.edges,2) .== closeVertices[2]))[1]
        enveUnitVector = nuc.enve.edgeNormalUnitVectors[edge]
    else
        enveUnitVector = nuc.enve.triangleNormalUnitVectors[tri]
    end

    if dot(unitVector,enveUnitVector) <= 0
        return false
    elseif closePointDistance < spar.repulsionDistance
        return false
    else
        return true
    end

end

function get_other_vertices_adh_init(chro,spar,vert,startInd,currentInd,ladCenterIdx,k,envelopeTree,nuc)

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
            if norm(newPoint - nuc.chro.vert[j]) < spar.repulsionDistance
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

                if check_if_outside_envelope(nuc,envelopeTree,newPoint,spar)
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

function check_adhesion_file!(ex,initState,importFolder,simset)

    if initState == "load"
        println("blob")
        if importFolder[1] == '.'
            folderTemp = importFolder
        else
            folderTemp = "./results/"*importFolder
        end
        
        if isfile(folderTemp*"\\adh.txt")
            println("blob2")
            open(".\\results\\"*ex.folderName*"\\adh.txt", "w") do file
                write(file, "adh")
            end
            simset.adh = true
            importNumber = get_import_number(folderTemp)
            planes = readdlm(folderTemp*"\\planes_"*importNumber*".csv")
            return planes
        end
    end

    return []
end