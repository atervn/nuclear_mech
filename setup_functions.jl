function setup_simulation(initType::String,simType::String,importFolder::String,parameterFile::String)

    # read model parameters from file
    ipar = inputParametersType();
    ipar = read_parameters(ipar,parameterFile);

    # get import folder
    importFolder = get_import_folder(initType,importFolder)

    nuc = nucleusType()
    nuc = setup_shell(nuc,ipar,initType,importFolder)
    
    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);
    
    nuc = setup_chromatin(nuc,spar,initType,importFolder)

    simset = simulationSettingsType()
    simset.frictionMatrix = get_friction_matrix(nuc,spar)
    simset.iLU = ilu(simset.frictionMatrix, τ=spar.iLUCutoff)
    simset.simType = simType;

    # setup aspiration
    if cmp(simType,"MA") == 0
        pip = generate_pipette_mesh();
        export_pipette_mesh(folderName,pip)
        # vector to store the aspiration lengths
        maxX = []
        ext = (pip,maxX)
    # setup micromanipulation 
    elseif cmp(simType,"MM") == 0
        mm = setup_micromanipulation(nuc)
        nuclearLength = []
        ext = (mm,nuclearLength)
    elseif cmp(simType,"INIT") == 0
        ext = ()
    end

    return nuc, spar, simset, ext

end

function setup_shell(nuc,ipar,initType,importFolder)

    if cmp(initType,"new") == 0 # create a new nuclear envelope

        printstyled("Creating nuclear envelope..."; color = :blue)
        radius = ipar.freeNucleusRadius/ipar.scalingLength
        nuc.enve = get_icosaherdon!(nuc.enve,radius);
        nuc.enve = subdivide_mesh!(nuc.enve,radius,ipar.nSubdivisions)
    elseif cmp(initType,"load") == 0 # load nuclear envelope from previous simulation
        printstyled("Loading nuclear envelope..."; color = :blue)
        nuc = import_envelope(nuc,importFolder)
    end

    nuc.enve.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))
    nuc.enve.forces.area = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))
    nuc.enve.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))
    nuc.enve.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))
    nuc.enve.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))
    nuc.enve.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.enve.vert))

    nuc.enve = setup_shell_data(nuc.enve)
    printstyled("Done!\n"; color = :blue)

    return nuc
end

function setup_shell_data(shellStruct)

    shellStruct.edgeVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));
    shellStruct.edgeUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));
    shellStruct.edgeVectorNorms = Vector{Float64}(undef, length(shellStruct.edges));
    get_edge_vectors!(shellStruct);

    shellStruct.edges3Vertex = Vector{Vector{Int64}}(undef,length(shellStruct.edges));
    for i = eachindex(shellStruct.edges)
        thirdVertex1 = shellStruct.tri[shellStruct.edgesTri[i][1]][.!(shellStruct.tri[shellStruct.edgesTri[i][1]] .== shellStruct.edges[i][1] .|| shellStruct.tri[shellStruct.edgesTri[i][1]] .== shellStruct.edges[i][2])][1];
        thirdVertex2 = shellStruct.tri[shellStruct.edgesTri[i][2]][.!(shellStruct.tri[shellStruct.edgesTri[i][2]] .== shellStruct.edges[i][1] .|| shellStruct.tri[shellStruct.edgesTri[i][2]] .== shellStruct.edges[i][2])][1];

        shellStruct.edges3Vertex[i] = [thirdVertex1, thirdVertex2]
    end

    shellStruct.triEdge1 = Vector{Int64}(undef,length(shellStruct.tri))
    shellStruct.triEdge2 = Vector{Int64}(undef,length(shellStruct.tri))
    shellStruct.edgeThirdVertices = Vector{Vector{Int64}}(undef,length(shellStruct.edges))
    

    for i = eachindex(shellStruct.tri)

        shellStruct.triEdge1[i] = findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.tri[i][2])[1]

        shellStruct.triEdge2[i] = findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.tri[i][3])[1]

    end

    get_triangle_normals!(shellStruct);

    for i = eachindex(shellStruct.edges)

        firstNeighbor = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges3Vertex[i][1])[1]
        secondNeighbor = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges3Vertex[i][2])[1]

        shellStruct.edgeThirdVertices[i] = [firstNeighbor, secondNeighbor]

    end

    shellStruct.normalVolume = get_volume!(shellStruct);
    shellStruct.normalTriangleAreas = get_area!(shellStruct);
    shellStruct.normalArea = sum(shellStruct.normalTriangleAreas);
    shellStruct.normalAngle = mean(get_triangle_angles(shellStruct));
    lengths = zeros(Float64,length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)  
        lengths[i] = norm(shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges[i][1]]);
    end
    shellStruct.normalLengths = lengths;

    return shellStruct

end

function setup_chromatin(nuc,spar,initType,importFolder)

    if cmp(initType,"new") == 0
        ladCenterIdx = get_lad_centers(nuc,spar)
        nuc.enve.lads = get_lad_enve_vertices(ladCenterIdx,nuc,spar)

        printstyled("Creating chromatin..."; color = :blue)

        nuc = create_all_chromsomes(nuc,spar,nuc.enve.vert[ladCenterIdx])
        
        nuc = get_lad_chro_vertices(nuc,spar)

        nuc.chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

        for i = 1:spar.chromatinNumber

            nuc.chro.crosslinked[nuc.chro.strandIdx[i][nuc.chro.lads[i]]] .= -1

        end
    elseif cmp(initType,"load") == 0
        printstyled("Loading chromatin..."; color = :blue)

        nuc = import_chromatin(nuc,importFolder)

        nuc = import_lads(nuc,importFolder,spar)

        nuc = import_crosslinks(nuc,importFolder,spar)

        nuc.chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

        for i = 1:spar.chromatinNumber

            nuc.chro.crosslinked[nuc.chro.strandIdx[i][nuc.chro.lads[i]]] .= -1

        end
    end

    nuc.chro.neighbors = Vector{Vector{Int64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    ind = 1
    for _ = 1:spar.chromatinNumber
        for j = 1:spar.chromatinLength

            if j == 1
                nuc.chro.neighbors[ind] = [ind, ind+1, ind+2]
            elseif j == 2
                nuc.chro.neighbors[ind] = [ind-1, ind, ind+1, ind+2]
            elseif j == spar.chromatinLength-1
                nuc.chro.neighbors[ind] = [ind-2, ind-1, ind, ind+1]
            elseif j == spar.chromatinLength
                nuc.chro.neighbors[ind] = [ind-2, ind-1, ind]
            else
                nuc.chro.neighbors[ind] = [ind-2, ind-1, ind, ind+1, ind+2]
            end
            
            ind += 1
        end
    end

    printstyled("Done!\n"; color = :blue)

    return nuc
end

function setup_micromanipulation(nuc)

    mm = micromanipulationType();

    mm.leftmostVertex = argmin(getindex.(nuc.enve.vert,1));
    mm.rightmostVertex = argmax(getindex.(nuc.enve.vert,1));

    mm.leftNeighbors = nuc.enve.neighbors[mm.leftmostVertex];
    mm.rightNeighbors = nuc.enve.neighbors[mm.rightmostVertex];

    mm.leftmostVertexPosition = nuc.enve.vert[mm.leftmostVertex]
    mm.leftNeigborPositions = nuc.enve.vert[mm.leftNeighbors]

    return mm

end

function setup_export(folderName::String,nuc,spar,nameDate::Bool,exportData::Bool)

    ex = exportSettingsType()

    ex.exportData = exportData

    if ex.exportData
        if nameDate
            ex.folderName = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")*folderName
        else
            ex.folderName = folderName
        end

        try
            mkdir(".\\results\\"*ex.folderName)
        catch
            for i = 1:1000
                try
                    mkdir(".\\results\\"*ex.folderName*"_"*string(i))
                    ex.folderName = ex.folderName*"_"*string(i)
                    break
                catch
                end
            end
        end

        ex.enveCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(nuc.enve.tri))
        for i = eachindex(nuc.enve.tri)
            ex.enveCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.enve.tri[i]);
        end

        ex.chroCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,spar.chromatinNumber)
        for i = 1:spar.chromatinNumber
            ex.chroCells[i] = MeshCell(PolyData.Lines(), nuc.chro.strandIdx[i]);
        end


        totalNum = sum(length.(nuc.enve.lads));

        ex.ladCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,totalNum)

        for i = 1:totalNum
            ex.ladCells[i] = MeshCell(PolyData.Lines(), [i, totalNum+i]);
        end

        ex.ladIdx = []
        ex.ladEnveVertices = []
        for i = 1:spar.chromatinNumber
            for j = 1:length(nuc.enve.lads[i])
                push!(ex.ladIdx, i)
                push!(ex.ladEnveVertices, nuc.enve.lads[i][j])
            end
        end

        ex.ladIdx = Int64.(ex.ladIdx)

        ex.ladChroVertices = []
        for i = 1:spar.chromatinNumber
            for j = 1:length(nuc.chro.lads[i])
                push!(ex.ladChroVertices, nuc.chro.strandIdx[i][nuc.chro.lads[i][j]])
            end
        end

        # export lad indices
        ladExport = Matrix{Int64}[];
        for i = 1:spar.chromatinNumber
            for j = 1:length(nuc.enve.lads[i])
                newLine = [i nuc.enve.lads[i][j] nuc.chro.lads[i][j]]
                push!(ladExport,newLine)
            end
        end

        writedlm(".\\results\\"*ex.folderName*"\\lads.csv", ladExport,',')

        ex.step = spar.exportStep
    end
    
    return ex
end

function get_friction_matrix(nuc,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.enve.vert)+spar.chromatinLength*spar.chromatinNumber,length(nuc.enve.vert)+spar.chromatinLength*spar.chromatinNumber);

    for j = 1:length(nuc.enve.vert)
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(nuc.enve.neighbors[j]);
        frictionMatrix[j,nuc.enve.neighbors[j]] .= -spar.laminaFriction;
    end

    chroStart = length(nuc.enve.vert)+1
    for j = chroStart:chroStart+spar.chromatinLength*spar.chromatinNumber-1
        frictionMatrix[j,j] = 1
    end

    

    for j = 1:spar.chromatinNumber
        startIdx = chroStart + spar.chromatinLength*(j-1)
        endIdx = chroStart + spar.chromatinLength*(j)-1

        constant = 0.0001;
        for k = startIdx:endIdx
            if k == startIdx
                frictionMatrix[k,k+1] = -spar.laminaFriction*constant
                frictionMatrix[k,k] += spar.laminaFriction*constant
            elseif k == endIdx
                frictionMatrix[k,k-1] = -spar.laminaFriction*constant
                frictionMatrix[k,k] += spar.laminaFriction*constant
            else
                frictionMatrix[k,k+1] = -spar.laminaFriction*constant
                frictionMatrix[k,k-1] = -spar.laminaFriction*constant
                frictionMatrix[k,k] += 2*spar.laminaFriction*constant
            end
        end

        constant = 0.0001;
        for k = 1:length(nuc.chro.lads[j])

            nucVertex = nuc.enve.lads[j][k]
            chroVertex = startIdx + nuc.chro.lads[j][k] - 1

            frictionMatrix[nucVertex,chroVertex] = -spar.laminaFriction*constant
            frictionMatrix[chroVertex,nucVertex] = -spar.laminaFriction*constant
            frictionMatrix[chroVertex,chroVertex] += spar.laminaFriction*constant
            frictionMatrix[nucVertex,nucVertex] += spar.laminaFriction*constant

        end

    end

    return frictionMatrix
end

function read_parameters(ipar,filePath)

    f = open(filePath)

    while !eof(f)
        line = split(readline(f), ',')
        setproperty!(ipar,Symbol(line[1]),parse(Float64,line[2]))
    end

    return ipar
end

function get_model_parameters(ipar,spar,nuc)

    spar.laminaStiffness = ipar.laminaStiffness/ipar.viscosity*ipar.scalingTime;

    spar.laminaFriction = ipar.laminaFriction/ipar.viscosity;

    spar.areaCompressionStiffness = ipar.areaCompressionModulus/(mean(nuc.enve.normalLengths).*ipar.scalingLength);
    spar.areaCompressionStiffness = spar.areaCompressionStiffness/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    # spar.bendingStiffness = ipar.laminaYoung*ipar.laminaThickness^3/(12*(1-ipar.poissonsRatio^2));
    spar.bendingStiffness = 3e-19/ipar.viscosity*ipar.scalingTime/ipar.scalingLength^2;

    spar.bulkModulus = ipar.bulkModulus/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    spar.repulsionConstant = ipar.repulsionConstant/ipar.viscosity*ipar.scalingTime#/ipar.scalingLength;

    spar.repulsionDistance = ipar.repulsionDistance/ipar.scalingLength;

    spar.freeNucleusRadius = ipar.freeNucleusRadius/ipar.scalingLength;

    spar.chroVertexDistance = ipar.chroVertexDistance/ipar.scalingLength;

    spar.chromatinBendingModulus = ipar.chromatinBendingModulus/ipar.viscosity/ipar.scalingLength^2*ipar.scalingTime
    spar.chromatinStiffness = ipar.chromatinStiffness/ipar.viscosity*ipar.scalingTime;

    spar.ladStrength = ipar.ladStrength/ipar.viscosity*ipar.scalingTime;

    spar.chromatinNormalAngle = ipar.chromatinNormalAngle*pi/180

    spar.scalingLength = ipar.scalingLength
    spar.scalingTime = ipar.scalingTime
    spar.chromatinNumber = ipar.chromatinNumber
    spar.chromatinLength = ipar.chromatinLength
    spar.viscosity = ipar.viscosity
    spar.maxDt = ipar.maxDt/ipar.scalingTime

    spar.boltzmannConst = ipar.boltzmannConst/ipar.viscosity/ipar.scalingLength^2*ipar.scalingTime
    spar.temperature = ipar.temperature

    spar.crosslinkingBindingProbability = ipar.crosslinkingBindingProbability
    spar.crosslinkingUnbindingProbability = ipar.crosslinkingUnbindingProbability

    spar.meanLaminaLength = mean(nuc.enve.normalLengths)

    spar.pullingForce = ipar.pullingForce/ipar.viscosity/ipar.scalingLength*ipar.scalingTime;

    spar.iLUCutoff = ipar.iLUCutoff
    spar.exportStep = ipar.exportStep

    return spar
end

function get_repl_comp_friction_matrix(nuc,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.repl.vert),length(nuc.repl.vert));

    
    for j = 1:length(nuc.repl.vert)
        frictionMatrix[j,j] = 1 + 20*spar.laminaFriction*length(nuc.repl.neighbors[j]);
        frictionMatrix[j,nuc.repl.neighbors[j]] .= -20*spar.laminaFriction;
    end
    
    return frictionMatrix
end

function setup_simulation_adh_init(initType::String,simType::String,importFolder::String,parameterFile::String)

    # read model parameters from file
    ipar = inputParametersType();
    ipar = read_parameters(ipar,parameterFile);

    # get import folder
    importFolder = get_import_folder(initType,importFolder)

    nuc.enve = setup_envelope(ipar,initType,importFolder)
    
    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);
    
    simset = simulationSettingsType()
    simset.frictionMatrix = get_friction_matrix_adh_init(nuc,spar)
    simset.iLU = ilu(simset.frictionMatrix, τ=spar.iLUCutoff)
    simset.simType = simType;

    topPlane = spar.freeNucleusRadius + spar.repulsionDistance;
    bottomPlane = -spar.freeNucleusRadius - spar.repulsionDistance;
    ext = [topPlane,bottomPlane,zeros(Bool,length(nuc.enve.vert)),0]

    return nuc, spar, simset, ext

end

function get_friction_matrix_adh_init(nuc,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.enve.vert),length(nuc.enve.vert));

    for j = 1:length(nuc.enve.vert)
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(nuc.enve.neighbors[j]);
        frictionMatrix[j,nuc.enve.neighbors[j]] .= -spar.laminaFriction;
    end

    return frictionMatrix
end

function setup_export_adh_init(folderName::String,nuc,spar,nameDate::Bool,exportData::Bool)

    ex = exportSettingsType()

    ex.exportData = exportData

    if ex.exportData
        if nameDate
            ex.folderName = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")*folderName
        else
            ex.folderName = folderName
        end

        try
            mkdir(".\\results\\"*ex.folderName)
        catch
            for i = 1:1000
                try
                    mkdir(".\\results\\"*ex.folderName*"_"*string(i))
                    ex.folderName = ex.folderName*"_"*string(i)
                    break
                catch
                end
            end
        end

        open(".\\results\\"*ex.folderName*"\\adh.txt", "w") do file
            write(file, "adh")
        end

        ex.enveCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(nuc.enve.tri))
        for i = eachindex(nuc.enve.tri)
            ex.enveCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.enve.tri[i]);
        end

        ex.step = spar.exportStep
    end
    
    return ex

end