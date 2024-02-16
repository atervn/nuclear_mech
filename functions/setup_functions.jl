function setup_simulation(
    initType::String,
    simType::String,
    importFolder::String,
    parameterFiles::Tuple,
    noChromatin::Bool,
    noEnveSolve::Bool,
    adherent::Bool,
    adherentStatic::Bool,
    maxT::Number,
    newEnvelopeMultipliers::Bool,
    importTime::Int,
    newTargetVolume::Number,
    stickyBottom::Bool,
    restLengthRemodelling::Bool,
    laminaDisintegration::Float64)

    # read model parameters from file
    ipar = read_parameters(parameterFiles);

    check_export_number(ipar,maxT,parameterFiles)

    # get import folder
    importFolder = get_import_folder(initType,importFolder)

    if importFolder == "" && initType == "load"
        return [],[],[],[],[],[]
    end

    # setup the envelope
    enve = setup_envelope(ipar,initType,importFolder,importTime,newEnvelopeMultipliers)


    # scale parameters
    spar = get_model_parameters(ipar,enve);
    
    # setup chromatin
    chro = setup_chromatin(enve,spar,initType,importFolder,importTime,noChromatin)

    # setup simulation settings
    simset,ext = setup_simulation_settings(enve,chro,spar,noChromatin,noEnveSolve,simType,importFolder,adherent,adherentStatic,initType,maxT,importTime,stickyBottom,restLengthRemodelling)
      
    if newTargetVolume != 0

        simset.newVolumeSimulation = true
        simset.exportNormalLengths = true

        enve.targetVolume = newTargetVolume*1e-18 / ipar.scalingLength^3

    end

    simset.laminaRemodel = "tension_remodeling" 

    enve = setup_lamina_disintegration(enve, laminaDisintegration)

    return enve, chro, spar, simset, ext, ipar, importFolder

end

function setup_envelope(ipar,initType,importFolder,importTime,newEnvelopeMultipliers)

    # create new `envelopeType` object
    enve = envelopeType()

    # if simulation type is `new`, create new nuclear envelope
    if cmp(initType,"new") == 0 # create a new nuclear envelope

        # print message to console
        printstyled("Creating nuclear envelope..."; color = :blue)

        # get radius of nuclear envelope
        radius = ipar.freeNucleusRadius/ipar.scalingLength

        # get icosahedron of nuclear envelope
        enve = get_icosaherdon!(enve,radius);

        # subdivide icosahedron
        enve = subdivide_mesh!(enve,radius,ipar.nSubdivisions)

    # if simulation type is `load`, load nuclear envelope from previous simulation
    elseif cmp(initType,"load") == 0

        # print message to console
        printstyled("Loading nuclear envelope..."; color = :blue)

        # load nuclear envelope from specified import folder
        enve = import_envelope(enve,importFolder,importTime,ipar)
    end

    # init forces on envelope
    enve.forces.volume = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.area = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.bending = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.elastic = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.pipetteRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.aspiration = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.micromanipulation = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.planeRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.replRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))
    enve.forces.afmRepulsion = Vector{Vec{3,Float64}}(undef, length(enve.vert))

    # initialize the lad forces
    for i = 1:length(enve.vert)
        enve.forces.ladEnveForces[i] = Vec(0.,0.,0.) 
    end

    # init envelope multipliers
    if initType == "new" || newEnvelopeMultipliers
        enve.envelopeMultipliers = 10 .^ (ipar.laminaVariabilityMultiplier .* randn(length(enve.edges)));
    end

    # set up shell data
    enve = setup_shell_data(enve,initType,"enve")

    # print message to console
    printstyled("Done!\n"; color = :blue)

    return enve

end

function setup_repl(initType,enve,spar,ex,importFolder,importTime,vrcGrowth)

    # initialize a replication compartment
    repl = replicationCompartmentType()

    replInitType = ""
    if initType == "new"
        replInitType = "new"
    elseif initType == "load"

        if isfile(importFolder*"\\inf.txt")
            replInitType = "load"
        else
            replInitType = "new"
        end
    end

    # if simulation type is `new`, create new nuclear envelope
    if replInitType == "new" # create a new vrc

        repl = create_replication_compartment(repl,enve,spar)

        repl.initVolume = get_volume!(repl);

    # if simulation type is `load`, load nuclear envelope from previous simulation
    elseif replInitType == "load"

        repl = import_replication_compartment(repl,spar,importFolder,importTime)

    end

    # initialize various force vectors for each vertex in the replication compartment
    repl.forces.volume = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.area = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.bending = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.elastic = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))

    
    tempVolume = repl.initVolume;

    # set up shell data for the replication compartment
    repl = setup_shell_data(repl,initType,"repl")

    repl.initVolume = tempVolume;

    # calculate various properties and vectors for the replication compartment
    get_edge_vectors!(repl);
    repl.triangleAreas = get_area!(repl)
    get_voronoi_areas!(repl);
    get_shell_normals!(repl);
    get_area_unit_vectors!(repl);

    # calculate the friction matrix for the replication compartment
    repl.frictionMatrix = get_repl_friction_matrix(repl,spar)
    
    # perform iLU (incomplete LU) factorization of the friction matrix with a specified cutoff value
    repl.iLU = ilu(repl.frictionMatrix, τ = spar.iLUCutoff)

    # calculate various properties and vectors for the envelope
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    repl.baseArea = mean(enve.triangleAreas)
    repl.normalVolume = get_volume!(repl)

    # export normal volume (initial volume)
    writedlm(".\\results\\"*ex.folderName*"\\replInitVolume.csv", repl.normalVolume.*spar.scalingLength^3,',')

    open(".\\results\\"*ex.folderName*"\\inf.txt", "w") do file
        write(file, "infected")
    end

    repl.growth = vrcGrowth

    return repl

end

function setup_shell_data(shellStruct,initType,shellType)

    # initialize the edge vectors, unit vectors, and norms
    shellStruct.edgeVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));
    shellStruct.edgeUnitVectors = Vector{Vec{3,Float64}}(undef, length(shellStruct.edges));
    shellStruct.edgeVectorNorms = Vector{Float64}(undef, length(shellStruct.edges));

    # get edge vectors
    get_edge_vectors!(shellStruct);

    # initialize the edges3Vertex array and get third verticex for the the neighboring triangles for each edge
    shellStruct.edges3Vertex = Vector{Vector{Int64}}(undef,length(shellStruct.edges));
    for i = eachindex(shellStruct.edges)
        thirdVertex1 = shellStruct.tri[shellStruct.edgesTri[i][1]][.!(shellStruct.tri[shellStruct.edgesTri[i][1]] .== shellStruct.edges[i][1] .|| shellStruct.tri[shellStruct.edgesTri[i][1]] .== shellStruct.edges[i][2])][1];
        thirdVertex2 = shellStruct.tri[shellStruct.edgesTri[i][2]][.!(shellStruct.tri[shellStruct.edgesTri[i][2]] .== shellStruct.edges[i][1] .|| shellStruct.tri[shellStruct.edgesTri[i][2]] .== shellStruct.edges[i][2])][1];

        shellStruct.edges3Vertex[i] = [thirdVertex1, thirdVertex2]
    end

    # initialize the triEdge1, triEdge2, and edgeThirdVertices arrays
    shellStruct.triEdge1 = Vector{Int64}(undef,length(shellStruct.tri))
    shellStruct.triEdge2 = Vector{Int64}(undef,length(shellStruct.tri))
    shellStruct.edgeThirdVertices = Vector{Vector{Int64}}(undef,length(shellStruct.edges))
    
    # get the edge indices between the first and second as well as the first and third vertices of each triangle
    for i = eachindex(shellStruct.tri)

        shellStruct.triEdge1[i] = findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.tri[i][2])[1]
        shellStruct.triEdge2[i] = findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.tri[i][3])[1]

    end

    # calculate the triangle normals
    get_shell_normals!(shellStruct);

    # for each edge,get the edge indices between the first vertex of the edge and the third vertices of the neighboring triangles
    for i = eachindex(shellStruct.edges)

        firstNeighbor = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges3Vertex[i][1])[1]
        secondNeighbor = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][1] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges3Vertex[i][2])[1]
        shellStruct.edgeThirdVertices[i] = [firstNeighbor, secondNeighbor]

    end

    # if the simulation is new or if the shell is not enve, calculate the normal properties
    if initType == "new" || !(shellType == "enve")
        shellStruct.normalVolume = get_volume!(shellStruct);
        shellStruct.normalTriangleAreas = get_area!(shellStruct);
        shellStruct.normalArea = sum(shellStruct.normalTriangleAreas);
        shellStruct.normalAngle = mean(get_triangle_angles(shellStruct)[shellStruct.firstEdges .== 1]);

        lengths = zeros(Float64,length(shellStruct.edges));
        for i = eachindex(shellStruct.edges)  
            lengths[i] = norm(shellStruct.vert[shellStruct.edges[i][2]] - shellStruct.vert[shellStruct.edges[i][1]]);
        end
        shellStruct.normalLengths = lengths;
    end

    return shellStruct

end

function setup_chromatin(enve,spar,initType,importFolder,importTime,noChromatin)

    # init new chromatin object
    chro = chromatinType()

    # if chromatin not disabled, create or load chromatin
    if !noChromatin

        # if `initType` is "load", check if the chromatin file exists. If it does not exist, create a new chromatin
        if initType == "load"
            importNumber = get_import_number(importFolder,importTime)
            if !isfile(importFolder *"/chro_"*importNumber*".vtp")
                createNew = true
            else
                createNew = false
            end
        end

        # if initType is new or createNew is true, create new chromatin
        if initType == "new" || createNew

            # get lad center indices
            ladCenterIdx = get_lad_centers(enve,spar)

            # get lad envelope vertices
            enve.lads = get_lad_enve_vertices(ladCenterIdx,enve,spar)

            # print a message to console
            printstyled("Creating chromatin..."; color = :blue)

            # create all chromosomes
            chro = create_all_chromsomes(enve,chro,spar,enve.vert[ladCenterIdx])

            # get chromatin lad vertices
            chro = get_lad_chro_vertices(enve,chro,spar)

            # init crosslinks array
            chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

            # set the crosslink state of all lad vertices to -1 (unable to crosslink)
            for i = 1:spar.chromatinNumber
                chro.crosslinked[chro.strandIdx[i][chro.lads[i]]] .= -1
            end
        elseif cmp(initType,"load") == 0

            # print a message to console
            printstyled("Loading chromatin..."; color = :blue)

            # load chromatin object from file
            chro = import_chromatin(chro,spar,importFolder,importTime)

            # import lads from file
            enve,chro = import_lads(enve,chro,spar,importFolder,)

            # import crosslinks from file
            chro = import_crosslinks(chro,spar,importFolder,importTime)

        end

        # create list of neighbors for each vertex in chromatin
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

        # print done message to console
        printstyled("Done!\n"; color = :blue)
    end

    return chro

end

function setup_micromanipulation(enve,spar)

    # init object
    mm = micromanipulationType();

    # find the leftmost and rightmost vertices based on their x-coordinate
    mm.leftmostVertex = argmin(getindex.(enve.vert,1))
    mm.rightmostVertex = argmax(getindex.(enve.vert,1))

    # store the position of the leftmost vertex
    mm.leftmostVertexPosition = enve.vert[mm.leftmostVertex]

    # set the pipette position
    mm.pipettePosition = enve.vert[mm.rightmostVertex]

    # load the pipette movements
    tempMovements = readdlm("./parameters/pipette_movement.csv",',')
    tempMovements[:,1] ./= spar.scalingTime
    tempMovements[:,2] ./= spar.scalingLength
    mm.pipetteMovements = tempMovements

    return mm

end

function setup_export(simType,folderName::String,enve,chro,ext,spar,simset,nameDate::Bool,exportData::Bool,noChromatin::Bool,ipar,newTargetVolume,importFolder,simulationDate)

    # init object
    ex = exportSettingsType()
    ex.exportData = exportData

    # if data is exported
    if ex.exportData

        # create a folder for exporting results with a given name and current date/time
        if nameDate
            if simulationDate != ""
                ex.folderName = simulationDate*folderName
            else
                ex.folderName = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")*folderName
            end
        else
            ex.folderName = folderName
        end

        # If the folder already exists, append a numerical suffix to create a new folder
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

        # create MeshCell objects for envelope
        ex.enveCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(enve.tri))
        for i = eachindex(enve.tri)
            ex.enveCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, enve.tri[i]);
        end

        # if chromatin is included
        if !noChromatin

            # create MeshCell objects for chromatin
            ex.chroCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,spar.chromatinNumber)
            for i = 1:spar.chromatinNumber
                ex.chroCells[i] = MeshCell(PolyData.Lines(), chro.strandIdx[i]);
            end

            # create MeshCell objects for lad cells
            totalNum = sum(length.(enve.lads));
            ex.ladCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,totalNum)
            for i = 1:totalNum
                ex.ladCells[i] = MeshCell(PolyData.Lines(), [i, totalNum+i]);
            end

            # store lad indices and vertices for export
            ex.ladIdx = []
            ex.ladEnveVertices = []
            for i = 1:spar.chromatinNumber
                for j = 1:length(enve.lads[i])
                    push!(ex.ladIdx, i)
                    push!(ex.ladEnveVertices, enve.lads[i][j])
                end
            end

            ex.ladIdx = Int64.(ex.ladIdx)

            ex.ladChroVertices = []
            for i = 1:spar.chromatinNumber
                for j = 1:length(chro.lads[i])
                    push!(ex.ladChroVertices, chro.strandIdx[i][chro.lads[i][j]])
                end
            end

            # export lad indices to a CSV file
            ladExport = Matrix{Int64}[];
            for i = 1:spar.chromatinNumber
                for j = 1:length(enve.lads[i])
                    newLine = [i enve.lads[i][j] chro.lads[i][j]]
                    push!(ladExport,newLine)
                end
            end

            writedlm(".\\results\\"*ex.folderName*"\\lads.csv", ladExport,',')

        end

        # Create a file to indicate adhesion
        if simset.adh.adherent
            open(".\\results\\"*ex.folderName*"\\adh.txt", "w") do file
                write(file, "adh")
            end
        end
        
        ex.step = spar.exportStep
    end
    
    # export the pipette mesh
    if cmp(simType,"MA") == 0
        export_pipette_mesh(ex.folderName,ext[1])
    end

    # export normal values
    export_normal_values(enve,ex,spar,simset)

    # export parameters
    export_parameters(ipar,ex)

    return ex

end

function get_friction_matrix(enve,chro,spar,noChromatin)

    # if chromatin is present, create a friction matrix with additional space for chromatin vertices
    if !noChromatin

        frictionMatrix = sparse(Int64[],Int64[],Float64[],length(enve.vert)+spar.chromatinLength*spar.chromatinNumber,length(enve.vert)+spar.chromatinLength*spar.chromatinNumber);
    
    # if no chromatin is present, create a friction matrix for envelope vertices only
    else

        frictionMatrix = sparse(Int64[],Int64[],Float64[],length(enve.vert),length(enve.vert));

    end

    # set diagonal elements and off-diagonal elements for envelope vertices
    for j = 1:length(enve.vert)

        # diagonal elements represent the self-friction of each envelope vertex
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(enve.neighbors[j]);

        # off-diagonal elements represent the friction between neighboring envelope vertices
        frictionMatrix[j,enve.neighbors[j]] .= -spar.laminaFriction;

    end


    if !noChromatin

        # get the first chromatin index
        chroStart = length(enve.vert)+1

        # set diagonal elements for chromatin vertices
        for j = chroStart:chroStart+spar.chromatinLength*spar.chromatinNumber-1

            frictionMatrix[j,j] = 1

        end

        # go through chorosomes
        for j = 1:spar.chromatinNumber

            # start and end indices for the chromosomes
            startIdx = chroStart + spar.chromatinLength*(j-1)
            endIdx = chroStart + spar.chromatinLength*(j)-1

            # go through the chromosome vertices
            for k = startIdx:endIdx

                # if first index
                if k == startIdx

                    # off-diagonal element for the first chromatin vertex in the segment
                    frictionMatrix[k,k+1] = -spar.chromatinFriction

                    # diagonal element for the first chromatin vertex in the segment
                    frictionMatrix[k,k] += spar.chromatinFriction

                # if last index
                elseif k == endIdx

                    # off-diagonal element for the last chromatin vertex in the segment
                    frictionMatrix[k,k-1] = -spar.chromatinFriction

                    # diagonal element for the last chromatin vertex in the segment
                    frictionMatrix[k,k] += spar.chromatinFriction

                # otherwise
                else

                    # off-diagonal elements for the chromatin vertices in the middle of the segment
                    frictionMatrix[k,k+1] = -spar.chromatinFriction
                    frictionMatrix[k,k-1] = -spar.chromatinFriction

                    # diagonal element for the chromatin vertices in the middle of the segment
                    frictionMatrix[k,k] += 2*spar.chromatinFriction

                end
            end

            # go through the lads
            for k = 1:length(chro.lads[j])

                # get the enve and chro indices
                enveVertex = enve.lads[j][k]
                chroVertex = startIdx + chro.lads[j][k] - 1

                # off-diagonal elements between envelope and chromatin vertices
                frictionMatrix[enveVertex,chroVertex] = -spar.ladFriction
                frictionMatrix[chroVertex,enveVertex] = -spar.ladFriction

                # diagonal elements for the chromatin vertices
                frictionMatrix[chroVertex,chroVertex] += spar.ladFriction
                frictionMatrix[enveVertex,enveVertex] += spar.ladFriction

            end
        end

        tempIdx = 1

        tempCrosslinked = chro.crosslinked

        for i = 1:length(tempCrosslinked)
            if chro.crosslinked[i] == 1
                tempCrosslinked[chro.crosslinks[tempIdx][1]] = 0
                tempCrosslinked[chro.crosslinks[tempIdx][2]] = 0
                frictionMatrix[length(enve.vert) + chro.crosslinks[tempIdx][1], length(enve.vert) + chro.crosslinks[tempIdx][2]] -= spar.crosslinkFriction
                frictionMatrix[length(enve.vert) + chro.crosslinks[tempIdx][2], length(enve.vert) + chro.crosslinks[tempIdx][1]] -= spar.crosslinkFriction
                frictionMatrix[length(enve.vert) + chro.crosslinks[tempIdx][1], length(enve.vert) + chro.crosslinks[tempIdx][1]] += spar.crosslinkFriction
                frictionMatrix[length(enve.vert) + chro.crosslinks[tempIdx][2], length(enve.vert) + chro.crosslinks[tempIdx][2]] += spar.crosslinkFriction
                tempIdx += 1
            end
        end
    end

    return frictionMatrix

end

function read_parameters(parameterFiles)

    # create a new inputParametersType object
    ipar = inputParametersType();

    for i = 1:length(parameterFiles)

        # open the file
        f = open(parameterFiles[i])

        # loop over the lines in the file
        while !eof(f)

            # split the line into a list of strings
            line = split(readline(f), ',')

            # get the name of the parameter
            name = Symbol(line[1])

            # get the value of the parameter
            value = parse(Float64, line[2])

            # set the property of the ipar object with the given name to the given value
            setproperty!(ipar, name, value)
        end

        # close the file
        close(f)
    end
    return ipar
end

function get_model_parameters(ipar,enve)

    # create new `scaledParametersType` object
    spar = scaledParametersType()

    # set lamina stiffness
    spar.laminaStiffness = ipar.laminaStiffness / ipar.viscosity * ipar.scalingTime

    # set lamina friction
    spar.laminaFriction = ipar.laminaFriction / ipar.viscosity

    # set area compression stiffness
    spar.areaCompressionStiffness = ipar.areaCompressionModulus / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    # set bending stiffness
    spar.laminaBendingStiffness = ipar.laminaBendingStiffness / ipar.viscosity * ipar.scalingTime / ipar.scalingLength^2

    spar.osmoticPressure = ipar.osmoticPressure / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    # set bulk modulus
    spar.bulkModulus = ipar.bulkModulus / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    # set repulsion constant
    spar.repulsionConstant = ipar.repulsionConstant / ipar.viscosity * ipar.scalingTime

    # set repulsion distance
    spar.repulsionDistance = ipar.repulsionDistance / ipar.scalingLength

    # set free nucleus radius
    spar.freeNucleusRadius = ipar.freeNucleusRadius / ipar.scalingLength

    # set chromatin vertex distance
    spar.chroVertexDistance = ipar.chroVertexDistance / ipar.scalingLength

    # set chromatin bending modulus
    spar.chromatinBendingModulus = ipar.chromatinBendingModulus / ipar.viscosity / ipar.scalingLength^2 * ipar.scalingTime

    # set chromatin stiffness
    spar.chromatinStiffness = ipar.chromatinStiffness / ipar.viscosity * ipar.scalingTime

    # set lad strength
    spar.ladStrength = ipar.ladStrength / ipar.viscosity * ipar.scalingTime


    spar.crosslinkStrength = ipar.crosslinkStrength / ipar.viscosity * ipar.scalingTime

    # set chromatin normal angle
    spar.chromatinNormalAngle = ipar.chromatinNormalAngle * pi / 180

    # set pipette radius
    spar.pipetteRadius = ipar.pipetteRadius / ipar.scalingLength

    # set aspiration pressure
    spar.aspirationPressure = ipar.aspirationPressure / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    # set Boltzmann constant
    spar.boltzmannConst = ipar.boltzmannConst / ipar.viscosity / ipar.scalingLength^2 * ipar.scalingTime

    # set pulling force
    spar.pullingForce = ipar.pullingForce / ipar.viscosity / ipar.scalingLength * ipar.scalingTime

    # set crosslinking binding probability
    spar.crosslinkingBindingProbability = ipar.crosslinkingBindingProbability * ipar.scalingTime

    # set crosslinking unbinding probability
    spar.crosslinkingUnbindingProbability = ipar.crosslinkingUnbindingProbability * ipar.scalingTime

    # set chromatin number
    spar.chromatinNumber = ipar.chromatinNumber

    # set chromatin length
    spar.chromatinLength = ipar.chromatinLength

    # set scaling length
    spar.scalingLength = ipar.scalingLength

    # set scaling time
    spar.scalingTime = ipar.scalingTime

    # set temperature
    spar.temperature = ipar.temperature

    # set viscosity
    spar.viscosity = ipar.viscosity

    # set maximum time step
    spar.dt = ipar.dt / ipar.scalingTime

    # set mean lamina length
    spar.meanLaminaLength = mean(enve.normalLengths)

    # set ILU cutoff
    spar.iLUCutoff = ipar.iLUCutoff

    # set export step
    spar.exportStep = ipar.exportStep

    # set maximum vertex movement
    spar.maxMovement = ipar.maxMovement / ipar.scalingLength

    # set minimum vertex movement (if movement below this, time step may be increased)
    spar.minMovement = ipar.minMovement / ipar.scalingLength

    # set adherens plane force
    spar.planeForce = ipar.planeForce / ipar.viscosity * ipar.scalingTime / ipar.scalingLength

    # set multiplier for the REPL compartment size
    spar.replSizeMultiplier = ipar.replSizeMultiplier

    # set maximum crosslink distance
    spar.maxCrosslinkDistance = ipar.maxCrosslinkDistance / ipar.scalingLength

    # get chromatin friction
    spar.chromatinFriction = ipar.chromatinFriction / ipar.viscosity

    # get LAD friction
    spar.ladFriction = ipar.ladFriction / ipar.viscosity

    # get crosslink  friction
    spar.crosslinkFriction = ipar.crosslinkFriction / ipar.viscosity

    # get REPL friction
    spar.replFriction = ipar.replFriction / ipar.viscosity

    # get maximum overlap between chromosome regions
    spar.maxChromosomeOverlap = ipar.maxChromosomeOverlap / ipar.scalingLength

    # get maximum distance multiplier to get the maximum size of LAD regions compared to the nucleus diameter
    spar.maxLadDistanceMultiplier = ipar.maxLadDistanceMultiplier

    # get minimum distance in chromsomes between two LAD vertices
    spar.minLadVertexDistance = ipar.minLadVertexDistance

    # get minimum number of LADs per chromsome
    spar.minLadNumber = ipar.minLadNumber

    # get maximum number of LADs per chromsome
    spar.maxLadNumber = ipar.maxLadNumber

    spar.staticPullingForceMultiplier = ipar.staticPullingForceMultiplier

    spar.ladLength = ipar.ladLength / ipar.scalingLength
    
    spar.crosslinkLength = ipar.crosslinkLength / ipar.scalingLength
    
    spar.cytoskeletonPlaneRadius = ipar.cytoskeletonPlaneRadius / ipar.scalingLength 

    spar.planeRepulsionConstant = ipar.planeRepulsionConstant / ipar.viscosity * ipar.scalingTime

    spar.replStiffness = ipar.replStiffness  / ipar.viscosity * ipar.scalingTime

    spar.replBendingStiffness = ipar.replBendingStiffness / ipar.viscosity * ipar.scalingTime / ipar.scalingLength^2

    spar.outsideRepulsionMultiplier = ipar.outsideRepulsionMultiplier

    spar.replPressure = ipar.replPressure / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    spar.laminaVariabilityMultiplier = ipar.laminaVariabilityMultiplier

    spar.replBulkModulus = ipar.replBulkModulus / ipar.viscosity * ipar.scalingTime * ipar.scalingLength

    spar.cantileverSpeed = ipar.cantileverSpeed / ipar.scalingLength * ipar.scalingTime 

    spar.cantileverFriction = ipar.cantileverFriction / ipar.viscosity

    spar.cantileverSpring = ipar.cantileverSpring / spar.viscosity * spar.scalingTime

    spar.cantileverMaxForce = ipar.cantileverMaxForce / spar.viscosity / spar.scalingLength *  spar.scalingTime

    spar.nucleusHeight = ipar.nucleusHeight / ipar.scalingLength

    spar.replTargetVolume = ipar.replTargetVolume / ipar.scalingLength^3

    return spar

end

function get_repl_friction_matrix(repl,spar)

    # initialize the matrix
    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(repl.vert),length(repl.vert));

    # Set elements
    for j = 1:length(repl.vert)
        frictionMatrix[j,j] = 1 + spar.replFriction*length(repl.neighbors[j]);
        frictionMatrix[j,repl.neighbors[j]] .= -spar.replFriction;
    end
    
    return frictionMatrix

end

function check_adhesion!(initType,spar,enve,importFolder,simset,adherent,adherentStatic,importTime,stickyBottom)

    # check adhesion based on the initialization type and adhesion settings
    if initType == "load"

        # adhesion information is loaded from a file
        if !(importFolder[1] == '.')
            folderTemp = importFolder
        else
            folderTemp = "./results/"*importFolder
        end

        # adherens file exists
        if isfile(folderTemp*"\\adh.txt")

            # adherent to true
            simset.adh.adherent = true

            # get import file number
            importNumber = get_import_number(folderTemp,importTime)

            # load the plane date
            planes = DataFrame(CSV.File(folderTemp*"\\planes_"*importNumber*".csv"))

            # load top and bottom planes
            simset.adh.topPlane = planes[1,3]
            simset.adh.bottomPlane = planes[2,3]

            # initialize touching vector
            simset.adh.touchingTop = zeros(Bool,length(enve.vert))

        elseif adherent

            # adherent to true
            simset.adh.adherent = true

            # define top and bottom planes
            simset.adh.topPlane = maximum(getindex.(enve.vert,3)) + spar.repulsionDistance - spar.cytoskeletonPlaneRadius
            simset.adh.bottomPlane = minimum(getindex.(enve.vert,3)) - spar.repulsionDistance

            # initialize touching vector
            simset.adh.touchingTop = zeros(Bool,length(enve.vert))

        end

    # if new simulation is adherent
    elseif adherent
        
        # adherent to true
        simset.adh.adherent = true

        # define top and bottom planes
        simset.adh.topPlane = spar.freeNucleusRadius + spar.repulsionDistance - spar.cytoskeletonPlaneRadius
        simset.adh.bottomPlane = -spar.freeNucleusRadius - spar.repulsionDistance

        # initialize touching vector
        simset.adh.touchingTop = zeros(Bool,length(enve.vert))

    end

    if stickyBottom
        simset.adh.stickyBottom = stickyBottom
        simset.adh.originalCoordinates = copy(enve.vert);
    end

    simset.adh.static = adherentStatic;

    return simset

end

function export_normal_values(enve,ex,spar,simset)

    # export normal envelope properties
    writedlm(".\\results\\"*ex.folderName*"\\normalArea.csv", [enve.normalArea].*spar.scalingLength^2,',')
    writedlm(".\\results\\"*ex.folderName*"\\normalTriangleAreas.csv", enve.normalTriangleAreas.*spar.scalingLength^2,',')
    writedlm(".\\results\\"*ex.folderName*"\\normalAngle.csv", [enve.normalAngle],',')
    writedlm(".\\results\\"*ex.folderName*"\\normalVolume.csv", enve.normalVolume.*spar.scalingLength^3,',')
    writedlm(".\\results\\"*ex.folderName*"\\envelope_multipliers.csv", enve.envelopeMultipliers,',')
    if !simset.newVolumeSimulation
        writedlm(".\\results\\"*ex.folderName*"\\normalLengths.csv", enve.normalLengths.*spar.scalingLength,',')
    end

end

function setup_simulation_settings(enve,chro,spar,noChromatin,noEnveSolve,simType,importFolder,adherent,adherentStatic,initType,maxT,importTime,stickyBottom,restLengthRemodelling)

    # create a new simulation settings object
    simset = simulationSettingsType()

    # get the friction matrix for the linear system
    simset.frictionMatrix = get_friction_matrix(enve,chro,spar,noChromatin)

    # get the incomplete LU factorization
    simset.iLU = ilu(simset.frictionMatrix, τ=spar.iLUCutoff)
    
    # set the simulation type
    simset.simType = simType;

    # set the enve solutions
    simset.noEnveSolve = noEnveSolve;

    # set the chromatin inclusion
    simset.noChromatin = noChromatin;

    # check for adhesion, if necessary
    simset = check_adhesion!(initType,spar,enve,importFolder,simset,adherent,adherentStatic,importTime,stickyBottom)

    # setup experimental aspiration
    if simType == "MA"
        
        # load the pipette mesh
        pip = get_pipette_mesh(spar,enve);
        
        # vector to store the aspiration lengths
        maxX = []

        # combine into a tuple
        exp = (pip,maxX)

    # setup experimental micromanipulation
    elseif simType == "MM"

        # setup micromanipulation
        mm = setup_micromanipulation(enve,spar)

        # initialize vectors for nuclear length and force data
        nuclearLength = []
        force = []

        # combine into a tuple
        exp = (mm,nuclearLength,force)

        # modify the friction matrix to add friction in the micromanipulation attachment vertices
        simset.frictionMatrix[mm.rightmostVertex,mm.rightmostVertex] += spar.laminaFriction
        simset.frictionMatrix[mm.leftmostVertex,mm.leftmostVertex] += spar.laminaFriction

    # if no experimental effect
    elseif simType == "INIT"
        exp = ()
    elseif simType == "AFM"

        afm = setup_afm(enve,spar)
 
        exp = (afm)

    end

    # create a progress bar
    simset.prog = Progress(Int64(round(maxT/(spar.scalingTime*spar.dt))), 0.1, "Simulating...", 100)

    return simset,exp

end

function get_pipette_mesh(spar,enve)

    # init a pipette object
    pip = pipetteType();

    # Load the STL file
    mesh = load("./parameters/pipette_v2_binary.stl")

    # init vector for the mesh triangles
    pip.tri = Vector{Vector{Int64}}(undef,length(mesh))

    # iterate through the mesh
    for i = eachindex(mesh)

        # init a temp triangle
        tempTri = zeros(Int64,3)

        # iterate through the triangle vertices
        for j = 1:3

            # get vertex coordinates
            coords = mesh[i][j]
            
            # check if the vertex already exists in the `pip.vert` vector
            pointComparison = zeros(Bool,length(pip.vert));
            for k = 1:length(pip.vert)
                if isapprox(pip.vert[k],Vec(coords[1],coords[2],coords[3]))
                    pointComparison[k] = true;
                    break
                end
            end

            # if the vertex does not exist, add it to the `pip.vert` vector and assign it a unique index
            if any(pointComparison)
                tempTri[j] = findall(pointComparison)[1]
            else
                tempTri[j] = length(pip.vert) + 1;
                push!(pip.vert,Vec(coords[1], coords[2], coords[3]))
            end
        end

        # save the new triangle
        pip.tri[i] = tempTri;

    end    

    # get approximate radius of the nucleus
    nucleusRadius = mean(norm.(enve.vert));

    # calculate the offset of the pipette tip from the nucleus, based on the pipette radius and the nucleus radius
    xOffset = spar.pipetteRadius*tan(acos(spar.pipetteRadius/(nucleusRadius + -0.5*spar.repulsionDistance)));
    
    # scale the pipette mesh to the desired size and offset it
    for i = eachindex(pip.vert)
        pip.vert[i] = pip.vert[i] .* Vec(spar.pipetteRadius,spar.pipetteRadius,spar.pipetteRadius)
        pip.vert[i] = pip.vert[i] + Vec(xOffset,0.,0.)
    end

    # create a vector of pairs of vertices that are connected by an edge
    pip.edges = Vector{Vector{Int64}}(undef,0);
    
    # create a vector of neighboring vertices for each vertex
    pip.neighbors = fill(Int[], length(pip.vert));

    # create a matrix that stores the indices of the triangles indices
    triMatrix = [getindex.(pip.tri,1) getindex.(pip.tri,2) getindex.(pip.tri,3)];

    # for each vertex, store the indices of the triangles that it belongs to
    pip.vertexTri = fill(Int[], length(pip.vert), 1);
    for i = 1:length(pip.vert)
        triangles = findall(triMatrix  .== i);
        pip.vertexTri[i] = [i[1] for i in triangles];
    end

    # for each vertex, find the indices of its neighboring vertices and the edges
    for i = 1:length(pip.vert)
        
        # find the indices of the triangles that the vertex belongs to
        hasVertex = findall(triMatrix .== i)

        # initialize a vector for the neighbors
        neighbors = Array{Int64}(undef, 0)
        
        # for each triangle that the vertex belongs to, add the other two vertices to the neighbors vector
        for j = eachindex(hasVertex)
            append!(neighbors,pip.tri[hasVertex[j][1]][pip.tri[hasVertex[j][1]] .!= i ])
        end
        
        # remove duplicates from the neighbors vector and save
        unique!(neighbors)
        pip.neighbors[i] = neighbors

        # for each neighbor, add an edge from the current vertex to the neighbor
        for j = eachindex(neighbors)
            push!(pip.edges, [i, neighbors[j]])
        end
    end

    # initialize mirrorEdges and firstEdges arrays
    pip.mirrorEdges = zeros(Int64,length(pip.edges));
    pip.firstEdges = zeros(Int64,length(pip.edges));

    # assign mirrorEdges and firstEdges values
    for i = eachindex(pip.edges)
            
        # find the mirror edge
        mirrorIdx = findall(getindex.(pip.edges,1) .== pip.edges[i][2] .&& getindex.(pip.edges,2) .== pip.edges[i][1]);
        
        # assign indices
        pip.mirrorEdges[i] = mirrorIdx[1];
        pip.mirrorEdges[mirrorIdx[1]] = i;

        # if first edge has not been assigned, assign it
        if pip.firstEdges[i] == 0 && pip.firstEdges[mirrorIdx[1]] == 0 
            pip.firstEdges[i] = 1;
        end
    end

    # create vector to store triangle indices for each edge
    pip.edgesTri = Vector{Vector{Int64}}(undef,length(pip.edges));

    for i = eachindex(pip.edges)

        # find indices of triangles that contain the edge
        triangleIdx = findall(any(triMatrix .== pip.edges[i][1], dims = 2) .&& any(triMatrix .== pip.edges[i][2], dims = 2))
    
        # add the triangles to the edge's list of triangles
        pip.edgesTri[i] = [triangleIdx[1] for triangleIdx in triangleIdx]

    end

    # vector for the edgeVectors
    pip.edgeVectors = Vector{Vec{3,Float64}}(undef, length(pip.edges));

    # calculate the edge vectors
    for i = eachindex(pip.edges)
        if pip.firstEdges[i] == 1
            pip.edgeVectors[i] = pip.vert[pip.edges[i][2]] - pip.vert[pip.edges[i][1]];
            pip.edgeVectors[pip.mirrorEdges[i]] = -pip.edgeVectors[i]
        end
    end

    # vectors for triangle edges
    pip.triEdge1 = Vector{Int64}(undef,length(pip.tri))
    pip.triEdge2 = Vector{Int64}(undef,length(pip.tri))
    
    # iterate throught the triangles and get the edges
    for i = eachindex(pip.tri)
        pip.triEdge1[i] = findall(getindex.(pip.edges,1) .== pip.tri[i][1] .&& getindex.(pip.edges,2) .== pip.tri[i][2])[1]
        pip.triEdge2[i] = findall(getindex.(pip.edges,1) .== pip.tri[i][1] .&& getindex.(pip.edges,2) .== pip.tri[i][3])[1]
    end

    # get triangle normal unit vectors
    pip.triangleNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.tri));
    for i = eachindex(pip.tri)
        normalVector = cross(pip.edgeVectors[pip.triEdge1[i]],pip.edgeVectors[pip.triEdge2[i]])
        pip.triangleNormalUnitVectors[i] = -normalVector/norm(normalVector);
    end

    # get edge normal unit vectors
    pip.edgeNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.edges));
    for i = eachindex(pip.edges)
        vector = pip.triangleNormalUnitVectors[pip.edgesTri[i][1]] + pip.triangleNormalUnitVectors[pip.edgesTri[i][2]]
        pip.edgeNormalUnitVectors[i] = vector./norm(vector);
    end

    # get vertex normal unit vectors
    pip.vertexNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.vert));

    for i = 1:length(pip.vert)
        vector = sum(pip.triangleNormalUnitVectors[pip.vertexTri[i]]);
        pip.vertexNormalUnitVectors[i] = vector./norm(vector);
    end

    # get pipette KD tree
    pip.pipetteTree = KDTree(pip.vert)

    return pip

end

function setup_afm(enve,spar)

    afm = afmType();

    afm.beadPosition = Vec(0.,0.,maximum(getindex.(enve.vert,3)) + spar.repulsionDistance + 3.31)

    afm.topPosition = afm.beadPosition + Vec(1e-6,1e-6,10.);

    afm.normDistance = norm(afm.topPosition - afm.beadPosition);

    return afm

end

function get_parameter_files(simType,nuclearMechPars,nuclearPropPars,expPars,simPars,sysPars,replPars)

    if simPars == "./parameters/simulation_parameters.txt"

        if simType == "MM"
            simPars = "./parameters/simulation_parameters_mm.txt"
        elseif simType == "MA"
            simPars = "./parameters/simulation_parameters_ma.txt"
        elseif simType == "AFM"
            simPars = "./parameters/simulation_parameters_afm.txt"
        end

    end

    return (nuclearMechPars,nuclearPropPars,expPars,simPars,sysPars,replPars)

end

function check_export_number(ipar,maxT,parameterFiles)

    nExports = (round(maxT/ipar.dt))/ipar.exportStep+1

    if nExports > 1000
        
        println("There are $nExports exported time points in the simulation. The export step can be changed in file \""*parameterFiles[4]* "\".")

    end
end

function setup_lamina_disintegration(enve, laminaDisintegration)
   
    firstEdgeIdx = findall(enve.firstEdges .== 1)
    
    nEdges = length(firstEdgeIdx);

    nDisintegrationSites = round(Int64,nEdges*laminaDisintegration)

    disintegrationSites = firstEdgeIdx[randperm(nEdges)[1:nDisintegrationSites]]

    laminaDisintegrationMultipliers = zeros(size(enve.firstEdges))
    laminaDisintegrationMultipliers[firstEdgeIdx] .= 1
    laminaDisintegrationMultipliers[disintegrationSites] .= 0.0001

    enve.laminaDisintegrationMultipliers = laminaDisintegrationMultipliers;

    return enve

end