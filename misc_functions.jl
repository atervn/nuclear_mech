function get_friction_matrix(nuc,chro,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.vert)+spar.chromatinLength*spar.chromatinNumber,length(nuc.vert)+spar.chromatinLength*spar.chromatinNumber);

    for j = 1:length(nuc.vert)
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(nuc.neighbors[j]);
        frictionMatrix[j,nuc.neighbors[j]] .= -spar.laminaFriction;
    end

    chroStart = length(nuc.vert)+1
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
        for k = 1:length(chro.lads[j])

            nucVertex = nuc.lads[j][k]
            chroVertex = startIdx + chro.lads[j][k] - 1

            frictionMatrix[nucVertex,chroVertex] = -spar.laminaFriction*constant
            frictionMatrix[chroVertex,nucVertex] = -spar.laminaFriction*constant
            frictionMatrix[chroVertex,chroVertex] += spar.laminaFriction*constant
            frictionMatrix[nucVertex,nucVertex] += spar.laminaFriction*constant

        end

    end

    return frictionMatrix
end

#########################################################################################################

function read_parameters(ipar,filePath)

    f = open(filePath)

    while !eof(f)
        line = split(readline(f), ',')
        setproperty!(ipar,Symbol(line[1]),parse(Float64,line[2]))
    end

    return ipar
end

#########################################################################################################

function get_model_parameters(ipar,spar,nuc)

    spar.laminaStiffness = ipar.laminaStiffness/ipar.viscosity*ipar.scalingTime;

    spar.laminaFriction = ipar.laminaFriction/ipar.viscosity;

    spar.areaCompressionStiffness = ipar.areaCompressionModulus/(mean(nuc.normalLengths).*ipar.scalingLength);
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

    spar.meanLaminaLength = mean(nuc.normalLengths)

    spar.pullingForce = ipar.pullingForce/ipar.viscosity/ipar.scalingLength*ipar.scalingTime;

    spar.iLUCutoff = ipar.iLUCutoff
    spar.exportStep = ipar.exportStep

    return spar

end

#########################################################################################################

function setup_export(folderName::String,nuc,chro,spar,nameDate::Bool,exportData::Bool)

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

        ex.enveCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(nuc.tri))
        for i = eachindex(nuc.tri)
            ex.enveCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.tri[i]);
        end

        ex.chroCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,spar.chromatinNumber)
        for i = 1:spar.chromatinNumber
            ex.chroCells[i] = MeshCell(PolyData.Lines(), chro.strandIdx[i]);
        end


        totalNum = sum(length.(nuc.lads));

        ex.ladCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,totalNum)

        for i = 1:totalNum
            ex.ladCells[i] = MeshCell(PolyData.Lines(), [i, totalNum+i]);
        end

        ex.ladIdx = []
        ex.ladEnveVertices = []
        for i = 1:spar.chromatinNumber
            for j = 1:length(nuc.lads[i])
                push!(ex.ladIdx, i)
                push!(ex.ladEnveVertices, nuc.lads[i][j])
            end
        end

        ex.ladIdx = Int64.(ex.ladIdx)

        ex.ladChroVertices = []
        for i = 1:spar.chromatinNumber
            for j = 1:length(chro.lads[i])
                push!(ex.ladChroVertices, chro.strandIdx[i][chro.lads[i][j]])
            end
        end

        # export lad indices
        ladExport = Matrix{Int64}[];
        for i = 1:spar.chromatinNumber
            for j = 1:length(nuc.lads[i])
                newLine = [i nuc.lads[i][j] chro.lads[i][j]]
                push!(ladExport,newLine)
            end
        end

        writedlm(".\\results\\"*ex.folderName*"\\lads.csv", ladExport,',')

        ex.step = spar.exportStep
    end
    
    return ex

end

#########################################################################################################

function export_pipette_mesh(folderName, pip)

    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(pip.tri,1))
    for i = 1:size(pip.tri,1)
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, pip.tri[i,:])
    end

    vtk_save(vtk_grid(".\\results\\"*folderName*"\\pipette", [getindex.(pip.vert,1) getindex.(pip.vert,2) getindex.(pip.vert,3)]', triCells))

end

#########################################################################################################

function get_strand_vectors!(chro,spar)

    for i = 1:spar.chromatinNumber
        chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];
        chro.vectorNorms[i] = norm.(chro.vectors[i])
    end

end

#########################################################################################################

function setup_micromanipulation(nuc)

    mm = micromanipulationType();

    mm.leftmostVertex = argmin(getindex.(nuc.vert,1));
    mm.rightmostVertex = argmax(getindex.(nuc.vert,1));

    mm.leftNeighbors = nuc.neighbors[mm.leftmostVertex];
    mm.rightNeighbors = nuc.neighbors[mm.rightmostVertex];

    mm.leftmostVertexPosition = nuc.vert[mm.leftmostVertex]
    mm.leftNeigborPositions = nuc.vert[mm.leftNeighbors]

    return mm

end

#########################################################################################################

function export_data(nuc,chro,spar,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(nuc.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],nuc.vert[ext[1].rightmostVertex][1] - nuc.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));

            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(nuc.vert,1) getindex.(nuc.vert,2) getindex.(nuc.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.triangleNormalUnitVectors,1) getindex.(nuc.triangleNormalUnitVectors,2) getindex.(nuc.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(nuc.forces.volume,1) getindex.(nuc.forces.volume,2) getindex.(nuc.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(nuc.forces.area,1) getindex.(nuc.forces.area,2) getindex.(nuc.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.forces.elastic,1) getindex.(nuc.forces.elastic,2) getindex.(nuc.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.envelopeRepulsion,1) getindex.(nuc.forces.envelopeRepulsion,2) getindex.(nuc.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.chromationRepulsion,1) getindex.(nuc.forces.chromationRepulsion,2) getindex.(nuc.forces.chromationRepulsion,3)]'

            end
            vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(chro.forces.ladChroForces,1) getindex.(chro.forces.ladChroForces,2) getindex.(chro.forces.ladChroForces,3)]'
                # vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                # vtk["Movement", VTKPointData()] = [dt*vX[length(nuc.vert)+1:end] dt*vY[length(nuc.vert)+1:end] dt*vZ[length(nuc.vert)+1:end]]'
            end

            # export lads
        
            tempVert = [chro.vert[ex.ladChroVertices] ; nuc.vert[ex.ladEnveVertices]]
            vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
                vtk["LAD ID"] = ex.ladIdx
            end

            crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(chro.crosslinks))
            for i = 1:length(chro.crosslinks)
                crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(chro.crosslinks)+i]);
            end
            tempVert = [chro.vert[getindex.(chro.crosslinks,1)] ; chro.vert[getindex.(chro.crosslinks,2)]];
            vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
            end

            writedlm(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0") * ".csv", [getindex.(chro.crosslinks,1) getindex.(chro.crosslinks,2)],',')

        end
    end
end

function solve_system!(nuc,chro,spar,simset,dt,ext)

    
    movements = Vector{Vec{3,Float64}}(undef,length(nuc.vert)+length(chro.vert))

    maxMovement::Float64 = 0;
    maxMovInd::Int64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = nuc.forces.total .+ enveFlucs;
        totalChro = chro.forces.total .+ fluctuationForces;

        solX = cg(simset.frictionMatrix,[getindex.(totalEnve,1);getindex.(totalChro,1)],Pl = simset.iLU)
        solY = cg(simset.frictionMatrix,[getindex.(totalEnve,2);getindex.(totalChro,2)],Pl = simset.iLU)
        solZ = cg(simset.frictionMatrix,[getindex.(totalEnve,3);getindex.(totalChro,3)],Pl = simset.iLU)

        movements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

        maxMovement = 0;

        for i = eachindex(movements)
            
            movementNorm = norm(movements[i]);
            if movementNorm >= maxMovement
                maxMovement = movementNorm
                maxMovInd = i
            end

            if movementNorm >= 0.5
                everythingIsFine = false
                simset.timeStepMultiplier = simset.timeStepMultiplier/2
                break
            end
        end

        if everythingIsFine
            break
        end

    end

    nuc.vert .+= movements[1:length(nuc.vert)]
    chro.vert .+= movements[length(nuc.vert)+1:end]

    if simset.simType == "PC"
        planeMovement = (9e5 - ext[4])*dt.*simset.timeStepMultiplier
        if planeMovement > 0.001
            ext[1] -= 10*dt.*simset.timeStepMultiplier
        else
            println((9e5 - ext[4]))
            ext[1] -= planeMovement
        end
    end


    multiplier::Rational = 1;
    if simset.timeStepMultiplier != 1 && maxMovement <= 0.2
        multiplier = 2;
    end
    
    while true
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier
            simset.timeStepProgress += simset.timeStepMultiplier
            if simset.timeStepProgress == 1
                simset.timeStepProgress = 0
            end
            break
        else
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;
        end
    end
end

function solve_system_init!(nuc,chro,spar,simset,dt,ext,noEnve)

    
    movements = Vector{Vec{3,Float64}}(undef,length(nuc.vert)+length(chro.vert))

    maxMovement::Float64 = 0;
    maxMovInd::Int64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = nuc.forces.total .+ enveFlucs;
        totalChro = chro.forces.total .+ fluctuationForces;

        solX = cg(simset.frictionMatrix,[getindex.(totalEnve,1);getindex.(totalChro,1)],Pl = simset.iLU)
        solY = cg(simset.frictionMatrix,[getindex.(totalEnve,2);getindex.(totalChro,2)],Pl = simset.iLU)
        solZ = cg(simset.frictionMatrix,[getindex.(totalEnve,3);getindex.(totalChro,3)],Pl = simset.iLU)

        movements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

        maxMovement = 0;

        for i = eachindex(movements)
            
            movementNorm = norm(movements[i]);
            if movementNorm >= maxMovement
                maxMovement = movementNorm
                maxMovInd = i
            end

            if movementNorm >= 0.5
                everythingIsFine = false
                simset.timeStepMultiplier = simset.timeStepMultiplier/2
                break
            end
        end

        if everythingIsFine
            break
        end

    end
    if !noEnve
        nuc.vert .+= movements[1:length(nuc.vert)]
    end
    chro.vert .+= movements[length(nuc.vert)+1:end]

    if simset.simType == "PC"
        planeMovement = (9e5 - ext[4])*dt.*simset.timeStepMultiplier
        if planeMovement > 0.001
            ext[1] -= 10*dt.*simset.timeStepMultiplier
        else
            println((9e5 - ext[4]))
            ext[1] -= planeMovement
        end
    end

    multiplier::Rational = 1;
    if simset.timeStepMultiplier != 1 && maxMovement <= 0.2
        multiplier = 2;
    end
    
    while true
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier
            simset.timeStepProgress += simset.timeStepMultiplier
            if simset.timeStepProgress == 1
                simset.timeStepProgress = 0
            end
            break
        else
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;
        end
    end
end

function get_nuclear_properties!(nuc,chro,simset,spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(nuc.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(nuc);
    nuc.triangleAreas = get_area!(nuc)
    get_voronoi_areas!(nuc);
    get_triangle_normals!(nuc);
    get_area_unit_vectors!(nuc);
    # get_local_curvatures!(nuc);
end

function save_specific_data!(nuc,ext,simset)
   
    # save data for analysis
    if simset.timeStepProgress == 0

        if cmp(simset.simType,"MA") == 0 # for MA, save maximum x-coordinate
            push!(ext[2],maximum(getindex.(nuc.vert,1)));
        elseif cmp(simset.simType,"MM") == 0 # for MM, save distance between the manipulated vertices
            push!(ext[2],nuc.vert[ext[1].rightmostVertex][1] - nuc.vert[ext[1].leftmostVertex][1]);
        end
    end
end

function check_simulation_type(simType)

    # check that the simulation type is correct
    if cmp(simType,"MA") != 0 && cmp(simType,"MM") != 0 && cmp(simType,"PC") != 0 && cmp(simType,"INIT") != 0
        printstyled("Unknown simulation type"; color = :blue)
        return true
    end

    return false
end

function setup_simulation(initType::String,simType::String,importFolder::String,parameterFile::String)

    # read model parameters from file
    ipar = inputParametersType();
    ipar = read_parameters(ipar,parameterFile);

    # get import folder
    importFolder = get_import_folder(initType,importFolder)

    # create nucleus
    nuc = nucleusType();
    if cmp(initType,"new") == 0 # create a new nuclear envelope

        printstyled("Creating nuclear envelope..."; color = :blue)
        nuc = create_icosahedron!(nuc,ipar);
        nuc = subdivide_mesh!(nuc,ipar)
    elseif cmp(initType,"load") == 0 # load nuclear envelope from previous simulation
        
        printstyled("Loading nuclear envelope..."; color = :blue)
        nuc,importNumber = import_envelope(nuc,importFolder)
    end

    nuc.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.vert))
    nuc.forces.area = Vector{Vec{3,Float64}}(undef, length(nuc.vert))
    nuc.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.vert))
    nuc.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.vert))
    nuc.forces.ladEnveForces = Vector{Vec{3,Float64}}(undef, length(nuc.vert))
    nuc.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert))

    nuc = setup_nucleus_data(nuc)
    printstyled("Done!\n"; color = :blue)
    
    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);
    
    chro = chromatinType();
    if cmp(initType,"new") == 0
        ladCenterIdx = get_lad_centers(nuc,spar)
        nuc.lads = get_lad_enve_vertices(ladCenterIdx,nuc,spar)

        printstyled("Creating chromatin..."; color = :blue)

        chro = create_all_chromsomes(chro,spar,nuc.vert[ladCenterIdx])
        
        chro.lads = get_lad_chro_vertices(nuc,spar)

        chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

        for i = 1:spar.chromatinNumber

            chro.crosslinked[chro.strandIdx[i][chro.lads[i]]] .= -1

        end
    elseif cmp(initType,"load") == 0
        printstyled("Loading chromatin..."; color = :blue)
        chro = import_chromatin(chro,importFolder,importNumber)

        nuc,chro = import_lads(nuc,chro,importFolder,spar)

        chro = import_crosslinks(chro,importFolder,importNumber,spar)

        chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

        for i = 1:spar.chromatinNumber

            chro.crosslinked[chro.strandIdx[i][chro.lads[i]]] .= -1

        end
    end

    chro.neighbors = Vector{Vector{Int64}}(undef,spar.chromatinNumber*spar.chromatinLength)
    ind = 1
    for i = 1:spar.chromatinNumber
        for j = 1:spar.chromatinLength

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

    

    printstyled("Done!\n"; color = :blue)

    simset = simulationSettingsType()
    simset.frictionMatrix = get_friction_matrix(nuc,chro,spar)
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
    elseif cmp(simType,"PC") == 0
        topPlane = spar.freeNucleusRadius + spar.repulsionDistance;
        bottomPlane = -spar.freeNucleusRadius - spar.repulsionDistance;
        ext = [topPlane,bottomPlane,zeros(Bool,length(nuc.vert)),0]
    elseif cmp(simType,"INIT") == 0
        ext = ()
    end

    return nuc, chro, spar, simset, ext

end

function import_envelope(nuc,importFolder)

    files = readdir(importFolder)
    ifNucFile = zeros(Bool,length(files))
    for i = eachindex(ifNucFile)
        ifNucFile[i] = cmp(files[i][1:4],"nucl") == 0
    end

    nucFileIdx = findall(ifNucFile)

    numTimePoints = length(nucFileIdx)

    numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

    timePointNumbers = zeros(Int64,numTimePoints)
    for i = eachindex(timePointNumbers)

        tempNum = [filter(isdigit, collect(s)) for s in files[nucFileIdx[i]]][end-(numOfDigitsInName+3):end-4]

        numString = ""
        for j = 1:numOfDigitsInName
            numString = numString*string(tempNum[j][1])
        end

        timePointNumbers[i] = parse(Int64,numString)
    end

    importNumber = lpad(maximum(timePointNumbers),numOfDigitsInName,"0")

    importName = "nucl_"*importNumber

    vtk = VTKFile(importFolder*"\\"*importName*".vtu")

    vert = get_points(vtk)

    nuc.vert = Vector{Vec{Float64,3}}(undef,size(vert)[2])
    for i = eachindex(vert[1,:])
        nuc.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    VTKCelldata = get_cells(vtk)
    tri = VTKCelldata.connectivity

    # convert data to the required format
    tri = reshape(tri,(3,:))
    tri = tri' .+ 1
    nuc.tri = Vector{Vector{Int64}}(undef,size(tri,1))
    for i = eachindex(tri[:,1])
        nuc.tri[i] = tri[i,:]
    end

    nuc = get_edges(nuc)
    nuc = get_vertex_triangles(nuc)
    
    return nuc,importNumber
end

function import_chromatin(chro,importFolder,importNumber)

    vtk = VTKFile(importFolder*"\\chro_"*importNumber*".vtp")
    vert = get_points(vtk)

    chro.vert = Vector{Vec{Float64,3}}(undef,size(vert)[2])
    for i = eachindex(vert[1,:])
        chro.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    endPoints = get_primitives(vtk,"Lines").offsets

    chro.strandIdx = Vector{Vector{Int64}}(undef,length(endPoints))
    for i = eachindex(endPoints)
        if i == 1
            chro.strandIdx[i] = collect(1:endPoints[i])
        else
            chro.strandIdx[i] = collect(endPoints[i-1] + 1:endPoints[i])
        end
    end

    chromatinLength = endPoints[1]
    chromatinNumber = length(endPoints);

    chro.forces.linear = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    chro.forces.bending = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,chromatinNumber)
    chro.vectorNorms = Vector{Vector{Float64}}(undef,chromatinNumber)
    chro.forces.strandLinear = Vector{Any}(undef,chromatinNumber)
    chro.forces.strandBending = Vector{Any}(undef,chromatinNumber)
    chro.forces.strandChroRepulsion = Vector{Any}(undef,chromatinNumber)
    chro.forces.strandEnveRepulsion = Vector{Any}(undef,chromatinNumber)
    chro.forces.strandLadChroForces = Vector{Any}(undef,chromatinNumber)
    initialize_chromatin_forces!(chro);
    
    chro.strandVert = Vector{Any}(undef,chromatinNumber)

    for i = 1:chromatinNumber
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

end

function import_lads(nuc,chro,importFolder,spar)

    nuc.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    chro.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    tempLads = readdlm(importFolder*"\\lads.csv", ',', Int64, '\n')

    for i = 1:spar.chromatinNumber

        tempIdx = findall(tempLads[:,1] .== i)
        nuc.lads[i] = tempLads[tempIdx,2]
        chro.lads[i] = tempLads[tempIdx,3]

    end

    return nuc,chro
end

function get_crosslinks!(nuc,chro,simset,spar)

    constant = 0.001;

    changesDone = false

    # remove crosslinks
    nLinked = length(chro.crosslinks)
    probs = rand(nLinked)
    for i = nLinked:-1:1

        if probs[i] < spar.crosslinkingUnbindingProbability*spar.maxDt*simset.timeStepMultiplier

            chro.crosslinked[chro.crosslinks[i][1]] = 0
            chro.crosslinked[chro.crosslinks[i][2]] = 0

            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][1], length(nuc.vert) + chro.crosslinks[i][2]] = 0
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][2], length(nuc.vert) + chro.crosslinks[i][1]] = 0
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][1], length(nuc.vert) + chro.crosslinks[i][1]] -= spar.laminaFriction*constant
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][2], length(nuc.vert) + chro.crosslinks[i][2]] -= spar.laminaFriction*constant

            chro.crosslinks = chro.crosslinks[1:end .!= i]

            changesDone = true

        end
    end

    dropzeros!(simset.frictionMatrix)

    # form crosslinks
    notCrosslinked = findall(chro.crosslinked .== 0)
    closestVerts = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)
    possiblyLinking =  zeros(Bool,spar.chromatinLength*spar.chromatinNumber)

    for i = notCrosslinked
        closest,distance = knn(simset.chromatinTree, chro.vert[i],1,true,j -> any(j .== chro.neighbors[i]))
        if distance[1] <= 0.5
            closestVerts[i] = closest[1]
            possiblyLinking[i] = true
        end
    end
    
    possibleLinkingIdx = findall(possiblyLinking)

    
    for i = possibleLinkingIdx

        if chro.crosslinked[i] == 0 && chro.crosslinked[closestVerts[i]] == 0

            if rand() < spar.crosslinkingBindingProbability*spar.maxDt*simset.timeStepMultiplier

                push!(chro.crosslinks, [i, closestVerts[i]])
                chro.crosslinked[i] = 1
                chro.crosslinked[closestVerts[i]] = 1

                simset.frictionMatrix[length(nuc.vert) + i, length(nuc.vert) + closestVerts[i]] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + closestVerts[i], length(nuc.vert) + i] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + i, length(nuc.vert) + i] += spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + closestVerts[i], length(nuc.vert) + closestVerts[i]] += spar.laminaFriction*constant

                changesDone = true
            end
        end
    end
    if changesDone
        iLU = ilu(simset.frictionMatrix, τ = spar.iLUCutoff)
    end

end

function import_crosslinks(chro,importFolder,folderNumber,spar)

    tempCrosslinks = readdlm(importFolder*"\\crosslinks_" * folderNumber * ".csv", ',', Int64, '\n')

    chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

    for i = eachindex(tempCrosslinks[:,1])

        push!(chro.crosslinks, tempCrosslinks[i,:])
        chro.crosslinked[tempCrosslinks[i,:]] .= 1

    end

    for i = 1:spar.chromatinNumber

        chro.crosslinked[chro.strandIdx[i][chro.lads[i]]] .= -1

    end

    return chro

end

function progress_time!(simset,intTime)

    if simset.timeStepProgress == 0
        next!(simset.prog)
        intTime += 1
    end

    return intTime

end

function post_export(ex,simset,ext)

    if ex.exportData
        if cmp(simset.simType, "MA") == 0
            writedlm(".\\results\\" * ex.folderName * "\\maxX.csv", maxX, ',')
            dL = ext[2] .- minimum(ext[2])

            J = 2 * pi .* dL ./ (3 * 2.1 * 3 * 1)

            plot(10*dt:dt:maxT*dt, J[11:end], yaxis=:log, xaxis=:log, xlim=(0.1, 200), ylim=(0.01, 10))
        elseif cmp(simset.simType, "MM") == 0
            writedlm(".\\results\\" * ex.folderName * "\\nuclearLength.csv", ext[2], ',')
        end
    end

end

function get_import_folder(initType,importFolder)
    
    if initType == "load"
        if importFolder == "" # open dialog if no folder was provided
            importFolder = pick_folder(pwd()*"\\results")
        else # get the path if folder was provided
            importFolder = pwd()*"\\results\\"*importFolder;
        end
    else
        importFolder = ""
    end

    return importFolder
end