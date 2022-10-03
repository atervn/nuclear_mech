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

    constant = 0.05;

    for j = 1:spar.chromatinNumber
        startIdx = chroStart + spar.chromatinLength*(j-1)
        endIdx = chroStart + spar.chromatinLength*(j)-1

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
    println(spar.laminaStiffness)
    spar.laminaFriction = ipar.laminaFriction/ipar.viscosity;

    spar.areaCompressionStiffness = ipar.areaCompressionModulus/(mean(nuc.normalLengths).*ipar.scalingLength);
    spar.areaCompressionStiffness = spar.areaCompressionStiffness/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    spar.bendingStiffness = ipar.laminaYoung*ipar.laminaThickness^3/(12*(1-ipar.poissonsRatio^2));
    spar.bendingStiffness = spar.bendingStiffness/ipar.viscosity*ipar.scalingTime/ipar.scalingLength^2;

    spar.bulkModulus = ipar.bulkModulus/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    spar.repulsionConstant = ipar.repulsionConstant/ipar.viscosity*ipar.scalingTime#/ipar.scalingLength;
    println(spar.repulsionConstant)
    spar.repulsionDistance = ipar.repulsionDistance/ipar.scalingLength;

    spar.freeNucleusRadius = ipar.freeNucleusRadius/ipar.scalingLength;

    spar.chroVertexDistance = ipar.chroVertexDistance/ipar.scalingLength;

    spar.chromatinBendingModulus = ipar.chromatinBendingModulus/ipar.viscosity/ipar.scalingLength*ipar.scalingTime

    spar.chromatinStiffness = ipar.chromatinStiffness/ipar.viscosity*ipar.scalingTime;
    println(spar.chromatinStiffness)
    spar.ladStrenght = ipar.ladStrenght/ipar.viscosity*ipar.scalingTime;

    spar.chromatinNormalAngle = ipar.chromatinNormalAngle*pi/180

    spar.scalingLength = ipar.scalingLength
    spar.scalingTime = ipar.scalingTime
    spar.chromatinNumber = ipar.chromatinNumber
    spar.chromatinLength = ipar.chromatinLength
    spar.viscosity = ipar.viscosity
    spar.dt = ipar.dt/ipar.scalingTime

    spar.boltzmannConst = ipar.boltzmannConst/ipar.viscosity/ipar.scalingLength^2*ipar.scalingTime
    spar.temperature = ipar.temperature

    return spar

end

#########################################################################################################

function setup_export(folderName,nuc,chro,spar,nameDate)


    ex = exportSettingsType()
    
    if nameDate == "yes"
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

    ex.enveCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(nuc.tri,1))
    for i = 1:size(nuc.tri,1)
        ex.enveCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.tri[i,:]);
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

function export_data(nuc,chro,spar,ex,t,simType)

if mod(t,ex.step) == 0
    vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(t,4,"0"), [getindex.(nuc.vert,1) getindex.(nuc.vert,2) getindex.(nuc.vert,3)]', ex.enveCells) do vtk
        if cmp(simType,"MA") == 0 
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
    vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(t,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
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
    vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(t,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
        vtk["LAD ID"] = ex.ladIdx
    end

    crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(chro.crosslinks))
    for i = 1:length(chro.crosslinks)
        crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(chro.crosslinks)+i]);
    end
    tempVert = [chro.vert[getindex.(chro.crosslinks,1)] ; chro.vert[getindex.(chro.crosslinks,2)]];
    vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(t,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
    end

end

end

function solve_system!(nuc,chro,spar,simset)

    vX = cg(simset.frictionMatrix,[getindex.(nuc.forces.total,1);getindex.(chro.forces.total,1)]);
    vY = cg(simset.frictionMatrix,[getindex.(nuc.forces.total,2);getindex.(chro.forces.total,2)]);
    vZ = cg(simset.frictionMatrix,[getindex.(nuc.forces.total,3);getindex.(chro.forces.total,3)]);

    for i = 1:length(nuc.vert)
        nuc.vert[i] += Vec(vX[i],vY[i],vZ[i])*spar.dt
    end
   
    for k = 1:spar.chromatinLength*spar.chromatinNumber
        chro.vert[k] += Vec(vX[length(nuc.vert)+k],vY[length(nuc.vert)+k],vZ[length(nuc.vert)+k])*spar.dt
    end

end

function get_iteration_properties!(nuc,chro,simset,spar)
   
    simset.envelopeTree = KDTree(nuc.vert);
    simset.chromatinTree = KDTree(chro.vert);
    get_strand_vectors!(chro,spar)
    initialize_chromatin_forces!(chro)

    get_edge_vectors!(nuc);
    get_voronoi_areas!(nuc);
    get_area_unit_vectors!(nuc);
    # get_local_curvatures!(nuc);
    get_triangle_normals!(nuc);

end

function save_specific_data!(nuc,ext,simset,t)
   
    if cmp(simset.simType,"MA") == 0
        ext[2][t+1] = maximum(getindex.(nuc.vert,1));
    elseif cmp(simset.simType,"MM") == 0
        ext[2][t+1] = nuc.vert[ext[1].rightmostVertex][1] - nuc.vert[ext[1].leftmostVertex][1];
    end

end

function check_simulation_type(simType)

    if cmp(simType,"MA") != 0 && cmp(simType,"MM") != 0 && cmp(simType,"PC") != 0 && cmp(simType,"INIT") != 0
        printstyled("Unknown simulation type"; color = :blue)
        return true
    end

    return false

end

function setup_simulation(initType,simType,maxT,importFolder)

    # model parameters
    ipar = inputParametersType();
    ipar = read_parameters(ipar,"./parameters.txt");

    if cmp(initType,"load") == 0
        if importFolder == ""
            importFolder = pick_folder(pwd()*"\\results")
        else
            importFolder = pwd()*"\\results\\"*importFolder
        end
    end
    # create nucleus
    nuc = nucleusType();
    if cmp(initType,"new") == 0
        printstyled("Creating nuclear envelope..."; color = :blue)
        nuc = create_icosahedron!(nuc,ipar);
        nuc = subdivide_mesh!(nuc,ipar)
    elseif cmp(initType,"load") == 0
        printstyled("Loading nuclear envelope..."; color = :blue)
        nuc,importNumber = import_envelope(nuc,importFolder)
    end

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

    elseif cmp(initType,"load") == 0
        printstyled("Loading chromatin..."; color = :blue)
        chro = import_chromatin(chro,importFolder,importNumber)

        nuc,chro = import_lads(nuc,chro,importFolder,spar)

    end

    chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

    printstyled("Done!\n"; color = :blue)

    simset = simulationSettingsType()
    simset.frictionMatrix = get_friction_matrix(nuc,chro,spar)
    simset.simType = simType;

    # setup aspiration
    if cmp(simType,"MA") == 0
        pip = generate_pipette_mesh();
        export_pipette_mesh(folderName,pip)
        # vector to store the aspiration lengths
        maxX = zeros(Float64,maxT+1)
        ext = (pip,maxX)
    # setup micromanipulation 
    elseif cmp(simType,"MM") == 0
        mm = setup_micromanipulation(nuc)
        nuclearLength = zeros(Float64,maxT+1)
        ext = (mm,nuclearLength)
    elseif cmp(simType,"PC") == 0
        plane = spar.freeNucleusRadius + spar.repulsionDistance;
        ext = (plane)
    elseif cmp(simType,"INIT") == 0
        ext = ()
    end

    return nuc, chro, spar, simset,ext

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
    for i = 1:size(vert,2)
        nuc.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    VTKCelldata = get_cells(vtk)
    tri = VTKCelldata.connectivity

    # convert data to the required format
    tri = reshape(tri,(3,:))
    tri = tri' .+ 1
    nuc.tri = tri

    nuc = get_edges(nuc)
    nuc = get_vertex_triangles(nuc)
    
    return nuc,importNumber
end

function import_chromatin(chro,importFolder,importNumber)

    vtk = VTKFile(importFolder*"\\chro_"*importNumber*".vtp")
    vert = get_points(vtk)

    chro.vert = Vector{Vec{Float64,3}}(undef,size(vert)[2])
    for i = 1:size(vert,2)
        chro.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    endPoints = get_primitives(vtk,"Lines").offsets

    chro.strandIdx = Vector{Vector{Int64}}(undef,length(endPoints))
    for i = 1:length(endPoints)
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

function get_crosslinks!(chro,simset,spar)

    # remove crosslinks
    nLinked = length(chro.crosslinks)
    probs = rand(nLinked)
    for i = nLinked:-1:1
        if probs[i] < 0.001

            chro.crosslinked[chro.crosslinks[i][1]] = 0
            chro.crosslinked[chro.crosslinks[i][2]] = 0

            chro.crosslinks = chro.crosslinks[1:end .!= i]

        end
    end
    # form crosslinks
    notCrosslinked = findall(chro.crosslinked .== 0)
    closestVerts = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)
    possiblyLinking =  zeros(Bool,spar.chromatinLength*spar.chromatinNumber)

    for i = notCrosslinked
        # closestVertsTemp = inrange(simset.chromatinTree, chro.vert[i], 0.5, true)
        closest,distance = knn(simset.chromatinTree, chro.vert[i],1,true,j -> any(j .== i))
        if distance[1] <= 0.5 
            closestVerts[i] = closest[1]
            possiblyLinking[i] = true
        end
    end
    
    possibleLinkingIdx = findall(possiblyLinking)
    for i = possibleLinkingIdx

        if chro.crosslinked[i] == 0 && chro.crosslinked[closestVerts[i]] == 0

            if rand() < 0.0001

                push!(chro.crosslinks, [i, closestVerts[i]])
                chro.crosslinked[i] = 1
                chro.crosslinked[closestVerts[i]] = 1

            end
        end
    end



end