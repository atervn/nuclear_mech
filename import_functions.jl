function import_envelope(nuc,importFolder)

    importNumber = get_import_number(importFolder)
    
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
    
    return nuc
end

function import_chromatin(chro,importFolder)

    importNumber = get_import_number(importFolder)

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

function import_crosslinks(chro,importFolder,spar)

    importNumber = get_import_number(importFolder)

    tempCrosslinks = []
    try
        tempCrosslinks = readdlm(importFolder*"\\crosslinks_" * importNumber * ".csv", ',', Int64, '\n')
        println("blob1")
    catch
    end
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

function get_import_number(importFolder)

    files = readdir(importFolder)
    ifNucFile = zeros(Bool,length(files))
    for i = eachindex(ifNucFile)
        ifNucFile[i] = cmp(files[i][1:5],"nucl_") == 0
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

    return lpad(maximum(timePointNumbers),numOfDigitsInName,"0")

end