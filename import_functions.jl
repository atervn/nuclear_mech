function import_envelope(nuc,importFolder)

    importNumber = get_import_number(importFolder)
    
    importName = "nucl_"*importNumber

    vtk = VTKFile(importFolder*"\\"*importName*".vtu")

    vert = get_points(vtk)

    nuc.enve.vert = Vector{Vec{Float64,3}}(undef,size(vert)[2])
    for i = eachindex(vert[1,:])
        nuc.enve.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    VTKCelldata = get_cells(vtk)
    tri = VTKCelldata.connectivity

    # convert data to the required format
    tri = reshape(tri,(3,:))
    tri = tri' .+ 1
    nuc.enve.tri = Vector{Vector{Int64}}(undef,size(tri,1))
    for i = eachindex(tri[:,1])
        nuc.enve.tri[i] = tri[i,:]
    end

    nuc.enve = get_edges(nuc.enve)
    nuc.enve = get_vertex_triangles(nuc.enve)
    
    return nuc
end

function import_chromatin(nuc,importFolder)

    importNumber = get_import_number(importFolder)

    vtk = VTKFile(importFolder*"\\chro_"*importNumber*".vtp")
    vert = get_points(vtk)

    nuc.chro.vert = Vector{Vec{Float64,3}}(undef,size(vert)[2])
    for i = eachindex(vert[1,:])
        nuc.chro.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    endPoints = get_primitives(vtk,"Lines").offsets

    nuc.chro.strandIdx = Vector{Vector{Int64}}(undef,length(endPoints))
    for i = eachindex(endPoints)
        if i == 1
            nuc.chro.strandIdx[i] = collect(1:endPoints[i])
        else
            nuc.chro.strandIdx[i] = collect(endPoints[i-1] + 1:endPoints[i])
        end
    end

    chromatinLength = endPoints[1]
    chromatinNumber = length(endPoints);

    nuc.chro.forces.linear = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    nuc.chro.forces.bending = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    nuc.chro.forces.chroRepulsion = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    nuc.chro.forces.enveRepulsion = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    nuc.chro.forces.ladChroForces = Vector{Vec{3,Float64}}(undef,chromatinNumber*chromatinLength)
    nuc.chro.vectors = Vector{Vector{Vec{3,Float64}}}(undef,chromatinNumber)
    nuc.chro.vectorNorms = Vector{Vector{Float64}}(undef,chromatinNumber)
    nuc.chro.forces.strandLinear = Vector{Any}(undef,chromatinNumber)
    nuc.chro.forces.strandBending = Vector{Any}(undef,chromatinNumber)
    nuc.chro.forces.strandChroRepulsion = Vector{Any}(undef,chromatinNumber)
    nuc.chro.forces.strandEnveRepulsion = Vector{Any}(undef,chromatinNumber)
    nuc.chro.forces.strandLadChroForces = Vector{Any}(undef,chromatinNumber)
    initialize_chromatin_forces!(chro);
    
    nuc.chro.strandVert = Vector{Any}(undef,chromatinNumber)

    for i = 1:chromatinNumber
        nuc.chro.strandVert[i] = @view nuc.chro.vert[nuc.chro.strandIdx[i]];
        nuc.chro.vectors[i] = nuc.chro.strandVert[i][2:end] .- nuc.chro.strandVert[i][1:end-1];
        nuc.chro.vectorNorms[i] = norm.(nuc.chro.vectors[i])
        nuc.chro.forces.strandLinear[i] = @view nuc.chro.forces.linear[nuc.chro.strandIdx[i]];
        nuc.chro.forces.strandCrosslink[i] = @view nuc.chro.forces.crosslink[nuc.chro.strandIdx[i]];
        nuc.chro.forces.strandBending[i] = @view nuc.chro.forces.bending[nuc.chro.strandIdx[i]];
        nuc.chro.forces.strandChroRepulsion[i] = @view nuc.chro.forces.chroRepulsion[nuc.chro.strandIdx[i]];
        nuc.chro.forces.strandEnveRepulsion[i] = @view nuc.chro.forces.enveRepulsion[nuc.chro.strandIdx[i]];
        nuc.chro.forces.strandLadChroForces[i] = @view nuc.chro.forces.ladChroForces[nuc.chro.strandIdx[i]];
    end

    return nuc

end

function import_lads(nuc,importFolder,spar)

    nuc.enve.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    nuc.chro.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    tempLads = readdlm(importFolder*"\\lads.csv", ',', Int64, '\n')

    for i = 1:spar.chromatinNumber

        tempIdx = findall(tempLads[:,1] .== i)
        nuc.enve.lads[i] = tempLads[tempIdx,2]
        nuc.chro.lads[i] = tempLads[tempIdx,3]

    end

    return nuc
end

function import_crosslinks(nuc,importFolder,spar)

    importNumber = get_import_number(importFolder)

    tempCrosslinks = []
    try
        tempCrosslinks = readdlm(importFolder*"\\crosslinks_" * importNumber * ".csv", ',', Int64, '\n')
        println("blob1")
    catch
    end
    nuc.chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

    for i = eachindex(tempCrosslinks[:,1])

        push!(nuc.chro.crosslinks, tempCrosslinks[i,:])
        nuc.chro.crosslinked[tempCrosslinks[i,:]] .= 1

    end

    for i = 1:spar.chromatinNumber

        nuc.chro.crosslinked[nuc.chro.strandIdx[i][nuc.chro.lads[i]]] .= -1

    end

    return nuc

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