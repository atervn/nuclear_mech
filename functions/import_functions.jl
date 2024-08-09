function import_envelope(enve,importFolder,importTime,ipar)

    # get the import number
    importNumber = get_import_number(importFolder,importTime)
    
    # get the import name
    importName = "enve_"*importNumber

    # load a VTK file object
    vtk = VTKFile(importFolder*"\\"*importName*".vtu")

    # get the points
    vert = get_points(vtk)

    # create a new array to store the vertices
    enve.vert = Vector{Vec{3,Float64}}(undef,size(vert)[2])

    # iterate over the vertices
    for i = eachindex(vert[1,:])

        # set the vertex coordinates
        enve.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])

    end

    # get the cells
    VTKCelldata = get_cells(vtk)

    # get the connectivity
    tri = VTKCelldata.connectivity

    # convert data to the required format
    tri = reshape(tri,(3,:))
    if minimum(tri) == 0
        tri = tri' .+ 1
    else
        tri = tri'
    end

    # create a new array to store the triangles
    enve.tri = Vector{Vector{Int64}}(undef,size(tri,1))

    # iterate over the triangles
    for i = eachindex(tri[:,1])

        # set the triangle vertices
        enve.tri[i] = tri[i,:]

    end

    # get the edges
    enve = get_edges(enve)

    # get the vertex triangles
    enve = get_vertex_triangles(enve)
    
    # set the normal properties
    enve.normalArea = readdlm(importFolder*"\\normalArea.csv")[1]/ipar.scalingLength^2
    enve.normalAngle = readdlm(importFolder*"\\normalAngle.csv")[1]
    enve.normalTriangleAreas = readdlm(importFolder*"\\normalTriangleAreas.csv")[:,1]./ipar.scalingLength^2
    enve.normalVolume = readdlm(importFolder*"\\normalVolume.csv")[1]/ipar.scalingLength^3
    try
        enve.normalLengths = readdlm(importFolder*"\\normalLengths_"* importNumber *".csv")[:,1]./ipar.scalingLength
        enve.normalTriangleAreas = readdlm(importFolder*"\\normalTriangleAreas_"* importNumber *".csv")[:,1]./ipar.scalingLength^2
        enve.normalArea = readdlm(importFolder*"\\normalArea_"* importNumber *".csv")[1]/ipar.scalingLength^2
    catch
        enve.normalLengths = readdlm(importFolder*"\\normalLengths.csv")[:,1]./ipar.scalingLength
        enve.normalTriangleAreas = readdlm(importFolder*"\\normalTriangleAreas.csv")[:,1]./ipar.scalingLength^2
        enve.normalArea = readdlm(importFolder*"\\normalArea.csv")[1]/ipar.scalingLength^2
    end

    enve.envelopeMultipliers = readdlm(importFolder*"\\envelope_multipliers.csv")[:]

    return enve
end

function import_chromatin(chro,spar,importFolder,importTime)

    # get the import number from the importFolder
    importNumber = get_import_number(importFolder,importTime)

    # load the VTK file
    vtk = VTKFile(importFolder*"\\chro_"*importNumber*".vtp")

    # get the vertex coordinates
    vert = get_points(vtk)

    # create a vector to store the vertex coordinates
    chro.vert = Vector{Vec{3,Float64}}(undef,size(vert)[2])

    # iterate over the vertex coordinates and add them to the chro.vert vector
    for i = eachindex(vert[1,:])
        chro.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
    end

    # get the strand end indices
    endPoints = get_primitives(vtk,"Lines").offsets

    # create a vector to store the strand indices
    chro.strandIdx = Vector{Vector{Int64}}(undef,length(endPoints))

    # iterate over the strands and add them to the chro.strandIdx vector
    for i = eachindex(endPoints)
        if i == 1
            chro.strandIdx[i] = collect(1:endPoints[i])
        else
            chro.strandIdx[i] = collect(endPoints[i-1] + 1:endPoints[i])
        end
    end

    # initialize the chromatin properties
    chro = initialize_chromatin_properties(chro,spar)

    return chro

end

function import_lads(enve,chro,spar,importFolder)

    # initialize the LADs vectors
    enve.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);
    chro.lads = Vector{Vector{Int64}}(undef, spar.chromatinNumber);

    # read the LADs CSV file
    tempLads = readdlm(importFolder*"\\lads.csv", ',', Int64, '\n')

    # iterate over the LADs and add them to the LADs vectors
    for i = 1:spar.chromatinNumber

        # find the indices of the LADs for chromosome i
        tempIdx = findall(tempLads[:,1] .== i)
        
        # add the LADs to the enve.lads and chro.lads vectors
        enve.lads[i] = tempLads[tempIdx,2]
        chro.lads[i] = tempLads[tempIdx,3]

    end

    return enve,chro

end

function import_crosslinks(chro,spar,importFolder,importTime)

    # get the import number from the importFolder
    importNumber = get_import_number(importFolder,importTime)

    # initialize the crosslinks vector
    chro.crosslinked = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)

    # try to read the crosslinks CSV file
    tempCrosslinks = []

    try
        vtk = VTKFile(importFolder*"\\crosslinks_" * importNumber * ".vtp")
        a = get_point_data(vtk)
        aa = get_data(a["crosslink ID"])
        tempCrosslinks = [aa[1:Int(length(aa)/2)] aa[Int(length(aa)/2)+1:end]]
    catch
    end
    
    # iterate over the crosslinks and add them to the crosslinks vector
    for i = eachindex(tempCrosslinks[:,1])

        # add the crosslink to the crosslinks vector
        push!(chro.crosslinks, tempCrosslinks[i,:])

        # set the state of the crosslink to 1
        chro.crosslinked[tempCrosslinks[i,:]] .= 1

    end

    # set the crosslinks in the LADs to -1
    for i = 1:spar.chromatinNumber
        chro.crosslinked[chro.strandIdx[i][chro.lads[i]]] .= -1
    end

    return chro

end

function get_import_folder(initType,importFolder)
    
    # get the type of simulation
    if initType == "load"

        # if no import folder was provided, open a dialog to get one
        if importFolder == "" # open dialog if no folder was provided
            importFolder = pick_folder(pwd()*"\\results")
        else

            # if an import folder was provided, get its path
            importFolder = pwd()*"\\results\\"*importFolder;
        end
    else

        # if the simulation type is not "load", set the import folder to an empty string
        importFolder = ""

    end

    return importFolder

end

function get_import_number(importFolder,importTime)

    # get all the files in the import folder
    files = readdir(importFolder)

    # create a boolean array to indicate whether each file is a nuclear file
    isNuclFile = zeros(Bool,length(files))
    for i = eachindex(isNuclFile)
        isNuclFile[i] = cmp(files[i][1:5],"enve_") == 0
    end

    # get the indices of the nuclear files
    nucFileIdx = findall(isNuclFile)

    # get the number of time points
    numTimePoints = length(nucFileIdx)

    # get the number of digits in the name of the first nuclear file
    numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

    # create an array to store the time point numbers
    timePointNumbers = zeros(Int64,numTimePoints)

    # for each nuclear file, get the time point number string
    for i = eachindex(timePointNumbers)

        tempNum = [filter(isdigit, collect(s)) for s in files[nucFileIdx[i]]][end-(numOfDigitsInName+3):end-4]

        numString = ""
        for j = 1:numOfDigitsInName
            numString = numString*string(tempNum[j][1])
        end

        timePointNumbers[i] = parse(Int64,numString)
    end

    # get the maximum time point number
    if importTime == 0
        timePointNumber = maximum(timePointNumbers)
    else
        timePointNumber = importTime
    end

    return lpad(timePointNumber,numOfDigitsInName,"0")

end

function import_replication_compartment(repl,spar,importFolder,importTime)

     # get the import number
     importNumber = get_import_number(importFolder,importTime)
    
     # get the import name
     importName = "repl_"*importNumber
 
     # load a VTK file object
     vtk = VTKFile(importFolder*"\\"*importName*".vtu")
 
     # get the points
     vert = get_points(vtk)
 
     # create a new array to store the vertices
     repl.vert = Vector{Vec{3,Float64}}(undef,size(vert)[2])
 
     # iterate over the vertices
     for i = eachindex(vert[1,:])
 
         # set the vertex coordinates
         repl.vert[i] = Vec(vert[1,i],vert[2,i],vert[3,i])
 
     end
 
     # get the cells
     VTKCelldata = get_cells(vtk)
 
     # get the connectivity
     tri = VTKCelldata.connectivity
 
     # convert data to the required format
     tri = reshape(tri,(3,:))
     if tri[1,1] == 0
        tri = tri' .+ 1
     else
        tri = tri'
     end
 
     # create a new array to store the triangles
     repl.tri = Vector{Vector{Int64}}(undef,size(tri,1))
 
     # iterate over the triangles
     for i = eachindex(tri[:,1])
 
         # set the triangle vertices
         repl.tri[i] = tri[i,:]
 
     end
 
     # get the edges
     repl = get_edges(repl)
 
     # get the vertex triangles
     repl = get_vertex_triangles(repl)
    
     repl.initVolume = readdlm(importFolder*"\\replInitVolume.csv")[1]/spar.scalingLength^3

     return repl

end