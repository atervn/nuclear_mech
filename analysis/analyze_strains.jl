using NativeFileDialog, Plots, DelimitedFiles, ReadVTK, Meshes, LinearAlgebra, ProgressMeter, Statistics


function analyze_strains()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    files = readdir(folder)
    ifNucFile = zeros(Bool,length(files))
    for i = eachindex(ifNucFile)
        ifNucFile[i] = cmp(files[i][1:5],"enve_") == 0
    end

    nucFileIdx = findall(ifNucFile)

    numTimePoints = length(nucFileIdx)

    numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

    means = zeros(Float64,numTimePoints)

    stds = zeros(Float64,numTimePoints)

    mainStrains = Vector{Vector{Float64}}(undef,4);

    normalLengths = readdlm(folder*"\\normalLengths.csv")[:,1]./1e-6

    prog = Progress(length(means), 0.01, "Progress:", 100)

    mainIdx = [1, floor(length(means)/3), floor(2*length(means)/3), length(means)]

    mainIdxTemp = 1;

    for i = eachindex(means)

        tempNum = [filter(isdigit, collect(s)) for s in files[nucFileIdx[i]]][end-(numOfDigitsInName+3):end-4]

        numString = ""
        for j = 1:numOfDigitsInName
            numString = numString*string(tempNum[j][1])
        end

        timePointNumber = parse(Int64,numString)
        importNumber = lpad(timePointNumber,numOfDigitsInName,"0")

        vtk = VTKFile(folder*"\\enve_" * importNumber * ".vtu")

        vert2 = get_points(vtk)

        vert = Vector{Any}(undef,size(vert2)[2])
        for i = eachindex(vert2[1,:])
            vert[i] = Vec(vert2[1,i],vert2[2,i],vert2[3,i])
        end

        VTKCelldata = get_cells(vtk)
        tri2 = VTKCelldata.connectivity

        tri2 = reshape(tri2,(3,:))
        tri2 = tri2' .+ 1

        tri = Vector{Any}(undef,size(tri2,1))
        for i = eachindex(tri2[:,1])
            tri[i] = tri2[i,:]
        end

        
        # init edges vector
        edges = Vector{Vector{Int64}}(undef,0);
        
        # init neighbors vector
        neighbors = fill(Int[], length(vert));

        # create matrix of triangles
        triMatrix = [getindex.(tri,1) getindex.(tri,2) getindex.(tri,3)];

        # iterate through vertices
        for i = 1:length(vert)
            
            # find triangles that contain vertex i
            hasVertex = findall(triMatrix .== i);
            
            # init neighbors vector
            neighborsTemp = Array{Int64}(undef, 0);
            
            # iterate through triangles that contain vertex i
            for j = eachindex(hasVertex)

                # add other 2 vertices in triangle to neighbors vector
                append!(neighborsTemp,tri[hasVertex[j][1]][tri[hasVertex[j][1]] .!= i ]);

            end
            
            # get unique neighbors and save them
            unique!(neighborsTemp)
            neighbors[i] = neighborsTemp;

            # add connections to edges vector
            for j = eachindex(neighborsTemp)
                push!(edges, [i, neighborsTemp[j]]);
            end  
        end

        next!(prog, showvalues = [(:index,i)])

        edgeLengths = zeros(Float64,length(edges))

        for i = 1:length(edges)
            edgeLengths[i] = norm(vert[edges[i][1]] - vert[edges[i][2]])
        end

        strain = (edgeLengths .- normalLengths)./normalLengths

        if any(i .== mainIdx)
            mainStrains[mainIdxTemp] = strain
            mainIdxTemp += 1
        end

        means[i] = mean(strain)
        stds[i] = std(strain)

    end

    return means,stds,mainStrains

end