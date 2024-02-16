using NativeFileDialog, PlotlyJS, DelimitedFiles, ReadVTK, Meshes, LinearAlgebra, ProgressMeter, Statistics


"""
Plot nuclear and area during a simulation.

syntax: plot = analyze_volume()
returns a PlotlyJS plot object
"""
function analyze_volume(;n = 1, sameName = true)

    traces1 = GenericTrace[]
    traces2 = GenericTrace[]
    traces3 = GenericTrace[]

    isRepl = false

    folder = "";
    folderBase = "";

    finalValues = []

    volumes = zeros(Float64,0)

    areas = zeros(Float64,0)


    padding = 0;

    for i = 1:n


        if i == 1 || !sameName

            folder = pick_folder(pwd()*"\\results")

            if folder == ""
                return
            end

            if sameName

                subStrings = split(folder, '_')

                padding = length(subStrings[end])


                folderBase = folder[1:end-padding]

            end
        else

            folder = folderBase*lpad(i,padding,'0')

        end

        files = readdir(folder)
        ifNucFile = zeros(Bool,length(files))
        for i = eachindex(ifNucFile)
            ifNucFile[i] = cmp(files[i][1:5],"enve_") == 0
        end

        nucFileIdx = findall(ifNucFile)

        numTimePoints = length(nucFileIdx)

        numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

        volumes = zeros(Float64,numTimePoints)

        areas = zeros(Float64,numTimePoints)

        if isfile(folder*"\\repl_0001.vtu")
            replicationVolumes = zeros(Float64,numTimePoints)
            repl = true
        else
            repl = false
        end

        prog = Progress(length(volumes), 0.01, "Progress:", 100)

        for i = eachindex(volumes)

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
            if tri2[1,1] == 0
                tri2 = tri2' .+ 1
            else
                tri2 = tri2'
            end

            tri = Vector{Any}(undef,size(tri2,1))
            for i = eachindex(tri2[:,1])
                tri[i] = tri2[i,:]
            end

            volumes2 = zeros(length(tri));
            areas2 = zeros(length(tri));

            for j = eachindex(tri)
                volumes2[j] = 1/6*dot(vert[tri[j][1]],cross(vert[tri[j][2]],vert[tri[j][3]]))

                vect1 = vert[tri[j][2]] - vert[tri[j][1]]
                vect2 = vert[tri[j][3]] - vert[tri[j][1]]

                areas2[j] = 0.5*norm(cross(vect1,vect2))  

            end

            
            volumes[i] = sum(volumes2)
            areas[i] = sum(areas2)

            if repl
                
                vtk = VTKFile(folder*"\\repl_" * importNumber * ".vtu")

                vert2 = get_points(vtk)

                vert = Vector{Any}(undef,size(vert2)[2])
                for i = eachindex(vert2[1,:])
                    vert[i] = Vec(vert2[1,i],vert2[2,i],vert2[3,i])
                end

                VTKCelldata = get_cells(vtk)
                tri2 = VTKCelldata.connectivity

                tri2 = reshape(tri2,(3,:))
                if tri2[1,1] == 0
                    tri2 = tri2' .+ 1
                else
                    tri2 = tri2'
                end

                tri = Vector{Any}(undef,size(tri2,1))
                for i = eachindex(tri2[:,1])
                    tri[i] = tri2[i,:]
                end

                volumes2 = zeros(length(tri));

                for j = eachindex(tri)
                    volumes2[j] = 1/6*dot(vert[tri[j][1]],cross(vert[tri[j][2]],vert[tri[j][3]]))
                end

                
                replicationVolumes[i] = sum(volumes2)

                isRepl = true

            end

            next!(prog, showvalues = [(:index,i)])

            x = 1:length(volumes)

            trace1 = scatter(x=x, y=volumes, mode="lines", showlegend=false,line = attr(color = "black"))
            trace2 = scatter(x=x, y=areas, mode="lines", showlegend=false,line = attr(color = "black"))

            if isRepl
                trace3 = scatter(x=x, y=replicationVolumes, mode="lines", showlegend=false,line = attr(color = "black"))
                push!(traces3,trace3)
            end

            push!(traces1,trace1)
            push!(traces2,trace2)

        end
        
        push!(finalValues,volumes[end])

    end

    if !isRepl
        subplot1 = plot(traces1, Layout(title="Nuclear volume (µm³)"))
        subplot2 = plot(traces2, Layout(title="Nuclear area (µm²)"))

        p = [subplot1; subplot2]
    else
        subplot1 = plot(traces1, Layout(title="Nuclear volume (µm³)"))
        subplot2 = plot(traces3, Layout(title="VRC volume (µm³)"))
        subplot3 = plot(traces2, Layout(title="Nuclear area (µm²)"))
    
        p = [subplot1; subplot2; subplot3]
    end

    if n == 1

        println("Final volume: "*string(round(finalValues[1];digits = 2))*" µm³")

    else
        meanV = mean(finalValues)
        stdV = std(finalValues)

        println("Final volume: "*string(round(meanV;digits = 2))*" ± "*string(round(stdV;digits = 2))*" µm³")
        
    end


    return p

end

"""
Analyze AFM results.

syntax: force, depth = analyze_afm(folder)
if no folder is given, it will ask for one
"""
function analyze_afm(;folder = "")

    if folder == ""
        folder = pick_folder(pwd()*"\\results")
    end

    if folder == ""
        return
    end

    allFiles = readdir(folder)
    
    temp = zeros(Int64,length(allFiles))

    for i = eachindex(temp)
        if allFiles[i][1:4] == "afm_"
            temp[i] = 1
        end

    end

    afmFiles = allFiles[findall(temp .== 1)]

    beadPositions = zeros(Float64,length(afmFiles))
    topPositions = zeros(Float64,length(afmFiles))
    enveTop = zeros(Float64,length(afmFiles))
    forces2 = zeros(Float64,length(afmFiles))
    

    vtk = VTKFile(folder*"\\enve_0001.vtu")

    # get the points
    vert = get_points(vtk)

    distances = sqrt.(vert[1,:].^2 .+ vert[2,:].^2)

    sortedIdx = sortperm(distances)


    maxIdx = argmax(vert[3,sortedIdx[1:10]])

    maxIdx = sortedIdx[maxIdx]

    prog = Progress(length(afmFiles), 0.01, "Progress:", 100)

    for i = 1:length(afmFiles)


        importNumber = lpad(i,4,"0")

        # load the VTK file
        vtk = VTKFile(folder*"\\afm_"*importNumber*".vtu")

        # get the vertex coordinates
        vert1 = get_points(vtk)

        beadPositions[i] = vert1[3,1]
        topPositions[i] = vert1[3,2]

        try
            vtk2 = VTKFile(folder*"\\enve_"*importNumber*".vtu")

            # get the points
            vert2 = get_points(vtk2)

            enveTop[i] = vert2[3,maxIdx]
        catch
            enveTop[i] = enveTop[i-1]
        end

        next!(prog, showvalues = [(:index,i)])

    end

    # enveTop .-= enveTop[1]

    deltaX = ((topPositions[1] - beadPositions[1]) .- (topPositions .- beadPositions)).*1e-6

    force = 0.05.*abs.(deltaX)

    enveTop .-= enveTop[1]

    depth = abs.(enveTop).*1e-6

    return force, depth

end