using NativeFileDialog, PlotlyJS, DelimitedFiles, ReadVTK, Meshes, LinearAlgebra, ProgressMeter

function analyze_volume_infection()

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

    volumes = zeros(Float64,numTimePoints)

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

        try
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

            volumes2 = zeros(length(tri));
        
            for j = eachindex(tri)
                volumes2[j] = 1/6*dot(vert[tri[j][1]],cross(vert[tri[j][2]],vert[tri[j][3]]))
            end

            volumes[i] = sum(volumes2)

        catch
            volumes[i] = volumes[i-1]
        end

        if repl
            
            try
                vtk = VTKFile(folder*"\\repl_" * importNumber * ".vtu")

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

                volumes2 = zeros(length(tri));

                for j = eachindex(tri)
                    volumes2[j] = 1/6*dot(vert[tri[j][1]],cross(vert[tri[j][2]],vert[tri[j][3]]))
                end

                
                replicationVolumes[i] = sum(volumes2)
            catch
                replicationVolumes[i] = replicationVolumes[i-1]
            end
            
        end

        next!(prog, showvalues = [(:index,i)])

    end

    chromatinVolume = volumes - replicationVolumes;

    times = 1:length(chromatinVolume)

    p = PlotlyJS.plot([PlotlyJS.scatter(x = times,y=volumes,name="Nuclear volume (µm³)"),PlotlyJS.scatter(x = times,y=replicationVolumes,name="VRC volume (µm³)"),PlotlyJS.scatter(x=times,y=chromatinVolume,name = "Difference (µm³)")], Layout(xaxis_title_text="Time",yaxis_title_text="Volume",xaxis_titlefont_size = 16,yaxis_titlefont_size = 16,legend_font_size = 16,legend=attr(x=0.5,y=1.0,xanchor = "center",yanchor="bottom",orientation="h")))

    return p
end







function analyze_afm()

    folder = pick_folder(pwd()*"\\results")

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
        # println(i)
        # point_data = get_point_data(vtk)
        # println(point_data)
        # println(point_data["Force on bead"])
        # temp = get_data(point_data["Force on bead"])
        # println(temp)
        # forces2[i] = temp[1]


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

    F = 0.05.*abs.(deltaX)

    enveTop .-= enveTop[1]

    depth = abs.(enveTop).*1e-6

    return F,depth,enveTop


end