using NativeFileDialog, Plots, DelimitedFiles, ReadVTK, Meshes, LinearAlgebra, ProgressMeter


function analyze_volume()

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
        tri2 = tri2' .+ 1

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
        end

        next!(prog, showvalues = [(:index,i)])

    end


    p = Plots.plot([volumes areas],layout = grid(2,1), title = ["Volume (µm³)"     "Area (µm²)"],legend = false, linewidth = 5)

    return volumes,areas,p
end