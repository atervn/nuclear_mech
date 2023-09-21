using PlotlyJS,FileIO,NativeFileDialog,ReadVTK,DelimitedFiles

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
    vtk = VTKFile(folder*"\\enve_0001.vtu")

    # get the points
    vert = get_points(vtk)

    distances = sqrt.(vert[1,:].^2 .+ vert[2,:].^2)

    sortedIdx = sortperm(distances)


    maxIdx = argmax(vert[3,sortedIdx[1:10]])

    maxIdx = sortedIdx[maxIdx]

    println(maxIdx)

    for i = 1:length(afmFiles)


        importNumber = lpad(i,4,"0")

        # load the VTK file
        vtk = VTKFile(folder*"\\afm_"*importNumber*".vtu")

        # get the vertex coordinates
        vert1 = get_points(vtk)

        beadPositions[i] = vert1[3,1]
        topPositions[i] = vert1[3,2]

        println(importNumber)
        vtk2 = VTKFile(folder*"\\enve_"*importNumber*".vtu")

        # get the points
        vert2 = get_points(vtk2)

        enveTop[i] = vert2[3,maxIdx]


    end

    # enveTop .-= enveTop[1]

    deltaX = ((topPositions[1] - beadPositions[1]) .- (topPositions .- beadPositions)).*1e-6

    F = 0.05.*abs.(deltaX)

    enveTop .-= enveTop[1]

    depth = abs.(enveTop).*1e-6

    # depth = abs.(enveTop).*1e-6

    # plotIdx = minimum(findall((depth) .> 1e-7))

    # E = F.*3/4*(1-0.5^2)./(sqrt(3.31e-6).*depth.^(3/2))

    # p = plot(E[plotIdx:end])
    # twinx()
    # plot!(p,depth[plotIdx:end]./3.31e-6.*100)

    # times = Float64.(collect(1:length(afmFiles)))

    # p = PlotlyJS.plot(

    #     [

    #         PlotlyJS.scatter(x=times,y=E[plotIdx:end], name="Stiffness"),

    #         PlotlyJS.scatter(x=times,y=depth[plotIdx:end]./3.31e-6, name="Intendation", yaxis="y2")

    #     ],

    #     Layout(

    #         xaxis_title_text="Time",

    #         yaxis_title_text="Stiffness (Pa)",

    #         yaxis2=attr(

    #             title="Intendation depth (µm)",

    #             overlaying="y",

    #             side="right"

    #         )

    #     )

    # )


    return F,depth


end