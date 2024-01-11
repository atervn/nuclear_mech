using DelimitedFiles, Statistics, NativeFileDialog, ReadVTK
# theme(:ggplot2)

include("../functions/setup_functions.jl")

if !(@isdefined envelopeType)
    include("../functions/NuclearMechTypes.jl")
    using .NuclearMechTypes
end

function analyze_tension()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    allFiles = readdir(folder)
    
    temp = zeros(Int64,length(allFiles))

    for i = eachindex(temp)

        if allFiles[i][1:5] == "enve_"
            temp[i] = 1
        end

    end

    enveFiles = allFiles[findall(temp .== 1)]

    vtk = VTKFile(folder*"\\"*enveFiles[1])
    point_data = get_point_data(vtk)
    elasticForces = point_data["Elastic forces"]
    data = get_data(elasticForces)

    means = Vector{Float64}(undef,length(enveFiles))
    stds = Vector{Float64}(undef,length(enveFiles))

    mainIdx = [1, floor(length(means)/3), floor(2*length(means)/3), length(means)]
    mainTensions = Vector{Vector{Float64}}(undef,4);

    mainIdxTemp = 1;

    for i = 1:length(means)

        vtk = VTKFile(folder*"\\"*enveFiles[i])
        
        point_data = get_point_data(vtk)
        data =  get_data(point_data["Volume forces"])

        values = zeros(size(data)[2])

        for j = 1:size(data)[2]
            values[j] = norm(data[:,j])
        end

        means[i] = mean(values)
        stds[i] = std(values)

        if any(i .== mainIdx)
            mainTensions[mainIdxTemp] = values
            mainIdxTemp += 1
        end

    end


    times = 1:length(means)

traces=[
    scatter(
        x=times,
        y=means,
        line=attr(color="rgb(255,0,0)"),
        mode="lines"
    ),

    scatter(
        x=vcat(times, reverse(times)), # x, then x reversed
        y=vcat(means.+stds, reverse(means.-stds)), # upper, then lower reversed
        fill="toself",
        fillcolor="rgba(255,0,0,0.2)",
        line=attr(color="rgba(255,255,255,0)"),
        hoverinfo="skip",
        showlegend=false
    )
]


    p = plot(traces)

    return means,stds,mainTensions,p

end