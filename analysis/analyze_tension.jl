using DelimitedFiles, Plots, Statistics, NativeFileDialog, ReadVTK
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

    ipar = read_parameters(folder*"\\parameters.txt")

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
        data =  get_data(point_data["Elastic forces"])

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


    times = 0:ipar.dt*ipar.exportStep:(length(means)-1)*ipar.dt*ipar.exportStep

    p = plot(times,means,ribbon=stds,alpha=0.8;linewidth = 4, legend=false)
    # plot!(p,times,histoneMeans,ribbon=telomereStds,alpha=0.8;linewidth = 4)
    # annotate!(p,[(1,0.48, ("Telomeres: "*string(round(mean(telomereMeans[end-5:end]);digits=3)), 12, :black, :left))])
    # annotate!(p,[(1,0.44, ("Histones: "*string(round(mean(histoneMeans[end-5:end]);digits=3)), 12, :black, :left))])
    
    return means,stds,mainTensions,p

end