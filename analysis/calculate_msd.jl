using DelimitedFiles, Plots, Statistics, NativeFileDialog, ReadVTK
# theme(:ggplot2)

include("../functions/setup_functions.jl")

if !(@isdefined envelopeType)
    include("../functions/NuclearMechTypes.jl")
    using .NuclearMechTypes
end

function calculate_msd()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    ipar = read_parameters(folder*"\\parameters.txt")

    allFiles = readdir(folder)
    
    temp = zeros(Int64,length(allFiles))

    for i = eachindex(temp)

        if allFiles[i][1:4] == "chro"
            temp[i] = 1
        end

    end

    chroFiles = allFiles[findall(temp .== 1)]

    vtk = VTKFile(folder*"\\"*chroFiles[1])
    nVerts = size(get_points(vtk))[2]
    endPoints = get_primitives(vtk,"Lines").offsets
    startPoints = cat([1],endPoints[1:end-1] .+ 1, dims = 1)

    histoneInds = collect(1:endPoints[end])

    telomerInds = sort(cat(endPoints,startPoints,dims=1))
    deleteat!(histoneInds,telomerInds)

    telomeres = Vector{Matrix{Float64}}(undef,length(telomerInds))
    histones = Vector{Matrix{Float64}}(undef,length(histoneInds))


    for i = eachindex(chroFiles)

        vtk = VTKFile(folder*"\\"*chroFiles[i])
        vert = get_points(vtk)
    
        histInd = 1
        for ii = histoneInds

            if i == 1
                histones[histInd] = vert[1:2,ii]'
            else
                histones[histInd] = vcat(histones[histInd], vert[1:2,ii]')
            end
            histInd += 1

        end

        teloInd = 1
        for ii = telomerInds

            if i == 1
                telomeres[teloInd] = vert[1:2,ii]'
            else
                telomeres[teloInd] = vcat(telomeres[teloInd], vert[1:2,ii]')
            end
            teloInd += 1

        end
    end

    histoneMSD = zeros(Float64,length(chroFiles),length(histoneInds))
    for i = 1:length(histoneInds)

        histoneMSD[:,i] = (histones[i][1,1] .- histones[i][:,1]).^2 + (histones[i][1,2] .- histones[i][:,2]).^2

    end

    telomereMSD = zeros(Float64,length(chroFiles),length(telomerInds))
    for i = 1:length(telomerInds)

        telomereMSD[:,i] = (telomeres[i][1,1] .- telomeres[i][:,1]).^2 + (telomeres[i][1,2] .- telomeres[i][:,2]).^2

    end

    telomereMeans = mean(telomereMSD,dims=2) 
    telomereStds = std(telomereMSD,dims=2) 
    histoneMeans = mean(histoneMSD,dims=2) 
    histoneStds = std(histoneMSD,dims=2) 

    times = 0:ipar.dt*ipar.exportStep:(length(telomereMeans)-1)*ipar.dt*ipar.exportStep

    p = plot(times,telomereMeans,ribbon=telomereStds,alpha=0.8;linewidth = 4,ylims=(-0.2, 0.5), legend=false)
    plot!(p,times,histoneMeans,ribbon=telomereStds,alpha=0.8;linewidth = 4)
    annotate!(p,[(1,0.48, ("Telomeres: "*string(round(mean(telomereMeans[end-5:end]);digits=3)), 12, :black, :left))])
    annotate!(p,[(1,0.44, ("Histones: "*string(round(mean(histoneMeans[end-5:end]);digits=3)), 12, :black, :left))])
    
    return p

end