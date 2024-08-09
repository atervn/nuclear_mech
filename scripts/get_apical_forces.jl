using NativeFileDialog, ZipFile, Statistics, ReadVTK, PlotlyJS

include("../functions/utils.jl")


file = pick_file(pwd()*".\\results";filterlist="zip")

if file == ""
    return
end

file = splitdir(file)[end][1:end-4]

numberOfNumbers = 0;
while true
    try
        parse(Int,file[end-numberOfNumbers:end])
        global numberOfNumbers += 1
    catch
        break
    end


end

if numberOfNumbers == 0
    println("No number at the end of the directory name.")
    return
end

namePart = file[1:end-numberOfNumbers];

lName = length(namePart)

allContent = readdir(".//results")

temp = zeros(Int64,length(allContent))

for i = eachindex(temp)

    if length(allContent[i]) >= lName && allContent[i][1:lName] == namePart
        temp[i] = 1
    end

end

numberOfSims = sum(temp)

triangleAreasMeans = zeros(Float64,numberOfSims)
triangleAreasStds = zeros(Float64,numberOfSims)

justVoronois = zeros(Float64,8)

allForces = Vector{Vector{Float64}}(undef,12)

offset = 8

#Threads.@threads 
for i = 1:8

    unzip("./results/"*namePart*lpad(i+offset,numberOfNumbers,'0')*".zip")

    vtk = VTKFile("./results/"*namePart*lpad(i+offset,numberOfNumbers,'0')*"/enve_0001.vtu")

    # get the points
    vert = get_points(vtk)'
    
    middleZ = mean(vert[:,3])
    
    thePoints = middleZ .< vert[:,3]
    
    thePoints = thePoints .&& sqrt.(vert[:,1].^2 + vert[:,2].^2) .< 6

    println(sum(thePoints))

    pointData = get_point_data(vtk)
    elasticForces =  get_data(pointData["Elastic forces"])[:,thePoints]

    allForces[i] = sqrt.(elasticForces[1,:].^2 + elasticForces[2,:].^2 + elasticForces[3,:].^2)
        
    rm("./results/"*namePart*lpad(i+offset,numberOfNumbers,'0'); recursive = true)

end