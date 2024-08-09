using NativeFileDialog, ZipFile

include("../analysis/analysis_functions.jl")
include("../functions/utils.jl")

# allowedForces = ["Volume forces", "Area forces", "Elastic forces", "chroRepulsion forces", "Bending forces", "LAD forces", "AFM forces", "Cytoskeleton forces", "Total forces"];

forceName = "Bending forces";

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

outputName = "./afm_analysis/"*namePart[19:end]

means = zeros(Float64,length(allContent))
stds = zeros(Float64,length(allContent))

#Threads.@threads 
for i = 1:numberOfSims
    println(i)
    unzip("./results/"*namePart*lpad(i,numberOfNumbers,'0')*".zip")

    means[i],stds[i] = analyze_forces(forceName; folder = ".\\results\\"*namePart*lpad(i,numberOfNumbers,'0'), output = "first")

    rm("./results/"*namePart*lpad(i,numberOfNumbers,'0'); recursive = true)

end