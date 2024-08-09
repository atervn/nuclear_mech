
include("../analysis/analysis_functions.jl")


folder = pick_folder(pwd()*"\\results")

if folder == ""
    return
end

folder = splitdir(folder)[end]

numberOfNumbers = 0;
while true
    try
        parse(Int,folder[end-numberOfNumbers:end])
        global numberOfNumbers += 1
    catch
        break
    end


end

if numberOfNumbers == 0
    println("No number at the end of the directory name.")
    return
end

namePart = folder[1:end-numberOfNumbers];

lName = length(namePart)

allDirs = readdir(".//results")

temp = zeros(Int64,length(allDirs))

for i = eachindex(temp)

    if length(allDirs[i]) >= lName && allDirs[i][1:lName] == namePart
        temp[i] = 1
    end

end

numberOfSims = sum(temp)

outputName = "./afm_analysis/"*namePart[19:end]

Threads.@threads for i = 1:numberOfSims

    f,d = analyze_afm(;folder = ".\\results\\"*namePart*lpad(i,2,'0'))

    writedlm(outputName*lpad(i,2,'0')*".csv",[d[2:end] f[2:end]])

end