using NativeFileDialog, Plots, DelimitedFiles, ReadVTK

function analyze_strain(yMax = 0)

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    nuclearLengths = []
    try
        nuclearLengths = readdlm(folder * "\\nuclearLength.csv");
    catch


        files = readdir(folder)
        ifNucFile = zeros(Bool,length(files))
        for i = eachindex(ifNucFile)
            ifNucFile[i] = cmp(files[i][1:4],"nucl") == 0
        end

        nucFileIdx = findall(ifNucFile)

        numTimePoints = length(nucFileIdx)

        numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

        nuclearLengths = zeros(Float64,numTimePoints)

        for i = eachindex(nuclearLengths)

            tempNum = [filter(isdigit, collect(s)) for s in files[nucFileIdx[i]]][end-(numOfDigitsInName+3):end-4]

            numString = ""
            for j = 1:numOfDigitsInName
                numString = numString*string(tempNum[j][1])
            end

            timePointNumber = parse(Int64,numString)
            importNumber = lpad(timePointNumber,numOfDigitsInName,"0")

            vtk = VTKFile(folder*"\\nucl_" * importNumber * ".vtu")

            vert = get_points(vtk)

            nuclearLengths[i] = maximum(vert[1,:]) - minimum(vert[1,:])

        end
    end
    strain = (nuclearLengths .- nuclearLengths[1])./nuclearLengths[1];

    if yMax == 0
        p = plot(strain,xaxis = "Time step", yaxis = "Strain",legend = false, linewidth = 5)
    else
        p = plot(strain,xaxis = "Time step", yaxis = "Strain",ylim = (0,yMax), legend = false, linewidth = 5)
    end

    println("Final strain: " * string(strain[end]))

    return p
end