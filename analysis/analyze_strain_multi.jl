function analyze_strain_multi()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    folder = splitdir(folder)[end]

    numberOfNumbers = 0;
    while true
        try
            parse(Int,folder[end-numberOfNumbers:end])
            numberOfNumbers += 1
        catch
            break
        end


    end

    if numberOfNumbers == 0
        println("No number at the end of the directory name.")
        return
    end

    namePart = folder[1:end-numberOfNumbers];

    allDirs = readdir(".//results")

    temp = zeros(Int64,length(allDirs))

    for i = eachindex(temp)

        if allDirs[i][1:end-numberOfNumbers] == namePart
            temp[i] = 1
        end

    end

    numberOfSims = sum(temp)

    strains = zeros(Float64,numberOfSims)

    for i = eachindex(strains)

        folder = pwd()*"\\results\\"*namePart*string(i,pad = numberOfNumbers)

        if !isfile(folder*"\\failed.txt")

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
        
            strains[i] = ((nuclearLengths .- nuclearLengths[1])./nuclearLengths[1])[end];
        end
    end

    return strains

end