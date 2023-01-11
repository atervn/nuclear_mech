using NativeFileDialog, Plots, DelimitedFiles, ReadVTK, LaTeXStrings

function analyze_aspiration()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    maxX = []
    try
        maxX = readdlm(folder * "\\maxX.csv");
    catch


        files = readdir(folder)
        ifNucFile = zeros(Bool,length(files))
        for i = eachindex(ifNucFile)
            ifNucFile[i] = cmp(files[i][1:4],"nucl") == 0
        end

        nucFileIdx = findall(ifNucFile)

        numTimePoints = length(nucFileIdx)

        numOfDigitsInName = sum(.!isempty.([filter(isdigit, collect(s)) for s in files[nucFileIdx[1]]]))

        maxX = zeros(Float64,numTimePoints)

        for i = eachindex(maxX)

            tempNum = [filter(isdigit, collect(s)) for s in files[nucFileIdx[i]]][end-(numOfDigitsInName+3):end-4]

            numString = ""
            for j = 1:numOfDigitsInName
                numString = numString*string(tempNum[j][1])
            end

            timePointNumber = parse(Int64,numString)
            importNumber = lpad(timePointNumber,numOfDigitsInName,"0")

            vtk = VTKFile(folder*"\\nucl_" * importNumber * ".vtu")

            vert = get_points(vtk)

            maxX[i] = maximum(vert[1,:])

        end
    end

    # prompt to input
    print("Aspiration pressure (in Pa): ") 
      
    # Calling rdeadline() function
    pressure = parse(Float64,readline())

    print("Pipette radius (in µm): ") 
      
    # Calling rdeadline() function
    pipetteRadius = parse(Float64,readline())

    print("Time step (in s): ") 
      
    # Calling rdeadline() function
    dt = parse(Float64,readline())

    dL = maxX .- minimum(maxX);

    J = 2*pi*2.1.*dL./(3*pressure*1e-3*pipetteRadius);

    p = plot(10*dt:dt:length(maxX)*dt-dt,J[11:end],
                yaxis=:log,
                xaxis=:log,
                xlabel="Time (s)",
                ylabel="Creep compliance kPa" * L"^{-1}",
                xlim = (0.1, 500),
                ylim = (0.01, 10))

    return p
end