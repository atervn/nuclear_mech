using NativeFileDialog, Plots, DelimitedFiles, ReadVTK, LaTeXStrings

include("setup_functions.jl")

if !(@isdefined envelopeType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

function analyze_aspiration(case)

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

    ipar = inputParametersType();
    read_parameters(ipar,folder*"\\parameters.txt")

    pipettePosition = ipar.pipetteRadius*tan(acos(ipar.pipetteRadius/(ipar.freeNucleusRadius + ipar.repulsionDistance)));

    if case == 1
        dL = maxX.*ipar.scalingLength .- pipettePosition;
    else
        dL = 0.001e-6.+(maxX .- minimum(maxX)).*ipar.scalingLength;
    end

    J = 2*pi*2.1.*dL./(3*ipar.aspirationPressure*1e-3*ipar.pipetteRadius);

    dt = ipar.maxDt

    p = plot(2*dt:dt:length(maxX)*dt-dt,J[3:end],
                yaxis=:log,
                xaxis=:log,
                xlabel="Time (s)",
                ylabel="Creep compliance kPa" * L"^{-1}",
    )
    
    dahl2005 = readdlm("dahl2005.txt",',');
    dahl2006 = readdlm("dahl2006.txt",',');
    witner2020 = readdlm("wintner2020.txt",',');
    guilak2000 = readdlm("guilak2000.txt",',');

    scatter!(p,dahl2005[:,1],dahl2005[:,2]; legend = false)
    scatter!(p,dahl2006[:,1],dahl2006[:,2]; legend = false)
    scatter!(p,witner2020[:,1],witner2020[:,2]; legend = false)

    return p
end