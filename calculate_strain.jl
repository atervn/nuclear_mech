using NativeFileDialog, Plots, DelimitedFiles

function calculate_strain(yMax = 0)

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    nuclearLengths = []
    try
        nuclearLengths = readdlm(folder * "\\nuclearLength.csv");
    catch
        println("no nuclearLength.csv found.")
        return
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