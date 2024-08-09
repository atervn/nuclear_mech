using DelimitedFiles, Statistics, NativeFileDialog, ReadVTK, PlotlyJS
# theme(:ggplot2)

include("../functions/setup_functions.jl")


function analyze_afm_forces()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    folder = splitdir(folder)[end]

    global numberOfNumbers = 0
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

    allDirs = readdir(".//results")

    temp = zeros(Int64,length(allDirs))

    for i = eachindex(temp)

        if allDirs[i][1:end-numberOfNumbers] == namePart
            temp[i] = 1
        end

    end

    numberOfSims = sum(temp)

    nTimeSteps = zeros(Int,numberOfSims)
    maxTimeSteps = 0;

    for i = 1:numberOfSims

        allFiles = readdir(".\\results\\"*namePart*lpad(i,numberOfNumbers,'0'))

        temp = zeros(Int64,length(allFiles))

        for i = eachindex(temp)
            if allFiles[i][1:4] == "afm_"
                temp[i] = 1
            end
        end

        nTimeSteps[i] = sum(temp)

        if maxTimeSteps < sum(temp)
            maxTimeSteps = sum(temp)
        end


    end

    forceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    areaForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    volumeForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    bendingForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    elasticForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    chroRepulsionForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    ladForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    cytoskeletonForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negareaForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negvolumeForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negbendingForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negelasticForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negchroRepulsionForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negladForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)
    negcytoskeletonForceOnBead = zeros(Float64,numberOfSims,maxTimeSteps)

    for i = 1:numberOfSims

        for j = 1:nTimeSteps[i]

            pointData = get_point_data(VTKFile(".\\results\\"*namePart*lpad(i,numberOfNumbers,'0')*"\\afm_"*lpad(j,4,'0')*".vtu"))

            forceOnBead[i,j] = get_data(pointData["Force on bead"])[1]
            areaForceOnBead[i,j] = get_data(pointData["Area force on bed"])[1]
            volumeForceOnBead[i,j] = get_data(pointData["Volume force on bed"])[1]
            bendingForceOnBead[i,j] = get_data(pointData["Bending force on bed"])[1]
            elasticForceOnBead[i,j] = get_data(pointData["Elastic force on bed"])[1]
            chroRepulsionForceOnBead[i,j] = get_data(pointData["chroRepulsion force on bed"])[1]
            ladForceOnBead[i,j] = get_data(pointData["LAD force on bed"])[1]
            cytoskeletonForceOnBead[i,j] = get_data(pointData["afm force on cell"])[1]
            negareaForceOnBead[i,j] = get_data(pointData["Negative Area force on bed"])[1]
            negvolumeForceOnBead[i,j] = get_data(pointData["Negative Volume force on bed"])[1]
            negbendingForceOnBead[i,j] = get_data(pointData["Negative Bending force on bed"])[1]
            negelasticForceOnBead[i,j] = get_data(pointData["Negative Elastic force on bed"])[1]
            negchroRepulsionForceOnBead[i,j] = get_data(pointData["Negative chroRepulsion force on bed"])[1]
            negladForceOnBead[i,j] = get_data(pointData["Negative LAD force on bed"])[1]
            negcytoskeletonForceOnBead[i,j] = get_data(pointData["Negative afm force on cell"])[1]
        end
    end

    forceOnBeadMean = zeros(maxTimeSteps)
    forceOnBeadStd = zeros(maxTimeSteps)

    areaForceOnBeadMean = zeros(maxTimeSteps)
    areaForceOnBeadStd = zeros(maxTimeSteps)

    volumeForceOnBeadMean = zeros(maxTimeSteps)
    volumeForceOnBeadStd = zeros(maxTimeSteps)

    bendingForceOnBeadMean = zeros(maxTimeSteps)
    bendingForceOnBeadStd = zeros(maxTimeSteps)

    elasticForceOnBeadMean = zeros(maxTimeSteps)
    elasticForceOnBeadStd = zeros(maxTimeSteps)

    chroRepulsionForceOnBeadMean = zeros(maxTimeSteps)
    chroRepulsionForceOnBeadStd = zeros(maxTimeSteps)

    ladForceOnBeadMean = zeros(maxTimeSteps)
    ladForceOnBeadStd = zeros(maxTimeSteps)

    cytoskeletonForceOnBeadMean = zeros(maxTimeSteps)
    cytoskeletonForceOnBeadStd = zeros(maxTimeSteps)

    negareaForceOnBeadMean = zeros(maxTimeSteps)
    negareaForceOnBeadStd = zeros(maxTimeSteps)

    negvolumeForceOnBeadMean = zeros(maxTimeSteps)
    negvolumeForceOnBeadStd = zeros(maxTimeSteps)

    negbendingForceOnBeadMean = zeros(maxTimeSteps)
    negbendingForceOnBeadStd = zeros(maxTimeSteps)

    negelasticForceOnBeadMean = zeros(maxTimeSteps)
    negelasticForceOnBeadStd = zeros(maxTimeSteps)

    negchroRepulsionForceOnBeadMean = zeros(maxTimeSteps)
    negchroRepulsionForceOnBeadStd = zeros(maxTimeSteps)

    negladForceOnBeadMean = zeros(maxTimeSteps)
    negladForceOnBeadStd = zeros(maxTimeSteps)

    negcytoskeletonForceOnBeadMean = zeros(maxTimeSteps)
    negcytoskeletonForceOnBeadStd = zeros(maxTimeSteps)

    for i = 1:maxTimeSteps

        forceOnBeadMean, forceOnBeadStd = get_force_means(forceOnBead[:,i],forceOnBeadMean,forceOnBeadStd,i)
        areaForceOnBeadMean, areaForceOnBeadStd = get_force_means(areaForceOnBead[:,i],areaForceOnBeadMean,areaForceOnBeadStd,i)
        volumeForceOnBeadMean, volumeForceOnBeadStd = get_force_means(volumeForceOnBead[:,i],volumeForceOnBeadMean,volumeForceOnBeadStd,i)
        bendingForceOnBeadMean, bendingForceOnBeadStd = get_force_means(bendingForceOnBead[:,i],bendingForceOnBeadMean,bendingForceOnBeadStd,i)
        elasticForceOnBeadMean, elasticForceOnBeadStd = get_force_means(elasticForceOnBead[:,i],elasticForceOnBeadMean,elasticForceOnBeadStd,i)
        chroRepulsionForceOnBeadMean, chroRepulsionForceOnBeadStd = get_force_means(chroRepulsionForceOnBead[:,i],chroRepulsionForceOnBeadMean,chroRepulsionForceOnBeadStd,i)
        ladForceOnBeadMean, ladForceOnBeadStd = get_force_means(ladForceOnBead[:,i],ladForceOnBeadMean,ladForceOnBeadStd,i)
        cytoskeletonForceOnBeadMean, cytoskeletonForceOnBeadStd = get_force_means(cytoskeletonForceOnBead[:,i],cytoskeletonForceOnBeadMean,cytoskeletonForceOnBeadStd,i)

        negareaForceOnBeadMean, negareaForceOnBeadStd = get_force_means(negareaForceOnBead[:,i],negareaForceOnBeadMean,negareaForceOnBeadStd,i)
        negvolumeForceOnBeadMean, negvolumeForceOnBeadStd = get_force_means(negvolumeForceOnBead[:,i],negvolumeForceOnBeadMean,negvolumeForceOnBeadStd,i)
        negbendingForceOnBeadMean, negbendingForceOnBeadStd = get_force_means(negbendingForceOnBead[:,i],negbendingForceOnBeadMean,negbendingForceOnBeadStd,i)
        negelasticForceOnBeadMean, negelasticForceOnBeadStd = get_force_means(negelasticForceOnBead[:,i],negelasticForceOnBeadMean,negelasticForceOnBeadStd,i)
        negchroRepulsionForceOnBeadMean, negchroRepulsionForceOnBeadStd = get_force_means(negchroRepulsionForceOnBead[:,i],negchroRepulsionForceOnBeadMean,negchroRepulsionForceOnBeadStd,i)
        negladForceOnBeadMean, negladForceOnBeadStd = get_force_means(negladForceOnBead[:,i],negladForceOnBeadMean,negladForceOnBeadStd,i)
        negcytoskeletonForceOnBeadMean, negcytoskeletonForceOnBeadStd = get_force_means(negcytoskeletonForceOnBead[:,i],negcytoskeletonForceOnBeadMean,negcytoskeletonForceOnBeadStd,i)

    end

    if true

        folder = pick_folder(pwd()*"\\results")

        if folder == ""
            return
        end
    
        folder = splitdir(folder)[end]
    
        global numberOfNumbers = 0
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
    
        allDirs = readdir(".//results")
    
        temp = zeros(Int64,length(allDirs))
    
        for i = eachindex(temp)
    
            if allDirs[i][1:end-numberOfNumbers] == namePart
                temp[i] = 1
            end
    
        end
    
        numberOfSims = sum(temp)
    
        nTimeSteps = zeros(Int,numberOfSims)
        maxTimeSteps = 0;
    
        for i = 1:numberOfSims
    
            allFiles = readdir(".\\results\\"*namePart*lpad(i,numberOfNumbers,'0'))
    
            temp = zeros(Int64,length(allFiles))
    
            for i = eachindex(temp)
                if allFiles[i][1:4] == "afm_"
                    temp[i] = 1
                end
            end
    
            nTimeSteps[i] = sum(temp)
    
            if maxTimeSteps < sum(temp)
                maxTimeSteps = sum(temp)
            end
    
    
        end
    
        forceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        areaForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        volumeForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        bendingForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        elasticForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        chroRepulsionForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        ladForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        cytoskeletonForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negareaForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negvolumeForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negbendingForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negelasticForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negchroRepulsionForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negladForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
        negcytoskeletonForceOnBead2 = zeros(Float64,numberOfSims,maxTimeSteps)
    
        for i = 1:numberOfSims
    
            for j = 1:nTimeSteps[i]
    
                pointData = get_point_data(VTKFile(".\\results\\"*namePart*lpad(i,numberOfNumbers,'0')*"\\afm_"*lpad(j,4,'0')*".vtu"))
    
                forceOnBead2[i,j] = get_data(pointData["Force on bead"])[1]
                areaForceOnBead2[i,j] = get_data(pointData["Area force on bed"])[1]
                volumeForceOnBead2[i,j] = get_data(pointData["Volume force on bed"])[1]
                bendingForceOnBead2[i,j] = get_data(pointData["Bending force on bed"])[1]
                elasticForceOnBead2[i,j] = get_data(pointData["Elastic force on bed"])[1]
                chroRepulsionForceOnBead2[i,j] = get_data(pointData["chroRepulsion force on bed"])[1]
                ladForceOnBead2[i,j] = get_data(pointData["LAD force on bed"])[1]
                cytoskeletonForceOnBead2[i,j] = get_data(pointData["afm force on cell"])[1]
                negareaForceOnBead2[i,j] = get_data(pointData["Negative Area force on bed"])[1]
                negvolumeForceOnBead2[i,j] = get_data(pointData["Negative Volume force on bed"])[1]
                negbendingForceOnBead2[i,j] = get_data(pointData["Negative Bending force on bed"])[1]
                negelasticForceOnBead2[i,j] = get_data(pointData["Negative Elastic force on bed"])[1]
                negchroRepulsionForceOnBead2[i,j] = get_data(pointData["Negative chroRepulsion force on bed"])[1]
                negladForceOnBead2[i,j] = get_data(pointData["Negative LAD force on bed"])[1]
                negcytoskeletonForceOnBead2[i,j] = get_data(pointData["Negative afm force on cell"])[1]
            end
        end
    
        forceOnBeadMean2 = zeros(maxTimeSteps)
        forceOnBeadStd2 = zeros(maxTimeSteps)
    
        areaForceOnBeadMean2 = zeros(maxTimeSteps)
        areaForceOnBeadStd2 = zeros(maxTimeSteps)
    
        volumeForceOnBeadMean2 = zeros(maxTimeSteps)
        volumeForceOnBeadStd2 = zeros(maxTimeSteps)
    
        bendingForceOnBeadMean2 = zeros(maxTimeSteps)
        bendingForceOnBeadStd2 = zeros(maxTimeSteps)
    
        elasticForceOnBeadMean2 = zeros(maxTimeSteps)
        elasticForceOnBeadStd2 = zeros(maxTimeSteps)
    
        chroRepulsionForceOnBeadMean2 = zeros(maxTimeSteps)
        chroRepulsionForceOnBeadStd2 = zeros(maxTimeSteps)
    
        ladForceOnBeadMean2 = zeros(maxTimeSteps)
        ladForceOnBeadStd2 = zeros(maxTimeSteps)
    
        cytoskeletonForceOnBeadMean2 = zeros(maxTimeSteps)
        cytoskeletonForceOnBeadStd2 = zeros(maxTimeSteps)

        negareaForceOnBeadMean2 = zeros(maxTimeSteps)
        negareaForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negvolumeForceOnBeadMean2 = zeros(maxTimeSteps)
        negvolumeForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negbendingForceOnBeadMean2 = zeros(maxTimeSteps)
        negbendingForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negelasticForceOnBeadMean2 = zeros(maxTimeSteps)
        negelasticForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negchroRepulsionForceOnBeadMean2 = zeros(maxTimeSteps)
        negchroRepulsionForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negladForceOnBeadMean2 = zeros(maxTimeSteps)
        negladForceOnBeadStd2 = zeros(maxTimeSteps)
    
        negcytoskeletonForceOnBeadMean2 = zeros(maxTimeSteps)
        negcytoskeletonForceOnBeadStd2 = zeros(maxTimeSteps)
    
    
        for i = 1:maxTimeSteps
    
            forceOnBeadMean2, forceOnBeadStd2 = get_force_means(forceOnBead2[:,i],forceOnBeadMean2,forceOnBeadStd2,i)
            areaForceOnBeadMean2, areaForceOnBeadStd2 = get_force_means(areaForceOnBead2[:,i],areaForceOnBeadMean2,areaForceOnBeadStd2,i)
            volumeForceOnBeadMean2, volumeForceOnBeadStd2 = get_force_means(volumeForceOnBead2[:,i],volumeForceOnBeadMean2,volumeForceOnBeadStd2,i)
            bendingForceOnBeadMean2, bendingForceOnBeadStd2 = get_force_means(bendingForceOnBead2[:,i],bendingForceOnBeadMean2,bendingForceOnBeadStd2,i)
            elasticForceOnBeadMean2, elasticForceOnBeadStd2 = get_force_means(elasticForceOnBead2[:,i],elasticForceOnBeadMean2,elasticForceOnBeadStd2,i)
            chroRepulsionForceOnBeadMean2, chroRepulsionForceOnBeadStd2 = get_force_means(chroRepulsionForceOnBead2[:,i],chroRepulsionForceOnBeadMean2,chroRepulsionForceOnBeadStd2,i)
            ladForceOnBeadMean2, ladForceOnBeadStd2 = get_force_means(ladForceOnBead2[:,i],ladForceOnBeadMean2,ladForceOnBeadStd2,i)
            cytoskeletonForceOnBeadMean2, cytoskeletonForceOnBeadStd2 = get_force_means(cytoskeletonForceOnBead2[:,i],cytoskeletonForceOnBeadMean2,cytoskeletonForceOnBeadStd2,i)
    
            negareaForceOnBeadMean2, negareaForceOnBeadStd2 = get_force_means(negareaForceOnBead2[:,i],negareaForceOnBeadMean2,negareaForceOnBeadStd2,i)
            negvolumeForceOnBeadMean2, negvolumeForceOnBeadStd2 = get_force_means(negvolumeForceOnBead2[:,i],negvolumeForceOnBeadMean2,negvolumeForceOnBeadStd2,i)
            negbendingForceOnBeadMean2, negbendingForceOnBeadStd2 = get_force_means(negbendingForceOnBead2[:,i],negbendingForceOnBeadMean2,negbendingForceOnBeadStd2,i)
            negelasticForceOnBeadMean2, negelasticForceOnBeadStd2 = get_force_means(negelasticForceOnBead2[:,i],negelasticForceOnBeadMean2,negelasticForceOnBeadStd2,i)
            negchroRepulsionForceOnBeadMean2, negchroRepulsionForceOnBeadStd2 = get_force_means(negchroRepulsionForceOnBead2[:,i],negchroRepulsionForceOnBeadMean2,negchroRepulsionForceOnBeadStd2,i)
            negladForceOnBeadMean2, negladForceOnBeadStd2 = get_force_means(negladForceOnBead2[:,i],negladForceOnBeadMean2,negladForceOnBeadStd2,i)
            negcytoskeletonForceOnBeadMean2, negcytoskeletonForceOnBeadStd2 = get_force_means(negcytoskeletonForceOnBead2[:,i],negcytoskeletonForceOnBeadMean2,negcytoskeletonForceOnBeadStd2,i)

        end

    end
    


    times = 1:length(forceOnBeadMean)
    times2 = 1:length(forceOnBeadMean2)

    traces=[
        scatter(
            x=times,
            y=forceOnBeadMean,
            line=attr(color="rgb(255,0,0)"),
            mode="lines",
            name="Force on bead",
            error_y=attr(type="data", array=forceOnBeadStd, visible=true),
        ),
    
        scatter(
            x=times,
            y=areaForceOnBeadMean,
            line=attr(color="rgb(0,255,0)"),
            mode="lines",
            name="Area force on bead",
            error_y=attr(type="data", array=areaForceOnBeadStd, visible=true),
        ),
            
        scatter(
            x=times,
            y=volumeForceOnBeadMean,
            line=attr(color="rgb(0,0,255)"),
            mode="lines",
            name="Volume force on bead",
            error_y=attr(type="data", array=volumeForceOnBeadStd, visible=true),
        ),
         
        scatter(
            x=times,
            y=bendingForceOnBeadMean,
            line=attr(color="rgb(255,0,255)"),
            mode="lines",
            name="Bending force on bead",
            error_y=attr(type="data", array=bendingForceOnBeadStd, visible=true),
        ),
                
        scatter(
            x=times,
            y=elasticForceOnBeadMean,
            line=attr(color="rgb(0,0,128)"),
            mode="lines",
            name="Elastic force on bead",
            error_y=attr(type="data", array=elasticForceOnBeadStd, visible=true),
        ),
                  
        scatter(
            x=times,
            y=chroRepulsionForceOnBeadMean,
            line=attr(color="rgb(255,102,0)"),
            mode="lines",
            name="Chro repulsion force on bead",
            error_y=attr(type="data", array=chroRepulsionForceOnBeadStd, visible=true),
        ),
                
        scatter(
            x=times,
            y=ladForceOnBeadMean,
            line=attr(color="rgb(0,255,255)"),
            mode="lines",
            name="LAD force on bead",
            error_y=attr(type="data", array=ladForceOnBeadStd, visible=true),
        ),

        scatter(
            x=times,
            y=forceOnBeadMean2,
            line=attr(color="rgb(255,0,0)", dash="dot"),
            mode="lines",
            name="Force on bead stiff",
            error_y=attr(type="data", array=forceOnBeadStd2, visible=true),
        ),
    
        scatter(
            x=times,
            y=areaForceOnBeadMean2,
            line=attr(color="rgb(0,255,0)", dash="dot"),
            mode="lines",
            name="Area force on bead stiff",
            error_y=attr(type="data", array=areaForceOnBeadStd2, visible=true),
        ),
            
        scatter(
            x=times,
            y=volumeForceOnBeadMean2,
            line=attr(color="rgb(0,0,255)", dash="dot"),
            mode="lines",
            name="Volume force on bead stiff",
            error_y=attr(type="data", array=volumeForceOnBeadStd2, visible=true),
        ),
         
        scatter(
            x=times,
            y=bendingForceOnBeadMean2,
            line=attr(color="rgb(255,0,255)", dash="dot"),
            mode="lines",
            name="Bending force on bead stiff",
            error_y=attr(type="data", array=bendingForceOnBeadStd2, visible=true),
        ),
                
        scatter(
            x=times,
            y=elasticForceOnBeadMean2,
            line=attr(color="rgb(0,0,128)", dash="dot"),
            mode="lines",
            name="Elastic force on bead stiff",
            error_y=attr(type="data", array=elasticForceOnBeadStd2, visible=true),
        ),
                  
        scatter(
            x=times,
            y=chroRepulsionForceOnBeadMean2,
            line=attr(color="rgb(255,102,0)", dash="dot"),
            mode="lines",
            name="Chro repulsion force on bead stiff",
            error_y=attr(type="data", array=chroRepulsionForceOnBeadStd2, visible=true),
        ),
                
        scatter(
            x=times,
            y=ladForceOnBeadMean2,
            line=attr(color="rgb(0,255,255)", dash="dot"),
            mode="lines",
            name="LAD force on bead stiff",
            error_y=attr(type="data", array=ladForceOnBeadStd2, visible=true),
        ),
    
        scatter(
            x=times,
            y=negareaForceOnBeadMean,
            line=attr(color="rgb(0,255,0)"),
            mode="lines",
            name="Area force on bead neg",
            error_y=attr(type="data", array=negareaForceOnBeadStd, visible=true),
        ),
            
        scatter(
            x=times,
            y=negvolumeForceOnBeadMean,
            line=attr(color="rgb(0,0,255)"),
            mode="lines",
            name="Volume force on bead neg",
            error_y=attr(type="data", array=negvolumeForceOnBeadStd, visible=true),
        ),
         
        scatter(
            x=times,
            y=negbendingForceOnBeadMean,
            line=attr(color="rgb(255,0,255)"),
            mode="lines",
            name="Bending force on bead neg",
            error_y=attr(type="data", array=negbendingForceOnBeadStd, visible=true),
        ),
                
        scatter(
            x=times,
            y=negelasticForceOnBeadMean,
            line=attr(color="rgb(0,0,128)"),
            mode="lines",
            name="Elastic force on bead neg",
            error_y=attr(type="data", array=negelasticForceOnBeadStd, visible=true),
        ),
                  
        scatter(
            x=times,
            y=negchroRepulsionForceOnBeadMean,
            line=attr(color="rgb(255,102,0)"),
            mode="lines",
            name="Chro repulsion force on bead neg",
            error_y=attr(type="data", array=negchroRepulsionForceOnBeadStd, visible=true),
        ),
                
        scatter(
            x=times,
            y=negladForceOnBeadMean,
            line=attr(color="rgb(0,255,255)"),
            mode="lines",
            name="LAD force on bead neg",
            error_y=attr(type="data", array=negladForceOnBeadStd, visible=true),
        ),
    
        scatter(
            x=times,
            y=negareaForceOnBeadMean2,
            line=attr(color="rgb(0,255,0)", dash="dot"),
            mode="lines",
            name="Area force on bead stiff neg",
            error_y=attr(type="data", array=negareaForceOnBeadStd2, visible=true),
        ),
            
        scatter(
            x=times,
            y=negvolumeForceOnBeadMean2,
            line=attr(color="rgb(0,0,255)", dash="dot"),
            mode="lines",
            name="Volume force on bead stiff neg",
            error_y=attr(type="data", array=negvolumeForceOnBeadStd2, visible=true),
        ),
         
        scatter(
            x=times,
            y=negbendingForceOnBeadMean2,
            line=attr(color="rgb(255,0,255)", dash="dot"),
            mode="lines",
            name="Bending force on bead stiff neg",
            error_y=attr(type="data", array=negbendingForceOnBeadStd2, visible=true),
        ),
                
        scatter(
            x=times,
            y=negelasticForceOnBeadMean2,
            line=attr(color="rgb(0,0,128)", dash="dot"),
            mode="lines",
            name="Elastic force on bead stiff neg",
            error_y=attr(type="data", array=negelasticForceOnBeadStd2, visible=true),
        ),
                  
        scatter(
            x=times,
            y=negchroRepulsionForceOnBeadMean2,
            line=attr(color="rgb(255,102,0)", dash="dot"),
            mode="lines",
            name="Chro repulsion force on bead stiff neg",
            error_y=attr(type="data", array=negchroRepulsionForceOnBeadStd2, visible=true),
        ),
                
        scatter(
            x=times,
            y=negladForceOnBeadMean2,
            line=attr(color="rgb(0,255,255)", dash="dot"),
            mode="lines",
            name="LAD force on bead stiff neg",
            error_y=attr(type="data", array=negladForceOnBeadStd2, visible=true),
        ),

    ]

    p = plot(traces)



end

function get_force_means(temp,means,stds,i)

    if sum(temp) > 0
        temp = [x == 0 ? missing : x for x in temp]
        means[i] = mean(skipmissing(temp))
        stds[i] = std(skipmissing(temp))
    else
        means[i] = mean(temp)
        stds[i] = std(temp)
    end

    stds[isnan.(stds)] .= 0

    return means,stds
end