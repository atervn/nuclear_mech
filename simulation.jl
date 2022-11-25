function simulation(simType::String, maxT, folderName::String, initState::String; importFolder::String="", nameDate::Bool=true, parameterFile::String="./parameters.txt", exportData::Bool=true)

    if check_simulation_type(simType)
        return
    end

    # setup
     nuc, spar, simset, ext = setup_simulation(initState, simType, importFolder, parameterFile)
    ex = setup_export(folderName, nuc, spar, nameDate,exportData)
    simset.prog = Progress(Int64(round(maxT/(spar.scalingTime*spar.maxDt))), 0.1, "Simulating...", 100)

    # timestepping variables
    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    printstyled("Starting simulation (" * Dates.format(now(), "YYYY-mm-dd HH:MM") * ")\n"; color=:blue)

    ####################################################################################################
    # run the simulation

    nCrosslinks = []

    try
        while intTime <= intMaxTime

            
            get_nuclear_properties!(nuc, simset, spar)
            
            get_crosslinks!(nuc, simset, spar)
            
            get_forces!(nuc, spar, ext, simset)

            export_data(nuc, spar, ex, ext, intTime, simset)

            solve_system!(nuc, spar, simset, dt, ext)

            intTime = progress_time!(simset,intTime);
            
            push!(nCrosslinks,length(nuc.chro.crosslinks[:,1]))

        end
    catch error
        println("\nSimulation failed\n")
        rethrow(error)
    end

    ##################################################################################################
    
    post_export(ex,simset,ext)
end