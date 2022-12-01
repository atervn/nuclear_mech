function simulation_init(simType::String, maxT, folderName::String, initState::String, noEnve::Bool;
    importFolder::String="", nameDate::Bool=true, parameterFile::String="./parameters.txt",
    exportData::Bool=true; replComp = false,)

    if check_simulation_type(simType)
        return
    end

    # setup
    enve, chro, spar, simset, ext = setup_simulation(initState, simType, importFolder, parameterFile)
    ex = setup_export(folderName, enve, chro, spar, nameDate,exportData)
    simset = check_adhesion_file!(ex,initState,importFolder,simset)

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

            get_nuclear_properties!(enve, chro, simset, spar)

            get_crosslinks!(enve, chro, simset, spar)

            get_forces!(enve, chro, spar, ext, simset)

            export_data(enve, chro, spar,ex ,ext, intTime, simset,adh)

            solve_system_init!(enve, chro, spar, simset, dt, ext, noEnve)

            intTime = progress_time!(simset,intTime);
            
            push!(nCrosslinks,length(chro.crosslinks[:,1]))

        end
    catch error
        println("\nSimulation failed\n")
        rethrow(error)
    end

    ##################################################################################################
    
    post_export(ex,simset,ext)

    return ex.folderName
end