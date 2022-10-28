function simulation(simType::String, maxT, folderName::String, initState::String; importFolder::String="", nameDate::Bool=true, parameterFile::String="./parameters.txt", exportData::Bool=true)

    if check_simulation_type(simType)
        return
    end

    # setup
    nuc, chro, spar, simset, ext = setup_simulation(initState, simType, importFolder, parameterFile)
    ex = setup_export(folderName, nuc, chro, spar, nameDate,exportData)
    simset.prog = Progress(Int64(round(maxT/(spar.scalingTime*spar.maxDt))), 0.1, "Simulating...", 100)

    # timestepping variables
    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    printstyled("Starting simulation (" * Dates.format(now(), "YYYY-mm-dd HH:MM") * ")\n"; color=:blue)

    ####################################################################################################
    # run the simulation

    try
        while intTime <= intMaxTime

            get_nuclear_properties!(nuc, chro, simset, spar)

            get_crosslinks!(nuc, chro, simset, spar)

            get_forces!(nuc, chro, spar, ext, simset)

            export_data(nuc, chro, spar,ex ,ext, intTime, simset)

            solve_system!(nuc, chro, spar, simset, dt)

            # if cmp(simType, "PC") == 0
            #     if ext > -spar.freeNucleusRadius / 3
            #         ext = ext - 0.2 * dt * simset.timeStepMultiplier
            #     end
            # end
            
            intTime = progress_time!(simset,intTime);
            
        end
    catch error
        println("\nSimulation failed\n")
        rethrow(error)
    end

    ##################################################################################################
    
    post_export(ex,ext)
    
end