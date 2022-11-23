function simulation_adh_init(simType::String, maxT, folderName::String, initState::String; importFolder::String="", nameDate::Bool=true, parameterFile::String="./parameters_adh_init.txt", exportData::Bool=true)

    if check_simulation_type(simType)
        return
    end

    # setup
    nuc, spar, simset, ext = setup_simulation_adh_init(initState, simType, importFolder, parameterFile)
    ex = setup_export_adh_init(folderName, nuc, spar, nameDate,exportData)
    simset.prog = Progress(Int64(round(maxT/(spar.scalingTime*spar.maxDt))), 0.1, "Simulating...", 100)

    # timestepping variables
    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    printstyled("Starting simulation (" * Dates.format(now(), "YYYY-mm-dd HH:MM") * ")\n"; color=:blue)

    ####################################################################################################
    # run the simulation

    topForces = [];

    try
        while intTime <= intMaxTime

            get_nuclear_properties_adh_init!(nuc, simset)

            planeRepulsion,touchTop = get_forces_adh_init!(nuc, spar, ext)

            push!(topForces,sum(getindex.(nuc.forces.volume[touchTop],3)))

            export_data_adh_init(nuc, ex, intTime, simset, planeRepulsion, ext)

            solve_system_adh_init!(nuc, spar, simset, dt, ext, sum(getindex.(nuc.forces.volume[touchTop],3)))

            intTime = progress_time!(simset,intTime);

        end
    catch error
        printstyled("\nSimulation FAILED ----------- shitty code\n\n"; color=:red)
        rethrow(error)
    end

    writedlm("test.csv",topForces,',')
    ##################################################################################################
    
    post_export(ex,simset,ext)

    return ".\\results\\"*ex.folderName

end