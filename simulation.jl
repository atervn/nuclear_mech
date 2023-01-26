function simulation(simType::String, maxT::Number, folderName::String, initState::String;
    importFolder::String="",
    nameDate::Bool=true,
    parameterFile::String="./parameters.txt",
    exportData::Bool=true,
    replComp::Bool = false,
    adherent::Bool = false,
    noEnveSolve::Bool = false,
    noChromatin::Bool = false,
    returnFoldername::Bool = false)

    if check_simulation_type(simType)
        return
    end

    # setup
    enve, chro, spar, simset, ext, ipar = setup_simulation(initState, simType, importFolder, parameterFile,noChromatin,adherent)
    
    ex = setup_export(simType,folderName, enve, chro, ext, spar, simset, nameDate,exportData,noChromatin,ipar)
    
    simset.prog = Progress(Int64(round(maxT/(spar.scalingTime*spar.maxDt))), 0.1, "Simulating...", 100)
    
    printstyled("Starting simulation (" * Dates.format(now(), "YYYY-mm-dd HH:MM") * ")\n"; color=:blue)

    if noChromatin
        run_simulation(enve, spar, ex, ext, simset, maxT)
    elseif replComp
        repl = create_replication_compartment(enve,spar)
        run_simulation(enve, chro, repl, spar, ex, ext, simset, maxT)
    else
        run_simulation(enve, chro, spar, ex, ext, simset, maxT, noEnveSolve)
    end
    
    
    if returnFoldername
        return ex.folderName
    end

end

function run_simulation(enve, chro, spar, ex, ext, simset, maxT, noEnveSolve)

    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    simset.timeStepTiming = now()
    try
        while intTime <= intMaxTime

            get_nuclear_properties!(enve, chro, simset, spar)

            get_crosslinks!(enve, chro, simset, spar)

            get_forces!(enve, chro, spar, ext, simset, intTime)

            export_data(enve, chro, spar, ex, ext, intTime, simset)

            solve_system!(enve, chro, spar, simset, dt, ext, noEnveSolve)

            intTime = progress_time!(simset,intTime);

        end
    catch error
        print_error(error)
    end

    post_export(ex,simset,ext)
end

function run_simulation(enve::envelopeType, chro::chromatinType, repl::replicationCompartmentType, spar::scaledParametersType, ex::exportSettingsType, ext, simset::simulationSettingsType, maxT::Number)

    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    simset.timeStepTiming = now()
    try
        while intTime <= intMaxTime

            get_nuclear_properties!(enve, chro, repl, simset, spar)

            get_crosslinks!(enve, chro, simset, spar)

            get_forces!(enve, chro, repl, spar, ext, simset)

            export_data(enve, chro, repl, spar, ex, ext, intTime, simset)

            solve_system!(enve, chro, repl, spar, simset, dt, ext)

            intTime = progress_time!(simset,intTime);

        end
    catch error
        print_error(error)
    end

    post_export(ex,simset,ext)
end

function run_simulation(enve, spar, ex, ext, simset, maxT)

    intMaxTime = Int64(floor(maxT/(spar.maxDt*spar.scalingTime)))
    intTime = Int64(0)
    dt = spar.maxDt

    simset.timeStepTiming = now()
    try
        while intTime <= intMaxTime

            get_nuclear_properties!(enve, simset, spar)

            get_forces!(enve, spar, ext, simset)

            export_data(enve, ex, ext, intTime, simset)

            solve_system!(enve, spar, simset, dt, ext)

            intTime = progress_time!(simset,intTime);

        end
    catch error
        print_error(error)
    end

    post_export(ex,simset,ext)
end