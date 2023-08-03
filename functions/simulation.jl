function simulation(
    simType::String,
    maxT::Number,
    folderName::String,
    initState::String;
    importFolder::String="",
    nameDate::Bool=true,
    parameterFile::String="./parameters/parameters.txt",
    exportData::Bool=true,
    replComp::Bool = false,
    adherent::Bool = false,
    noEnveSolve::Bool = false,
    noChromatin::Bool = false,
    returnFoldername::Bool = false)

    # check if the simulation type is valid. If not, return
    if check_simulation_type(simType)
        return
    end

    # setup the simulation environment
    enve, chro, spar, simset, ext, ipar = setup_simulation(initState, simType, importFolder, parameterFile,noChromatin,noEnveSolve,adherent,maxT)

    # setup the export settings
    ex = setup_export(simType,folderName, enve, chro, ext, spar, simset, nameDate,exportData,noChromatin,ipar)
    
    # print a message to indicate that the simulation is starting
    printstyled("Starting simulation (" * Dates.format(now(), "YYYY-mm-dd HH:MM") * ")\n"; color=:blue)

    # run the simulation
    if noChromatin
        run_simulation(enve, spar, ex, ext, simset, maxT)
    elseif replComp
        repl = create_replication_compartment(enve,spar,simset.initType)
        run_simulation(enve, chro, repl, spar, ex, ext, simset, maxT)
    else
        run_simulation(enve, chro, spar, ex, ext, simset, maxT)
    end
    
    # if the user wants to return the folder name, do so
    if returnFoldername
        return ex.folderName
    end
end

function run_simulation(enve, chro, spar, ex, ext, simset, maxT)

    # maximum number of time steps
    intMaxTime = Int64(floor(maxT/(spar.dt*spar.scalingTime)))

    # current full time step
    intTime = Int64(0)

    # set the timing for the time step
    simset.timeStepTiming = now()

    try

        # perform simulation steps for each time step
        while intTime <= intMaxTime

            # get nuclear properties
            get_nuclear_properties!(enve, chro, simset, spar)

            # get crosslinks
            get_crosslinks!(enve, chro, simset, spar)

            # calculate forces
            get_forces!(enve, chro, spar, ext, simset, intTime)

            # export data
            export_data(enve, chro, spar, ex, ext, intTime, simset)

            # solve the system
            solve_system!(enve, chro, spar, simset, ext, intTime)

            # update the time step
            intTime = progress_time!(simset,intTime);

        end

    catch error

        # handle any errors that occur during simulation
        print_error(error)

    end

    # perform post-simulation export operations
    post_export(ex,simset,ext)

end

function run_simulation(enve::envelopeType, chro::chromatinType, repl::replicationCompartmentType, spar::scaledParametersType, ex::exportSettingsType, ext, simset::simulationSettingsType, maxT::Number)

    # maximum number of time steps
    intMaxTime = Int64(floor(maxT/(spar.dt*spar.scalingTime)))

    # current full time step
    intTime = Int64(0)

    # set the timing for the time step
    simset.timeStepTiming = now()

    try

        # perform simulation steps for each time step
        while intTime <= intMaxTime

            # get nuclear properties
            get_nuclear_properties!(enve, chro, repl, simset, spar)

            # get crosslinks
            get_crosslinks!(enve, chro, simset, spar)

            # calculate forces
            get_forces!(enve, chro, repl, spar, ext, simset)

            # export data
            export_data(enve, chro, repl, spar, ex, ext, intTime, simset)

            # solve the system
            solve_system!(enve, chro, repl, spar, simset, ext)

            # update the time step
            intTime = progress_time!(simset,intTime);

        end

    catch error

        # handle any errors that occur during simulation
        print_error(error)

    end

    # perform post-simulation export operations
    post_export(ex,simset,ext)

end

function run_simulation(enve, spar, ex, ext, simset, maxT)

    # maximum number of time steps
    intMaxTime = Int64(floor(maxT/(spar.dt*spar.scalingTime)))

    # current full time step
    intTime = Int64(0)

    # set the timing for the time step
    simset.timeStepTiming = now()

    try

        # perform simulation steps for each time step
        while intTime <= intMaxTime

            # get nuclear properties
            get_nuclear_properties!(enve, simset)

            # calculate forces
            get_forces!(enve, spar, ext, simset, intTime)

            # export data
            export_data(enve, ex, ext, intTime, simset)

            # solve the system
            solve_system!(enve, spar, simset, ext, intTime)

            # update the time step
            intTime = progress_time!(simset,intTime);

        end

    catch error

        # handle any errors that occur during simulation  
        print_error(error)

    end

    # perform post-simulation export operations
    post_export(ex,simset,ext)
    
end