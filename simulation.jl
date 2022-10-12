function simulation(simType,maxT,folderName, initState; importFolder ="", nameDate = "yes", parameterFile = "./parameters.txt")

    if check_simulation_type(simType); return; end

    # setup simulation
    nuc,chro,spar,simset,ext = setup_simulation(initState,simType,importFolder,parameterFile)

    # setup nucleus export
    ex = setup_export(folderName,nuc,chro,spar,nameDate)
    ex.step = 2*spar.dt;

    maxT = maxT/spar.scalingTime

    # run_simulation!(nuc,chro,spar,maxT,frictionMatrix)
    printstyled("Starting simulation ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :blue)

    # setup progress bar
    simset.prog = Progress(Int64(round(maxT/spar.dt)), 0.1, "Simulating...", 100)

    time = 0;
    intTime = Int64(0)

    intMaxTime = Int64(floor(maxT/spar.dt))

    dt = spar.dt;

    nuc.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
    nuc.forces.area = Vector{Vec{3,Float64}}(undef, length(nuc.vert));

    nCrosslinks = []

    ####################################################################################################
    # run the simulation

    try
        while intTime <= intMaxTime
            
            save_specific_data!(nuc,ext,simset)
    
            get_iteration_properties!(nuc,chro,simset,spar)
        
            get_crosslinks!(nuc,chro,simset,spar)

            get_forces!(nuc,chro,spar,ext,simset)

            export_data(nuc,chro,spar,ex,time,simType,simset)
            
            solve_system!(nuc,chro,spar,simset,dt)
            
            if cmp(simType,"PC") == 0
                if ext > -spar.freeNucleusRadius/3
                    ext = ext - 0.2*dt*simset.timeStepMultiplier;
                end
            end

            if simset.timeStepProgress == 0
                next!(simset.prog)
                intTime += 1
            end

            time += dt*simset.timeStepMultiplier

            push!(nCrosslinks,length(chro.crosslinks[:,1]))

            # println(simset.timeStepMultiplier)
        end
    catch e
        println("Simulation failed")
        rethrow(e)
    end

    ##################################################################################################

    if cmp(simType,"MA") == 0
        writedlm(".\\results\\"*ex.folderName*"\\maxX.csv", maxX,',')
        dL = ext[2] .- minimum(ext[2]);
    
        J = 2*pi.*dL./(3*2.1*3*1);
        
        plot(10*dt:dt:maxT*dt,J[11:end],yaxis=:log,xaxis=:log,xlim = (0.1, 200),ylim = (0.01, 10))
    elseif cmp(simType,"MM") == 0
        writedlm(".\\results\\"*ex.folderName*"\\nuclearLength.csv", ext[2],',')
    end
    
    return nCrosslinks

end