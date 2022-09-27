function simulation(simType,maxT,folderName)

    if check_simulation_type(simType); return; end

    # setup simulation
    nuc,chro,spar,simset,ext = setup_simulation("load",simType)

    # setup nucleus export
    ex = setup_export(folderName,nuc,chro,spar)
    ex.step = 10;

    # run_simulation!(nuc,chro,spar,maxT,frictionMatrix)
    printstyled("Starting simulation ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :blue)

    # setup progress bar
    simset.prog = Progress(maxT, 0.1, "Simulating...", 100)
    
    ####################################################################################################
    # run the simulation
    for t = 0:maxT
        
        save_specific_data!(nuc,ext,simset)
  
        get_iteration_properties!(nuc,chro,simset,spar)
    
        get_forces!(nuc,chro,spar,ext,simset)

        export_data(nuc,chro,spar,ex,t,simType)
        
        solve_system!(nuc,chro,spar,simset)
        
        if cmp(simType,"PC") == 0
            ext = ext - 0.1*spar.dt;
        end

        next!(simset.prog)
    
    end
    
    ##################################################################################################

    if cmp(simType,"MA") == 0
        writedlm(".\\results\\"*folderName*"\\maxX.csv", maxX,',')
        dL = ext(2) .- minimum(maxX);
    
        J = 2*pi.*dL./(3*2.1*3*1);
        
        plot(10*dt:dt:maxT*dt,J[11:end],yaxis=:log,xaxis=:log,xlim = (0.1, 200),ylim = (0.01, 10))
    elseif cmp(simType,"MM") == 0
        writedlm(".\\results\\"*folderName*"\\nuclearLength.csv", ext(2),',')
    end
    
    
end