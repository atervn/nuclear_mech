function simulation(simType,maxT,folderName)

    if cmp(simType,"MA") != 0 && cmp(simType,"MM") != 0 && cmp(simType,"PC") != 0 && cmp(simType,"GR") != 0
        printstyled("Unknown simulation type"; color = :blue)
        return
    end
    
    # model parameters
    ipar = inputParametersType();
    ipar = read_parameters(ipar,"./parameters.txt");

    # create nucleus
    printstyled("Creating nuclear envelope..."; color = :blue)
    nuc = nucleusType();
    nuc = create_icosahedron!(nuc,ipar);
    nuc = subdivide_mesh!(nuc,ipar)
    nuc = setup_nucleus_data(nuc)
    printstyled("Done!\n"; color = :blue)
    
    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);
    
    ladCenterIdx = get_lad_centers(nuc,spar)
    nuc.lads = get_lad_enve_vertices(ladCenterIdx,nuc,spar)


    printstyled("Creating chromatin..."; color = :blue)
    chro = chromatinType();
    chro = create_all_chromsomes(chro,spar,nuc.vert[ladCenterIdx])
    
    chro.lads = get_lad_chro_vertices(nuc,spar)


    # nuc,chro = create_lads(chro,nuc,spar)
    printstyled("Done!\n"; color = :blue)

    # setup nucleus export
    ex = setup_export(folderName,nuc,chro,spar)
    ex.step = 10;

    # setup aspiration
    if cmp(simType,"MA") == 0
        pip = generate_pipette_mesh();
        export_pipette_mesh(folderName,pip)
        # vector to store the aspiration lengths
        maxX = zeros(Float64,maxT+1)
        ext = (pip,maxX)
    # setup micromanipulation 
    elseif cmp(simType,"MM") == 0
        mm = setup_micromanipulation(nuc)
        nuclearLength = zeros(Float64,maxT+1)
        ext = (mm,nuclearLength)
    elseif cmp(simType,"PC") == 0
        plane = spar.freeNucleusRadius + spar.repulsionDistance;
        ext = (plane)
    elseif cmp(simType,"GR") == 0
        ext = ()
    end

    simset = simulationSettingsType()
    simset.frictionMatrix = get_friction_matrix(nuc,spar);
    simset.simType = simType;
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
        
        plane = plane - 0.1*spar.dt;

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