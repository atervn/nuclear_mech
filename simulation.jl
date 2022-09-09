function simulation(simType,maxT,folderName)

    if cmp(simType,"MA") != 0 && cmp(simType,"MM") != 0
        printstyled("Unknown simulation type"; color = :blue)
        return
    end

    # number of subdivisions when creating the geometry
    nSubdivisions = 4;
    
    # time step
    dt = 0.05;
    
    # simulation time
    # maxT = 10000;
    
    # model parameters
    ipar = inputParametersType();
    ipar = read_parameters(ipar,"./parameters.txt");
    
    # create nucleus
    printstyled("Creating nuclear envelope..."; color = :blue)
    nuc = nucleusType();
    nuc = create_icosahedron!(nuc,ipar);
    nuc = subdivide_mesh!(nuc,ipar,nSubdivisions)
    nuc = setup_nucleus_data(nuc)
    printstyled("Done!\n"; color = :blue)
    
    # scale parameters
    spar = scaledParametersType();
    spar = get_model_parameters(ipar,spar,nuc);
    
    printstyled("Creating chromatin..."; color = :blue)
    chro = chromatinType();
    chro = create_all_chromatins(chro,spar)
    printstyled("Done!\n"; color = :blue)

    # setup nucleus export
    triCells, lineCells,folderName = setup_export(folderName,nuc,chro,spar)

    # setup aspiration
    if cmp(simType,"MA") == 0
        pip = generate_pipette_mesh();
        export_pipette_mesh(folderName,pip)
        # vector to store the aspiration lengths
        maxX = zeros(Float64,maxT+1)
    
    # setup micromanipulation 
    elseif cmp(simType,"MM") == 0
        mm = setup_micromanipulation(nuc)
        nuclearLength = zeros(Float64,maxT+1)
    end

    # create the friction matrix
    frictionMatrix = get_friction_matrix(nuc,spar);
   
    printstyled("Starting simulation ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :blue)

    # setup progress bar
    p = Progress(maxT, 1, "Simulating...", 100)
    
    # run the simulation
    for t = 0:maxT
        
        if cmp(simType,"MA") == 0
            maxX[t+1] = maximum(getindex.(nuc.vert,1));
        elseif cmp(simType,"MM") == 0
            nuclearLength[t+1] = nuc.vert[mm.rightmostVertex][1] - nuc.vert[mm.leftmostVertex][1];
        end

        envelopeTree = KDTree(nuc.vert);
        chromatinTree = KDTree(chro.vert);
        get_strand_vectors!(chro,spar)
        initialize_chromatin_forces!(chro)

        get_edge_vectors!(nuc);
        get_voronoi_areas!(nuc);
        get_area_unit_vectors!(nuc);
        # get_local_curvatures!(nuc);
        get_triangle_normals!(nuc);
    
        get_volume_forces!(nuc,spar);
        get_area_forces!(nuc,spar);
        get_bending_forces!(nuc,spar);
        get_elastic_forces!(nuc,spar);
        get_repulsion_forces!(nuc,spar,envelopeTree);
        if cmp(simType,"MA") == 0 
            repulsion = get_aspiration_repulsion_forces(nuc,pip,spar);
            aspiration = get_aspiration_forces(nuc,pip,spar);
        elseif cmp(simType,"MM") == 0
            micromanipulation = get_micromanipulation_forces(nuc,mm,spar)
        end
        get_linear_chromatin_forces!(chro,spar);
        get_bending_chromatin_forces!(chro,spar,130)
        get_chromation_chromation_repulsion_forces!(chro,spar,chromatinTree)
        fluctuationForces = get_random_fluctuations(spar)
        get_envelope_chromatin_repulsion_forces!(nuc,chro,spar,envelopeTree)


        local totalForces = nuc.forces.volume .+ nuc.forces.area .+ nuc.forces.elastic .+ nuc.forces.envelopeRepulsion .+ nuc.forces.chromationRepulsion;
        if cmp(simType,"MA") == 0 
            totalForces .+=  repulsion .+ aspiration
        elseif cmp(simType,"MM") == 0
            totalForces .+= micromanipulation
        end

        totalChromatinForces = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ fluctuationForces .+ chro.forces.enveRepulsion;
        

        if mod(t,20) == 0
            vtk_grid(".\\results\\"*folderName*"\\nucl_" * lpad(t,4,"0"), [getindex.(nuc.vert,1) getindex.(nuc.vert,2) getindex.(nuc.vert,3)]', triCells) do vtk
                if cmp(simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.triangleNormalUnitVectors,1) getindex.(nuc.triangleNormalUnitVectors,2) getindex.(nuc.triangleNormalUnitVectors,3)]'
            end
            vtk_grid(".\\results\\"*folderName*"\\chro_" * lpad(t,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', lineCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
            end
        end
    
        local vX = cg(frictionMatrix,[getindex.(totalForces,1);getindex.(totalChromatinForces,1)]);
        local vY = cg(frictionMatrix,[getindex.(totalForces,2);getindex.(totalChromatinForces,2)]);
        local vZ = cg(frictionMatrix,[getindex.(totalForces,3);getindex.(totalChromatinForces,3)]);
    
        for i = 1:length(nuc.vert)
            nuc.vert[i] += Vec(vX[i],vY[i],vZ[i])*dt
        end
        for k = 1:spar.chromatinLength*spar.chromatinNumber
            chro.vert[k] += Vec(vX[length(nuc.vert)+k],vY[length(nuc.vert)+k],vZ[length(nuc.vert)+k])*0.01
        end

        next!(p)
    
    end
    
    if cmp(simType,"MA") == 0
        writedlm(".\\results\\"*folderName*"\\maxX.csv", maxX,',')
        dL = maxX .- minimum(maxX);
    
        J = 2*pi.*dL./(3*2.1*3*1);
        
        plot(10*dt:dt:maxT*dt,J[11:end],yaxis=:log,xaxis=:log,xlim = (0.1, 200),ylim = (0.01, 10))
    elseif cmp(simType,"MM") == 0
        writedlm(".\\results\\"*folderName*"\\nuclearLength.csv", nuclearLength,',')
    end
    
    
end