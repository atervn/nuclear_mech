function simulation(simType,maxT,folderName)

    if cmp(simType,"MA") != 0 && cmp(simType,"MM") != 0 && cmp(simType,"PC") != 0
        printstyled("Unknown simulation type"; color = :blue)
        return
    end

    # number of subdivisions when creating the geometry
    nSubdivisions = 4;
    
    # time step
    dt = 0.02;
    
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
    chro = create_all_chromsomes(chro,spar)
    nuc,chro = create_lads(chro,nuc,spar)
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
    p = Progress(maxT, 0.1, "Simulating...", 100)
    
    plane = spar.freeNucleusRadius + spar.repulsionDistance;
  
    # run the simulation
    for t = 0:maxT
        
        # dt = 0.05;

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
        elseif cmp(simType,"PC") == 0
            planeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.vert));
            for i = eachindex(nuc.vert)
                planeRepulsion[i] = Vec(0.,0.,0.);
            end
            for i = eachindex(nuc.vert)
                if plane - nuc.vert[i][3] < spar.repulsionDistance
                    distance = plane - nuc.vert[i][3];
                    planeRepulsion[i] += Vec(0.,0.,-spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
                end
            end
            for i = eachindex(nuc.vert)
                if abs(-spar.freeNucleusRadius - nuc.vert[i][3]) < spar.repulsionDistance
                    distance = abs(-spar.freeNucleusRadius - nuc.vert[i][3]);
                    planeRepulsion[i] += Vec(0.,0.,spar.repulsionConstant*(spar.repulsionDistance - distance)^(3/2))
                end
            end
        end
        get_linear_chromatin_forces!(chro,spar);
        get_bending_chromatin_forces!(chro,spar,130)
        get_chromation_chromation_repulsion_forces!(chro,spar,chromatinTree)
        fluctuationForces = get_random_fluctuations(spar)
        get_envelope_chromatin_repulsion_forces!(nuc,chro,spar,envelopeTree)
        ladEnveForces, ladChroForces = get_lad_forces(nuc,chro,spar)


        local totalForces = nuc.forces.volume .+ nuc.forces.area .+ nuc.forces.elastic .+ nuc.forces.envelopeRepulsion .+ nuc.forces.chromationRepulsion .+ ladEnveForces;
        
        if cmp(simType,"MA") == 0 
            totalForces .+=  repulsion .+ aspiration
        elseif cmp(simType,"MM") == 0
            totalForces .+= micromanipulation
        else
            #totalForces .+= planeRepulsion
        end

        totalChromatinForces = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ fluctuationForces .+ chro.forces.enveRepulsion .+ ladChroForces;
    
        local vX = cg(frictionMatrix,[getindex.(totalForces,1);getindex.(totalChromatinForces,1)]);
        local vY = cg(frictionMatrix,[getindex.(totalForces,2);getindex.(totalChromatinForces,2)]);
        local vZ = cg(frictionMatrix,[getindex.(totalForces,3);getindex.(totalChromatinForces,3)]);
    
        if mod(t,50) == 0
            vtk_grid(".\\results\\"*folderName*"\\nucl_" * lpad(t,4,"0"), [getindex.(nuc.vert,1) getindex.(nuc.vert,2) getindex.(nuc.vert,3)]', triCells) do vtk
                if cmp(simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.triangleNormalUnitVectors,1) getindex.(nuc.triangleNormalUnitVectors,2) getindex.(nuc.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(nuc.forces.volume,1) getindex.(nuc.forces.volume,2) getindex.(nuc.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(nuc.forces.area,1) getindex.(nuc.forces.area,2) getindex.(nuc.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.forces.elastic,1) getindex.(nuc.forces.elastic,2) getindex.(nuc.forces.elastic,3)]'
                vtk["enveRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.envelopeRepulsion,1) getindex.(nuc.forces.envelopeRepulsion,2) getindex.(nuc.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.chromationRepulsion,1) getindex.(nuc.forces.chromationRepulsion,2) getindex.(nuc.forces.chromationRepulsion,3)]'

            end
            vtk_grid(".\\results\\"*folderName*"\\chro_" * lpad(t,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', lineCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
                vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                vtk["Movement", VTKPointData()] = [dt*vX[length(nuc.vert)+1:end] dt*vY[length(nuc.vert)+1:end] dt*vZ[length(nuc.vert)+1:end]]'
            end

            # export lads
            totalNum = sum(length.(nuc.lads));

            vertices = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,totalNum)
            for i = 1:totalNum
                vertices[i] = MeshCell(VTKCellTypes.VTK_VERTEX, Int64.([i]));
            end

            vertexIDs = []
            allVertices = []
            for i = 1:spar.chromatinNumber
                for j = 1:length(nuc.lads[i])
                    push!(vertexIDs, i)
                    push!(allVertices, nuc.lads[i][j])
                end
            end
            vertexIDs = Int64.(vertexIDs)

            vtk_grid(".\\results\\"*folderName*"\\enve_lads_" * lpad(t,4,"0"), [getindex.(nuc.vert[allVertices],1) getindex.(nuc.vert[allVertices],2) getindex.(nuc.vert[allVertices],3)]', vertices) do vtk
                vtk["vertex_id"] = vertexIDs
            end

            # export chrolads
            totalNum = sum(length.(chro.lads));

            vertices = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,totalNum)
            for i = 1:totalNum
                vertices[i] = MeshCell(VTKCellTypes.VTK_VERTEX, Int64.([i]));
            end

            vertexIDs = []
            allVertices = []
            for i = 1:spar.chromatinNumber
                for j = 1:length(chro.lads[i])
                    push!(vertexIDs, i)
                    push!(allVertices, chro.strandIdx[i][chro.lads[i][j]])
                end
            end
            vertexIDs = Int64.(vertexIDs)

            vtk_grid(".\\results\\"*folderName*"\\chro_lads_" * lpad(t,4,"0"), [getindex.(chro.vert[allVertices],1) getindex.(chro.vert[allVertices],2) getindex.(chro.vert[allVertices],3)]', vertices) do vtk
                vtk["vertex_id"] = vertexIDs
            end   





        end

        # enveVelocities = Vector{Vec{3,Float64}}(undef,length(nuc.vert))
        for i = 1:length(nuc.vert)
            # enveVelocities[i] = Vec(vX[i],vY[i],vZ[i])
            nuc.vert[i] += Vec(vX[i],vY[i],vZ[i])*dt
        end
        # chroVelocities = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
        for k = 1:spar.chromatinLength*spar.chromatinNumber
            # chroVelocities[k] = Vec(vX[length(nuc.vert)+k],vY[length(nuc.vert)+k],vZ[length(nuc.vert)+k])
            chro.vert[k] += Vec(vX[length(nuc.vert)+k],vY[length(nuc.vert)+k],vZ[length(nuc.vert)+k])*dt
        end

        # enveMovements = Vector{Vec{3,Float64}}(undef,length(nuc.vert))
        # chroMovements = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)

        # while true

        #     enveMovements = enveVelocities.*dt;
        #     if all(norm.(enveMovements) .< 0.05)
        #         chroMovements = chroVelocities.*dt;
        #         if all(norm.(chroMovements) .< 0.05)
        #             break
        #         else
        #             dt = dt/2;
        #         end
        #     else
        #         dt = dt/2;
        #     end
        # end

        # nuc.vert .+= enveMovements
        # chro.vert .+= chroMovements

        plane = plane - 0.1*dt;

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