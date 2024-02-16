function export_pipette_mesh(folderName, pip)

    # create an array to store triangular mesh cells
    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(pip.tri,1))
    
    # iterate over each row of the pip.tri array
    for i = 1:size(pip.tri,1)

        # create a MeshCell object of type VTK_TRIANGLE with the current triangle indices
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, pip.tri[i])
    end

    # generate a vtk_grid object with vertex and cell information
    vtkGridObject = vtk_grid(".\\results\\" * folderName * "\\pipette", [getindex.(pip.vert, 1) getindex.(pip.vert, 2) getindex.(pip.vert, 3)]', triCells)
    
    # save the vtk_grid_object as a VTK file
    vtk_save(vtkGridObject)

end

function get_strand_vectors!(chro,spar)

    # iterate over each chromatin strand
    for i = 1:spar.chromatinNumber

        # calculate the strand vectors by subtracting consecutive vertices
        chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];

        # calculate the norms of the strand vectors
        chro.vectorNorms[i] = norm.(chro.vectors[i])
    end

end

function export_data(enve::envelopeType,chro::chromatinType,spar::scaledParametersType,ex,ext,intTime,simset)

     # check if data export is enabled
    if ex.exportData

        # save data for experiment analysis
        save_analysis_data(simset,ext,enve)

        # export data at full time and export time steps 
        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            # determine the export number for file naming
            exportNumber = string(Int64(intTime/ex.step+1));
            
            # export envelope data
            export_envelope_data(enve,ex,simset,exportNumber,spar)

            # export chromatin data
            export_chromatin_data(enve,chro,spar,ex,exportNumber)

            if simset.simType == "AFM"

                cell = [MeshCell(VTKCellTypes.VTK_VERTEX,[1]),MeshCell(VTKCellTypes.VTK_VERTEX,[2])]

                points = [ext.beadPosition[1] ext.beadPosition[2] ext.beadPosition[3] ; ext.topPosition[1] ext.topPosition[2] ext.topPosition[3]]'

                distance = norm(ext.topPosition - ext.beadPosition);

                springConstant = 0.05/spar.viscosity*spar.scalingTime;

                cantileverForce = -springConstant*(distance - ext.normDistance)

                vtk_grid(".\\results\\"*ex.folderName*"\\afm_" * lpad(exportNumber,4,"0"), points, cell) do vtk
                    # assign line IDs
                    vtk["point_id"] = 1:2
                    vtk["Force on bead"] = [ext.forceOnBead, ext.forceOnBead]
                    vtk["Spring force"] = [cantileverForce, cantileverForce]
                end
            end

        end
    end
end

function export_data(enve,chro,repl,spar,ex,ext,intTime,simset)

     # check if data export is enabled
    if ex.exportData

        # save data for experiment analysis
        save_analysis_data(simset,ext,enve)

        # export data at full time and export time steps 
        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            # determine the export number for file naming
            exportNumber = string(Int64(intTime/ex.step+1));
            
            # export envelope data
            export_envelope_data(enve,ex,simset,exportNumber,spar)

            # export chromatin data
            export_chromatin_data(enve,chro,spar,ex,exportNumber)

            # export replication compartment data
            export_replication_compartment_data(repl,ex,exportNumber)

            if simset.simType == "AFM"

                cell = [MeshCell(VTKCellTypes.VTK_VERTEX,[1]),MeshCell(VTKCellTypes.VTK_VERTEX,[2])]

                points = [ext.beadPosition[1] ext.beadPosition[2] ext.beadPosition[3] ; ext.topPosition[1] ext.topPosition[2] ext.topPosition[3]]'

                distance = norm(ext.topPosition - ext.beadPosition);

                springConstant = 0.05/spar.viscosity*spar.scalingTime;

                cantileverForce = -springConstant*(distance - ext.normDistance)

                vtk_grid(".\\results\\"*ex.folderName*"\\afm_" * lpad(exportNumber,4,"0"), points, cell) do vtk
                    # assign line IDs
                    vtk["point_id"] = 1:2
                    vtk["Force on bead"] = [ext.forceOnBead, ext.forceOnBead]
                    vtk["Spring force"] = [cantileverForce, cantileverForce]
                end
            end

        end
    end
end

function export_data(enve,spar,ex,ext,intTime,simset)

    # check if data export is enabled
    if ex.exportData

        # save data for experiment analysis
        save_analysis_data(simset,ext,enve)

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            # determine the export number for file naming
            exportNumber = string(Int64(intTime/ex.step+1));
            
            # export envelope data
            export_envelope_data(enve,ex,simset,exportNumber,spar)

            if simset.simType == "AFM"

                cell = [MeshCell(VTKCellTypes.VTK_VERTEX,[1]),MeshCell(VTKCellTypes.VTK_VERTEX,[2])]

                points = [ext.beadPosition[1] ext.beadPosition[2] ext.beadPosition[3] ; ext.topPosition[1] ext.topPosition[2] ext.topPosition[3]]'

                vtk_grid(".\\results\\"*ex.folderName*"\\afm_" * lpad(exportNumber,4,"0"), points, cell) do vtk
                    # assign line IDs
                    vtk["point_id"] = 1:2
                end
            end
        end
    end
end

function get_nuclear_properties!(enve, chro, simset, spar, intTime, intMaxTime)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_shell_normals!(enve);
    get_area_unit_vectors!(enve);

    if simset.newVolumeSimulation

        volume = get_volume!(enve);

        if volume >= enve.targetVolume
            direction = -1
        else
            direction = 1
        end

        enve.normalLengths = enve.normalLengths.*(1 .+ direction*0.01*spar.dt*simset.timeStepMultiplier)

        enve.normalTriangleAreas = enve.normalTriangleAreas.*(1 .+ direction*0.01*2.1*spar.dt*simset.timeStepMultiplier)

    end

    enve.volume = get_volume!(enve)

end

function get_nuclear_properties!(enve, chro, repl::replicationCompartmentType, simset, spar, intTime)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_shell_normals!(enve);
    get_area_unit_vectors!(enve);

    # get properties of the replication compartment
    add_repl_comp_triangles!(repl,spar,simset)
    repl.tree = KDTree(repl.vert);

    if simset.newVolumeSimulation

        volume = get_volume!(enve);

        if volume >= enve.targetVolume
            direction = -1
        else
            direction = 1
        end

        if abs(volume - enve.targetVolume) < 1
            enve.targetVolumeReached = true
        end

        enve.normalLengths = enve.normalLengths.*(1 .+ direction*0.01*spar.dt*simset.timeStepMultiplier)

        enve.normalTriangleAreas = enve.normalTriangleAreas.*(1 .+ direction*0.01*2.1*spar.dt*simset.timeStepMultiplier)

        replVolume = get_volume!(repl)

        if replVolume - spar.replTargetVolume > -1
            repl.growthDone = true
        end
        
        if  enve.targetVolumeReached && repl.growthDone
            simset.stopSimulation = true
        end

    end

    enve.volume = get_volume!(enve)


end

function get_nuclear_properties!(enve, simset, spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);

    # get various nuclear properties needed in the solution
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_shell_normals!(enve);
    get_area_unit_vectors!(enve);

    if simset.newVolumeSimulation

        volume = get_volume!(enve);

        if volume >= enve.targetVolume
            direction = -1
        else
            direction = 1
        end

        enve.normalLengths = enve.normalLengths.*(1 .+ direction*0.01*spar.dt*simset.timeStepMultiplier)

        enve.normalTriangleAreas = enve.normalTriangleAreas.*(1 .+ direction*0.01*2.1*spar.dt*simset.timeStepMultiplier)

    end

end

function check_simulation_type(simType)

    # check if the given simulation type does not match any of the valid simulation types: "MA", "MM", "INIT"
    if !any(simType .== ["MA", "MM", "INIT", "AFM"])

        # print a message indicating that the simulation type is unknown
        printstyled("Unknown simulation type"; color = :blue)

        # return true to indicate that the simulation type is invalid
        return true
    end

    # if the simulation type matches one of the valid simulation types, return false to indicate that it is valid
    return false
end


function get_crosslinks!(enve, chro, simset, spar)

    # initialize a variable to track if any changes are made
    changesDone = false

    # get number of crosslinks
    nLinked = length(chro.crosslinks)

    # get random numbers for each crosslink
    probs = rand(nLinked)

    # iterate through the crosslinks
    for i = nLinked:-1:1

        # check if the crosslink should be removed based on probability
        if probs[i] < spar.crosslinkingUnbindingProbability*spar.dt*simset.timeStepMultiplier

            # set crosslinked values to 0
            chro.crosslinked[chro.crosslinks[i][1]] = 0
            chro.crosslinked[chro.crosslinks[i][2]] = 0

            # update frictionMatrix values
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][1], length(enve.vert) + chro.crosslinks[i][2]] = 0
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][2], length(enve.vert) + chro.crosslinks[i][1]] = 0
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][1], length(enve.vert) + chro.crosslinks[i][1]] -= spar.crosslinkFriction
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][2], length(enve.vert) + chro.crosslinks[i][2]] -= spar.crosslinkFriction

            # remove the crosslink
            chro.crosslinks = chro.crosslinks[1:end .!= i]

            # set changesDone flag to true
            changesDone = true

        end
    end

    # remove zero values from frictionMatrix
    dropzeros!(simset.frictionMatrix)

    # find all vertices not crosslinked or lads
    notCrosslinked = findall(chro.crosslinked .== 0)

    # init closest vertices and possibly linking vertices
    closestVerts = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)
    possiblyLinking =  zeros(Bool,spar.chromatinLength*spar.chromatinNumber)

    # iterate through the notCrosslinked
    for i = notCrosslinked

        # find the closest vertex in the chromatinTree
        closest,distance = knn(simset.chromatinTree, chro.vert[i],1,true,j -> any(j .== chro.neighbors[i]))
        if length(distance) > 0

            # check if the distance is within a threshold
            if distance[1] <= spar.maxCrosslinkDistance
                closestVerts[i] = closest[1]
                possiblyLinking[i] = true
            end
        end
    end
    
    # get the indices of the possibly linking vertices
    possibleLinkingIdx = findall(possiblyLinking)
    
    # iterate through the possibly linking vertices
    for i = possibleLinkingIdx

        # check if the vertices are not already crosslinked
        if chro.crosslinked[i] == 0 && chro.crosslinked[closestVerts[i]] == 0

            # check if a crosslink should be formed based on probability
            if rand() < spar.crosslinkingBindingProbability*spar.dt*simset.timeStepMultiplier

                # add the crosslink and update values
                push!(chro.crosslinks, [i, closestVerts[i]])
                chro.crosslinked[i] = 1
                chro.crosslinked[closestVerts[i]] = 1

                # update frictionMatrix values
                simset.frictionMatrix[length(enve.vert) + i, length(enve.vert) + closestVerts[i]] -= spar.crosslinkFriction
                simset.frictionMatrix[length(enve.vert) + closestVerts[i], length(enve.vert) + i] -= spar.crosslinkFriction
                simset.frictionMatrix[length(enve.vert) + i, length(enve.vert) + i] += spar.crosslinkFriction
                simset.frictionMatrix[length(enve.vert) + closestVerts[i], length(enve.vert) + closestVerts[i]] += spar.crosslinkFriction

                # set changesDone flag to true
                changesDone = true
            end
        end
    end

    # update iLU if changes are made
    # if changesDone
    #     simset.iLU = ilu(simset.frictionMatrix, τ = spar.iLUCutoff)
    # end

end

function progress_time!(simset,intTime,enve)
    
    # checks if the time step progress is 0, indicating the beginning of a new time step
    if simset.timeStepProgress == 0

        # calculates the duration of the previous time step by subtracting the current time from the time at the start of the time step
        stepTime = now() - simset.timeStepTiming
        
        nucleusVolume = get_volume!(enve)

        # updates the progress bar by displaying the duration of the previous time step
        next!(simset.prog, showvalues = [(:stepTime,stepTime), (:volume,round(nucleusVolume,digits=2))])
        
        # increments the internal time counter by 1
        intTime += 1
        simset.timeStepTiming = now()
    end
    
    # updates the timing reference to the current time, indicating the start of the next time step
    return intTime

end

function post_export(ex,simset,ext)

    # checks if data export is enabled in the export settings
    if ex.exportData

        # checks if the simulation type is "MA"
        if cmp(simset.simType, "MA") == 0
            
            # writes the maximum X-coordinate values from the `ext` array to a CSV file
            writedlm(".\\results\\" * ex.folderName * "\\maxX.csv", ext[2], ',')

        # checks if the simulation type is "MM"
        elseif cmp(simset.simType, "MM") == 0

            # writes the nuclear length values from the `ext` array to a CSV file
            writedlm(".\\results\\" * ex.folderName * "\\nuclearLength.csv", ext[2], ',')
            
            # writes the nuclear forces from the `ext` array to a CSV file
            writedlm(".\\results\\" * ex.folderName * "\\nuclearForces.csv", ext[3], ',')
        end
    end
end

function create_replication_compartment(repl,enve,spar)

    # calculate the radius as 0.1 times the free nucleus radius
    radius = spar.replSizeMultiplier*spar.freeNucleusRadius

    # generate an icosahedral mesh for the replication compartment
    repl = get_icosaherdon!(repl,radius)

    # subdivide the mesh of the replication compartment
    repl = subdivide_mesh!(repl,radius,2)

    # calculate the maximum and minimum values of envelope z-coordinates
    enveMax = maximum(getindex.(enve.vert,3))
    enveMin = minimum(getindex.(enve.vert,3))

    # calculate the center position along the Z-axis
    centerZ = (enveMax + enveMin)/2

    # translate the vertices of the replication compartment to align with the center position
    for i = 1:length(repl.vert)
        repl.vert[i] += Vec(0.,0.,centerZ)
    end
    
    return repl

end

function add_repl_comp_triangles!(repl,spar,simset)
    
    # calculate various properties and vectors for the replication compartment
    get_edge_vectors!(repl);
    repl.triangleAreas = get_area!(repl)

    # identify triangles with areas larger than the base area
    tooLargeAreas = repl.triangleAreas .> repl.baseArea

    # if more than 50 % of the triangles have areas larger than the base area
    if sum(tooLargeAreas)/length(repl.triangleAreas) > 0.5
        
        # subdivide the mesh of the replication compartment
        repl = subdivide_mesh!(repl,0,1)

        # initialize various force vectors for each vertex in the replication compartment
        repl.forces.volume = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.area = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.bending = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.elastic = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))


        # set up shell data for the replication compartment
        
        volumeTemp = repl.normalVolume
        repl = setup_shell_data(repl,simset.initType,"repl")
        repl.normalVolume = volumeTemp 
        
        get_edge_vectors!(repl);
        repl.triangleAreas = get_area!(repl)

        # calculate the friction matrix for the replication compartment
        repl.frictionMatrix = get_repl_friction_matrix(repl,spar)

        # perform iLU (incomplete LU) factorization of the friction matrix with a specified cutoff value
        repl.iLU = ilu(repl.frictionMatrix, τ = spar.iLUCutoff)

    end

    # calculate various properties and vectors for the replication compartment
    get_voronoi_areas!(repl);
    get_shell_normals!(repl);
    get_area_unit_vectors!(repl);

end

function export_parameters(ipar,ex)

    # open the file for writing
    f = open(".\\results\\"*ex.folderName*"\\parameters.txt", "w")

    # get the parameter names of the ipar object
    fields = propertynames(ipar)

    # iterate over the indices of the fields
    for i = eachindex(fields)

        # write the field name and value to the file
        write(f,string(fields[i]) *","* string(getfield(ipar,fields[i])) * "\n" )

    end

    # close the file
    close(f)

end

function export_envelope_data(enve,ex,simset,exportNumber,spar)

    # export envelope data
    vtk_grid(".\\results\\"*ex.folderName*"\\enve_" * lpad(exportNumber,4,"0"), [getindex.(enve.vert,1) getindex.(enve.vert,2) getindex.(enve.vert,3)]', ex.enveCells) do vtk
                
        # if micropipette aspiration
        if cmp(simset.simType,"MA") == 0
            # export aspiration forces
            vtk["Aspiration forces", VTKPointData()] = [getindex.(enve.forces.aspiration,1) getindex.(enve.forces.aspiration,2) getindex.(enve.forces.aspiration,3)]'
            # export pipette repulsion forces
            vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(enve.forces.pipetteRepulsion,1) getindex.(enve.forces.pipetteRepulsion,2) getindex.(enve.forces.pipetteRepulsion,3)]'
        end

        # export element areas
        vtk["Element areas", VTKCellData()] = enve.triangleAreas
        # export element normals
        vtk["Element normals", VTKCellData()] = [getindex.(enve.triangleNormalUnitVectors,1) getindex.(enve.triangleNormalUnitVectors,2) getindex.(enve.triangleNormalUnitVectors,3)]'
        # export volume forces
        vtk["Volume forces", VTKPointData()] = [getindex.(enve.forces.volume,1) getindex.(enve.forces.volume,2) getindex.(enve.forces.volume,3)]'
        # export area forces
        vtk["Area forces", VTKPointData()] = [getindex.(enve.forces.area,1) getindex.(enve.forces.area,2) getindex.(enve.forces.area,3)]'
        # export elastic forces
        vtk["Elastic forces", VTKPointData()] = [getindex.(enve.forces.elastic,1) getindex.(enve.forces.elastic,2) getindex.(enve.forces.elastic,3)]'
        # export chromatin repulsion forces
        vtk["chroRepulsion forces", VTKPointData()] = [getindex.(enve.forces.chromationRepulsion,1) getindex.(enve.forces.chromationRepulsion,2) getindex.(enve.forces.chromationRepulsion,3)]'
        # export bending forces
        vtk["Bending forces", VTKPointData()] = [getindex.(enve.forces.bending,1) getindex.(enve.forces.bending,2) getindex.(enve.forces.bending,3)]'
        # export LAD forces
        vtk["LAD forces", VTKPointData()] = [getindex.(enve.forces.ladEnveForces,1) getindex.(enve.forces.ladEnveForces,2) getindex.(enve.forces.ladEnveForces,3)]'
        vtk["AFM forces", VTKPointData()] = [getindex.(enve.forces.afmRepulsion,1) getindex.(enve.forces.afmRepulsion,2) getindex.(enve.forces.afmRepulsion,3)]'
        vtk["Cytoskeleton forces", VTKPointData()] = [getindex.(enve.forces.planeRepulsion,1) getindex.(enve.forces.planeRepulsion,2) getindex.(enve.forces.planeRepulsion,3)]'
        # export total forces
        vtk["Total forces", VTKPointData()] = [getindex.(enve.forces.total,1) getindex.(enve.forces.total,2) getindex.(enve.forces.total,3)]'
    end

    # export planes if adhesion is enabled
    if simset.adh.adherent
        df = DataFrame(x = [0, 0], y = [0, 0], z = [simset.adh.topPlane, simset.adh.bottomPlane], rad = [spar.cytoskeletonPlaneRadius,0])
        CSV.write(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",df)
    end

    if simset.exportNormalLengths
        writedlm(".\\results\\"*ex.folderName*"\\normalLengths_" * lpad(exportNumber,4,"0") *".csv", enve.normalLengths.*spar.scalingLength,',')
    end

end

function export_chromatin_data(enve,chro,spar,ex,exportNumber)

    # export chromatin data
    vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
        # assign line IDs
        vtk["line_id"] = 1:spar.chromatinNumber
        # export linear forces
        vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
        # export bending forces
        vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
        # export chromatin repulsion forces
        vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
        # export envelope repulsion forces
        vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
        # export LAD forces
        vtk["LAD forces", VTKPointData()] = [getindex.(chro.forces.ladChroForces,1) getindex.(chro.forces.ladChroForces,2) getindex.(chro.forces.ladChroForces,3)]'
        # crosslink forces
        vtk["Crosslink forces", VTKPointData()] = [getindex.(chro.forces.crosslink,1) getindex.(chro.forces.crosslink,2) getindex.(chro.forces.crosslink,3)]'
        # crosslink forces
        vtk["Total forces", VTKPointData()] = [getindex.(chro.forces.total,1) getindex.(chro.forces.total,2) getindex.(chro.forces.total,3)]'
    end

    # export LADs
    tempVert = [chro.vert[ex.ladChroVertices] ; enve.vert[ex.ladEnveVertices]]
    vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
        vtk["LAD ID"] = ex.ladIdx
    end

    # export crosslinks
    crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(chro.crosslinks))
    for i = 1:length(chro.crosslinks)
        crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(chro.crosslinks)+i]);
    end
    tempVert = [chro.vert[getindex.(chro.crosslinks,1)] ; chro.vert[getindex.(chro.crosslinks,2)]];
    vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
    end

    # write crosslinks to CSV
    writedlm(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0") * ".csv", [getindex.(chro.crosslinks,1) getindex.(chro.crosslinks,2)],',')

end

function export_replication_compartment_data(repl,ex,exportNumber)

    # create MeshCell data for the repicalition compartment
    vrcCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(repl.tri))
    for i = eachindex(repl.tri)
        vrcCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, repl.tri[i]);
    end

    # export replication compartment data
    vtk_grid(".\\results\\"*ex.folderName*"\\repl_" * lpad(exportNumber,4,"0"), [getindex.(repl.vert,1) getindex.(repl.vert,2) getindex.(repl.vert,3)]', vrcCells) do vtk
        # export total forces
        vtk["Total forces", VTKPointData()] = [getindex.(repl.forces.total,1) getindex.(repl.forces.total,2) getindex.(repl.forces.total,3)]'
        # export elastic forces
        vtk["Elastic forces", VTKPointData()] = [getindex.(repl.forces.elastic,1) getindex.(repl.forces.elastic,2) getindex.(repl.forces.elastic,3)]'
        # export area forces
        vtk["Area forces", VTKPointData()] = [getindex.(repl.forces.area,1) getindex.(repl.forces.area,2) getindex.(repl.forces.area,3)]'
        # export volume forces
        vtk["Volume forces", VTKPointData()] = [getindex.(repl.forces.volume,1) getindex.(repl.forces.volume,2) getindex.(repl.forces.volume,3)]'
        # exoirt bending forces
        vtk["Bending forces", VTKPointData()] = [getindex.(repl.forces.bending,1) getindex.(repl.forces.bending,2) getindex.(repl.forces.bending,3)]'
        # export chromatin repulsion forces
        vtk["Chromatin repulsion forces", VTKPointData()] = [getindex.(repl.forces.chromationRepulsion,1) getindex.(repl.forces.chromationRepulsion,2) getindex.(repl.forces.chromationRepulsion,3)]'
        # export emvelope repulsion forces
        vtk["Envelope repulsion forces", VTKPointData()] = [getindex.(repl.forces.envelopeRepulsion,1) getindex.(repl.forces.envelopeRepulsion,2) getindex.(repl.forces.envelopeRepulsion,3)]'
    end
end

function save_analysis_data(simset,ext,enve)

    # perform export only at each full time step
    if simset.timeStepProgress == 0

        # check simulation type
        if cmp(simset.simType,"MA") == 0

            # record the maximum x-coordinate of envelope vertices
            push!(ext[2],maximum(getindex.(enve.vert,1)));

        elseif cmp(simset.simType,"MM") == 0
        
            # calculate the distance between the rightmost and leftmost envelope vertices
            push!(ext[2],enve.vert[ext[1].rightmostVertex][1] - enve.vert[ext[1].leftmostVertex][1]);
        
            # calculate the 
            push!(ext[3],norm(enve.forces.micromanipulation))

        end
    end

end

function print_error(error)

    # check if the error is not an instance of InterruptException
    if !(error isa InterruptException)

        # print "Simulation failed..." message in red color
        printstyled("Simulation failed...\n"; color=:red)

        # rethrow the error to propagate it further
        rethrow(error)

    else

        # print "Simulation stopped" message in red color
        printstyled("Simulation stopped\n"; color=:red)

    end
end

function move_micromanipulator!(ext, spar, simset, intTime)
    
    if simset.simType == "MM"

        mm = ext[1]

        time = Float64(intTime + simset.timeStepProgress)*spar.dt
        
        prevTime = maximum(ext[1].pipetteMovements[ext[1].pipetteMovements[:,1] .<= time,1])
        nextTime = minimum(ext[1].pipetteMovements[ext[1].pipetteMovements[:,1] .> time,1])

        prevInd = findall(ext[1].pipetteMovements[:,1] .== prevTime)[1]
        nextInd = findall(ext[1].pipetteMovements[:,1] .== nextTime)[1]

        if prevTime == time

            # get the displacement from the initial position
            ext[1].pipettePosition = Vec(ext[1].pipetteMovements[prevInd,2],0.,0.)
        
        # otherwise
        else

            # interpolate the displacement from the initial position
            intermediatePosition = mm.pipetteMovements[prevInd,2] + (mm.pipetteMovements[nextInd,2] - mm.pipetteMovements[prevInd,2])*(time - mm.pipetteMovements[prevInd,1])/(mm.pipetteMovements[nextInd,1] - mm.pipetteMovements[prevInd,1])
            ext[1].pipettePosition = Vec(intermediatePosition,0.,0.)
        end
    end
end