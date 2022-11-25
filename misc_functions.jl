function export_pipette_mesh(folderName, pip)

    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(pip.tri,1))
    for i = 1:size(pip.tri,1)
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, pip.tri[i,:])
    end

    vtk_save(vtk_grid(".\\results\\"*folderName*"\\pipette", [getindex.(pip.vert,1) getindex.(pip.vert,2) getindex.(pip.vert,3)]', triCells))

end

function get_strand_vectors!(nuc,spar)

    for i = 1:spar.chromatinNumber
        nuc.chro.vectors[i] = nuc.chro.strandVert[i][2:end] .- nuc.chro.strandVert[i][1:end-1];
        nuc.chro.vectorNorms[i] = norm.(nuc.chro.vectors[i])
    end

end

function export_data(nuc,spar,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(nuc.enve.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],nuc.enve.vert[ext[1].rightmostVertex][1] - nuc.enve.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));
            
            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(nuc.enve.vert,1) getindex.(nuc.enve.vert,2) getindex.(nuc.enve.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.enve.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.enve.triangleNormalUnitVectors,1) getindex.(nuc.enve.triangleNormalUnitVectors,2) getindex.(nuc.enve.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(nuc.enve.forces.volume,1) getindex.(nuc.enve.forces.volume,2) getindex.(nuc.enve.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(nuc.enve.forces.area,1) getindex.(nuc.enve.forces.area,2) getindex.(nuc.enve.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.enve.forces.elastic,1) getindex.(nuc.enve.forces.elastic,2) getindex.(nuc.enve.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(nuc.enve.forces.envelopeRepulsion,1) getindex.(nuc.enve.forces.envelopeRepulsion,2) getindex.(nuc.enve.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(nuc.enve.forces.chromationRepulsion,1) getindex.(nuc.enve.forces.chromationRepulsion,2) getindex.(nuc.enve.forces.chromationRepulsion,3)]'

            end
            vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(nuc.chro.vert,1) getindex.(nuc.chro.vert,2) getindex.(nuc.chro.vert,3)]', ex.chroCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(nuc.chro.forces.linear,1) getindex.(nuc.chro.forces.linear,2) getindex.(nuc.chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(nuc.chro.forces.bending,1) getindex.(nuc.chro.forces.bending,2) getindex.(nuc.chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(nuc.chro.forces.chroRepulsion,1) getindex.(nuc.chro.forces.chroRepulsion,2) getindex.(nuc.chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(nuc.chro.forces.enveRepulsion,1) getindex.(nuc.chro.forces.enveRepulsion,2) getindex.(nuc.chro.forces.enveRepulsion,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(nuc.chro.forces.ladChroForces,1) getindex.(nuc.chro.forces.ladChroForces,2) getindex.(nuc.chro.forces.ladChroForces,3)]'
                # vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                # vtk["Movement", VTKPointData()] = [dt*vX[length(nuc.enve.vert)+1:end] dt*vY[length(nuc.enve.vert)+1:end] dt*vZ[length(nuc.enve.vert)+1:end]]'
            end

            if simset.simType == "VRC"

                vrcCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(nuc.repl.tri))
                for i = eachindex(nuc.repl.tri)
                    vrcCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.repl.tri[i]);
                end


                vtk_grid(".\\results\\"*ex.folderName*"\\replcomp_" * lpad(exportNumber,4,"0"), [getindex.(nuc.repl.vert,1) getindex.(nuc.repl.vert,2) getindex.(nuc.repl.vert,3)]', vrcCells) do vtk
                    vtk["Total forces", VTKPointData()] = [getindex.(nuc.repl.forces.total,1) getindex.(nuc.repl.forces.total,2) getindex.(nuc.repl.forces.total,3)]'
                    vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.repl.forces.elastic,1) getindex.(nuc.repl.forces.elastic,2) getindex.(nuc.repl.forces.elastic,3)]'
                    vtk["Area forces", VTKPointData()] = [getindex.(nuc.repl.forces.area,1) getindex.(nuc.repl.forces.area,2) getindex.(nuc.repl.forces.area,3)]'
                    vtk["Volume forces", VTKPointData()] = [getindex.(nuc.repl.forces.volume,1) getindex.(nuc.repl.forces.volume,2) getindex.(nuc.repl.forces.volume,3)]'
                    vtk["Bending forces", VTKPointData()] = [getindex.(nuc.repl.forces.bending,1) getindex.(nuc.repl.forces.bending,2) getindex.(nuc.repl.forces.bending,3)]'
                end
            end

            # export lads
        
            tempVert = [nuc.chro.vert[ex.ladChroVertices] ; nuc.enve.vert[ex.ladEnveVertices]]
            vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
                vtk["LAD ID"] = ex.ladIdx
            end

            crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(nuc.chro.crosslinks))
            for i = 1:length(nuc.chro.crosslinks)
                crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(nuc.chro.crosslinks)+i]);
            end
            tempVert = [nuc.chro.vert[getindex.(nuc.chro.crosslinks,1)] ; nuc.chro.vert[getindex.(nuc.chro.crosslinks,2)]];
            vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
            end

            writedlm(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0") * ".csv", [getindex.(nuc.chro.crosslinks,1) getindex.(nuc.chro.crosslinks,2)],',')

            if nuc.adh.adherent
                writedlm(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",[nuc.adh.topPlane, nuc.adh.bottomPlane],',')
            end
        end
    end
end

function get_nuclear_properties!(nuc, simset, spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(nuc.enve.vert);
    simset.chromatinTree = KDTree(nuc.chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(nuc,spar)
    get_edge_vectors!(nuc.enve);
    nuc.enve.triangleAreas = get_area!(nuc.enve)
    get_voronoi_areas!(nuc.enve);
    get_triangle_normals!(nuc.enve);
    get_area_unit_vectors!(nuc.enve);
    # get_local_curvatures!(nuc.enve);

    if simset.simType == "VRC"
        add_repl_comp_triangles!(nuc,spar)
        nuc.repl.tree = KDTree(nuc.repl.vert);
    end

end

function check_simulation_type(simType)

    # check that the simulation type is correct
    if !any(simType .== ["MA", "MM", "PC", "INIT", "VRC"])
        printstyled("Unknown simulation type"; color = :blue)
        return true
    end

    return false
end

function get_crosslinks!(nuc,simset,spar)

    constant = 0.001;

    changesDone = false

    # remove crosslinks
    nLinked = length(nuc.chro.crosslinks)
    probs = rand(nLinked)
    for i = nLinked:-1:1

        if probs[i] < spar.crosslinkingUnbindingProbability*spar.maxDt*simset.timeStepMultiplier

            nuc.chro.crosslinked[nuc.chro.crosslinks[i][1]] = 0
            nuc.chro.crosslinked[nuc.chro.crosslinks[i][2]] = 0

            simset.frictionMatrix[length(nuc.enve.vert) + nuc.chro.crosslinks[i][1], length(nuc.enve.vert) + nuc.chro.crosslinks[i][2]] = 0
            simset.frictionMatrix[length(nuc.enve.vert) + nuc.chro.crosslinks[i][2], length(nuc.enve.vert) + nuc.chro.crosslinks[i][1]] = 0
            simset.frictionMatrix[length(nuc.enve.vert) + nuc.chro.crosslinks[i][1], length(nuc.enve.vert) + nuc.chro.crosslinks[i][1]] -= spar.laminaFriction*constant
            simset.frictionMatrix[length(nuc.enve.vert) + nuc.chro.crosslinks[i][2], length(nuc.enve.vert) + nuc.chro.crosslinks[i][2]] -= spar.laminaFriction*constant

            nuc.chro.crosslinks = nuc.chro.crosslinks[1:end .!= i]

            changesDone = true

        end
    end

    dropzeros!(simset.frictionMatrix)

    # form crosslinks
    notCrosslinked = findall(nuc.chro.crosslinked .== 0)
    closestVerts = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)
    possiblyLinking =  zeros(Bool,spar.chromatinLength*spar.chromatinNumber)

    for i = notCrosslinked
        closest,distance = knn(simset.chromatinTree, nuc.chro.vert[i],1,true,j -> any(j .== nuc.chro.neighbors[i]))
        if distance[1] <= 0.5
            closestVerts[i] = closest[1]
            possiblyLinking[i] = true
        end
    end
    
    possibleLinkingIdx = findall(possiblyLinking)

    
    for i = possibleLinkingIdx

        if nuc.chro.crosslinked[i] == 0 && nuc.chro.crosslinked[closestVerts[i]] == 0

            if rand() < spar.crosslinkingBindingProbability*spar.maxDt*simset.timeStepMultiplier

                push!(nuc.chro.crosslinks, [i, closestVerts[i]])
                nuc.chro.crosslinked[i] = 1
                nuc.chro.crosslinked[closestVerts[i]] = 1

                simset.frictionMatrix[length(nuc.enve.vert) + i, length(nuc.enve.vert) + closestVerts[i]] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.enve.vert) + closestVerts[i], length(nuc.enve.vert) + i] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.enve.vert) + i, length(nuc.enve.vert) + i] += spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.enve.vert) + closestVerts[i], length(nuc.enve.vert) + closestVerts[i]] += spar.laminaFriction*constant

                changesDone = true
            end
        end
    end
    if changesDone
        simset.iLU = ilu(simset.frictionMatrix, τ = spar.iLUCutoff)
    end

end

function progress_time!(simset,intTime)

    if simset.timeStepProgress == 0
        next!(simset.prog)
        intTime += 1
    end

    return intTime

end

function post_export(ex,simset,ext)

    if ex.exportData
        if cmp(simset.simType, "MA") == 0
            writedlm(".\\results\\" * ex.folderName * "\\maxX.csv", maxX, ',')
            dL = ext[2] .- minimum(ext[2])

            J = 2 * pi .* dL ./ (3 * 2.1 * 3 * 1)

            plot(10*dt:dt:maxT*dt, J[11:end], yaxis=:log, xaxis=:log, xlim=(0.1, 200), ylim=(0.01, 10))
        elseif cmp(simset.simType, "MM") == 0
            writedlm(".\\results\\" * ex.folderName * "\\nuclearLength.csv", ext[2], ',')
        end
    end

end

function create_replication_compartment(nuc,spar)

    replComp = replicationCompartmentType()

    radius = 0.1*spar.freeNucleusRadius
    replComp = get_icosaherdon!(replComp,radius)
    replComp = subdivide_mesh!(replComp,radius,2)

    replComp.forces.volume = Vector{Vec{3,Float64}}(undef, length(replComp.vert))
    replComp.forces.area = Vector{Vec{3,Float64}}(undef, length(replComp.vert))
    replComp.forces.bending = Vector{Vec{3,Float64}}(undef, length(replComp.vert))
    replComp.forces.elastic = Vector{Vec{3,Float64}}(undef, length(replComp.vert))
    replComp.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(replComp.vert))
    replComp.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(replComp.vert))

    replComp = setup_shell_data(replComp)

    get_edge_vectors!(replComp);
    replComp.triangleAreas = get_area!(replComp)
    get_voronoi_areas!(replComp);
    get_triangle_normals!(replComp);
    get_area_unit_vectors!(replComp);

    replComp.frictionMatrix = get_repl_comp_friction_matrix(replComp,spar)
    replComp.iLU = ilu(replComp.frictionMatrix, τ = spar.iLUCutoff)

    get_edge_vectors!(nuc.enve);
    nuc.enve.triangleAreas = get_area!(nuc.enve)
    replComp.baseArea = mean(nuc.enve.triangleAreas)

    return replComp

end

function add_repl_comp_triangles!(nuc,spar)
    
    get_edge_vectors!(nuc.repl);
    nuc.repl.triangleAreas = get_area!(nuc.repl)

    tooLargeAreas = zeros(length(nuc.repl.triangleAreas))
    for i = eachindex(nuc.repl.triangleAreas)

        tooLargeAreas[i] = nuc.repl.triangleAreas[i] > nuc.repl.baseArea

    end

    if sum(tooLargeAreas)/length(nuc.repl.triangleAreas) > 0.5
        
        println("blob")
        nuc.repl = subdivide_mesh!(nuc.repl,0,1)

        nuc.repl.forces.volume = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))
        nuc.repl.forces.area = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))
        nuc.repl.forces.bending = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))
        nuc.repl.forces.elastic = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))
        nuc.repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))
        nuc.repl.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(nuc.repl.vert))

        nuc.repl = setup_shell_data(nuc.repl)
        get_edge_vectors!(nuc.repl);
        nuc.repl.triangleAreas = get_area!(nuc.repl)

        nuc.repl.frictionMatrix = get_repl_comp_friction_matrix(nuc,spar)
        nuc.repl.iLU = ilu(nuc.repl.frictionMatrix, τ = spar.iLUCutoff)

    end
    get_voronoi_areas!(nuc.repl);
    get_triangle_normals!(nuc.repl);
    get_area_unit_vectors!(nuc.repl);
end

function get_nuclear_properties_adh_init!(nuc, simset)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(nuc.enve.vert);

    # get various nuclear properties needed in the solution
    get_edge_vectors!(nuc.enve);
    nuc.enve.triangleAreas = get_area!(nuc.enve)
    get_voronoi_areas!(nuc.enve);
    get_triangle_normals!(nuc.enve);
    get_area_unit_vectors!(nuc.enve);
    # get_local_curvatures!(nuc.enve);

end

function export_data_adh_init(nuc,ex,intTime,simset,planeRepulsion,ext)

    if ex.exportData

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));
            
            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(nuc.enve.vert,1) getindex.(nuc.enve.vert,2) getindex.(nuc.enve.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.enve.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.enve.triangleNormalUnitVectors,1) getindex.(nuc.enve.triangleNormalUnitVectors,2) getindex.(nuc.enve.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(nuc.enve.forces.volume,1) getindex.(nuc.enve.forces.volume,2) getindex.(nuc.enve.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(nuc.enve.forces.area,1) getindex.(nuc.enve.forces.area,2) getindex.(nuc.enve.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.enve.forces.elastic,1) getindex.(nuc.enve.forces.elastic,2) getindex.(nuc.enve.forces.elastic,3)]'
                vtk["Bending forces", VTKPointData()] = [getindex.(nuc.enve.forces.bending,1) getindex.(nuc.enve.forces.bending,2) getindex.(nuc.enve.forces.bending,3)]'
                vtk["Plane Compression", VTKPointData()] = [getindex.(planeRepulsion,1) getindex.(planeRepulsion,2) getindex.(planeRepulsion,3)]'
                planeRepulsion
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(nuc.enve.forces.envelopeRepulsion,1) getindex.(nuc.enve.forces.envelopeRepulsion,2) getindex.(nuc.enve.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(nuc.enve.forces.chromationRepulsion,1) getindex.(nuc.enve.forces.chromationRepulsion,2) getindex.(nuc.enve.forces.chromationRepulsion,3)]'
            end

            writedlm(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",[ext[1], ext[2]],',')

        end
    end
end