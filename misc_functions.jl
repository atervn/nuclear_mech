function export_pipette_mesh(folderName, pip)

    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(pip.tri,1))
    for i = 1:size(pip.tri,1)
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, pip.tri[i])
    end

    vtk_save(vtk_grid(".\\results\\"*folderName*"\\pipette", [getindex.(pip.vert,1) getindex.(pip.vert,2) getindex.(pip.vert,3)]', triCells))

end

function get_strand_vectors!(chro,spar)

    for i = 1:spar.chromatinNumber
        chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];
        chro.vectorNorms[i] = norm.(chro.vectors[i])
    end

end

function export_data(enve,chro,spar,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(enve.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],enve.vert[ext[1].rightmostVertex][1] - enve.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));
            
            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(enve.vert,1) getindex.(enve.vert,2) getindex.(enve.vert,3)]', ex.enveCells) do vtk
                # if cmp(simset.simType,"MA") == 0 
                #     vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                #     vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                # end
                # vtk["Curvature"] = enve.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(enve.triangleNormalUnitVectors,1) getindex.(enve.triangleNormalUnitVectors,2) getindex.(enve.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(enve.forces.volume,1) getindex.(enve.forces.volume,2) getindex.(enve.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(enve.forces.area,1) getindex.(enve.forces.area,2) getindex.(enve.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(enve.forces.elastic,1) getindex.(enve.forces.elastic,2) getindex.(enve.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(enve.forces.envelopeRepulsion,1) getindex.(enve.forces.envelopeRepulsion,2) getindex.(enve.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(enve.forces.chromationRepulsion,1) getindex.(enve.forces.chromationRepulsion,2) getindex.(enve.forces.chromationRepulsion,3)]'
                vtk["Bending forces", VTKPointData()] = [getindex.(enve.forces.bending,1) getindex.(enve.forces.bending,2) getindex.(enve.forces.bending,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(enve.forces.ladEnveForces,1) getindex.(enve.forces.ladEnveForces,2) getindex.(enve.forces.ladEnveForces,3)]'

            end
            vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(chro.forces.ladChroForces,1) getindex.(chro.forces.ladChroForces,2) getindex.(chro.forces.ladChroForces,3)]'
                # vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                # vtk["Movement", VTKPointData()] = [dt*vX[length(enve.vert)+1:end] dt*vY[length(enve.vert)+1:end] dt*vZ[length(enve.vert)+1:end]]'
            end

            # export lads
        
            tempVert = [chro.vert[ex.ladChroVertices] ; enve.vert[ex.ladEnveVertices]]
            vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
                vtk["LAD ID"] = ex.ladIdx
            end

            crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(chro.crosslinks))
            for i = 1:length(chro.crosslinks)
                crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(chro.crosslinks)+i]);
            end
            tempVert = [chro.vert[getindex.(chro.crosslinks,1)] ; chro.vert[getindex.(chro.crosslinks,2)]];
            vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
            end

            writedlm(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0") * ".csv", [getindex.(chro.crosslinks,1) getindex.(chro.crosslinks,2)],',')

            if simset.adh.adherent
                writedlm(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",[simset.adh.topPlane, simset.adh.bottomPlane],',')
            end
        end
    end
end

function export_data(enve,chro,repl,spar,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(enve.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],enve.vert[ext[1].rightmostVertex][1] - enve.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));
            
            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(enve.vert,1) getindex.(enve.vert,2) getindex.(enve.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = enve.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(enve.triangleNormalUnitVectors,1) getindex.(enve.triangleNormalUnitVectors,2) getindex.(enve.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(enve.forces.volume,1) getindex.(enve.forces.volume,2) getindex.(enve.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(enve.forces.area,1) getindex.(enve.forces.area,2) getindex.(enve.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(enve.forces.elastic,1) getindex.(enve.forces.elastic,2) getindex.(enve.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(enve.forces.envelopeRepulsion,1) getindex.(enve.forces.envelopeRepulsion,2) getindex.(enve.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(enve.forces.chromationRepulsion,1) getindex.(enve.forces.chromationRepulsion,2) getindex.(enve.forces.chromationRepulsion,3)]'

            end
            vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(chro.forces.ladChroForces,1) getindex.(chro.forces.ladChroForces,2) getindex.(chro.forces.ladChroForces,3)]'
                # vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                # vtk["Movement", VTKPointData()] = [dt*vX[length(enve.vert)+1:end] dt*vY[length(enve.vert)+1:end] dt*vZ[length(enve.vert)+1:end]]'
            end

            vrcCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(repl.tri))
            for i = eachindex(repl.tri)
                vrcCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, repl.tri[i]);
            end

            vtk_grid(".\\results\\"*ex.folderName*"\\replcomp_" * lpad(exportNumber,4,"0"), [getindex.(repl.vert,1) getindex.(repl.vert,2) getindex.(repl.vert,3)]', vrcCells) do vtk
                vtk["Total forces", VTKPointData()] = [getindex.(repl.forces.total,1) getindex.(repl.forces.total,2) getindex.(repl.forces.total,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(repl.forces.elastic,1) getindex.(repl.forces.elastic,2) getindex.(repl.forces.elastic,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(repl.forces.area,1) getindex.(repl.forces.area,2) getindex.(repl.forces.area,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(repl.forces.volume,1) getindex.(repl.forces.volume,2) getindex.(repl.forces.volume,3)]'
                vtk["Bending forces", VTKPointData()] = [getindex.(repl.forces.bending,1) getindex.(repl.forces.bending,2) getindex.(repl.forces.bending,3)]'
            end

            # export lads
        
            tempVert = [chro.vert[ex.ladChroVertices] ; enve.vert[ex.ladEnveVertices]]
            vtk_grid(".\\results\\"*ex.folderName*"\\lads_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', ex.ladCells) do vtk
                vtk["LAD ID"] = ex.ladIdx
            end

            crossLinkCells =  Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,length(chro.crosslinks))
            for i = 1:length(chro.crosslinks)
                crossLinkCells[i] = MeshCell(PolyData.Lines(), [i, length(chro.crosslinks)+i]);
            end
            tempVert = [chro.vert[getindex.(chro.crosslinks,1)] ; chro.vert[getindex.(chro.crosslinks,2)]];
            vtk_grid(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0"), [getindex.(tempVert,1) getindex.(tempVert,2) getindex.(tempVert,3)]', crossLinkCells) do vtk
            end

            writedlm(".\\results\\"*ex.folderName*"\\crosslinks_" * lpad(exportNumber,4,"0") * ".csv", [getindex.(chro.crosslinks,1) getindex.(chro.crosslinks,2)],',')

            if simset.adh.adherent
                writedlm(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",[simset.adh.topPlane, simset.adh.bottomPlane],',')
            end
        end
    end
end

function export_data(enve,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(enve.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],enve.vert[ext[1].rightmostVertex][1] - enve.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));
            
            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(enve.vert,1) getindex.(enve.vert,2) getindex.(enve.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = enve.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(enve.triangleNormalUnitVectors,1) getindex.(enve.triangleNormalUnitVectors,2) getindex.(enve.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(enve.forces.volume,1) getindex.(enve.forces.volume,2) getindex.(enve.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(enve.forces.area,1) getindex.(enve.forces.area,2) getindex.(enve.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(enve.forces.elastic,1) getindex.(enve.forces.elastic,2) getindex.(enve.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(enve.forces.envelopeRepulsion,1) getindex.(enve.forces.envelopeRepulsion,2) getindex.(enve.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(enve.forces.chromationRepulsion,1) getindex.(enve.forces.chromationRepulsion,2) getindex.(enve.forces.chromationRepulsion,3)]'

            end
           
            if simset.adh.adherent
                writedlm(".\\results\\"*ex.folderName*"\\planes_" * lpad(exportNumber,4,"0")*".csv",[simset.adh.topPlane, simset.adh.bottomPlane],',')
            end
        end
    end
end

function get_nuclear_properties!(enve, chro, simset, spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_triangle_normals!(enve);
    get_area_unit_vectors!(enve);
    # get_local_curvatures!(nuc.enve);

end

function get_nuclear_properties!(enve, chro, repl, simset, spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_triangle_normals!(enve);
    get_area_unit_vectors!(enve);
    # get_local_curvatures!(nuc.enve);

    add_repl_comp_triangles!(repl,spar)
    repl.tree = KDTree(repl.vert);

end

function get_nuclear_properties!(enve, simset, spar)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(enve.vert);

    # get various nuclear properties needed in the solution
    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    get_voronoi_areas!(enve);
    get_triangle_normals!(enve);
    get_area_unit_vectors!(enve);
    # get_local_curvatures!(nuc.enve);

end

function check_simulation_type(simType)

    # check that the simulation type is correct
    if !any(simType .== ["MA", "MM", "INIT"])
        printstyled("Unknown simulation type"; color = :blue)
        return true
    end

    return false
end

function get_crosslinks!(enve, chro, simset,spar)

    constant = 0.001;

    changesDone = false

    # remove crosslinks
    nLinked = length(chro.crosslinks)
    probs = rand(nLinked)
    for i = nLinked:-1:1

        if probs[i] < spar.crosslinkingUnbindingProbability*spar.maxDt*simset.timeStepMultiplier

            chro.crosslinked[chro.crosslinks[i][1]] = 0
            chro.crosslinked[chro.crosslinks[i][2]] = 0

            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][1], length(enve.vert) + chro.crosslinks[i][2]] = 0
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][2], length(enve.vert) + chro.crosslinks[i][1]] = 0
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][1], length(enve.vert) + chro.crosslinks[i][1]] -= spar.laminaFriction*constant
            simset.frictionMatrix[length(enve.vert) + chro.crosslinks[i][2], length(enve.vert) + chro.crosslinks[i][2]] -= spar.laminaFriction*constant

            chro.crosslinks = chro.crosslinks[1:end .!= i]

            changesDone = true

        end
    end

    dropzeros!(simset.frictionMatrix)

    # form crosslinks
    notCrosslinked = findall(chro.crosslinked .== 0)

    closestVerts = zeros(Int64,spar.chromatinLength*spar.chromatinNumber)
    possiblyLinking =  zeros(Bool,spar.chromatinLength*spar.chromatinNumber)
    for i = notCrosslinked
        closest,distance = knn(simset.chromatinTree, chro.vert[i],1,true,j -> any(j .== chro.neighbors[i]))
        if length(distance) > 0
            if distance[1] <= 0.5
                closestVerts[i] = closest[1]
                possiblyLinking[i] = true
            end
        end
    end
    
    possibleLinkingIdx = findall(possiblyLinking)
    
    for i = possibleLinkingIdx

        if chro.crosslinked[i] == 0 && chro.crosslinked[closestVerts[i]] == 0

            if rand() < spar.crosslinkingBindingProbability*spar.maxDt*simset.timeStepMultiplier

                push!(chro.crosslinks, [i, closestVerts[i]])
                chro.crosslinked[i] = 1
                chro.crosslinked[closestVerts[i]] = 1

                simset.frictionMatrix[length(enve.vert) + i, length(enve.vert) + closestVerts[i]] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(enve.vert) + closestVerts[i], length(enve.vert) + i] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(enve.vert) + i, length(enve.vert) + i] += spar.laminaFriction*constant
                simset.frictionMatrix[length(enve.vert) + closestVerts[i], length(enve.vert) + closestVerts[i]] += spar.laminaFriction*constant

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
        stepTime = now() - simset.timeStepTiming
        next!(simset.prog, showvalues = [(:stepTime,stepTime)])
        intTime += 1
        simset.timeStepTiming = now()
    end
    
    return intTime

end

function post_export(ex,simset,ext)

    if ex.exportData
        if cmp(simset.simType, "MA") == 0
            writedlm(".\\results\\" * ex.folderName * "\\maxX.csv", ext[2], ',')
            dL = ext[2] .- minimum(ext[2])

       elseif cmp(simset.simType, "MM") == 0
            writedlm(".\\results\\" * ex.folderName * "\\nuclearLength.csv", ext[2], ',')
        end
    end

end

function create_replication_compartment(enve,spar)

    repl = replicationCompartmentType()

    radius = 0.1*spar.freeNucleusRadius
    repl = get_icosaherdon!(repl,radius)

    repl = subdivide_mesh!(repl,radius,2)

    enveMax = maximum(getindex.(enve.vert,3))
    enveMin = minimum(getindex.(enve.vert,3))
    centerZ = (enveMax + enveMin)/2

    for i = 1:length(repl.vert)
        repl.vert[i] += Vec(0.,0.,centerZ)
    end

    repl.forces.volume = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.area = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.bending = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.elastic = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))
    repl.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))

    repl = setup_shell_data(repl)

    get_edge_vectors!(repl);
    repl.triangleAreas = get_area!(repl)
    get_voronoi_areas!(repl);
    get_triangle_normals!(repl);
    get_area_unit_vectors!(repl);

    repl.frictionMatrix = get_repl_comp_friction_matrix(repl,spar)
    repl.iLU = ilu(repl.frictionMatrix, τ = spar.iLUCutoff)

    get_edge_vectors!(enve);
    enve.triangleAreas = get_area!(enve)
    repl.baseArea = mean(enve.triangleAreas)

    return repl

end

function add_repl_comp_triangles!(repl,spar)
    
    get_edge_vectors!(repl);
    repl.triangleAreas = get_area!(repl)

    tooLargeAreas = zeros(length(repl.triangleAreas))
    for i = eachindex(repl.triangleAreas)

        tooLargeAreas[i] = repl.triangleAreas[i] > repl.baseArea

    end

    if sum(tooLargeAreas)/length(repl.triangleAreas) > 0.5
        
        repl = subdivide_mesh!(repl,0,1)

        repl.forces.volume = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.area = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.bending = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.elastic = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.chromationRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))
        repl.forces.envelopeRepulsion = Vector{Vec{3,Float64}}(undef, length(repl.vert))

        repl = setup_shell_data(repl)
        get_edge_vectors!(repl);
        repl.triangleAreas = get_area!(repl)

        repl.frictionMatrix = get_repl_comp_friction_matrix(repl,spar)
        repl.iLU = ilu(repl.frictionMatrix, τ = spar.iLUCutoff)

    end
    get_voronoi_areas!(repl);
    get_triangle_normals!(repl);
    get_area_unit_vectors!(repl);
end