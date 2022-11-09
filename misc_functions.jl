function export_pipette_mesh(folderName, pip)

    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(pip.tri,1))
    for i = 1:size(pip.tri,1)
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, pip.tri[i,:])
    end

    vtk_save(vtk_grid(".\\results\\"*folderName*"\\pipette", [getindex.(pip.vert,1) getindex.(pip.vert,2) getindex.(pip.vert,3)]', triCells))

end

function get_strand_vectors!(chro,spar)

    for i = 1:spar.chromatinNumber
        chro.vectors[i] = chro.strandVert[i][2:end] .- chro.strandVert[i][1:end-1];
        chro.vectorNorms[i] = norm.(chro.vectors[i])
    end

end

function export_data(nuc,chro,spar,ex,ext,intTime,simset)

    if ex.exportData
        if simset.timeStepProgress == 0
            if cmp(simset.simType,"MA") == 0
                push!(ext[2],maximum(getindex.(nuc.vert,1)));
            elseif cmp(simset.simType,"MM") == 0
                push!(ext[2],nuc.vert[ext[1].rightmostVertex][1] - nuc.vert[ext[1].leftmostVertex][1]);
            end
        end

        if mod(intTime,ex.step) == 0 && simset.timeStepProgress == 0

            exportNumber = string(Int64(intTime/ex.step+1));

            vtk_grid(".\\results\\"*ex.folderName*"\\nucl_" * lpad(exportNumber,4,"0"), [getindex.(nuc.vert,1) getindex.(nuc.vert,2) getindex.(nuc.vert,3)]', ex.enveCells) do vtk
                if cmp(simset.simType,"MA") == 0 
                    vtk["Aspiration forces", VTKPointData()] = [getindex.(aspiration,1) getindex.(aspiration,2) getindex.(aspiration,3)]'
                    vtk["Pipette repulsion forces", VTKPointData()] = [getindex.(repulsion,1) getindex.(repulsion,2) getindex.(repulsion,3)]'
                end
                # vtk["Curvature"] = nuc.curvatures;
                vtk["Element normals", VTKCellData()] = [getindex.(nuc.triangleNormalUnitVectors,1) getindex.(nuc.triangleNormalUnitVectors,2) getindex.(nuc.triangleNormalUnitVectors,3)]'
                vtk["Volume forces", VTKPointData()] = [getindex.(nuc.forces.volume,1) getindex.(nuc.forces.volume,2) getindex.(nuc.forces.volume,3)]'
                vtk["Area forces", VTKPointData()] = [getindex.(nuc.forces.area,1) getindex.(nuc.forces.area,2) getindex.(nuc.forces.area,3)]'
                vtk["Elastic forces", VTKPointData()] = [getindex.(nuc.forces.elastic,1) getindex.(nuc.forces.elastic,2) getindex.(nuc.forces.elastic,3)]'
                # vtk["enveRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.envelopeRepulsion,1) getindex.(nuc.forces.envelopeRepulsion,2) getindex.(nuc.forces.envelopeRepulsion,3)]'
                vtk["chroRepulsion forces", VTKPointData()] = [getindex.(nuc.forces.chromationRepulsion,1) getindex.(nuc.forces.chromationRepulsion,2) getindex.(nuc.forces.chromationRepulsion,3)]'

            end
            vtk_grid(".\\results\\"*ex.folderName*"\\chro_" * lpad(exportNumber,4,"0"), [getindex.(chro.vert,1) getindex.(chro.vert,2) getindex.(chro.vert,3)]', ex.chroCells) do vtk
                vtk["line_id"] = 1:spar.chromatinNumber
                vtk["Linear Forces", VTKPointData()] = [getindex.(chro.forces.linear,1) getindex.(chro.forces.linear,2) getindex.(chro.forces.linear,3)]'
                vtk["Bending Forces", VTKPointData()] = [getindex.(chro.forces.bending,1) getindex.(chro.forces.bending,2) getindex.(chro.forces.bending,3)]'
                vtk["chroRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.chroRepulsion,1) getindex.(chro.forces.chroRepulsion,2) getindex.(chro.forces.chroRepulsion,3)]'
                vtk["enveRepulsion Forces", VTKPointData()] = [getindex.(chro.forces.enveRepulsion,1) getindex.(chro.forces.enveRepulsion,2) getindex.(chro.forces.enveRepulsion,3)]'
                vtk["LAD forces", VTKPointData()] = [getindex.(chro.forces.ladChroForces,1) getindex.(chro.forces.ladChroForces,2) getindex.(chro.forces.ladChroForces,3)]'
                # vtk["Fluc Forces", VTKPointData()] = [getindex.(fluctuationForces,1) getindex.(fluctuationForces,2) getindex.(fluctuationForces,3)]'
                # vtk["Movement", VTKPointData()] = [dt*vX[length(nuc.vert)+1:end] dt*vY[length(nuc.vert)+1:end] dt*vZ[length(nuc.vert)+1:end]]'
            end

            if simset.simType == "VRC"

                vrcCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,length(ext.tri))
                for i = eachindex(ext.tri)
                    vrcCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, ext.tri[i]);
                end


                vtk_grid(".\\results\\"*ex.folderName*"\\replcomp_" * lpad(exportNumber,4,"0"), [getindex.(ext.vert,1) getindex.(ext.vert,2) getindex.(ext.vert,3)]', vrcCells) do vtk
                    vtk["Index", VTKPointData()] = 1:length(ext.vert)
                end
            end

            # export lads
        
            tempVert = [chro.vert[ex.ladChroVertices] ; nuc.vert[ex.ladEnveVertices]]
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

        end
    end
end


function get_nuclear_properties!(nuc, chro, simset, spar, ext)
   
    # form the trees for the vertex distance search
    simset.envelopeTree = KDTree(nuc.vert);
    simset.chromatinTree = KDTree(chro.vert);

    # get various nuclear properties needed in the solution
    get_strand_vectors!(chro,spar)
    get_edge_vectors!(nuc);
    nuc.triangleAreas = get_area!(nuc)
    get_voronoi_areas!(nuc);
    get_triangle_normals!(nuc);
    get_area_unit_vectors!(nuc);
    # get_local_curvatures!(nuc);

    if simset.simType == "VRC"
        add_repl_comp_triangles!(ext)
    end

end

function save_specific_data!(nuc,ext,simset)
   
    # save data for analysis
    if simset.timeStepProgress == 0

        if cmp(simset.simType,"MA") == 0 # for MA, save maximum x-coordinate
            push!(ext[2],maximum(getindex.(nuc.vert,1)));
        elseif cmp(simset.simType,"MM") == 0 # for MM, save distance between the manipulated vertices
            push!(ext[2],nuc.vert[ext[1].rightmostVertex][1] - nuc.vert[ext[1].leftmostVertex][1]);
        end
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

function get_crosslinks!(nuc,chro,simset,spar)

    constant = 0.001;

    changesDone = false

    # remove crosslinks
    nLinked = length(chro.crosslinks)
    probs = rand(nLinked)
    for i = nLinked:-1:1

        if probs[i] < spar.crosslinkingUnbindingProbability*spar.maxDt*simset.timeStepMultiplier

            chro.crosslinked[chro.crosslinks[i][1]] = 0
            chro.crosslinked[chro.crosslinks[i][2]] = 0

            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][1], length(nuc.vert) + chro.crosslinks[i][2]] = 0
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][2], length(nuc.vert) + chro.crosslinks[i][1]] = 0
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][1], length(nuc.vert) + chro.crosslinks[i][1]] -= spar.laminaFriction*constant
            simset.frictionMatrix[length(nuc.vert) + chro.crosslinks[i][2], length(nuc.vert) + chro.crosslinks[i][2]] -= spar.laminaFriction*constant

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
        if distance[1] <= 0.5
            closestVerts[i] = closest[1]
            possiblyLinking[i] = true
        end
    end
    
    possibleLinkingIdx = findall(possiblyLinking)

    
    for i = possibleLinkingIdx

        if chro.crosslinked[i] == 0 && chro.crosslinked[closestVerts[i]] == 0

            if rand() < spar.crosslinkingBindingProbability*spar.maxDt*simset.timeStepMultiplier

                push!(chro.crosslinks, [i, closestVerts[i]])
                chro.crosslinked[i] = 1
                chro.crosslinked[closestVerts[i]] = 1

                simset.frictionMatrix[length(nuc.vert) + i, length(nuc.vert) + closestVerts[i]] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + closestVerts[i], length(nuc.vert) + i] -= spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + i, length(nuc.vert) + i] += spar.laminaFriction*constant
                simset.frictionMatrix[length(nuc.vert) + closestVerts[i], length(nuc.vert) + closestVerts[i]] += spar.laminaFriction*constant

                changesDone = true
            end
        end
    end
    if changesDone
        iLU = ilu(simset.frictionMatrix, τ = spar.iLUCutoff)
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

function create_replication_compartment(spar)

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

    return replComp

end

function add_repl_comp_triangles!(ext)
    get_edge_vectors!(nuc);
    nuc.triangleAreas = get_area!(nuc)
    get_voronoi_areas!(nuc);
    get_triangle_normals!(nuc);
    get_area_unit_vectors!(nuc);
end