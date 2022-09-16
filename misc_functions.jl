function get_friction_matrix(nuc,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.vert)+spar.chromatinLength*spar.chromatinNumber,length(nuc.vert)+spar.chromatinLength*spar.chromatinNumber);

    for j = 1:length(nuc.vert)
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(nuc.neighbors[j]);
        frictionMatrix[j,nuc.neighbors[j]] .= -spar.laminaFriction;
    end

    for j = length(nuc.vert)+1:length(nuc.vert)+spar.chromatinLength*spar.chromatinNumber
        frictionMatrix[j,j] = 1
    end

    return frictionMatrix
end

function read_parameters(ipar,filePath)

    f = open(filePath)

    while !eof(f)
        line = split(readline(f), ',')
        setproperty!(ipar,Symbol(line[1]),parse(Float64,line[2]))
    end

    return ipar
end

function get_model_parameters(ipar,spar,nuc)

    spar.laminaStiffness = 2/sqrt(3)*ipar.laminaYoung*ipar.laminaThickness;
    spar.laminaStiffness = spar.laminaStiffness/ipar.viscosity*ipar.scalingTime;

    spar.laminaFriction = ipar.laminaFriction/ipar.viscosity;

    spar.areaCompressionStiffness = ipar.areaCompressionModulus/(mean(nuc.normalLengths).*ipar.scalingLength);
    spar.areaCompressionStiffness = spar.areaCompressionStiffness/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    spar.bendingStiffness = ipar.laminaYoung*ipar.laminaThickness^3/(12*(1-ipar.poissonsRatio^2));
    spar.bendingStiffness = spar.bendingStiffness/ipar.viscosity*ipar.scalingTime/ipar.scalingLength^2;

    spar.bulkModulus = ipar.bulkModulus/ipar.viscosity*ipar.scalingTime*ipar.scalingLength;

    spar.repulsionConstant = ipar.repulsionConstant/ipar.viscosity*ipar.scalingTime#/ipar.scalingLength;

    spar.repulsionDistance = ipar.repulsionDistance/ipar.scalingLength;

    spar.freeNucleusRadius = ipar.freeNucleusRadius/ipar.scalingLength;

    spar.chroVertexDistance = ipar.chroVertexDistance/ipar.scalingLength;

    spar.chromatinNumber = ipar.chromatinNumber
    spar.chromatinLength = ipar.chromatinLength

    return spar

end

function setup_export(folderName,nuc,chro,spar)

    try
        mkdir(".\\results\\"*folderName)
    catch
        for i = 1:1000
            try
                mkdir(".\\results\\"*folderName*"_"*string(i))
                folderName = folderName*"_"*string(i)
                break
            catch
            end
        end
    end
    

    triCells = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef,size(nuc.tri,1))
    for i = 1:size(nuc.tri,1)
        triCells[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, nuc.tri[i,:]);
    end

    lineCells = Vector{MeshCell{PolyData.Lines, Vector{Int64}}}(undef,spar.chromatinNumber)
    for i = 1:spar.chromatinNumber
        lineCells[i] = MeshCell(PolyData.Lines(), chro.strandIdx[i]);
    end
    
    return triCells, lineCells, folderName

end

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

function setup_micromanipulation(nuc)

    mm = micromanipulationType();

    mm.leftmostVertex = argmin(getindex.(nuc.vert,1));
    mm.rightmostVertex = argmax(getindex.(nuc.vert,1));

    mm.leftNeighbors = nuc.neighbors[mm.leftmostVertex];
    mm.rightNeighbors = nuc.neighbors[mm.rightmostVertex];

    mm.leftmostVertexPosition = nuc.vert[mm.leftmostVertex]
    mm.leftNeigborPositions = nuc.vert[mm.leftNeighbors]

    return mm

end

export_data(nuc,chro,spar,folderName,step)

if mod(t,step) == 0
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