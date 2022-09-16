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