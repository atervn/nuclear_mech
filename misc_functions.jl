function get_friction_matrix(nuc,par)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.x),length(nuc.x));

    for j = 1:length(nuc.x)
        frictionMatrix[j,j] = par.viscosity + par.laminaFriction*length(nuc.neighbors[j]);
        frictionMatrix[j,nuc.neighbors[j]] .= -par.laminaFriction;
    end

    return frictionMatrix

end

function read_parameters(filePath)

    f = open(filePath)

    inputParameters = Dict{String,Float64}()

    while !eof(f)
        line = split(readline(f), ',')
        inputParameters[line[1]] = parse(Float64,line[2])
    end

    return inputParameters
end

function calculate_model_parameters(ipar)

    mpar.laminaStiffness
    laminaStiffness,
areaCompressionStiffness,
bendingStiffness,
nucleusBulkModulus,
repulsionConstant,
viscosity,
laminaFriction,
    return mpar

end
