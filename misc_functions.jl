function get_friction_matrix(nuc,spar)

    frictionMatrix = sparse(Int64[],Int64[],Float64[],length(nuc.x),length(nuc.x));

    for j = 1:length(nuc.x)
        frictionMatrix[j,j] = 1 + spar.laminaFriction*length(nuc.neighbors[j]);
        frictionMatrix[j,nuc.neighbors[j]] .= -spar.laminaFriction;
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

    spar.repulsionConstant = ipar.repulsionConstant/ipar.viscosity*ipar.scalingTime/ipar.scalingLength;

    spar.repulsionDistance = ipar.repulsionDistance/ipar.scalingLength;

    return spar

end
