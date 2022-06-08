function get_volume_forces(nucleus)

    nucleusVolume = get_volume!(nucleus);

    pressure = -log10(nucleusVolume/nucleus.normalVolume);

    volumeForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:length(nucleus.x)
        forceMagnitude = pressure*sum(nucleus.voronoiAreas[i]);
        volumeForces[i,:] = forceMagnitude.*nucleus.vertexNormalUnitVectors[i,:];
    end

    return volumeForces

end

function get_area_forces(nucleus)

    nucleusArea = get_area!(nucleus);

    forceMagnitude = (nucleusArea - nucleus.normalArea)/nucleus.normalArea;

    areaForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:length(nucleus.x)

        areaForces[i,:] = forceMagnitude.*sum(nucleus.areaUnitVectors[i],dims=1);

    end

    return areaForces

end

function get_bending_forces(nucleus)

    bendingForces = zeros(Float64,length(nucleus.x),3);
    angles = zeros(Float64,size(nucleus.edges,1),1);
    for i = 1:size(nucleus.edges,1)

        normalVector1 = nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,1],:];
        normalVector2 = nucleus.triangleNormalUnitVectors[nucleus.edgesTri[i,2],:];

        dotProduct = normalVector1[1]*normalVector2[1] + normalVector1[2]*normalVector2[2] + normalVector1[3]*normalVector2[3];
        angles[i] = acosd(dotProduct);
    end

    return angles

end