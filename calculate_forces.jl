function get_volume_forces(nucleus)

    nucleusVolume = get_volume!(nucleus);

    pressure = -log10(nucleusVolume/nucleus.normalVolume);

    volumeForces = zeros(Float64,length(nucleus.x),3);

    for i = 1:length(nucleus.x)
        forceMagnitude = pressure*sum(nucleus.voronoiAreas[i]);
        volumeForces[i,:] = forceMagnitude.*nucleus.surfaceNormalUnitVectors[i,:];
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