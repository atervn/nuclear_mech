function generate_pipette_mesh()

    pip = pipetteType();

    mesh = load("./test2.stl")
    pip.tri = zeros(Int64, length(mesh), 3);

    for i = eachindex(mesh)

        temptri = zeros(Int64,3)

        for j = 1:3
            coords = mesh[i][j]
            
            pointComparison = zeros(Bool,length(pip.vert));
            for k = 1:length(pip.vert)
                if isapprox(pip.vert[k],Vec(coords[3],coords[2],coords[1]))
                    pointComparison[k] = true;
                end
            end

            if any(pointComparison)
                temptri[j] = findall(pointComparison)[1]
            else
                temptri[j] = length(pip.vert) + 1;
                push!(pip.vert,Vec(coords[3], coords[2], coords[1]))
            end
        end

        pip.tri[i,:] = temptri;

    end

    pip.vertexTri = fill(Int[], length(pip.vert), 1);

    for i = 1:length(pip.vert)
        aa = findall(pip.tri  .== i);
        pip.vertexTri[i] = [i[1] for i in aa];
    end

    xOffset = 5.99;
    radius = 3; 

    for i = eachindex(pip.vert)
        pip.vert[i] = pip.vert[i] .* Vec(radius,radius,radius)
        pip.vert[i] = pip.vert[i] + Vec(xOffset,0.,0.)
    end

    return pip
end