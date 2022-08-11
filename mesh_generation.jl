function generate_pipette_mesh()

    pip = pipetteType();

    mesh = load("./test.stl")
    pip.tri = zeros(Int64, length(mesh), 3);

    for i = eachindex(mesh)

        temptri = zeros(Int64,3)

        for j = 1:3
            coords = mesh[i][j]
            pointComparison = isapprox.(pip.x,coords[1]) .&& isapprox.(pip.y,coords[2]) .&& isapprox.(pip.z,coords[3])
            if any(pointComparison)
                temptri[j] = findall(pointComparison)[1]
            else
                temptri[j] = length(pip.x) + 1;
                append!(pip.x,coords[1]);
                append!(pip.y,coords[2]);
                append!(pip.z,coords[3]);
            end
        end

        pip.tri[i,:] = temptri;

    end

    pip.vertexTri = fill(Int[], length(pip.x), 1);

    for i = 1:length(pip.x)
        aa = findall(pip.tri  .== i);
        pip.vertexTri[i] = [i[1] for i in aa];
    end

    pip.z,pip.x = pip.x,pip.z

    pip.x .+= 1.0583;

    pip.y = 0.3.*pip.y;
    pip.z = 0.3.*pip.z;

    return pip
end