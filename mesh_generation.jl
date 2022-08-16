function generate_pipette_mesh()

    pip = pipetteType();

    mesh = load("./test.stl")
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
                # append!(pip.x,coords[1]);
                # append!(pip.y,coords[2]);
                # append!(pip.z,coords[3]);
            end
        end

        pip.tri[i,:] = temptri;

    end

    pip.vertexTri = fill(Int[], length(pip.vert), 1);

    for i = 1:length(pip.vert)
        aa = findall(pip.tri  .== i);
        pip.vertexTri[i] = [i[1] for i in aa];
    end

    # pip.z,pip.x = pip.x,pip.z

    xOffset = 1.0583;
    radius = 0.3; 

    for i = eachindex(pip.vert)
        pip.vert[i] = pip.vert[i] + Vec(xOffset,0.,0.)
        pip.vert[i] = pip.vert[i] .* Vec(1.,radius,radius)
    end

# pip.x .+= 1.0583;

#     pip.y = 0.3.*pip.y;
#     pip.z = 0.3.*pip.z;

    return pip
end