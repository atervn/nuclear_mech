function generate_pipette_mesh()

    pip = pipetteType();

    mesh = load("./pipette_mesh.stl")
    pip.tri = Vector{Vector{Int64}}(undef,length(mesh))

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

        pip.tri[i] = temptri;

    end

    

    xOffset = 4.1;
    radius = 3; 

    for i = eachindex(pip.vert)
        pip.vert[i] = pip.vert[i] .* Vec(radius,radius,radius)
        pip.vert[i] = pip.vert[i] + Vec(xOffset,0.,0.)
    end

    # a vector for the pairs
    pip.edges = Vector{Vector{Int64}}(undef,0);
    
    # a vector for the connected neighboring vertices
    pip.neighbors = fill(Int[], length(pip.vert));

    # matrix form of the trignle vector
    triMatrix = [getindex.(pip.tri,1) getindex.(pip.tri,2) getindex.(pip.tri,3)];

    pip.vertexTri = fill(Int[], length(pip.vert), 1);

    for i = 1:length(pip.vert)
        aa = findall(triMatrix  .== i);
        pip.vertexTri[i] = [i[1] for i in aa];
    end

    # go through the vertices
    for i = 1:length(pip.vert)
        
        # find in which triangles vertex i is included
        hasVertex = findall(triMatrix .== i);
        # initialize a vector for the neighbors
        neighbors = Array{Int64}(undef, 0);
        
        # go through the triangles
        for j = eachindex(hasVertex)
            append!(neighbors,pip.tri[hasVertex[j][1]][pip.tri[hasVertex[j][1]] .!= i ]);
        end
        
        # get the unique neighbors and save them
        unique!(neighbors)
        pip.neighbors[i] = neighbors;

        # add the connections to the edges matrix
        for j = eachindex(neighbors)
            push!(pip.edges, [i, neighbors[j]]);
        end  
    end

    pip.mirrorEdges = zeros(Int64,length(pip.edges));
    pip.firstEdges = zeros(Int64,length(pip.edges));

    for i = eachindex(pip.edges)
            
        mirrorIdx = findall(getindex.(pip.edges,1) .== pip.edges[i][2] .&& getindex.(pip.edges,2) .== pip.edges[i][1]);
        
        pip.mirrorEdges[i] = mirrorIdx[1];
        pip.mirrorEdges[mirrorIdx[1]] = i;
        if pip.firstEdges[i] == 0 && pip.firstEdges[mirrorIdx[1]] == 0 
            pip.firstEdges[i] = 1;
        end
    end

    pip.edgesTri = Vector{Vector{Int64}}(undef,length(pip.edges));

    pip.edgesTri = Vector{Vector{Int64}}(undef,length(pip.edges));

    bufferedDoubleTriMatrix = [triMatrix[:,3] triMatrix triMatrix[:,1]]

    for i = eachindex(pip.edges)

        pip.edgesTri[i] = [0,0]

        firstEdges = findall(triMatrix .== pip.edges[i][1])

        for j = 1:length(firstEdges)

            if bufferedDoubleTriMatrix[firstEdges[j][1],firstEdges[j][2]+2] == pip.edges[i][2]
                pip.edgesTri[i][1] = firstEdges[j][1]
            elseif bufferedDoubleTriMatrix[firstEdges[j][1],firstEdges[j][2]] == pip.edges[i][2]
                pip.edgesTri[i][2] = firstEdges[j][1]
            end

        end

    
    end

    pip.edgeVectors = Vector{Vec{3,Float64}}(undef, length(pip.edges));

    for i = eachindex(pip.edges)
        if pip.firstEdges[i] == 1
            pip.edgeVectors[i] = pip.vert[pip.edges[i][2]] - pip.vert[pip.edges[i][1]];
            pip.edgeVectors[pip.mirrorEdges[i]] = -pip.edgeVectors[i]
        end
    end

    pip.triEdge1 = Vector{Int64}(undef,length(pip.tri))
    pip.triEdge2 = Vector{Int64}(undef,length(pip.tri))
    

    for i = eachindex(pip.tri)
        pip.triEdge1[i] = findall(getindex.(pip.edges,1) .== pip.tri[i][1] .&& getindex.(pip.edges,2) .== pip.tri[i][2])[1]
        pip.triEdge2[i] = findall(getindex.(pip.edges,1) .== pip.tri[i][1] .&& getindex.(pip.edges,2) .== pip.tri[i][3])[1]
    end

    triangleNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.tri));

    for i = eachindex(pip.tri)
        normalVector = cross(pip.edgeVectors[pip.triEdge1[i]],pip.edgeVectors[pip.triEdge2[i]])
        triangleNormalUnitVectors[i] = normalVector/norm(normalVector);
    end

    pip.triangleNormalUnitVectors = triangleNormalUnitVectors;

    edgeNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.edges));

    for i = eachindex(pip.edges)

        vector = pip.triangleNormalUnitVectors[pip.edgesTri[i][1]] + pip.triangleNormalUnitVectors[pip.edgesTri[i][2]]
        edgeNormalUnitVectors[i] = vector./norm(vector);

    end
    
    pip.edgeNormalUnitVectors = edgeNormalUnitVectors;


    vertexNormalUnitVectors = Vector{Vec{3,Float64}}(undef, length(pip.vert));
    
    for i = 1:length(pip.vert)

        vector = sum(pip.triangleNormalUnitVectors[pip.vertexTri[i]]);
        vertexNormalUnitVectors[i] = vector./norm(vector);

    end

    pip.vertexNormalUnitVectors = vertexNormalUnitVectors;

    return pip
end