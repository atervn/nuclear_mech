function get_icosaherdon!(shellStruct,radius)

    # define the icosahedron vertex coordinates
    a = (1 + sqrt(5))/2;

    # add the icosahedron vertices
    push!(shellStruct.vert,Vec(-1, a, 0))
    push!(shellStruct.vert,Vec( 1, a, 0))
    push!(shellStruct.vert,Vec(-1,-a, 0))
    push!(shellStruct.vert,Vec( 1,-a, 0))
    push!(shellStruct.vert,Vec( 0,-1, a))
    push!(shellStruct.vert,Vec( 0, 1, a))
    push!(shellStruct.vert,Vec( 0,-1,-a))
    push!(shellStruct.vert,Vec( 0, 1,-a))
    push!(shellStruct.vert,Vec( a, 0,-1))
    push!(shellStruct.vert,Vec( a, 0, 1))
    push!(shellStruct.vert,Vec(-a, 0,-1))
    push!(shellStruct.vert,Vec(-a, 0, 1))

    # normalize their distance from origo to the radius
    for i = eachindex(shellStruct.vert)
        shellStruct.vert[i] = shellStruct.vert[i]./norm(shellStruct.vert[i]).*radius
    end

    # create a vector for the triangles and add the icosahedron triangles (defined counterclockwise)
    shellStruct.tri = Vector{Vector{Int64}}(undef,0)
    push!(shellStruct.tri, [1, 12, 6])
    push!(shellStruct.tri, [1, 6, 2])
    push!(shellStruct.tri, [1, 2, 8])
    push!(shellStruct.tri, [1, 8, 11])
    push!(shellStruct.tri, [1, 11, 12])
    push!(shellStruct.tri, [2, 6, 10])
    push!(shellStruct.tri, [6, 12, 5])
    push!(shellStruct.tri, [12, 11, 3])
    push!(shellStruct.tri, [11, 8, 7])
    push!(shellStruct.tri, [8, 2, 9])
    push!(shellStruct.tri, [4, 10, 5])
    push!(shellStruct.tri, [4, 5, 3])
    push!(shellStruct.tri, [4, 3, 7])
    push!(shellStruct.tri, [4, 7, 9])
    push!(shellStruct.tri, [4, 9, 10])
    push!(shellStruct.tri, [5, 10, 6])
    push!(shellStruct.tri, [3, 5, 12])
    push!(shellStruct.tri, [7, 3, 11])
    push!(shellStruct.tri, [9, 7, 8])
    push!(shellStruct.tri, [10, 9, 2])

    return shellStruct
end

function get_edges(shellStruct)

    # a vector for the pairs
    shellStruct.edges = Vector{Vector{Int64}}(undef,0);
    
    # a vector for the connected neighboring vertices
    shellStruct.neighbors = fill(Int[], length(shellStruct.vert));

    # matrix form of the trignle vector
    triMatrix = [getindex.(shellStruct.tri,1) getindex.(shellStruct.tri,2) getindex.(shellStruct.tri,3)];

    # go through the vertices
    for i = 1:length(shellStruct.vert)
        
        # find in which triangles vertex i is included
        hasVertex = findall(triMatrix .== i);
        # initialize a vector for the neighbors
        neighbors = Array{Int64}(undef, 0);
        
        # go through the triangles
        for j = eachindex(hasVertex)
            append!(neighbors,shellStruct.tri[hasVertex[j][1]][shellStruct.tri[hasVertex[j][1]] .!= i ]);
        end
        
        # get the unique neighbors and save them
        unique!(neighbors)
        shellStruct.neighbors[i] = neighbors;

        # add the connections to the edges matrix
        for j = eachindex(neighbors)
            push!(shellStruct.edges, [i, neighbors[j]]);
        end  
    end

 
    # vector to store the first time the a pair vector is seen
    shellStruct.mirrorEdges = zeros(Int64,length(shellStruct.edges));
    
    # vector to store the corresponding edge pairs
    shellStruct.firstEdges = zeros(Int64,length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)
            
        mirrorIdx = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][2] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges[i][1]);
        
        shellStruct.mirrorEdges[i] = mirrorIdx[1];
        shellStruct.mirrorEdges[mirrorIdx[1]] = i;
        if shellStruct.firstEdges[i] == 0 && shellStruct.firstEdges[mirrorIdx[1]] == 0 
            shellStruct.firstEdges[i] = 1;
        end
    end
    
    shellStruct.vertexEdges = fill(Int[], length(shellStruct.vert));

    for i = 1:length(shellStruct.vert)
        shellStruct.vertexEdges[i] = findall(getindex(shellStruct.edges,1) .== i);
    end

    shellStruct.edgesTri = Vector{Vector{Int64}}(undef,length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)

        # neighboringTriangles = findall(any(triMatrix .== shellStruct.edges[i][1],dims=2) .&& any(triMatrix .== shellStruct.edges[i][2],dims=2));
        shellStruct.edgesTri[i] = [j[1] for j in findall(any(triMatrix .== shellStruct.edges[i][1],dims=2) .&& any(triMatrix .== shellStruct.edges[i][2],dims=2))];
    end

    return shellStruct
        
end

function add_middle_vertices!(shellStruct,i,radius)
    p1 = shellStruct.edges[i][1];
    p2 = shellStruct.edges[i][2];
    
    newVertex = (shellStruct.vert[p1] + shellStruct.vert[p2])./2;
    if radius != 0
        push!(shellStruct.vert,newVertex.*radius/norm(newVertex));
    else
        push!(shellStruct.vert,newVertex);
    end
    return shellStruct
end

function subdivide_triangles(shellStruct,radius)

    newVertexIdx = zeros(Int64,length(shellStruct.edges));

    for i = eachindex(shellStruct.edges)
        if shellStruct.firstEdges[i] == 1
            shellStruct = add_middle_vertices!(shellStruct,i,radius);
            newVertexIdx[i] = length(shellStruct.vert);
            newVertexIdx[shellStruct.mirrorEdges[i]] = length(shellStruct.vert);
        end
    end
    
    newTriangles = Vector{Vector{Int64}}(undef,0);
    
    for i = eachindex(shellStruct.tri)

        nextNeighbors = circshift(shellStruct.tri[i],-1);
        prevNeighbors = circshift(shellStruct.tri[i],1);

        middleTriangle = zeros(Int64,3);
        for j = 1:3
            firstVertex = shellStruct.tri[i][j];
            secondVertex = newVertexIdx[findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][j] .&& getindex.(shellStruct.edges,2) .== nextNeighbors[j])[1]];
            thirdVertex = newVertexIdx[findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][j] .&& getindex.(shellStruct.edges,2) .== prevNeighbors[j])[1]];
            middleTriangle[j] = secondVertex;
            
            push!(newTriangles, [firstVertex, secondVertex, thirdVertex]);
            
        end
        
        push!(newTriangles, middleTriangle)
    end

    shellStruct.tri = newTriangles;
    shellStruct = get_vertex_triangles(shellStruct) 
    return shellStruct
end

function subdivide_mesh!(shellStruct,radius,nSubdivisions)
    
    # get triangle edges
    shellStruct = get_edges(shellStruct);

    for i = 1:nSubdivisions
        shellStruct = subdivide_triangles(shellStruct,radius);
        shellStruct = get_edges(shellStruct);
    end

    return shellStruct
end

function get_vertex_triangles(shellStruct)
    
    shellStruct.vertexTri = Vector{Vector{Int64}}(undef, length(shellStruct.vert));

    triMatrix = [getindex.(shellStruct.tri,1) getindex.(shellStruct.tri,2) getindex.(shellStruct.tri,3)];

    for i = 1:length(shellStruct.vert)
        aa = findall(triMatrix  .== i);
        shellStruct.vertexTri[i] = [i[1] for i in aa];
    end

    return shellStruct

end