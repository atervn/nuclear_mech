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

    # init edges vector
    shellStruct.edges = Vector{Vector{Int64}}(undef,0);
    
    # init neighbors vector
    shellStruct.neighbors = fill(Int[], length(shellStruct.vert));

    # create matrix of triangles
    triMatrix = [getindex.(shellStruct.tri,1) getindex.(shellStruct.tri,2) getindex.(shellStruct.tri,3)];

    # iterate through vertices
    for i = 1:length(shellStruct.vert)
        
        # find triangles that contain vertex i
        hasVertex = findall(triMatrix .== i);
        
        # init neighbors vector
        neighbors = Array{Int64}(undef, 0);
        
        # iterate through triangles that contain vertex i
        for j = eachindex(hasVertex)

            # add other 2 vertices in triangle to neighbors vector
            append!(neighbors,shellStruct.tri[hasVertex[j][1]][shellStruct.tri[hasVertex[j][1]] .!= i ]);

        end
        
        # get unique neighbors and save them
        unique!(neighbors)
        shellStruct.neighbors[i] = neighbors;

        # add connections to edges vector
        for j = eachindex(neighbors)
            push!(shellStruct.edges, [i, neighbors[j]]);
        end  
    end

    # init mirrorEdges vector
    shellStruct.mirrorEdges = zeros(Int64,length(shellStruct.edges));
    
    # init firstEdges vector
    shellStruct.firstEdges = zeros(Int64,length(shellStruct.edges));

    # iterate through edges
    for i = eachindex(shellStruct.edges)
            
        # find index of edge that is reverse of current edge
        mirrorIdx = findall(getindex.(shellStruct.edges,1) .== shellStruct.edges[i][2] .&& getindex.(shellStruct.edges,2) .== shellStruct.edges[i][1]);
        
        # save index of reverse edge
        shellStruct.mirrorEdges[i] = mirrorIdx[1];
        shellStruct.mirrorEdges[mirrorIdx[1]] = i;

        # if this is first time we have seen either edge, mark it as such
        if shellStruct.firstEdges[i] == 0 && shellStruct.firstEdges[mirrorIdx[1]] == 0 
            shellStruct.firstEdges[i] = 1;
        end
    end
    
    # init vertexEdges vector
    shellStruct.vertexEdges = fill(Int[], length(shellStruct.vert));

    # iterate through vertices
    for i = 1:length(shellStruct.vert)

        # get indices of edges that are connected to vertex i
        shellStruct.vertexEdges[i] = findall(getindex(shellStruct.edges,1) .== i);

    end


    bufferedDoubleTriMatrix = [triMatrix[:,3] triMatrix triMatrix[:,1]]

    # init edgesTri vector
    shellStruct.edgesTri = Vector{Vector{Int64}}(undef,length(shellStruct.edges));

    # iterate through edges
    for i = eachindex(shellStruct.edges)

        # # find indices of triangles that contain the edge
        # triangleIdx = findall(any(triMatrix .== shellStruct.edges[i][1], dims = 2) .&& any(triMatrix .== shellStruct.edges[i][2], dims = 2))
    
        # # add the triangles to the edge's list of triangles
        # shellStruct.edgesTri[i] = [triangleIdx[1] for triangleIdx in triangleIdx]

        shellStruct.edgesTri[i] = [0,0]

        firstEdges = findall(triMatrix .== shellStruct.edges[i][1])

        for j = 1:length(firstEdges)

            if bufferedDoubleTriMatrix[firstEdges[j][1],firstEdges[j][2]+2] == shellStruct.edges[i][2]
                shellStruct.edgesTri[i][1] = firstEdges[j][1]
            elseif bufferedDoubleTriMatrix[firstEdges[j][1],firstEdges[j][2]] == shellStruct.edges[i][2]
                shellStruct.edgesTri[i][2] = firstEdges[j][1]
            end
        end
    end

    return shellStruct
        
end

function add_middle_vertices!(shellStruct,i,radius)

    # get the two endpoints of the edge
    p1 = shellStruct.edges[i][1];
    p2 = shellStruct.edges[i][2];
    
    # calculate the midpoint of the edge
    newVertex = (shellStruct.vert[p1] + shellStruct.vert[p2])./2;

    # scale the midpoint by radius and divide by its norm and save it to shell vertices
    if radius != 0
        push!(shellStruct.vert,newVertex.*radius/norm(newVertex));
    else
        push!(shellStruct.vert,newVertex);
    end

    return shellStruct

end

function subdivide_triangles(shellStruct,radius)

    # create a new array to store the indices of the new vertices
    newVertexIdx = zeros(Int64,length(shellStruct.edges));

    # iterate over all edges
    for i = eachindex(shellStruct.edges)

        # if the edge is a first edge, add a new vertex in the middle of the edge
        if shellStruct.firstEdges[i] == 1

            # add a new vertex to the shellstruct
            shellStruct = add_middle_vertices!(shellStruct,i,radius);

            # get the index of the new vertex
            newVertexIdx[i] = length(shellStruct.vert);

            # get the index of the mirror vertex
            newVertexIdx[shellStruct.mirrorEdges[i]] = length(shellStruct.vert);

        end
    end
    
    # create a new array to store the new triangles
    newTriangles = Vector{Vector{Int64}}(undef,0);
    
    # iterate over all triangles
    for i = eachindex(shellStruct.tri)

        # get the neighbors of the current triangle
        nextNeighbors = circshift(shellStruct.tri[i],-1);
        prevNeighbors = circshift(shellStruct.tri[i],1);

        # create a new triangle with the current vertex and its two neighbors
        middleTriangle = zeros(Int64,3);

        # iterate through the triangle vertices
        for j = 1:3

            # get the index of the current vertex
            firstVertex = shellStruct.tri[i][j];

            # get the index of the new vertex that is connected to the current vertex by the first edge
            secondVertex = newVertexIdx[findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][j] .&& getindex.(shellStruct.edges,2) .== nextNeighbors[j])[1]];
            
            # get the index of the new vertex that is connected to the current vertex by the second edge
            thirdVertex = newVertexIdx[findall(getindex.(shellStruct.edges,1) .== shellStruct.tri[i][j] .&& getindex.(shellStruct.edges,2) .== prevNeighbors[j])[1]];
            
            # add the second vertex to the new triangle
            middleTriangle[j] = secondVertex;
            
            # add the triangle to newTriangles
            push!(newTriangles, [firstVertex, secondVertex, thirdVertex]);
            
        end
        
        # add the middle triangle to newTriangles
        push!(newTriangles, middleTriangle)
    end

    # save the new triangles
    shellStruct.tri = newTriangles;

    # get the vertex triangles
    shellStruct = get_vertex_triangles(shellStruct)

    return shellStruct

end

function subdivide_mesh!(shellStruct,radius,nSubdivisions)
    
    # get triangle edges
    shellStruct = get_edges(shellStruct);

    # for each subdivision
    for i = 1:nSubdivisions

        # subdivide triangles
        shellStruct = subdivide_triangles(shellStruct,radius);

        # get new triangle edges
        shellStruct = get_edges(shellStruct);
    end

    return shellStruct

end

function get_vertex_triangles(shellStruct)
    
    # create a new array to store the triangles for each vertex
    shellStruct.vertexTri = Vector{Vector{Int64}}(undef, length(shellStruct.vert));

    # create a matrix of the triangle indices
    triMatrix = [getindex.(shellStruct.tri,1) getindex.(shellStruct.tri,2) getindex.(shellStruct.tri,3)];

    # iterate over all vertices
    for i = 1:length(shellStruct.vert)
        
        # find all triangles that contain the vertex
        triangleIdx = findall(triMatrix  .== i);
        
        # add the triangles to the vertex's list of triangles
        shellStruct.vertexTri[i] = [triangleIdx[1] for triangleIdx in triangleIdx];

    end

    return shellStruct

end