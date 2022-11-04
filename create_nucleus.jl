function create_icosahedron!(nuc,ipar)

    # get the scaled radius
    radius = ipar.freeNucleusRadius/ipar.scalingLength;

    # define the icosahedron vertex coordinates
    a = (1 + sqrt(5))/2;

    # add the icosahedron vertices
    push!(nuc.vert,Vec(-1, a, 0))
    push!(nuc.vert,Vec( 1, a, 0))
    push!(nuc.vert,Vec(-1,-a, 0))
    push!(nuc.vert,Vec( 1,-a, 0))
    push!(nuc.vert,Vec( 0,-1, a))
    push!(nuc.vert,Vec( 0, 1, a))
    push!(nuc.vert,Vec( 0,-1,-a))
    push!(nuc.vert,Vec( 0, 1,-a))
    push!(nuc.vert,Vec( a, 0,-1))
    push!(nuc.vert,Vec( a, 0, 1))
    push!(nuc.vert,Vec(-a, 0,-1))
    push!(nuc.vert,Vec(-a, 0, 1))

    # normalize their distance from origo to the radius
    for i = eachindex(nuc.vert)
        nuc.vert[i] = nuc.vert[i]./norm(nuc.vert[i]).*radius
    end

    # create a vector for the triangles and add the icosahedron triangles (defined counterclockwise)
    nuc.tri = Vector{Vector{Int64}}(undef,0)
    push!(nuc.tri, [1, 12, 6])
    push!(nuc.tri, [1, 6, 2])
    push!(nuc.tri, [1, 2, 8])
    push!(nuc.tri, [1, 8, 11])
    push!(nuc.tri, [1, 11, 12])
    push!(nuc.tri, [2, 6, 10])
    push!(nuc.tri, [6, 12, 5])
    push!(nuc.tri, [12, 11, 3])
    push!(nuc.tri, [11, 8, 7])
    push!(nuc.tri, [8, 2, 9])
    push!(nuc.tri, [4, 10, 5])
    push!(nuc.tri, [4, 5, 3])
    push!(nuc.tri, [4, 3, 7])
    push!(nuc.tri, [4, 7, 9])
    push!(nuc.tri, [4, 9, 10])
    push!(nuc.tri, [5, 10, 6])
    push!(nuc.tri, [3, 5, 12])
    push!(nuc.tri, [7, 3, 11])
    push!(nuc.tri, [9, 7, 8])
    push!(nuc.tri, [10, 9, 2])

    return nuc
end

function get_edges(nuc)

    # a vector for the pairs
    nuc.edges = Vector{Vector{Int64}}(undef,0);
    
    # a vector for the connected neighboring vertices
    nuc.neighbors = fill(Int[], length(nuc.vert));

    # matrix form of the trignle vector
    triMatrix = [getindex.(nuc.tri,1) getindex.(nuc.tri,2) getindex.(nuc.tri,3)];

    # go through the vertices
    @inbounds for i = 1:length(nuc.vert)
        
        # find in which triangles vertex i is included
        hasVertex = triMatrix .== i;
    
        # initialize a vector for the neighbors
        neighbors = [];
        
        # go through the triangles
        @inbounds for j = eachindex(nuc.tri)
            
            # check if vertex i is included in the triangle
            if any(hasVertex[j,:])
                
                # add the other vertices in the triangle to the list of neighboring vertices
                append!(neighbors,nuc.tri[j][.!hasVertex[j,:]]);
            end
        end
        
        # get the unique neighbors and save them
        unique!(neighbors);
        nuc.neighbors[i] = neighbors;

        # add the connections to the edges matrix
        for j = eachindex(neighbors)
            push!(nuc.edges, [i, neighbors[j]]);
        end  
    end
        
    # vector to store the first time the a pair vector is seen
    nuc.mirrorEdges = zeros(Int64,length(nuc.edges));
    
    # vector to store the corresponding edge pairs
    nuc.firstEdges = zeros(Int64,length(nuc.edges));
    
    @inbounds for i = eachindex(nuc.edges)
            
        mirrorIdx = findall(getindex.(nuc.edges,1) .== nuc.edges[i][2] .&& getindex.(nuc.edges,2) .== nuc.edges[i][1]);
        
        nuc.mirrorEdges[i] = mirrorIdx[1];
        nuc.mirrorEdges[mirrorIdx[1]] = i;
        if nuc.firstEdges[i] == 0 && nuc.firstEdges[mirrorIdx[1]] == 0 
            nuc.firstEdges[i] = 1;
        end
    end

    nuc.vertexEdges = fill(Int[], length(nuc.vert));
    
    for i = 1:length(nuc.vert)
        nuc.vertexEdges[i] = findall(getindex(nuc.edges,1) .== i);
    end

    nuc.edgesTri = Vector{Vector{Int64}}(undef,length(nuc.edges));

    for i = eachindex(nuc.edges)
        neighboringTriangles = findall(sum(triMatrix .== nuc.edges[i][1],dims=2) .> 0 .&& sum(triMatrix .== nuc.edges[i][2],dims=2) .> 0);
        nuc.edgesTri[i] = [j[1] for j in neighboringTriangles];
    end
    return nuc
        
end

function add_middle_vertices!(nuc,i,radius)
    p1 = nuc.edges[i][1];
    p2 = nuc.edges[i][2];
    
    newVertex = (nuc.vert[p1] + nuc.vert[p2])./2;

    push!(nuc.vert,newVertex.*radius/norm(newVertex));
    return nuc
end

function subdivide_triangles(nuc,radius)

    newVertexIdx = zeros(Int64,length(nuc.edges));

    @inbounds for i = eachindex(nuc.edges)
        if nuc.firstEdges[i] == 1
            nuc = add_middle_vertices!(nuc,i,radius);
            newVertexIdx[i] = length(nuc.vert);
            newVertexIdx[nuc.mirrorEdges[i]] = length(nuc.vert);
        end
    end
    
    newTriangles = Vector{Vector{Int64}}(undef,0);
    
    @inbounds for i = eachindex(nuc.tri)

        nextNeighbors = circshift(nuc.tri[i],-1);
        prevNeighbors = circshift(nuc.tri[i],1);

        middleTriangle = zeros(Int64,3);
        @inbounds for j = 1:3
            firstVertex = nuc.tri[i][j];
            secondVertex = newVertexIdx[findall(getindex.(nuc.edges,1) .== nuc.tri[i][j] .&& getindex.(nuc.edges,2) .== nextNeighbors[j])[1]];
            thirdVertex = newVertexIdx[findall(getindex.(nuc.edges,1) .== nuc.tri[i][j] .&& getindex.(nuc.edges,2) .== prevNeighbors[j])[1]];
            middleTriangle[j] = secondVertex;
            
            push!(newTriangles, [firstVertex, secondVertex, thirdVertex]);
            
        end
        
        push!(newTriangles, middleTriangle)
    end

    nuc.tri = newTriangles;
    nuc = get_vertex_triangles(nuc) 
    return nuc
end

function subdivide_mesh!(nuc,ipar)
    
    # get the scaled radius
    radius = ipar.freeNucleusRadius/ipar.scalingLength;

    # get triangle edges
    nuc = get_edges(nuc);

    for i = 1:ipar.nSubdivisions
        nuc = subdivide_triangles(nuc,radius);
        nuc = get_edges(nuc);
    end

    return nuc
end

function setup_nucleus_data(nuc)

    nuc.edgeVectors = Vector{Vec{3,Float64}}(undef, length(nuc.edges));
    nuc.edgeUnitVectors = Vector{Vec{3,Float64}}(undef, length(nuc.edges));
    nuc.edgeVectorNorms = Vector{Float64}(undef, length(nuc.edges));
    get_edge_vectors!(nuc);

    nuc.edges3Vertex = Vector{Vector{Int64}}(undef,length(nuc.edges));
    for i = eachindex(nuc.edges)
        thirdVertex1 = nuc.tri[nuc.edgesTri[i][1]][.!(nuc.tri[nuc.edgesTri[i][1]] .== nuc.edges[i][1] .|| nuc.tri[nuc.edgesTri[i][1]] .== nuc.edges[i][2])][1];
        thirdVertex2 = nuc.tri[nuc.edgesTri[i][2]][.!(nuc.tri[nuc.edgesTri[i][2]] .== nuc.edges[i][1] .|| nuc.tri[nuc.edgesTri[i][2]] .== nuc.edges[i][2])][1];

        nuc.edges3Vertex[i] = [thirdVertex1, thirdVertex2]
    end

    nuc.triEdge1 = Vector{Int64}(undef,length(nuc.tri))
    nuc.triEdge2 = Vector{Int64}(undef,length(nuc.tri))
    nuc.edgeThirdVertices = Vector{Vector{Int64}}(undef,length(nuc.edges))
    

    for i = eachindex(nuc.tri)

        nuc.triEdge1[i] = findall(getindex.(nuc.edges,1) .== nuc.tri[i][1] .&& getindex.(nuc.edges,2) .== nuc.tri[i][2])[1]

        nuc.triEdge2[i] = findall(getindex.(nuc.edges,1) .== nuc.tri[i][1] .&& getindex.(nuc.edges,2) .== nuc.tri[i][3])[1]

    end

    get_triangle_normals!(nuc);

    for i = eachindex(nuc.edges)

        firstNeighbor = findall(getindex.(nuc.edges,1) .== nuc.edges[i][1] .&& getindex.(nuc.edges,2) .== nuc.edges3Vertex[i][1])[1]
        secondNeighbor = findall(getindex.(nuc.edges,1) .== nuc.edges[i][1] .&& getindex.(nuc.edges,2) .== nuc.edges3Vertex[i][2])[1]

        nuc.edgeThirdVertices[i] = [firstNeighbor, secondNeighbor]

    end

    nuc.normalVolume = get_volume!(nuc);
    nuc.normalTriangleAreas = get_area!(nuc);
    nuc.normalArea = sum(nuc.normalTriangleAreas);
    nuc.normalAngle = mean(get_triangle_angles(nuc));
    lengths = zeros(Float64,length(nuc.edges));

    for i = eachindex(nuc.edges)  
        lengths[i] = norm(nuc.vert[nuc.edges[i][2]] - nuc.vert[nuc.edges[i][1]]);
    end
    nuc.normalLengths = lengths;

    return nuc

end

function get_vertex_triangles(nuc)
    
    nuc.vertexTri = Vector{Vector{Int64}}(undef, length(nuc.vert));

    triMatrix = [getindex.(nuc.tri,1) getindex.(nuc.tri,2) getindex.(nuc.tri,3)];

    for i = 1:length(nuc.vert)
        aa = findall(triMatrix  .== i);
        nuc.vertexTri[i] = [i[1] for i in aa];
    end

    return nuc

end