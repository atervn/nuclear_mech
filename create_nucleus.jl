function create_icosahedron!(nuc,ipar)

    radius = ipar.freeNucleusRadius/ipar.scalingLength;

    # define the icosahedron vertex coordinates
    a = (1 + sqrt(5))/2;

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

    for i = eachindex(nuc.vert)
        nuc.vert[i] = nuc.vert[i]./norm(nuc.vert[i]).*radius
    end

    # define the triangles so that they are defined 
    # in counterclockwise direction
    nuc.tri = [1 12 6; 1 6 2;   1 2 8;
                1 8 11;  1 11 12; 2 6 10;
                6 12 5;  12 11 3; 11 8 7;
                8 2 9;   4 10 5;  4 5 3;
                4 3 7;   4 7 9;   4 9 10;
                5 10 6;  3 5 12;  7 3 11;
                9 7 8;   10 9 2];

    return nuc
end

function get_edges(nuc)

    # initialize a vector for the pairs
    nuc.edges = zeros(0,2);
    
    nuc.neighbors = fill(Int[], length(nuc.vert));

    # go through the vertices
    @inbounds for i = 1:length(nuc.vert)
        
        # find in which triangles vertex i is included
        hasVertex = nuc.tri .== i;
    
        # initialize a vector for the neighbors
        neighbors = [];
        
        # go through the triangles
        @inbounds for j = 1:length(nuc.tri[:,1])
            
            # check if vertex i is included in the triangle
            if any(hasVertex[j,:])
                
                # add the other vertices in the triangle to the list of neighboring vertices
                append!(neighbors,nuc.tri[j,.!hasVertex[j,:]]);
            end
        end
        
        # get the unique neighbors
        unique!(neighbors);
        
        nuc.neighbors[i] = neighbors;

        # add the connections to the edges matrix
        nuc.edges = [nuc.edges ; [i.*Int.(ones(length(neighbors),1)) neighbors]]  ;  
    end
        
    # vector to store the first time the a pair vector is seen
    nuc.mirrorEdges = zeros(Int64,size(nuc.edges,1));
    
    # vector to store the corresponding edge pairs
    nuc.firstEdges = zeros(Int64,size(nuc.edges,1));
    
    @inbounds for i = 1:size(nuc.edges,1)
            
        b = findall(nuc.edges[:,1] .== nuc.edges[i,2] .&& nuc.edges[:,2] .== nuc.edges[i,1]);
        
        nuc.mirrorEdges[i] = b[1];
        nuc.mirrorEdges[b[1]] = i;
        if nuc.firstEdges[i] == 0 && nuc.firstEdges[b[1]] == 0 
            nuc.firstEdges[i] = 1;
        end
    end

    nuc.vertexEdges = fill(Int[], length(nuc.vert));
    
    for i = 1:length(nuc.vert)
        nuc.vertexEdges[i] = findall(nuc.edges[:,1] .== i);
    end

    nuc.edgesTri = zeros(Int64,size(nuc.edges,1),2);

    for i = 1:size(nuc.edges,1)
        neighboringTriangles = findall(sum(nuc.tri .== nuc.edges[i,1],dims=2) .> 0 .&& sum(nuc.tri .== nuc.edges[i,2],dims=2) .> 0);
        nuc.edgesTri[i,:] = [j[1] for j in neighboringTriangles];
    end
    return nuc
        
end

function add_middle_vertices!(nuc,i,radius)
    p1 = nuc.edges[i,1];
    p2 = nuc.edges[i,2];
    
    newVertex = (nuc.vert[p1] + nuc.vert[p2])./2;

    push!(nuc.vert,newVertex.*radius/norm(newVertex));
    return nuc
end

function subdivide_triangles(nuc,radius)

    newVertexIdx = zeros(Int64,size(nuc.edges,1));

    @inbounds for i = 1:size(nuc.edges,1)
        if nuc.firstEdges[i] == 1
            nuc = add_middle_vertices!(nuc,i,radius);
            newVertexIdx[i] = length(nuc.vert);
            newVertexIdx[nuc.mirrorEdges[i]] = length(nuc.vert);
        end
    end
    
    nextNeighbors = circshift(nuc.tri,(0,-1));
    prevNeighbors = circshift(nuc.tri,(0,1));
    
    newTriangles = zeros(Int64,0,3);
    
    @inbounds for i = 1:length(nuc.tri[:,1])
        
        middleTriangle = zeros(Int64,1,3);
        @inbounds for j = 1:3
            firstVertex = nuc.tri[i,j];
            secondVertex = newVertexIdx[findall(nuc.edges[:,1] .== nuc.tri[i,j] .&& nuc.edges[:,2] .== nextNeighbors[i,j])[1]];
            thirdVertex = newVertexIdx[findall(nuc.edges[:,1] .== nuc.tri[i,j] .&& nuc.edges[:,2] .== prevNeighbors[i,j])[1]];
            middleTriangle[j] = secondVertex;
            
            newTriangles = [newTriangles ; [firstVertex secondVertex thirdVertex]];
            
        end
            newTriangles = [newTriangles; middleTriangle];
    end

    nuc.tri = newTriangles;
    nuc = get_vertex_triangles(nuc) 
    return nuc
end

function subdivide_mesh!(nuc,ipar)
    
    radius = ipar.freeNucleusRadius/ipar.scalingLength;

    nuc = get_edges!(nuc);

    for i = 1:ipar.nSubdivisions
        nuc = subdivide_triangles(nuc,radius);
        nuc = get_edges(nuc);
    end

    return nuc
end

function setup_nucleus_data(nuc)

    get_edge_vectors!(nuc);

    nuc.neighboringTriangles = zeros(Int64,size(nuc.edges,1),2)
    nuc.edges3vertex = zeros(Int64,size(nuc.edges,1),2);
    for i = 1:size(nuc.edges,1)
        temp = findall(sum(nuc.tri .== nuc.edges[i,1],dims=2) .> 0 .&& sum(nuc.tri .== nuc.edges[i,2],dims=2) .> 0);
        nuc.neighboringTriangles[i,:] = [j[1] for j in temp];
        thirdVertex1 = nuc.tri[nuc.neighboringTriangles[i,1],.!(nuc.tri[nuc.neighboringTriangles[i,1],:] .== nuc.edges[i,1] .|| nuc.tri[nuc.neighboringTriangles[i,1],:] .== nuc.edges[i,2])][1];
        thirdVertex2 = nuc.tri[nuc.neighboringTriangles[i,2],.!(nuc.tri[nuc.neighboringTriangles[i,2],:] .== nuc.edges[i,1] .|| nuc.tri[nuc.neighboringTriangles[i,2],:] .== nuc.edges[i,2])][1];

        nuc.edges3vertex[i,1] = thirdVertex1;
        nuc.edges3vertex[i,2] = thirdVertex2; 
    end

    nuc.testview = Array{Any}(undef,length(nuc.edges))
    nuc.p1 = Array{Any}(undef,length(nuc.tri))
    nuc.p2 = Array{Any}(undef,length(nuc.tri))
    nuc.p3 = Array{Any}(undef,length(nuc.tri))
    nuc.trii = Array{Any}(undef,length(nuc.tri))
    nuc.tri21 = Array{Any}(undef,length(nuc.tri))
    nuc.tri31 = Array{Any}(undef,length(nuc.tri))
    nuc.tri12 = Array{Any}(undef,length(nuc.tri))
    nuc.tri13 = Array{Any}(undef,length(nuc.tri))
    nuc.ep1 = Array{Any}(undef,length(nuc.edges))
    nuc.ep2 = Array{Any}(undef,length(nuc.edges))
    nuc.ep31 = Array{Any}(undef,length(nuc.edges))
    nuc.ep32 = Array{Any}(undef,length(nuc.edges))

    for i = 1:size(nuc.tri,1)
        nuc.p1[i] = @view nuc.vert[nuc.tri[i,1]];
        nuc.p2[i] = @view nuc.vert[nuc.tri[i,2]];
        nuc.p3[i] = @view nuc.vert[nuc.tri[i,3]];
        nuc.trii[i] = @view nuc.vert[nuc.tri[i,:]];
        
        edgeIdx = findall(nuc.edges[:,1] .== nuc.tri[i,2] .&& nuc.edges[:,2] .== nuc.tri[i,1])
        nuc.tri21[i] = @view nuc.edgeVectors[edgeIdx];

        edgeIdx = findall(nuc.edges[:,1] .== nuc.tri[i,3] .&& nuc.edges[:,2] .== nuc.tri[i,1])
        nuc.tri31[i] = @view nuc.edgeVectors[edgeIdx];

        edgeIdx = findall(nuc.edges[:,1] .== nuc.tri[i,1] .&& nuc.edges[:,2] .== nuc.tri[i,2])
        nuc.tri12[i] = @view nuc.edgeVectors[edgeIdx];

        edgeIdx = findall(nuc.edges[:,1] .== nuc.tri[i,1] .&& nuc.edges[:,2] .== nuc.tri[i,3])
        nuc.tri13[i] = @view nuc.edgeVectors[edgeIdx];

    end

    get_triangle_normals!(nuc);

    for i = 1:size(nuc.edges,1)
        nuc.testview[i] = @view nuc.triangleNormalUnitVectors[nuc.edgesTri[i,:]];
        nuc.ep1[i] = @view nuc.vert[nuc.edges[i,1]];
        nuc.ep2[i] = @view nuc.vert[nuc.edges[i,2]];
        nuc.ep31[i] = @view nuc.vert[nuc.edges3vertex[i,1]];
        nuc.ep32[i] = @view nuc.vert[nuc.edges3vertex[i,2]];
    end

    nuc.normalVolume = get_volume!(nuc);
    nuc.normalTriangleAreas = get_area!(nuc);
    nuc.normalArea = sum(nuc.normalTriangleAreas);
    nuc.normalAngle = mean(get_triangle_angles(nuc));
    lengths = zeros(Float64,Int64(size(nuc.edges,1)));

    for i = 1:size(nuc.edges,1)  
        lengths[i] = norm(nuc.vert[nuc.edges[i,2]] - nuc.vert[nuc.edges[i,1]]);
    end
    nuc.normalLengths = lengths;

    return nuc

end

function get_vertex_triangles(nuc)
    
    nuc.vertexTri = fill(Int[], length(nuc.vert), 2);

    for i = 1:length(nuc.vert)
        aa = findall(nuc.tri  .== i);
        nuc.vertexTri[i,1] = [i[1] for i in aa];
        nuc.vertexTri[i,2] = [i[2]-1 for i in aa];
    end

    return nuc

end