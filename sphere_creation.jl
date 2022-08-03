function create_icosahedron!(nuc,radius)
    # define the icosahedron vertex coordinates
    a = (1 + sqrt(5))/2;
    nuc.x = [-1, 1, -1, 1, 0, 0, 0, 0, a, a, -a, -a];
    nuc.y = [a, a, -a, -a, -1, 1, -1, 1, 0, 0, 0, 0];
    nuc.z = [0, 0, 0, 0, a, a, -a, -a, -1, 1, -1, 1];

    currentRadius = sqrt(nuc.x[1]^2 + nuc.y[1]^2 + nuc.z[1]^2);

    nuc.x = nuc.x./currentRadius.*radius;
    nuc.y = nuc.y./currentRadius.*radius;
    nuc.z = nuc.z./currentRadius.*radius;

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

function get_edges!(nuc)

    # initialize a vector for the pairs
    nuc.edges = zeros(0,2);
    
    nuc.neighbors = fill(Int[], length(nuc.x));

    # go through the vertices
    @inbounds for i = 1:length(nuc.x)
        
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

    nuc.vertexEdges = fill(Int[], length(nuc.x));
    
    for i = 1:length(nuc.x)
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
    newX = (nuc.x[p1] + nuc.x[p2])/2;
    newY = (nuc.y[p1] + nuc.y[p2])/2;
    newZ = (nuc.z[p1] + nuc.z[p2])/2;
    
    vertexNorm = sqrt(newX^2 + newY^2 + newZ^2);

    nuc.x = cat(nuc.x, newX*radius/vertexNorm,dims=1);
    nuc.y = cat(nuc.y, newY*radius/vertexNorm,dims=1);
    nuc.z = cat(nuc.z, newZ*radius/vertexNorm,dims=1);
    return nuc
end

function subdivide_triangles!(nuc,radius)
    
    newVertexIdx = zeros(Int64,size(nuc.edges,1));
    
    @inbounds for i = 1:size(nuc.edges,1)
        if nuc.firstEdges[i] == 1
            nuc = add_middle_vertices!(nuc,i,radius);
            newVertexIdx[i] = length(nuc.x);
            newVertexIdx[nuc.mirrorEdges[i]] = length(nuc.x);
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

    nuc.vertexTri = fill(Int[], length(nuc.x), 2);

    for i = 1:length(nuc.x)
        aa = findall(nuc.tri  .== i);
        nuc.vertexTri[i,1] = [i[1] for i in aa];
        nuc.vertexTri[i,2] = [i[2]-1 for i in aa];
    end

    return nuc
        
end

function subdivide_mesh!(nuc,radius,nSubdivisions)
    
    @inbounds for i = 1:nSubdivisions
        nuc = subdivide_triangles!(nuc,radius);
        nuc = get_edges!(nuc);
    end

    return nuc
end