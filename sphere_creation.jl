function create_icosahedron!(nucleus,radius)
    # define the icosahedron vertex coordinates
    a = (1 + sqrt(5))/2;
    nucleus.x = [-1, 1, -1, 1, 0, 0, 0, 0, a, a, -a, -a];
    nucleus.y = [a, a, -a, -a, -1, 1, -1, 1, 0, 0, 0, 0];
    nucleus.z = [0, 0, 0, 0, a, a, -a, -a, -1, 1, -1, 1];

    currentRadius = sqrt(nucleus.x[1]^2 + nucleus.y[1]^2 + nucleus.z[1]^2);

    nucleus.x = nucleus.x./currentRadius.*radius;
    nucleus.y = nucleus.y./currentRadius.*radius;
    nucleus.z = nucleus.z./currentRadius.*radius;

    # define the triangles so that they are defined 
    # in counterclockwise direction
    nucleus.tri = [1 12 6; 1 6 2;   1 2 8;
                1 8 11;  1 11 12; 2 6 10;
                6 12 5;  12 11 3; 11 8 7;
                8 2 9;   4 10 5;  4 5 3;
                4 3 7;   4 7 9;   4 9 10;
                5 10 6;  3 5 12;  7 3 11;
                9 7 8;   10 9 2];

    return nucleus
end

function get_edges!(nucleus)

    # initialize a vector for the pairs
    nucleus.edges = zeros(0,2);
    
    # go through the vertices
    @inbounds for i = 1:length(nucleus.x)
        
        # find in which triangles vertex i is included
        hasVertex = nucleus.tri .== i;
    
        # initialize a vector for the neighbors
        neighbors = [];
        
        # go through the triangles
        @inbounds for j = 1:length(nucleus.tri[:,1])
            
            # check if vertex i is included in the triangle
            if any(hasVertex[j,:])
                
                # add the other vertices in the triangle to the list of neighboring vertices
                append!(neighbors,nucleus.tri[j,.!hasVertex[j,:]]);
            end
        end
        
        # get the unique neighbors
        unique!(neighbors);
        
        # add the connections to the edges matrix
        nucleus.edges = [nucleus.edges ; [i.*Int.(ones(length(neighbors),1)) neighbors]]  ;  
    end
        
    # vector to store the first time the a pair vector is seen
    nucleus.mirrorEdges = zeros(Int64,size(nucleus.edges,1));
    
    # vector to store the corresponding edge pairs
    nucleus.firstEdges = zeros(Int64,size(nucleus.edges,1));
    
    @inbounds for i = 1:size(nucleus.edges,1)
            
        b = findall(nucleus.edges[:,1] .== nucleus.edges[i,2] .&& nucleus.edges[:,2] .== nucleus.edges[i,1]);
        
        nucleus.mirrorEdges[i] = b[1];
        nucleus.mirrorEdges[b[1]] = i;
        if nucleus.firstEdges[i] == 0 && nucleus.firstEdges[b[1]] == 0 
            nucleus.firstEdges[i] = 1;
        end
    end

    nucleus.vertexEdges = fill(Int[], length(nucleus.x));
    
    for i = 1:length(nucleus.x)
        nucleus.vertexEdges[i] = findall(nucleus.edges[:,1] .== i);
    end

    nucleus.edgesTri = zeros(Int64,size(nucleus.edges,1),2);

    for i = 1:size(nucleus.edges,1)
        neighboringTriangles = findall(sum(nucleus.tri .== nucleus.edges[i,1],dims=2) .> 0 .&& sum(nucleus.tri .== nucleus.edges[i,2],dims=2) .> 0);
        nucleus.edgesTri[i,:] = [j[1] for j in neighboringTriangles];
    end
    return nucleus
        
end

function add_middle_vertices!(nucleus,i,radius)
    p1 = nucleus.edges[i,1];
    p2 = nucleus.edges[i,2];
    newX = (nucleus.x[p1] + nucleus.x[p2])/2;
    newY = (nucleus.y[p1] + nucleus.y[p2])/2;
    newZ = (nucleus.z[p1] + nucleus.z[p2])/2;
    
    vertexNorm = sqrt(newX^2 + newY^2 + newZ^2);

    nucleus.x = cat(nucleus.x, newX*radius/vertexNorm,dims=1);
    nucleus.y = cat(nucleus.y, newY*radius/vertexNorm,dims=1);
    nucleus.z = cat(nucleus.z, newZ*radius/vertexNorm,dims=1);
    return nucleus
end

function subdivide_triangles!(nucleus,radius)
    
    newVertexIdx = zeros(Int64,size(nucleus.edges,1));
    
    @inbounds for i = 1:size(nucleus.edges,1)
        if nucleus.firstEdges[i] == 1
            nucleus = add_middle_vertices!(nucleus,i,radius);
            newVertexIdx[i] = length(nucleus.x);
            newVertexIdx[nucleus.mirrorEdges[i]] = length(nucleus.x);
        end
    end
    
    nextNeighbors = circshift(nucleus.tri,(0,-1));
    prevNeighbors = circshift(nucleus.tri,(0,1));
    
    newTriangles = zeros(Int64,0,3);
    
    @inbounds for i = 1:length(nucleus.tri[:,1])
        
        middleTriangle = zeros(Int64,1,3);
        @inbounds for j = 1:3
            firstVertex = nucleus.tri[i,j];
            secondVertex = newVertexIdx[findall(nucleus.edges[:,1] .== nucleus.tri[i,j] .&& nucleus.edges[:,2] .== nextNeighbors[i,j])[1]];
            thirdVertex = newVertexIdx[findall(nucleus.edges[:,1] .== nucleus.tri[i,j] .&& nucleus.edges[:,2] .== prevNeighbors[i,j])[1]];
            middleTriangle[j] = secondVertex;
            
            newTriangles = [newTriangles ; [firstVertex secondVertex thirdVertex]];
            
        end
            newTriangles = [newTriangles; middleTriangle];
    end

    nucleus.tri = newTriangles;

    nucleus.vertexTri = fill(Int[], length(nucleus.x), 2);

    for i = 1:length(nucleus.x)
        aa = findall(newTriangles .== i);
        nucleus.vertexTri[i,1] = [i[1] for i in aa];
        nucleus.vertexTri[i,2] = [i[2]-1 for i in aa];
    end

    return nucleus
        
end

function subdivide_mesh!(nucleus,radius,nSubdivisions)
    
    @inbounds for i = 1:nSubdivisions
        nucleus = subdivide_triangles!(nucleus,radius);
        nucleus = get_edges!(nucleus);
    end

    return nucleus
end