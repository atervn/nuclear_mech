using Meshes, NearestNeighbors, LinearAlgebra, IterativeSolvers, Plots, Statistics

nSubdiv = 0

k = 0.41

maxT = 200000

dist = 2.0
if nSubdiv == 0
    vert = Vector{Vec{2,Float64}}(undef,3)
    tri = Vector{Vector{Int64}}(undef,1)
elseif nSubdiv == 1
    vert = Vector{Vec{2,Float64}}(undef,6)
    tri = Vector{Vector{Int64}}(undef,4)
elseif nSubdiv == 2
    vert = Vector{Vec{2,Float64}}(undef,15)
    tri = Vector{Vector{Int64}}(undef,16)
elseif nSubdiv == 3
    vert = Vector{Vec{2,Float64}}(undef,45)
    tri = Vector{Vector{Int64}}(undef,64)
end


vert[1] = Vec(0.,0.)
vert[2] = Vec(dist,0.)
vert[3] = Vec(dist/2,dist/2*tand(60))

if nSubdiv == 1

    vert[4] = (vert[1] + vert[2])/2
    vert[5] = (vert[1] + vert[3])/2
    vert[6] = (vert[2] + vert[3])/2

    tri[1] = [1,4,5]
    tri[2] = [4,2,6]
    tri[3] = [5,6,3]
    tri[4] = [4,6,5]
    
    bottomRow = [1, 4, 2]

elseif nSubdiv == 2

    vert[4] = (vert[1] + vert[2])/2
    vert[5] = (vert[1] + vert[3])/2
    vert[6] = (vert[2] + vert[3])/2
    vert[7] = (vert[1] + vert[4])/2
    vert[8] = (vert[1] + vert[5])/2
    vert[9] = (vert[4] + vert[5])/2
    vert[10] = (vert[4] + vert[2])/2
    vert[11] = (vert[4] + vert[6])/2
    vert[12] = (vert[6] + vert[2])/2
    vert[13] = (vert[5] + vert[6])/2
    vert[14] = (vert[5] + vert[3])/2
    vert[15] = (vert[6] + vert[3])/2
    
    tri[1] = [1,7,8]
    tri[2] = [7,8,9]
    tri[3] = [7,4,9]
    tri[4] = [4,9,11]
    tri[5] = [4,10,11]
    tri[6] = [10,11,12]
    tri[7] = [10,12,2]
    tri[8] = [8,9,5]
    tri[9] = [9,5,13]
    tri[10] = [9,11,13]
    tri[11] = [11,13,6]
    tri[12] = [11,12,6]
    tri[13] = [5,13,14]
    tri[14] = [13,14,15]
    tri[15] = [13,15,6]
    tri[16] = [14,15,3]

    bottomRow = [1, 4, 2, 7, 10]

elseif nSubdiv == 3

    vert[4] = (vert[1] + vert[2])/2
    vert[5] = (vert[1] + vert[3])/2
    vert[6] = (vert[2] + vert[3])/2
    vert[7] = (vert[1] + vert[4])/2
    vert[8] = (vert[1] + vert[5])/2
    vert[9] = (vert[4] + vert[5])/2
    vert[10] = (vert[4] + vert[2])/2
    vert[11] = (vert[4] + vert[6])/2
    vert[12] = (vert[6] + vert[2])/2
    vert[13] = (vert[5] + vert[6])/2
    vert[14] = (vert[5] + vert[3])/2
    vert[15] = (vert[6] + vert[3])/2
    vert[16] = (vert[1] + vert[7])/2
    vert[17] = (vert[1] + vert[8])/2
    vert[18] = (vert[7] + vert[8])/2
    vert[19] = (vert[7] + vert[4])/2
    vert[20] = (vert[7] + vert[9])/2
    vert[21] = (vert[4] + vert[9])/2
    vert[22] = (vert[4] + vert[10])/2
    vert[23] = (vert[4] + vert[11])/2
    vert[24] = (vert[10] + vert[11])/2
    vert[25] = (vert[10] + vert[2])/2
    vert[26] = (vert[10] + vert[12])/2
    vert[27] = (vert[12] + vert[2])/2
    vert[28] = (vert[8] + vert[9])/2
    vert[29] = (vert[9] + vert[11])/2
    vert[30] = (vert[11] + vert[12])/2
    vert[31] = (vert[5] + vert[8])/2
    vert[32] = (vert[5] + vert[9])/2
    vert[33] = (vert[9] + vert[13])/2
    vert[34] = (vert[11] + vert[13])/2
    vert[35] = (vert[6] + vert[11])/2
    vert[36] = (vert[6] + vert[12])/2
    vert[37] = (vert[5] + vert[13])/2
    vert[38] = (vert[6] + vert[13])/2
    vert[39] = (vert[5] + vert[14])/2
    vert[40] = (vert[13] + vert[14])/2
    vert[41] = (vert[13] + vert[15])/2
    vert[42] = (vert[6] + vert[15])/2
    vert[43] = (vert[14] + vert[15])/2
    vert[44] = (vert[3] + vert[14])/2
    vert[45] = (vert[15] + vert[3])/2
    
    tri[1] = [1,16,17]
    tri[2] = [16,18,17]
    tri[3] = [16,7,18]
    tri[4] = [7,20,18]
    tri[5] = [7,19,20]
    tri[6] = [19,21,20]
    tri[7] = [19,4,21]
    tri[8] = [4,23,21]
    tri[9] = [4,22,23]
    tri[10] = [22,24,23]
    tri[11] = [22,10,24]
    tri[12] = [10,26,24]
    tri[13] = [10,25,26]
    tri[14] = [25,27,26]
    tri[15] = [25,2,27]

    tri[16] = [17,18,8]
    tri[17] = [18,28,8]
    tri[18] = [18,20,28]
    tri[19] = [20,9,28]
    tri[20] = [20,21,9]
    tri[21] = [21,29,9]
    tri[22] = [21,23,29]
    tri[23] = [23,11,29]
    tri[24] = [23,24,11]
    tri[25] = [24,30,11]
    tri[26] = [24,26,30]
    tri[27] = [26,12,30]
    tri[28] = [26,27,12]

    tri[29] = [8,28,31]
    tri[30] = [28,32,31]
    tri[31] = [28,9,32]
    tri[32] = [9,33,32]
    tri[33] = [9,29,33]
    tri[34] = [29,34,33]
    tri[35] = [29,11,34]
    tri[36] = [11,35,34]
    tri[37] = [11,30,35]
    tri[38] = [30,36,35]
    tri[39] = [30,12,36]

    tri[40] = [31,32,5]
    tri[41] = [32,37,5]
    tri[42] = [32,33,37]
    tri[43] = [33,13,37]
    tri[44] = [33,34,13]
    tri[45] = [34,38,13]
    tri[46] = [34,35,38]
    tri[47] = [35,6,38]
    tri[48] = [35,36,6]

    tri[49] = [5,37,39]
    tri[50] = [37,40,39]
    tri[51] = [37,13,40]
    tri[52] = [13,41,40]
    tri[53] = [13,38,41]
    tri[54] = [38,42,41]
    tri[55] = [38,6,42]
    
    tri[56] = [39,40,14]
    tri[57] = [40,43,14]
    tri[58] = [40,41,43]
    tri[59] = [41,15,43]
    tri[60] = [41,42,15]
    
    tri[61] = [14,43,44]
    tri[62] = [43,45,44]
    tri[63] = [43,15,45]

    tri[64] = [44,45,3]

    bottomRow = [1, 4, 2, 7, 10, 16, 19, 22, 25]

else
    tri[1] = [1,2,3]
    bottomRow = [1, 2]
end

normalAreas = zeros(length(tri));

for i = eachindex(tri)
    
    vect1 = vert[tri[i][2]] - vert[tri[i][1]];
    vect2 = vert[tri[i][3]] - vert[tri[i][1]];
    dotProduct = dot(vect1,vect2)
    θ = acos(dotProduct/(norm(vect1)*norm(vect2)));
    normalAreas[i] = 1/2*norm(vect1)*norm(vect2)*sin(θ);
end

frictionMatrix = zeros(length(vert), length(vert))

for i = eachindex(vert)
    frictionMatrix[i,i] = 1.0;
end

dt = 0.0002;

for t = 1:maxT

    forces = Vector{Vec{2,Float64}}(undef,length(vert))

    areas = zeros(length(tri));

    for i = eachindex(tri)
        
        vect1 = vert[tri[i][2]] - vert[tri[i][1]];
        vect2 = vert[tri[i][3]] - vert[tri[i][1]];
        dotProduct = dot(vect1,vect2)
        θ = acos(dotProduct/(norm(vect1)*norm(vect2)));
        areas[i] = 1/2*norm(vect1)*norm(vect2)*sin(θ);
    end

    for i = eachindex(vert)
        forces[i] = Vec(0.,0.)
    end

    for i = eachindex(tri)

        baryocenter = mean(vert[tri[i]]);
        magnitude = k*(areas[i] - normalAreas[i]);
        
        for j = 1:3

            vector = vert[tri[i][j]] - baryocenter;
            unitVector = vector/norm(vector)

            forces[tri[i][j]] += -magnitude*unitVector
        end

    end

    forces[3] += Vec(0.,1.)

    solX = cg(frictionMatrix,getindex.(forces,1))
    solY = cg(frictionMatrix,getindex.(forces,2))

    movements = Vec.(solX,solY).*dt

    for i = eachindex(vert)
        if !(any(i .== bottomRow)) 
            vert[i] = vert[i] + movements[i];
        end
    end

end

println(getindex(vert[3],2))

plot(getindex.(vert,1),getindex.(vert,2),
    seriestype=:scatter,
    xlims = (-2,dist+2),
    ylims = (-2,getindex(vert[3],2)+2))
