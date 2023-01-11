using Meshes, NearestNeighbors, LinearAlgebra, IterativeSolvers, Plots, Statistics

# fine mesh

nRows = 3;
nCols = 6;

k = 1.5

dist = 2.0;

xIncrement = dist;
yIncrement = dist/2*tand(60)

vert = Vector{Vec{2,Float64}}(undef,nRows*nCols+floor(Int,nRows/2))

bottomRow = []
topRow = []

global ind = 1
for row = 1:nRows
    if mod(row,2) == 1
        for col = 1:nCols
            vert[ind] = Vec((col-1)*xIncrement,(row-1)*yIncrement)
            if row == 1
                push!(bottomRow,ind)
            elseif row == nRows
                push!(topRow,ind)
            end
            global ind += 1
        end
    else
        for col = 1:nCols+1
            vert[ind] = Vec((col-1)*xIncrement - dist/2,(row-1)*yIncrement)
            if row == 1
                push!(bottomRow,ind)
            elseif row == nRows
                push!(topRow,ind)
            end
            global ind += 1
        end
    end
end

tree = KDTree(vert)
neighbors = Vector{Vector{Int}}(undef,nRows*nCols+floor(Int,nRows/2))

for i = eachindex(vert)

    closeby = inrange(tree, vert[i], dist*1.1)
    deleteat!(closeby,closeby .== i)
    neighbors[i] = closeby;
end

frictionMatrix = zeros(length(vert), length(vert))

for i = eachindex(vert)

    frictionMatrix[i,i] = 1.0 + 0.5*length(neighbors[i]);
    for j = neighbors[i]

        frictionMatrix[i,j] = 0.5

    end
end

dt = 0.1;

bottom1 = mean(getindex.(vert,2)[bottomRow])
top1 = mean(getindex.(vert,2)[topRow])

for t = 1:2000

    forces = Vector{Vec{2,Float64}}(undef,nRows*nCols+floor(Int,nRows/2))

    for i = eachindex(vert)

        forces[i] = Vec(0.,0.)
        for j = neighbors[i]

            vector = vert[j] - vert[i]
            vectorLength = norm(vector)
            unitVector = vector./vectorLength

            forces[i] += -k*(dist - vectorLength)*unitVector;

        end

        if any(bottomRow .== i) 
            forces[i] += Vec(0.,-1.)
        elseif any(topRow .== i) 
            forces[i] += Vec(0.,1.)
        end

    end

    solX = cg(frictionMatrix,getindex.(forces,1))
    solY = cg(frictionMatrix,getindex.(forces,2))

    movements = Vec.(solX,solY).*dt

    vert .= vert .+ movements;
end


bottom2 = mean(getindex.(vert,2)[bottomRow])
top2 = mean(getindex.(vert,2)[topRow])

strain = (top2 - bottom2 - (top1 - bottom1))/(top1 - bottom1)

println("Strain: $strain" )

plot(getindex.(vert,1),getindex.(vert,2),
    seriestype=:scatter,
    xlims = (-2,(nCols-1)*dist+2),
    ylims = (-2,(nRows-1)*yIncrement + 2))
