using Plots
using Statistics
using LinearAlgebra
using IterativeSolvers
using SparseArrays
using ProgressMeter
# using PlotlyJS

include("sphere_creation.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")
include("misc_functions.jl")

Base.@kwdef mutable struct forcesType
    volumeX::Vector{Float64} = []
    volumeY::Vector{Float64} = []
    volumeZ::Vector{Float64} = []
    areaX::Vector{Float64} = []
    areaY::Vector{Float64} = []
    areaZ::Vector{Float64} = []
    bendingX::Vector{Float64} = []
    bendingY::Vector{Float64} = []
    bendingZ::Vector{Float64} = []
    elasticX::Vector{Float64} = []
    elasticY::Vector{Float64} = []
    elasticZ::Vector{Float64} = []
    repulsionX::Vector{Float64} = []
    repulsionY::Vector{Float64} = []
    repulsionZ::Vector{Float64} = []
end

Base.@kwdef mutable struct nucleusType
    x::Vector{Float64} = []
    y::Vector{Float64} = []
    z::Vector{Float64} = []
    neighbors::Vector{Vector{Int64}} = []
    tri::Array{Int64} = Array{Int64}[]
    edges::Array{Int64} = Array{Int64}[]
    mirrorEdges::Vector{Int64} = []
    firstEdges::Vector{Int64} = []
    vertexEdges::Vector{Vector{Int64}} = []
    vertexTri::Array{Vector{Int64}} = Array{Int64}[]
    edgesTri::Array{Int64} = Array{Int64}[]
    edges3vertex::Array{Int64} = Array{Int64}[]
    voronoiAreas::Vector{Vector{Float64}} = []
    curvatures::Vector{Float64} = []
    vertexNormalUnitVectors::Array{Float64} = Array{Int64}[]
    triangleNormalUnitVectors::Array{Float64} = Array{Int64}[]
    areaUnitVectors::Vector{Array{Float64}} = []
    normalVolume::Float64 = 0
    normalArea::Float64 = 0
    normalAngle::Float64 = 0
    normalLengths::Vector{Float64} = []
    forces = forcesType()
end

Base.@kwdef mutable struct parametersType
    bulkModulus::Float64 = 0;
    areaCompressionStiffness::Float64 = 0;
    bendingStiffness::Float64 = 0;
    laminaStiffness::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    laminaFriction::Float64 = 0;
    viscosity::Float64 = 0;
end

Base.@kwdef mutable struct scaledParametersType
    
    bulkModulus::Float64 = 0;
    areaCompressionStiffness::Float64 = 0;
    bendingStiffness::Float64 = 0;
    laminaStiffness::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    laminaFriction::Float64 = 0;
    viscosity::Float64 = 0;
end

nuc = nucleusType();

par = parametersType();

ipar = read_parameters("./parameters.txt");

mpar = calculate_model_parameters(ipar);

radius = 1;
nSubdivisions = 3;

nuc = create_icosahedron!(nuc,radius);
nuc = get_edges!(nuc);
nuc = subdivide_mesh!(nuc,radius,nSubdivisions)
get_triangle_normals!(nuc);


nuc.normalVolume = get_volume!(nuc);
nuc.normalArea = get_area!(nuc);
nuc.normalAngle = mean(get_triangle_angles(nuc));
lengths = zeros(Float64,Int64(size(nuc.edges,1)));

nuc.z .+= 0

for i = 1:size(nuc.edges,1)

        vector = [nuc.x[nuc.edges[i,2]] - nuc.x[nuc.edges[i,1]],
                  nuc.y[nuc.edges[i,2]] - nuc.y[nuc.edges[i,1]],
                  nuc.z[nuc.edges[i,2]] - nuc.z[nuc.edges[i,1]]];
                      
        lengths[i] = norm(vector);

end
nuc.normalLengths = lengths;

frictionMatrix = get_friction_matrix(nuc,par);

#println("the volume of the sphere is $nuc.normalVolume")
#println("the area of the sphere is $nuc.normalArea")

maxT = 666;

p = Progress(maxT)
anim = Animation();

# trace = PlotlyJS.mesh3d(
#     x=nuc.x,
#     y=nuc.z,
#     z=nuc.y,
#     i=nuc.tri[:,1] .- 1,
#     j=nuc.tri[:,2] .- 1,
#     k=nuc.tri[:,3] .- 1
# );
# data = [trace];
# pl = PlotlyJS.plot(data,layout)

for t = 1:maxT
    
    get_voronoi_areas!(nuc);
    get_local_curvatures!(nuc);
    get_area_unit_vectors!(nuc);
    get_triangle_normals!(nuc);

    get_volume_forces!(nuc,par);
    get_area_forces!(nuc,par);
    get_bending_forces!(nuc,par);
    get_elastic_forces!(nuc,par);
    get_repulsion_forces!(nuc,par);

    flatRepulsion = flat_repulsion_forces(nuc);
    if t >= 0
        flatRepulsion2 = flat_repulsion_forces2(nuc,t);
    else
        flatRepulsion2 = zeros(Float64,length(nuc.x),3);
    end

    local totalForcesX = nuc.forces.volumeX .+ nuc.forces.areaX .+ nuc.forces.bendingX .+ nuc.forces.elasticX .+ nuc.forces.repulsionX;
    local totalForcesY = nuc.forces.volumeY .+ nuc.forces.areaY .+ nuc.forces.bendingY .+ nuc.forces.elasticY .+ nuc.forces.repulsionY;
    local totalForcesZ = nuc.forces.volumeZ .+ nuc.forces.areaZ .+ nuc.forces.bendingZ .+ nuc.forces.elasticZ .+ nuc.forces.repulsionZ .+ flatRepulsion2[:,3];
    
    local vX = cg(frictionMatrix,totalForcesX);
    local vY = cg(frictionMatrix,totalForcesY);
    local vZ = cg(frictionMatrix,totalForcesZ);
    
    nuc.x = nuc.x .+ vX.*0.01;
    nuc.y = nuc.y .+ vY.*0.01;
    nuc.z = nuc.z .+ vZ.*0.01;

    if mod(t,5) == 0
        plot_sphere!(nuc,[vX vY vZ],t)
        Plots.frame(anim)
    end
    next!(p)
end

gif(anim, "test.gif", fps = 30)