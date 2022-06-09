using Plots
using Statistics
using LinearAlgebra
include("sphere_creation.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")

Base.@kwdef mutable struct nucleusType
    x::Vector{Float64} = []
    y::Vector{Float64} = []
    z::Vector{Float64} = []
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
    normalLength::Float64 = 0
end

struct parameters
    bulkModulus::Float64
    areaCompressionStiffness::Float64
end

nucleus = nucleusType();

radius = 1;
nSubdivisions = 2;

nucleus = create_icosahedron!(nucleus,radius);
nucleus = get_edges!(nucleus);
nucleus = subdivide_mesh!(nucleus,radius,nSubdivisions)
get_triangle_normals!(nucleus);

@time begin
nucleus.normalVolume = get_volume!(nucleus);
nucleus.normalArea = 0.5*get_area!(nucleus);
nucleus.normalAngle = mean(get_triangle_angles(nucleus));
lengths = zeros(Float64,Int64(size(nucleus.edges,1)/2),1);

local j = 1;
for i = 1:size(nucleus.edges,1)
    if nucleus.firstEdges[i] == 1

        vector = [nucleus.x[nucleus.edges[i,2]] - nucleus.x[nucleus.edges[i,1]],
                  nucleus.y[nucleus.edges[i,2]] - nucleus.y[nucleus.edges[i,1]],
                  nucleus.z[nucleus.edges[i,2]] - nucleus.z[nucleus.edges[i,1]]];
                      
        lengths[j] = norm(vector);
        j += 1;
    end
end
nucleus.normalLength = mean(lengths);


#println("the volume of the sphere is $nucleus.normalVolume")
#println("the area of the sphere is $nucleus.normalArea")

#for t = 1:10
    
    get_voronoi_areas!(nucleus);
    get_local_curvatures!(nucleus);
    get_area_unit_vectors!(nucleus);
    get_triangle_normals!(nucleus);

    volumeForces = get_volume_forces(nucleus);
    areaForces = get_area_forces(nucleus);
    bendingForces = get_bending_forces(nucleus);
    elasticForces = get_elastic_forces(nucleus);

#end

plot_sphere!(nucleus,areaForces)

end