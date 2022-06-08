using Plots
using Statistics

include("sphere_creation.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")

mutable struct nucleusType
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    tri::Array{Int64}
    edges::Array{Int64}
    mirrorEdges::Vector{Int64}
    firstEdges::Vector{Int64}
    vertexEdges::Vector{Vector{Int64}}
    vertexTri::Array{Vector{Int64}}
    edgesTri::Array{Int64}
    voronoiAreas::Vector{Vector{Float64}}
    curvatures::Vector{Float64}
    vertexNormalUnitVectors::Array{Float64}
    triangleNormalUnitVectors::Array{Float64}
    areaUnitVectors::Vector{Array{Float64}}
    normalVolume::Float64
    normalArea::Float64
    normalAngle::Float64
end

struct parameters
    bulkModulus::Float64
    areaCompressionStiffness::Float64
end

nucleus = nucleusType([],[],[],Array{Int64}[],Array{Int64}[],[],[],[],Array{Int64}[],Array{Int64}[],[],[],Array{Float64}[],Array{Float64}[],[],0,0,0);

radius = 1;
nSubdivisions = 2;

nucleus = create_icosahedron!(nucleus,radius);
nucleus = get_edges!(nucleus);
nucleus = subdivide_mesh!(nucleus,radius,nSubdivisions)
get_triangle_normals!(nucleus);

@time begin
nucleus.normalVolume = get_volume!(nucleus);
nucleus.normalArea = get_area!(nucleus);
nucleus.normalAngle = mean(get_triangle_angles(nucleus));
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


#end

#plot_sphere!(nucleus)

end

