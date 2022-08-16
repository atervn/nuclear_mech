using Plots
using Statistics
using LinearAlgebra
using IterativeSolvers
using SparseArrays
using ProgressMeter
using Meshes
using FileIO
using MeshIO
using NearestNeighbors
using ProfileView

include("sphere_creation.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")
include("misc_functions.jl")
include("mesh_generation.jl")

Base.@kwdef mutable struct forcesType
    volume::Vector{Vec{3,Float64}} = []
    area::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    elastic::Vector{Vec{3,Float64}} = []
    repulsion::Vector{Vec{3,Float64}} = []
end

Base.@kwdef mutable struct nucleusType
    vert::Vector{Vec{3,Float64}} = []
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
    vertexNormalUnitVectors::Vector{Vec{3,Float64}} = []
    triangleNormalUnitVectors::Vector{Vec{3,Float64}} = []
    areaUnitVectors::Vector{Vector{Vec{3,Float64}}} = []
    normalVolume::Float64 = 0
    normalArea::Float64 = 0
    normalAngle::Float64 = 0
    normalLengths::Vector{Float64} = []
    normalTriangleAreas::Vector{Float64} = []
    neighboringTriangles::Array{Int64} = []
    forces = forcesType()
    testview::Array{Any} = [];
    p1::Array{Any} = [];
    p2::Array{Any} = [];
    p3::Array{Any} = [];
    ep1::Array{Any} = [];
    ep2::Array{Any} = [];
    ep31::Array{Any} = [];
    ep32::Array{Any} = [];
    trii::Array{Any} = [];
end

Base.@kwdef mutable struct inputParametersType
    
    freeNucleusRadius::Float64 = 0;
    laminaYoung::Float64 = 0;
    laminaThickness::Float64 = 0;
    laminaFriction::Float64 = 0;
    areaCompressionModulus::Float64 = 0;
    poissonsRatio::Float64 = 0;
    bulkModulus::Float64 = 0;
    viscosity::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    scalingLength::Float64 = 0;
    scalingTime::Float64 = 0;

end

Base.@kwdef mutable struct scaledParametersType
    
    bulkModulus::Float64 = 0;
    areaCompressionStiffness::Float64 = 0;
    bendingStiffness::Float64 = 0;
    laminaStiffness::Float64 = 0;
    repulsionConstant::Float64 = 0;
    repulsionDistance::Float64 = 0;
    laminaFriction::Float64 = 0;
end

Base.@kwdef mutable struct pipetteType
    vert::Vector{Vec{3,Float64}} = []
    x::Vector{Float64} = []
    y::Vector{Float64} = []
    z::Vector{Float64} = []
    tri::Array{Int64} = Array{Int64}[]
    vertexTri::Array{Vector{Int64}} = Array{Int64}[]
end

include("main.jl")

# ProfileView.@profview main(100)

ProfileView.@profview main(1000)