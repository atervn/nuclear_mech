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
# using ProfileView
using WriteVTK
using DelimitedFiles
using Dates

include("create_nucleus.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")
include("misc_functions.jl")
include("mesh_generation.jl")
include("create_chromatin.jl")
include("simulation.jl")

Base.@kwdef mutable struct forcesType
    volume::Vector{Vec{3,Float64}} = []
    area::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    elastic::Vector{Vec{3,Float64}} = []
    envelopeRepulsion::Vector{Vec{3,Float64}} = []
    chromationRepulsion::Vector{Vec{3,Float64}} = []
end

Base.@kwdef mutable struct nucleusType
    vert::Vector{Vec{3,Float64}} = []
    neighbors::Vector{Vector{Int64}} = []
    tri::Array{Int64} = Array{Int64}[]
    edges::Array{Int64} = Array{Int64}[]
    edgeVectors::Vector{Vec{3,Float64}} = []
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
    tri21::Array{Any} = [];
    tri31::Array{Any} = [];
    tri12::Array{Any} = [];
    tri13::Array{Any} = [];
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
    chroVertexDistance::Float64 = 0;
    chromatinNumber::UInt64 = 0;
    chromatinLength::UInt64 = 0;
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
    freeNucleusRadius::Float64 = 0;
    chroVertexDistance::Float64 = 0;
    chromatinNumber::UInt64 = 0;
    chromatinLength::UInt64 = 0;
end

Base.@kwdef mutable struct chromatinForceType
    linear::Vector{Vec{3,Float64}} = []
    bending::Vector{Vec{3,Float64}} = []
    chroRepulsion::Vector{Vec{3,Float64}} = []
    enveRepulsion::Vector{Vec{3,Float64}} = []
    strandLinear::Vector{Any} = []
    strandBending::Vector{Any} = []
    strandChroRepulsion::Vector{Any} = []
    strandEnveRepulsion::Vector{Any} = []
end

Base.@kwdef mutable struct chromatinType
    vert::Vector{Vec{3,Float64}} = []
    strandIdx::Vector{Vector{Int64}} = []
    strandVert::Vector{Any} = []
    forces = chromatinForceType()
    vectors::Vector{Vector{Vec{3,Float64}}} = []
    vectorNorms::Vector{Vector{Float64}} = []
end

Base.@kwdef mutable struct pipetteType
    vert::Vector{Vec{3,Float64}} = []
    tri::Array{Int64} = Array{Int64}[]
    vertexTri::Array{Vector{Int64}} = Array{Int64}[]
end

Base.@kwdef mutable struct micromanipulationType
    rightmostVertex::Int64 = 0
    leftmostVertex::Int64 = 0
    leftmostVertexPosition::Vec{3,Float64} = Vec(0.,0.,0.)
    rightNeighbors::Vector{Int64} = []
    leftNeighbors::Vector{Int64} = []
    leftNeigborPositions::Vector{Vec{3,Float64}} = []
end

simulation("MM",5,"misc")

simulation("MM",2000,"mm_sim")