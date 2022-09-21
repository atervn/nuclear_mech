using Plots, Statistics, LinearAlgebra, IterativeSolvers, SparseArrays, ProgressMeter,
Meshes, FileIO, MeshIO, NearestNeighbors, ProfileView, WriteVTK, DelimitedFiles, Dates,
StatsBase

include("create_nucleus.jl")
include("plotting.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")
include("misc_functions.jl")
include("mesh_generation.jl")
include("create_chromatin.jl")
include("simulation.jl")
include("lad_creation.jl")

if !(@isdefined nucleusType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

# simulation("PC",5,"misc")

simulation("PC",1000,"blaa")