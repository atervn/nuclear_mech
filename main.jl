using Statistics, LinearAlgebra, IterativeSolvers, SparseArrays,
ProgressMeter, Meshes, FileIO, MeshIO, NearestNeighbors, WriteVTK,
DelimitedFiles, Dates, StatsBase, ReadVTK, NativeFileDialog, Random,
IncompleteLU

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


# simulation("INIT",0.1,"fggfg","new";exportData = false)
simulation("MM",5,"MM_TEST","load")
