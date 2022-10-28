using Plots, Statistics, LinearAlgebra, IterativeSolvers, SparseArrays, ProgressMeter,
Meshes, FileIO, MeshIO, NearestNeighbors, ProfileView, WriteVTK, DelimitedFiles, Dates,
StatsBase, ReadVTK, NativeFileDialog, Random, IncompleteLU, LinearSolve, Preconditioners, LimitedLDLFactorizations, LinearSolvePardiso

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

# simulation("MM",0.1,"mm_te","load";importFolder = "2022-10-18_141029_final",exportData = false)
Threads.@threads for i = 1:35
    printstyled("Starting simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :cyan)
    simulation("MM",10,"mm_sim_params_"*string(i),"load"; importFolder = "2022-10-18_141029_final", parameterFile = "./parameters/parameters_" * string(i) * ".txt", nameDate = false);
    printstyled("Finishing simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :cyan)
end