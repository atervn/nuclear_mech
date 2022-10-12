using Plots, Statistics, LinearAlgebra, IterativeSolvers, SparseArrays, ProgressMeter,
Meshes, FileIO, MeshIO, NearestNeighbors, ProfileView, WriteVTK, DelimitedFiles, Dates,
StatsBase, ReadVTK, NativeFileDialog, Random

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

Threads.@threads for i = 1:16
    printstyled("Starting simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :blue)
    simulation("MM",1000,"mm_sim_"*string(i),"load"; importFolder = "2022-10-10_235009_init_example_3", parameterFile = "./parameters/parameters_" * string(i) * ".txt", nameDate = "no");
    printstyled("Finishing simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :blue)
end
