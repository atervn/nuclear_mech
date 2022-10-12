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

# simulation("PC",5,"misc")
# @ProfileView.profview
# simulation("INIT",5,"crosslink_testing","load";importFolder = "blaaa_61")
# nCrosslinks = Vector{Vector{Int64}}(undef,20)
# for i = 1:20
# simulation("INIT",30,"sim_init_1","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_2","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_3","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_4","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_5","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_6","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_7","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_8","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_9","new"; nameDate = "no")
# simulation("INIT",30,"sim_init_10","new"; nameDate = "no")
# end

# nCrosslinks = Vector{Vector{Int64}}(undef,10)

test = simulation("INIT",60,"init_testing","load"; importFolder = "2022-10-10_235009_init_example_3")

# nCrosslinks[1] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_1")
# nCrosslinks[2] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_2")
# nCrosslinks[3] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_3")
# nCrosslinks[4] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_4")
# nCrosslinks[5] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_5")
# nCrosslinks[6] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_6")
# nCrosslinks[7] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_7")
# nCrosslinks[8] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_8")
# nCrosslinks[9] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_9")
# nCrosslinks[10] = simulation("INIT",3000,"crosslink_testing","load"; importFolder = "sim_init_10")

# l = length(nCrosslinks[1])
# a = getindex.(nCrosslinks,l)
# b = mean(a)
# println(b)
# b/46/128*2