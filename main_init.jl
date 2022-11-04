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
include("simulation_init.jl")
include("lad_creation.jl")

if !(@isdefined nucleusType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

# simulation("INIT",0.1,"fggfg","new";exportData = false)
Threads.@threads for i = 1:5
    fileName1 = simulation_init("INIT",10,"initP1","new",true; parameterFile = "parameters_init_1.txt")
    fileName2 = simulation_init("INIT",1000,"initP2","load",false; importFolder = fileName1, parameterFile = "parameters_init_2.txt")
    simulation_init("INIT",10,"init_final"*string(i),"load",false; importFolder = fileName2, parameterFile = "parameters_init_1.txt")
    rm(".\\results\\"*fileName1; recursive = true)
    rm(".\\results\\"*fileName2; recursive = true)
end