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
include("simulation_adh_init.jl")
include("lad_creation.jl")
include("setup_functions.jl")
include("import_functions.jl")
include("solve_system.jl")
include("import_functions.jl")
# include("simulation_init.jl")

if !(@isdefined nucleusType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end


envelopeFolderName = simulation_adh_init("PC",1000,"fggfg_0","new")
# chromatinFolderName = add_chromatin_adh_init(envelopeFolderName,"test");
# fileName1 = simulation_init("INIT",10,"initP1","load",true; importFolder = chromatinFolderName, parameterFile = "parameters_init_1.txt")
# fileName2 = simulation_init("INIT",1000,"initP2","load",false; importFolder = fileName1, parameterFile = "parameters_init_2.txt")
# simulation_init("INIT",10,"init_final_test","load",false; importFolder = fileName2, parameterFile = "parameters_init_1.txt")
# simulation("INIT",0.5,"MM_TEST","new")