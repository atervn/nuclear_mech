using Statistics, LinearAlgebra, IterativeSolvers, SparseArrays,
ProgressMeter, Meshes, FileIO, MeshIO, NearestNeighbors, WriteVTK,
DelimitedFiles, Dates, StatsBase, ReadVTK, NativeFileDialog, Random,
IncompleteLU

# include.(filter(contains(r".jl$"), readdir(dir; join=true)))
if !(@isdefined envelopeType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

include("create_nucleus.jl")
include("geometric_functions.jl")
include("calculate_forces.jl")
include("misc_functions.jl")
include("mesh_generation.jl")
include("create_chromatin.jl")
include("lad_creation.jl")
include("setup_functions.jl")
include("import_functions.jl")
include("solve_system.jl")
include("import_functions.jl")
include("print_error.jl")
include("get_forces.jl")
include("simulation.jl")

sim = 1

if sim == 1 # initialize a suspended nucleus

    # create the nucleus and let chromatin relax around the LADs
    fileName1 = simulation("INIT" ,10, "init_P1", "new"; noEnveSolve = true, parameterFile = "parameters_init_1.txt", returnFoldername = true)
    
    # create the crosslinks
    fileName2 = simulation("INIT" ,1000, "init_P2", "load"; importFolder = fileName1, parameterFile = "parameters_init_2.txt", returnFoldername = true)
    
    # relax the whole system
    simulation("INIT", 100, "init_chro_stiff_30", "load"; importFolder = fileName2, parameterFile = "parameters_init_1.txt")
    
    # remove the extra results
    rm(".\\results\\"*fileName1; recursive = true)
    rm(".\\results\\"*fileName2; recursive = true)

elseif sim == 2 # initialize adherent nucleus

    # load the squished nucleus, add chormatin and let chromatin relax around the LADs
    fileName1 = simulation("INIT" ,10, "init_P1", "load"; importFolder = "adherent_enve", noEnveSolve = true, parameterFile = "parameters_init_1.txt", returnFoldername = true)
    
    # create the crosslinks
    fileName2 = simulation("INIT" ,1000, "init_P2", "load"; importFolder = fileName1, parameterFile = "parameters_init_2.txt", returnFoldername = true)
    
    # relax the whole system
    simulation("INIT", 2, "adheret_nucleus", "load"; importFolder = fileName2, parameterFile = "parameters_init_1.txt")
    
    # remove the extra results
    rm(".\\results\\"*fileName1; recursive = true)
    rm(".\\results\\"*fileName2; recursive = true)

elseif sim == 3 # parallel example

    Threads.@threads for i = 1:35
        printstyled("Starting simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :cyan)
        simulation("MM",10,"mm_sim_params_" * string(i),"load"; importFolder = "2022-10-18_141029_final", parameterFile = "./parameters/parameters_" * string(i) * ".txt", nameDate = false);
        printstyled("Finishing simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n"; color = :cyan)
    end

elseif sim == 4 # mm simulation

    simulation("MM",100,"mm_test","load")

elseif sim == 5 # ma simulation

    simulation("MA",20,"ma_test","load"; importFolder = "2023-01-09_125914_init_chro_stiff_30")

elseif sim == 6

    simulation("INIT", 2, "init_final_test", "load"; replComp = true)

elseif sim == 7

    simulation("INIT", 20, "INIT_OSMO_300_Pa_25", "load"; parameterFile = "parameters_init_1.txt")

end