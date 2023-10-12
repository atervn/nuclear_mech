using Statistics, LinearAlgebra, IterativeSolvers, SparseArrays,
ProgressMeter, Meshes, FileIO, MeshIO, NearestNeighbors, WriteVTK,
DelimitedFiles, Dates, StatsBase, ReadVTK, NativeFileDialog, Random,
IncompleteLU,CSV,DataFrames

# include.(filter(contains(r".jl$"), readdir(dir; join=true)))
if !(@isdefined envelopeType)
    include("./functions/NuclearMechTypes.jl")
    using .NuclearMechTypes
end

include("./functions/create_nucleus.jl")
include("./functions/geometric_functions.jl")
include("./functions/calculate_forces.jl")
include("./functions/misc_functions.jl")
include("./functions/create_chromatin.jl")
include("./functions/lad_creation.jl")
include("./functions/setup_functions.jl")
include("./functions/import_functions.jl")
include("./functions/solve_system.jl")
include("./functions/import_functions.jl")
include("./functions/get_forces.jl")
include("./functions/simulation.jl")

sim = 9
if sim == 1 # initialize a suspended nucleus

    # create the nucleus and let chromatin relax around the LADs
    fileName1 = simulation("INIT" ,5, "init_P1", "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt", returnFoldername = true)
    
    # create the crosslinks
    fileName2 = simulation("INIT" ,1000, "init_P2", "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt", returnFoldername = true)
    
    # relax the whole system
    simulation("INIT", 200, "NEW_INIT_1x_cl", "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt")
    
    # remove the extra results
    rm(".\\results\\"*fileName1; recursive = true)
    rm(".\\results\\"*fileName2; recursive = true)

elseif sim == 2 # initialize adherent nucleus

    # simulation("INIT",400,"ADHERENT_INIT","load"; adherent = true, simPars = "./parameters/simulation_parameters_init_adh.txt")

    simulation("INIT",400,"ADHERENT_INIT_4_um","load"; adherent = true)

    # # load the squished nucleus, add chormatin and let chromatin relax around the LADs
    # fileName1 = simulation("INIT" ,10, "init_P1", "load"; importFolder = "2023-08-22_145359_adherend_shell", noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt", returnFoldername = true, newEnvelopeMultipliers = true)
    
    # # create the crosslinks
    # fileName2 = simulation("INIT" ,1000, "init_P2", "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt", returnFoldername = true)
    
    # # relax the whole system
    # simulation("INIT", 100, "adherent_nucleus_mod_1000", "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt")
    
    # # remove the extra results
    # rm(".\\results\\"*fileName1; recursive = true)
    # rm(".\\results\\"*fileName2; recursive = true)

elseif sim == 3 # parallel example

    # Threads.@threads for i = 1:5
    #     printstyled("Starting simulation " * string(i) *" ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n";
    #         color = :cyan)
    #     simulation("MM",150,"MM_SIM_PARAMS_" * string(i),"load"; importFolder = "2023-02-13_145413_INIT",
    #         parameterFile = "./parameters/parameters_" * string(i) * ".txt", nameDate = false);
    #     printstyled("Finishing simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n";
    #         color = :cyan)
    # end

elseif sim == 4 # mm simulation

    simulation("MM",100,"mm_test_5","load")

elseif sim == 5 # ma simulation

    simulation("MA",50,"ma_test_5000_Pa_Vol_10000_area_10_lam_50_chro_50_visc_100","load")

elseif sim == 6

    simulation("INIT", 100, "adherent_nucleus_mod_1k_op_0", "load")
    # simulation("AFM", 20  , "afm_test", "load"; importFolder = results)

elseif sim == 7

    simulation("INIT", 500, "ADH_INIT_900", "load")
    # simulation("INIT", 500, "INF_TEST_crosslink_slower_dynamics", "load"; vrc = true, adherentStatic = true)
elseif sim == 8

    simulation("INIT", 50  , "volume_test", "load")
elseif sim == 9
    simulation("AFM", 3  , "afm_test", "load"; adherentStatic = true)
elseif sim == 10

    simulation("INIT",50,"adherent_nucleus","load")
elseif sim == 11
    simulation("INIT" ,5000, "adherent_nucleus_mod_100_klam_20", "load")
elseif sim == 13

    notdone = Int.(readdlm("done_numbers.txt"))

    Threads.@threads for i = notdone

        printstyled("Starting simulation " * string(i) *" ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n";
            color = :cyan)
        simulation("MM",1000,"mm_test_chro_"*string(i,pad=3),"load"; importFolder = "2023-02-13_145413_INIT", parameterFile = "./parameters/parameters_" * string(i) * ".txt", nameDate = false)
        printstyled("Finishing simulation " * string(i) * " ("*Dates.format(now(), "YYYY-mm-dd HH:MM")*")\n";
            color = :cyan)
    end
end