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

simulationDate = "2023-12-13_101010_"

Threads.@threads for i = 1:64

    if any(i .== 1:4)
        mechFile = "nuclear_mechanics_1.txt"
    elseif any(i .== 5:8)
        mechFile = "nuclear_mechanics_2.txt"
    elseif any(i .== 9:12)
        mechFile = "nuclear_mechanics_3.txt"
    elseif any(i .== 13:16)
        mechFile = "nuclear_mechanics_4.txt"
    elseif any(i .== 17:20)
        mechFile = "nuclear_mechanics_5.txt"   
    elseif any(i .== 21:24)
        mechFile = "nuclear_mechanics_6.txt"
    elseif any(i .== 25:28)
        mechFile = "nuclear_mechanics_7.txt"
    elseif any(i .== 29:32)
        mechFile = "nuclear_mechanics_8.txt"
    elseif any(i .== 33:36)
        mechFile = "nuclear_mechanics_9.txt"
    elseif any(i .== 37:40)
        mechFile = "nuclear_mechanics_10.txt"
    elseif any(i .== 41:44)
        mechFile = "nuclear_mechanics_11.txt"
    elseif any(i .== 45:48)
        mechFile = "nuclear_mechanics_12.txt"
    elseif any(i .== 49:52)
        mechFile = "nuclear_mechanics_13.txt"
    elseif any(i .== 53:56)
        mechFile = "nuclear_mechanics_14.txt"
    elseif any(i .== 57:60)
        mechFile = "nuclear_mechanics_15.txt"
    elseif any(i .== 61:64)
        mechFile = "nuclear_mechanics_16.txt"
    end

    fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new";
                    noEnveSolve = true,
                    simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,
                    noChromatin = true,
                    simulationDate = simulationDate,
                    nuclearMechPars = mechFile)
    
    fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load";
                    importFolder = fileName1,
                    simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,
                    simulationDate = simulationDate,
                    noChromatin = true,
                    nuclearMechPars = mechFile)

    fileName3 = simulation("INIT",40,"ADHERENT_INIT_TEMP_"*string(i),"load";
                    adherent = true,
                    importFolder = fileName2,
                    returnFoldername = true,
                    simulationDate = simulationDate,
                    noChromatin = true,
                    nuclearMechPars = mechFile)

    fileName4 = simulation("INIT" ,5, "init_P1_2_"*string(i), "load";
                    adherent = true,
                    importFolder = fileName3,
                    noEnveSolve = true,
                    simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,
                    simulationDate = simulationDate,
                    nuclearMechPars = mechFile)

    # create the crosslinks
    fileName5 = simulation("INIT" ,1000, "init_P2_2_"*string(i), "load";
                    adherent = true,
                    noEnveSolve = true,
                    importFolder = fileName4,
                    simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true,
                    simulationDate = simulationDate,
                    nuclearMechPars = mechFile)

    # relax the whole system
    fileName6 = simulation("INIT", 20, "ADHERENT_INIT_"*string(i), "load";
                    adherent = true,
                    importFolder = fileName5,
                    simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,
                    simulationDate = simulationDate,
                    nuclearMechPars = mechFile)
    
    # remove the extra results
    rm(".\\results\\"*fileName1; recursive = true)
    rm(".\\results\\"*fileName2; recursive = true)
    rm(".\\results\\"*fileName3; recursive = true)
    rm(".\\results\\"*fileName4; recursive = true)
    rm(".\\results\\"*fileName5; recursive = true)

    simulation("AFM", 3  , "AFM_FITTING_NI_"*string(i), "load";
                    adherentStatic = true,
                    stickyBottom = true,
                    importFolder = fileName6,
                    simulationDate = simulationDate,
                    nuclearMechPars = mechFile)

end