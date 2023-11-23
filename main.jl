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

sim = 4

if sim == 1 # initialize a suspended nucleus

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:4

        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true)
        
        # create the crosslinks
        fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                        returnFoldername = true)
        
        # relax the whole system
        fileName3 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true,simulationDate = simulationDate)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)
        rm(".\\results\\"*fileName2; recursive = true)

        fileName4 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate)

        rm(".\\results\\"*fileName3; recursive = true)

    end

elseif sim == 2 # initialize adherent nucleus


    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")
        volumes = 1600

        simulation("INIT",250,"ADHERENT_INIT_FINAL_VI_1200","load"; adherent = true, nuclearPropPars = "nuclear_properties_temp.txt",simulationDate = simulationDate,newTargetVolume = volumes)

#     Threads.@threads for volumes = [1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000]
        # for i = 1
        #      name = simulation("INIT",20,"ADHERENT_INIT_ENLARGEMENT_VI_"*string(volumes)*"_"*string(i),"load"; importFolder = "2023-11-14_135558_ADHERENT_INIT_4_um_100_Pa_"*string(i), adherent = true, newTargetVolume = volumes, nuclearPropPars = "nuclear_properties_temp.txt",returnFoldername = true,simulationDate = simulationDate)
        #      simulation("INIT",40,"ADHERENT_INIT_FINAL_VI_"*string(volumes)*"_"*string(i),"load"; importFolder = name, adherent = true, nuclearPropPars = "nuclear_properties_temp.txt",simulationDate = simulationDate,restLengthRemodelling = true)
        # end
#     end

elseif sim == 3 # parallel example

    simulation("AFM", 1.5  , "afm_sim_test", "load"; adherentStatic = true, stickyBottom = true); 


elseif sim == 4 # mm simulation

    simulation("INIT",100,"inf_test_5","load"; vrc = true, nuclearPropPars = "nuclear_properties_temp.txt", importFolder = "2023-11-22_155526_ADHERENT_INIT_FINAL_VI_1050_1")

elseif sim == 5 # ma simulation

    simulation("MA",50,"ma_test_5000_Pa_Vol_10000_area_10_lam_50_chro_50_visc_100","load")

elseif sim == 6

    simulation("MM", 100, "adherent_nucleus_mod_1k_op_0", "load")
    # simulation("AFM", 20  , "afm_test", "load"; importFolder = results)

elseif sim == 7

    simulation("INIT", 500, "ADH_INIT_900", "load")
    # simulation("INIT", 500, "INF_TEST_crosslink_slower_dynamics", "load"; vrc = true, adherentStatic = true)
elseif sim == 8

    simulation("INIT", 50  , "volume_test", "load")
elseif sim == 9

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")
    # simulationDate = "2023-11-02_111635_"
    # Threads.@threads 

    Threads.@threads for i = 1:108

        if any(i .== 1:12)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_1.txt", nuclearPropPars = "nuclear_properties_1.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_1.txt", nuclearPropPars = "nuclear_properties_1.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_1.txt", nuclearPropPars = "nuclear_properties_1.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_1.txt", nuclearPropPars = "nuclear_properties_1.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_1.txt", nuclearPropPars = "nuclear_properties_1.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 13:24)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_2.txt", nuclearPropPars = "nuclear_properties_2.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_2.txt", nuclearPropPars = "nuclear_properties_2.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_2.txt", nuclearPropPars = "nuclear_properties_2.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_2.txt", nuclearPropPars = "nuclear_properties_2.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_2.txt", nuclearPropPars = "nuclear_properties_2.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 25:36)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_3.txt", nuclearPropPars = "nuclear_properties_3.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_3.txt", nuclearPropPars = "nuclear_properties_3.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_3.txt", nuclearPropPars = "nuclear_properties_3.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_3.txt", nuclearPropPars = "nuclear_properties_3.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_3.txt", nuclearPropPars = "nuclear_properties_3.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 37:48)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_4.txt", nuclearPropPars = "nuclear_properties_4.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_4.txt", nuclearPropPars = "nuclear_properties_4.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_4.txt", nuclearPropPars = "nuclear_properties_4.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_4.txt", nuclearPropPars = "nuclear_properties_4.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_4.txt", nuclearPropPars = "nuclear_properties_4.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 49:60)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_5.txt", nuclearPropPars = "nuclear_properties_5.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_5.txt", nuclearPropPars = "nuclear_properties_5.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_5.txt", nuclearPropPars = "nuclear_properties_5.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_5.txt", nuclearPropPars = "nuclear_properties_5.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_5.txt", nuclearPropPars = "nuclear_properties_5.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 61:72)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_6.txt", nuclearPropPars = "nuclear_properties_6.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_6.txt", nuclearPropPars = "nuclear_properties_6.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_6.txt", nuclearPropPars = "nuclear_properties_6.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_6.txt", nuclearPropPars = "nuclear_properties_6.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_6.txt", nuclearPropPars = "nuclear_properties_6.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 73:84)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_7.txt", nuclearPropPars = "nuclear_properties_7.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_7.txt", nuclearPropPars = "nuclear_properties_7.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_7.txt", nuclearPropPars = "nuclear_properties_7.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_7.txt", nuclearPropPars = "nuclear_properties_7.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_7.txt", nuclearPropPars = "nuclear_properties_7.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 85:96)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_8.txt", nuclearPropPars = "nuclear_properties_8.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_8.txt", nuclearPropPars = "nuclear_properties_8.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_8.txt", nuclearPropPars = "nuclear_properties_8.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_8.txt", nuclearPropPars = "nuclear_properties_8.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_8.txt", nuclearPropPars = "nuclear_properties_8.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        elseif any(i .== 97:108)
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true, nuclearMechPars = "nuclear_mechanics_9.txt", nuclearPropPars = "nuclear_properties_9.txt")

            # create the crosslinks
            fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; noEnveSolve = true, importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true, nuclearMechPars = "nuclear_mechanics_9.txt", nuclearPropPars = "nuclear_properties_9.txt")

            # relax the whole system
            fileName3 = simulation("INIT", 20, "NEW_INIT_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                    returnFoldername = true,simulationDate = simulationDate, nuclearMechPars = "nuclear_mechanics_9.txt", nuclearPropPars = "nuclear_properties_9.txt")

            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
            rm(".\\results\\"*fileName2; recursive = true)

            fileName4 = simulation("INIT",40,"ADHERENT_INIT_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_9.txt", nuclearPropPars = "nuclear_properties_9.txt")

            rm(".\\results\\"*fileName3; recursive = true)

            simulation("AFM", 1.5  , "afm_sim_final_set_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4, simulationDate = simulationDate,
                    nuclearMechPars = "nuclear_mechanics_9.txt", nuclearPropPars = "nuclear_properties_9.txt")

            rm(".\\results\\"*fileName4; recursive = true)
        end
    end
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
elseif sim == 14

    simulationDate = "2023-10-21_235643_"
    
    Threads.@threads for i = 13:18

        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_3.txt")
        
        # create the crosslinks
        fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
                        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_3.txt")
        
        # relax the whole system
        fileName3 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_3.txt",simulationDate = simulationDate)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)
        rm(".\\results\\"*fileName2; recursive = true)

        fileName4 = simulation("INIT",40,"ADHERENT_INIT_4_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true,
                        nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_3.txt",simulationDate = simulationDate)

        rm(".\\results\\"*fileName3; recursive = true)

        simulation("AFM", 1.5  , "afm_sim_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4,
                        nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_3.txt",simulationDate = simulationDate)

        rm(".\\results\\"*fileName4; recursive = true)

    end

elseif sim == 15

    simulationDate = "2023-10-21_235643_"
    
    Threads.@threads for i = 67:72

        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_12.txt")

        # create the crosslinks
        fileName2 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_2.txt",
        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_12.txt")

        # relax the whole system
        fileName3 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true, nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_12.txt",simulationDate = simulationDate)

        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)
        rm(".\\results\\"*fileName2; recursive = true)

        fileName4 = simulation("INIT",40,"ADHERENT_INIT_4_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName3, returnFoldername = true,
        nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_12.txt",simulationDate = simulationDate)

        rm(".\\results\\"*fileName3; recursive = true)

        simulation("AFM", 1.5  , "afm_sim_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName4,
        nuclearPropPars = "nuclear_properties_1.txt", nuclearMechPars = "nuclear_mechanics_12.txt",simulationDate = simulationDate)

        rm(".\\results\\"*fileName4; recursive = true)


    end
end