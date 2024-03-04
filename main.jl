using Statistics, LinearAlgebra, IterativeSolvers, SparseArrays,
ProgressMeter, Meshes, FileIO, MeshIO, NearestNeighbors, WriteVTK,
DelimitedFiles, Dates, StatsBase, ReadVTK, NativeFileDialog, Random,
IncompleteLU,CSV,DataFrames,ZipFile

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

sim = 31

if sim == 1 # initialize a suspended nucleus

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    for i = 1:2

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
        

        Threads.@threads for i = 1:12

            simulation("INIT",250,"ADHERENT_INIT_FINAL_VI_"*string(volumes)*"_"*string(i),"load"; importFolder = "2024-01-04_112804_INIT_NI_"*string(i), adherent = true, nuclearPropPars = "nuclear_properties_temp.txt",simulationDate = simulationDate,newTargetVolume = 875)

        end
 
#     Threads.@threads for volumes = [1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000]
        # for i = 1
        #      name = simulation("INIT",20,"ADHERENT_INIT_ENLARGEMENT_VI_"*string(volumes)*"_"*string(i),"load"; importFolder = "2023-11-14_135558_ADHERENT_INIT_4_um_100_Pa_"*string(i), adherent = true, newTargetVolume = volumes, nuclearPropPars = "nuclear_properties_temp.txt",returnFoldername = true,simulationDate = simulationDate)
        #      simulation("INIT",40,"ADHERENT_INIT_FINAL_VI_"*string(volumes)*"_"*string(i),"load"; importFolder = name, adherent = true, nuclearPropPars = "nuclear_properties_temp.txt",simulationDate = simulationDate,restLengthRemodelling = true)
        # end
#     end

elseif sim == 3 # parallel example

        Threads.@threads for i = 5:8
                simulation("AFM", 1.5  , "AFM_NI_"*string(i), "load"; adherentStatic = true, stickyBottom = true, importFolder = "2023-11-14_135558_ADHERENT_INIT_4_um_100_Pa_"*string(i)); 
        end

elseif sim == 4 # mm simulation

    simulation("INIT",150,"inf_test_5","load"; vrc = true, nuclearPropPars = "nuclear_properties_temp.txt", importFolder = "2023-11-23_122954_ADHERENT_INIT_FINAL_VI_1200_1")

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

elseif sim == 10

    simulation("INIT" ,5, "test", "new"; laminaDisintegration = 0.001)
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

elseif sim == 30 # initialize NI nucleus

        simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")
    
        #Threads.@threads for i = 1:4
            Threads.@threads for i = 1:12
            # create the nucleus and let chromatin relax around the LADs
            fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                            returnFoldername = true, noChromatin = true)
                        
            # relax the whole system
            fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_1.txt",
                            returnFoldername = true,simulationDate = simulationDate, noChromatin = true)
            
            # remove the extra results
            rm(".\\results\\"*fileName1; recursive = true)
    
            fileName3 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName2, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)
    
            fileName4 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName3, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
            returnFoldername = true)

            # create the crosslinks
            fileName5 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName4, simPars = "./parameters/simulation_parameters_init_2.txt",
                        returnFoldername = true)

            # relax the whole system
            fileName6 = simulation("INIT", 100, "NI_"*string(i), "load"; adherent = true, importFolder = fileName5,
                returnFoldername = true,simulationDate = simulationDate)


            rm(".\\results\\"*fileName2; recursive = true)
            rm(".\\results\\"*fileName3; recursive = true)
            rm(".\\results\\"*fileName4; recursive = true)
            rm(".\\results\\"*fileName5; recursive = true)

        end

elseif sim == 31 # initialize 8 hpi nucleus

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:4
        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, noChromatin = true)
                    
        # relax the whole system
        fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true,simulationDate = simulationDate, noChromatin = true)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)

        fileName3 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName2, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

        fileName4 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName3, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true)

        # create the crosslinks
        fileName5 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName4, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true)

        # relax the whole system
        fileName6 = simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, importFolder = fileName5,
            returnFoldername = true,simulationDate = simulationDate)


        rm(".\\results\\"*fileName2; recursive = true)
        rm(".\\results\\"*fileName3; recursive = true)
        rm(".\\results\\"*fileName4; recursive = true)
        rm(".\\results\\"*fileName5; recursive = true)

        fileName7 = simulation("INIT",250,"INF_12hpi_"*lpad(i,2,"0"),"load"; returnFoldername = true, importFolder = fileName6, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate,newTargetVolume = 1025, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt")
 
        rm(".\\results\\"*fileName6; recursive = true)

        fileName8 = simulation("INIT",250,"INF_12hpi_VRC_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = fileName7, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt",replPars = "./parameters/replication_compartment_mechanics_12hpi.txt")

        simulation("AFM", 5  , "INF_12hpi_AFM_lam_deg_0.2_0.0001_"*lpad(i,2,"0"), "load"; vrc = true, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",
                    replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", adherentStatic = true, stickyBottom = true, importFolder = fileName8, simulationDate = simulationDate, laminaDisintegration = 0.20);

    end

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:4
        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, noChromatin = true)
                    
        # relax the whole system
        fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true,simulationDate = simulationDate, noChromatin = true)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)

        fileName3 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName2, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

        fileName4 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName3, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true)

        # create the crosslinks
        fileName5 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName4, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true)

        # relax the whole system
        fileName6 = simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, importFolder = fileName5,
            returnFoldername = true,simulationDate = simulationDate)


        rm(".\\results\\"*fileName2; recursive = true)
        rm(".\\results\\"*fileName3; recursive = true)
        rm(".\\results\\"*fileName4; recursive = true)
        rm(".\\results\\"*fileName5; recursive = true)

        fileName7 = simulation("INIT",250,"INF_12hpi_"*lpad(i,2,"0"),"load"; returnFoldername = true, importFolder = fileName6, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate,newTargetVolume = 1025, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt")
 
        rm(".\\results\\"*fileName6; recursive = true)

        fileName8 = simulation("INIT",250,"INF_12hpi_VRC_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = fileName7, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt",replPars = "./parameters/replication_compartment_mechanics_12hpi.txt")

        simulation("AFM", 5  , "INF_12hpi_AFM_lam_deg_0.5_0.0001_"*lpad(i,2,"0"), "load"; vrc = true, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",
                    replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", adherentStatic = true, stickyBottom = true, importFolder = fileName8, simulationDate = simulationDate, laminaDisintegration = 0.50);

    end

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:4
        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, noChromatin = true)
                    
        # relax the whole system
        fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true,simulationDate = simulationDate, noChromatin = true)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)

        fileName3 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName2, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

        fileName4 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName3, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true)

        # create the crosslinks
        fileName5 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName4, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true)

        # relax the whole system
        fileName6 = simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, importFolder = fileName5,
            returnFoldername = true,simulationDate = simulationDate)


        rm(".\\results\\"*fileName2; recursive = true)
        rm(".\\results\\"*fileName3; recursive = true)
        rm(".\\results\\"*fileName4; recursive = true)
        rm(".\\results\\"*fileName5; recursive = true)

        fileName7 = simulation("INIT",250,"INF_12hpi_"*lpad(i,2,"0"),"load"; returnFoldername = true, importFolder = fileName6, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate,newTargetVolume = 1025, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt")
 
        rm(".\\results\\"*fileName6; recursive = true)

        fileName8 = simulation("INIT",250,"INF_12hpi_VRC_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = fileName7, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt",replPars = "./parameters/replication_compartment_mechanics_12hpi.txt")

        simulation("AFM", 5  , "INF_12hpi_AFM_lam_deg_0.8_0.0001_"*lpad(i,2,"0"), "load"; vrc = true, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",
                    replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", adherentStatic = true, stickyBottom = true, importFolder = fileName8, simulationDate = simulationDate, laminaDisintegration = 0.80);

    end

elseif sim == 32 # initialize 12 hpi nucleus

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    #Threads.@threads for i = 1:4
        Threads.@threads for i = 1:12
        # create the nucleus and let chromatin relax around the LADs
        fileName1 = simulation("INIT" ,5, "init_P1_"*string(i), "new"; noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true, noChromatin = true)
                    
        # relax the whole system
        fileName2 = simulation("INIT", 20, "NEW_INIT_100_Pa_"*string(i), "load"; importFolder = fileName1, simPars = "./parameters/simulation_parameters_init_1.txt",
                        returnFoldername = true,simulationDate = simulationDate, noChromatin = true)
        
        # remove the extra results
        rm(".\\results\\"*fileName1; recursive = true)

        fileName3 = simulation("INIT",40,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"load"; adherent = true,importFolder = fileName2, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

        fileName4 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName3, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
        returnFoldername = true)

        # create the crosslinks
        fileName5 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName4, simPars = "./parameters/simulation_parameters_init_2.txt",
                    returnFoldername = true)

        # relax the whole system
        fileName6 = simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, importFolder = fileName5,
            returnFoldername = true,simulationDate = simulationDate)


        rm(".\\results\\"*fileName2; recursive = true)
        rm(".\\results\\"*fileName3; recursive = true)
        rm(".\\results\\"*fileName4; recursive = true)
        rm(".\\results\\"*fileName5; recursive = true)

        simulation("INIT",250,"INF_12hpi_"*lpad(i,2,"0"),"load"; importFolder = fileName6, adherent = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt",simulationDate = simulationDate,newTargetVolume = 875)
            
        rm(".\\results\\"*fileName6; recursive = true)

    end

elseif sim == 33

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:12
        simulation("INIT",300,"INF_8hpi_"*lpad(i,2,"0"),"load"; vrc = true, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", importFolder = "2024-01-08_094354_INF_8hpi_"*lpad(i,2,"0"),simulationDate = simulationDate)
    end
elseif sim == 34
    
    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    Threads.@threads for i = 1:12
        simulation("INIT",300,"INF_12hpi_"*lpad(i,2,"0"),"load"; vrc = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", replPars = "./parameters/replication_compartment_mechanics_12hpi.txt",  importFolder = "2024-01-08_121324_INF_12hpi_"*lpad(i,2,"0"),simulationDate = simulationDate)
    end

elseif sim == 35

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

    # Threads.@threads for i = 1:4
    #     simulation("INIT", 50  , "NI_relax_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true,importFolder = "2024-01-08_105449_NI_"*string(i),simulationDate = simulationDate,laminaDisintegration = 0.6);
    # end

    # Threads.@threads for i = 1:4
    #     simulation("AFM", 5  , "NI_AFM_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true,importFolder = "2024-01-08_105449_NI_"*string(i),simulationDate = simulationDate,laminaDisintegration = 0.6);
    # end
    Threads.@threads for i = 1:4
        simulation("AFM", 5  , "INF_8hpi_AFM_lam_dis_0.95_0.0001_"*lpad(i,2,"0"), "load"; vrc = true, nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", adherentStatic = true, stickyBottom = true,importFolder = "2024-01-08_142634_INF_8hpi_"*lpad(i,2,"0"),simulationDate = simulationDate,laminaDisintegration = 0.95);
    end
    # Threads.@threads for i = 1:12
    #     simulation("AFM", 5  , "INF_12hpi_AFM_"*lpad(i,2,"0"), "load"; vrc = true, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", adherentStatic = true, stickyBottom = true,importFolder = "2024-01-08_152833_INF_12hpi_"*lpad(i,2,"0"),simulationDate = simulationDate);
    # end

elseif sim == 36
    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")
    Threads.@threads for i = 1:4
        simulation("INIT",250,"INF_8hpi_1.75xlamStiff"*lpad(i,2,"0"),"load"; importFolder = "2024-01-08_105449_NI_"*string(i), adherent = true, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", nuclearMechPars = "./parameters/nuclear_mechanics_temp.txt",simulationDate = simulationDate,newTargetVolume = 650)
    end
end