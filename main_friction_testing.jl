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
include("analysis_functions.jl")

sim = 32

if sim == 31 # initialize 8 hpi nucleus

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	ii = 0;
	# Threads.@threads 
	Threads.@threads for i = 1:4

		fileName2 = simulation("INIT",80,"ADHERENT_INIT_8_um_100_Pa_"*string(i),"new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

		fileName3 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, importFolder = fileName2, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true)

		# create the crosslinks
		fileName4 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, noEnveSolve = true, importFolder = fileName3, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true)

		# relax the whole system
		fileName5= simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, importFolder = fileName4,
			returnFoldername = true,simulationDate = simulationDate)

		rm(".\\results\\"*fileName2; recursive = true)
		rm(".\\results\\"*fileName3; recursive = true)
		rm(".\\results\\"*fileName4; recursive = true)

		# simulation("AFM", 5  , "NI_AFM_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName5, simulationDate = simulationDate)
	
		# rm(".\\results\\"*fileName5; recursive = true)

    end

elseif sim == 32 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")
	
		
	# Threads.@threads 
	for i = 1

		fileName6 = simulation("INIT",300,"INF_8hpi_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-02-14_210117_INIT_NI_"*lpad(i,1,"0"), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720)

		simulation("AFM", 5  , "NI_AFM_"*lpad(i,2,"0"), "load"; vrc = true, vrcGrowth = false, adherentStatic = true, stickyBottom = true, importFolder = fileName6, simulationDate = simulationDate)
	
		# rm(".\\results\\"*fileName5; recursive = true)

	end

elseif sim == 33

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName2 = simulation("INIT",80,"ADHERENT_INIT_8_um_100_Pa","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

	fileName6 = simulation("INIT",100,"INF_8hpi","load"; importFolder = fileName2,  returnFoldername = true, adherent = true, simulationDate = simulationDate, newTargetVolume = 720, noChromatin = true)

	fileName2 = simulation("INIT",80,"ADHERENT_test","load"; importFolder = fileName6, adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true)

end