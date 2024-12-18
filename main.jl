using Statistics, LinearAlgebra, IterativeSolvers, SparseArrays,
ProgressMeter, GeometryBasics, FileIO, MeshIO, NearestNeighbors, WriteVTK,
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
include("./functions/utils.jl")
include("./analysis/analysis_functions.jl")

sim = 2

#####################################################################################################
# initialize NI nuclei
#####################################################################################################

if sim == 1 

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	# get one adherent shell
	fileName1 = simulation("INIT",600,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720, resetVertexDistancesTime = 400.)

	# Threads.@threads 
	Threads.@threads for i = 1:1

		# insert chromatin
		fileName2 = simulation("INIT" ,5, "init_P1_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true,  importFolder = fileName1, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true, simulationDate = simulationDate)

		# create the crosslinks
		fileName3 = simulation("INIT" ,1000, "init_P2_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true,  noEnveSolve = true, importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true, simulationDate = simulationDate)

		# relax the whole system
		fileName4 = simulation("INIT", 20, "INIT_NI_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true, importFolder = fileName3, simulationDate = simulationDate,returnFoldername = true,)

		rm(".\\results\\"*fileName2; recursive = true)
		rm(".\\results\\"*fileName3; recursive = true)

		zip("./results/"*fileName4)

		rm(".\\results\\"*fileName4; recursive = true)

    end

	rm(".\\results\\"*fileName1; recursive = true)

#####################################################################################################
# 8 hpi test
#####################################################################################################

elseif sim == 2 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	# Threads.@threads 
	Threads.@threads for i = 1:1

		inportName = "2024-12-18_092428_INIT_NI_"*lpad(i,3,"0")

		unzip(".\\results\\"*inportName*".zip")

		# relax a 8 hpi nuclei into the required volume and height
		fileName1 = simulation("INIT",10000,"INIT_8hpi_"*lpad(i,3,"0"),"load"; returnFoldername = true, importFolder = inportName, adherent = true, simulationDate = simulationDate, newTargetVolume = 720)

		rm(".\\results\\"*inportName; recursive = true)

		# run the AFM simulation
		fileName2 = simulation("AFM", 5  , "AFM_8hpi_"*lpad(i,3,"0"), "load"; returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName1, simulationDate = simulationDate)
	
		rm(".\\results\\"*fileName1; recursive = true)

		zip(".\\results\\"*fileName2)

		rm(".\\results\\"*fileName2; recursive = true)
	end
end