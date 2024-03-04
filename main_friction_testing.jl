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

sim = 2

#####################################################################################################
# initialize NI nuclei
#####################################################################################################

if sim == 1 

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",80,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720)

	# Threads.@threads 
	Threads.@threads for i = 1:32


		fileName2 = simulation("INIT" ,5, "init_P1_"*string(i), "load"; adherent = true, adherentStatic = true,  importFolder = fileName1, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true, simulationDate = simulationDate)

		# create the crosslinks
		fileName3 = simulation("INIT" ,1000, "init_P2_"*string(i), "load"; adherent = true, adherentStatic = true,  noEnveSolve = true, importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true, simulationDate = simulationDate)

		# relax the whole system
		simulation("INIT", 20, "INIT_NI_"*string(i), "load"; adherent = true, adherentStatic = true, importFolder = fileName3, simulationDate = simulationDate)

		rm(".\\results\\"*fileName2; recursive = true)
		rm(".\\results\\"*fileName3; recursive = true)

    end

	rm(".\\results\\"*fileName1; recursive = true)


#####################################################################################################
# 8 hpi test
#####################################################################################################

elseif sim == 2 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 6

	simPerCase = 12

	numInitSims = 32;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 70:numCases*simPerCase

		if any(i .== 1:12)
			dis = 0.10
		elseif any(i .== 13:24)
			dis = 0.20
		elseif any(i .== 25:36)
			dis = 0.30
		elseif any(i .== 37:48)
			dis = 0.40
		elseif any(i .== 49:60)
			dis = 0.50
		elseif any(i .== 61:72)
			dis = 0.60
		end


		fileName1 = simulation("INIT",10000,"INIT_8hpi_lamina_disintegration_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-02-25_185915_INIT_NI_"*lpad(cases[i],2,"0"), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720, laminaDisintegration = dis)

		simulation("AFM", 5  , "AFM_8hpi_lamina_disintegration_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName1, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, laminaDisintegration = dis)
	
	end

end