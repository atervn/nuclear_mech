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
include("./analysis/analysis_functions.jl")

sim = 9

#####################################################################################################
# initialize NI nuclei
#####################################################################################################

if sim == 1 

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",600,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720, resetVertexDistancesTime = 400.)

	# Threads.@threads 
	Threads.@threads for i = 1:4


		fileName2 = simulation("INIT" ,5, "init_P1_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true,  importFolder = fileName1, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true, simulationDate = simulationDate)

		# create the crosslinks
		fileName3 = simulation("INIT" ,1000, "init_P2_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true,  noEnveSolve = true, importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true, simulationDate = simulationDate)

		# relax the whole system
		fileName4 = simulation("INIT", 20, "INIT_NI_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true, importFolder = fileName3, simulationDate = simulationDate,returnFoldername = true,)



    end


#####################################################################################################
# 8 hpi test
#####################################################################################################

elseif sim == 2 # initialize 8 hpi nucleus

	simulationDate = "2024-03-18_174556_" #Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 10

	simPerCase = 12

	numInitSims = 64;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 61:numCases*simPerCase

		if any(i .== 1:12)
			dis = 0.0	
		elseif any(i .== 13:24)
			dis = 0.1
		elseif any(i .== 25:36)
			dis = 0.2
		elseif any(i .== 37:48)
			dis = 0.3
		elseif any(i .== 49:60)
			dis = 0.4
		elseif any(i .== 61:72)
			dis = 0.0
		elseif any(i .== 73:84)
			dis = 0.1		
		elseif any(i .== 85:96)
			dis = 0.2
		elseif any(i .== 95:108)
			dis = 0.3
		elseif any(i .== 109:120)
			dis = 0.4
		end

		if any(i .== 1:60)
			name = "lam_disintegration_"
			offset = 0
			nucMec = "./temp_pars/nuclear_mechanics_dis_01.txt"
		else
			name = "lam_disintegration_lam_increase_"
			offset = 60
			nucMec = "./temp_pars/nuclear_mechanics_dis_02.txt"
		end

		fileName1 = simulation("INIT",10000,"INIT_8hpi_"*name*lpad(i-offset,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-03-15_073352_INIT_NI_"*lpad(cases[i],2,"0"), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720, nuclearMechPars = nucMec, laminaDisintegration = dis)

		simulation("AFM", 5  , "AFM_8hpi_"*name*lpad(i-offset,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName1, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = dis)
	
	end

elseif sim == 3 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 10

	simPerCase = 12

	numInitSims = 64;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 61:numCases*simPerCase

		if any(i .== 1:12)
			nucMec = ".\\temp_pars\\nuclear_mechanics_01.txt"
		elseif any(i .== 13:24)
			nucMec = ".\\temp_pars\\nuclear_mechanics_02.txt"
		elseif any(i .== 25:36)
			nucMec = ".\\temp_pars\\nuclear_mechanics_03.txt"
		elseif any(i .== 37:48)
			nucMec = ".\\temp_pars\\nuclear_mechanics_04.txt"
		elseif any(i .== 49:60)
			nucMec = ".\\temp_pars\\nuclear_mechanics_05.txt"
		elseif any(i .== 61:72)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_01.txt"
		elseif any(i .== 73:84)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_02.txt"
		elseif any(i .== 85:96)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_03.txt"
		elseif any(i .== 95:108)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_04.txt"
		elseif any(i .== 109:120)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_05.txt"
		end

		if any(i .== 1:60)
			name = "osmo_reduction_"
			offset = 0
		else
			name = "osmo_reduction_lam_increase_"
			offset = 60
		end

		fileName1 = simulation("INIT",10000,"INIT_8hpi_"*name*lpad(i-offset,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-03-15_073352_INIT_NI_"*lpad(cases[i],2,"0"), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720, nuclearMechPars = nucMec)

		simulation("AFM", 5  , "AFM_8hpi_"*name*lpad(i-offset,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName1, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
	end

elseif sim == 4 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 5

	simPerCase = 12

	numInitSims = 64;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 1:numCases*simPerCase

		if any(i .== 1:12)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_change_01.txt"
		elseif any(i .== 13:24)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_change_02.txt"
		elseif any(i .== 25:36)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_change_03.txt"
		elseif any(i .== 37:48)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_change_04.txt"
		elseif any(i .== 49:60)
			nucMec = ".\\temp_pars\\nuclear_mechanics_lam_change_05.txt"
		end


		fileName1 = simulation("INIT",10000,"INIT_8hpi_lam_change_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-03-15_073352_INIT_NI_"*lpad(cases[i],2,"0"), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720, nuclearMechPars = nucMec)

		simulation("AFM", 5  , "AFM_8hpi_lam_change_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName1, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
	end

elseif sim == 5 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 4

	numInitSims = 12;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 1:numCases*simPerCase

		simulation("AFM", 5  , "AFM_NI_increased_spring_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = "2024-03-25_105820_INIT_NI_"*lpad(cases[i],2,"0"), simulationDate = simulationDate)

	end

elseif sim == 6

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",80,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720)

	# Threads.@threads 
	Threads.@threads for i = 51:54

		if any(i .== 1:36)
			sysPars = "./temp_pars/system_parameters_01.txt"
			if any(i .== 1:12)
				nucMec = "./temp_pars/nuclear_mechanics_dis_01.txt"
			elseif any(i .== 13:24)
				nucMec = "./temp_pars/nuclear_mechanics_dis_02.txt"
			else
				nucMec = "./temp_pars/nuclear_mechanics_dis_03.txt"
			end
		elseif any(i .== 37:72)
			sysPars = "./temp_pars/system_parameters_02.txt"
			if any(i .== 37:48)
				nucMec = "./temp_pars/nuclear_mechanics_dis_01.txt"
			elseif any(i .== 49:60)
				nucMec = "./temp_pars/nuclear_mechanics_dis_02.txt"
			else
				nucMec = "./temp_pars/nuclear_mechanics_dis_03.txt"
			end
		end

		fileName2 = simulation("INIT" ,5, "init_P1_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true,  importFolder = fileName1, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true, simulationDate = simulationDate, sysPars = sysPars)

		# create the crosslinks
		fileName3 = simulation("INIT" ,1000, "init_P2_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true,  noEnveSolve = true, importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true, simulationDate = simulationDate)

		# relax the whole system
		fileName4 = simulation("INIT", 20, "INIT_NI_LAD_TEST_"*lpad(i,2,"0"), "load"; adherent = true, adherentStatic = true, importFolder = fileName3, simulationDate = simulationDate,returnFoldername = true,)

		rm(".\\results\\"*fileName2; recursive = true)
		rm(".\\results\\"*fileName3; recursive = true)

		fileName5 = simulation("INIT",10000,"INIT_8hpi_lad_test_"*lpad(i,2,"0"),"load"; vrc = true, returnFoldername = true, importFolder = fileName4, adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720, nuclearMechPars = nucMec)

		simulation("AFM", 5  , "AFM_8hpi_lad_test_"*lpad(i,2,"0"), "load"; adherentStatic = true, stickyBottom = true, importFolder = fileName5, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)

	end



	rm(".\\results\\"*fileName1; recursive = true)

elseif sim == 7

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",200,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720)

	simulation("AFM", 5  , "AFM_NI_no_chromatin_test_high_lam", "load"; adherentStatic = true, stickyBottom = true, noChromatin = true,  importFolder = fileName1, simulationDate = simulationDate)


elseif sim == 8

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",200,"ADHERENT_SHELL_RADIUS TEST","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720)


elseif sim == 9 

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 2

	numInitSims = 4;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 1:numCases*simPerCase
		fileName1 = simulation("INIT",10000,"INIT_8hpi_rounder_test_low_stiffness_"*lpad(i,2,"0"),"load"; returnFoldername = true, importFolder = "2024-04-05_122524_INIT_NI_"*lpad(cases[i],2,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720)

		fileName2 = simulation("AFM", 5  , "AFM_8hpi_rounder_test_low_stiffness_"*lpad(i,2,"0"), "load"; adherentStatic = true, returnFoldername = true, stickyBottom = true, importFolder = fileName1, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt",nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate)
		
		rm(".\\results\\"*fileName1; recursive = true)

		zip("./results/"*fileName2*".zip","./results/"*fileName2)

		rm(".\\results\\"*fileName2; recursive = true)

	end
end