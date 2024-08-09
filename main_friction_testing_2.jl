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

sim = 20


#####################################################################################################
# initialize NI nuclei
#####################################################################################################

if sim == 1 

    simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",600,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 720, resetVertexDistancesTime = 400.)

	# Threads.@threads 
	Threads.@threads for i = 54:64


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

	simulationDate = "2024-04-11_105335_" #Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 243

	simPerCase = 4

	numInitSims = 32;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = inds

		if any(i .== 1:4)
			nucMec = "./temp_pars/nuclear_mechanics_01.txt"
		elseif any(i .== 5:8)
			nucMec = "./temp_pars/nuclear_mechanics_02.txt"
		elseif any(i .== 9:12)
			nucMec = "./temp_pars/nuclear_mechanics_03.txt"
		elseif any(i .== 13:16)
			nucMec = "./temp_pars/nuclear_mechanics_04.txt"
		elseif any(i .== 17:20)
			nucMec = "./temp_pars/nuclear_mechanics_05.txt"
		elseif any(i .== 21:24)
			nucMec = "./temp_pars/nuclear_mechanics_06.txt"
		elseif any(i .== 25:28)
			nucMec = "./temp_pars/nuclear_mechanics_07.txt"
		elseif any(i .== 29:32)
			nucMec = "./temp_pars/nuclear_mechanics_08.txt"
		elseif any(i .== 33:36)
			nucMec = "./temp_pars/nuclear_mechanics_09.txt"
		elseif any(i .== 37:40)
			nucMec = "./temp_pars/nuclear_mechanics_10.txt"
		elseif any(i .== 41:44)
			nucMec = "./temp_pars/nuclear_mechanics_11.txt"
		elseif any(i .== 45:48)
			nucMec = "./temp_pars/nuclear_mechanics_12.txt"
		elseif any(i .== 49:52)
			nucMec = "./temp_pars/nuclear_mechanics_13.txt"
		elseif any(i .== 53:56)
			nucMec = "./temp_pars/nuclear_mechanics_14.txt"
		elseif any(i .== 57:60)
			nucMec = "./temp_pars/nuclear_mechanics_15.txt"
		elseif any(i .== 61:64)
			nucMec = "./temp_pars/nuclear_mechanics_16.txt"
		elseif any(i .== 65:68)
			nucMec = "./temp_pars/nuclear_mechanics_17.txt"
		elseif any(i .== 69:72)
			nucMec = "./temp_pars/nuclear_mechanics_18.txt"
		elseif any(i .== 73:76)
			nucMec = "./temp_pars/nuclear_mechanics_19.txt"
		elseif any(i .== 77:80)
			nucMec = "./temp_pars/nuclear_mechanics_20.txt"
		elseif any(i .== 81:84)
			nucMec = "./temp_pars/nuclear_mechanics_21.txt"
		elseif any(i .== 85:88)
			nucMec = "./temp_pars/nuclear_mechanics_22.txt"
		elseif any(i .== 89:92)
			nucMec = "./temp_pars/nuclear_mechanics_23.txt"
		elseif any(i .== 93:96)
			nucMec = "./temp_pars/nuclear_mechanics_24.txt"
		elseif any(i .== 97:100)
			nucMec = "./temp_pars/nuclear_mechanics_25.txt"
		elseif any(i .== 101:104)
			nucMec = "./temp_pars/nuclear_mechanics_26.txt"
		elseif any(i .== 105:108)
			nucMec = "./temp_pars/nuclear_mechanics_27.txt"
		elseif any(i .== 109:112)
			nucMec = "./temp_pars/nuclear_mechanics_28.txt"
		elseif any(i .== 113:116)
			nucMec = "./temp_pars/nuclear_mechanics_29.txt"
		elseif any(i .== 117:120)
			nucMec = "./temp_pars/nuclear_mechanics_30.txt"
		elseif any(i .== 121:124)
			nucMec = "./temp_pars/nuclear_mechanics_31.txt"
		elseif any(i .== 125:128)
			nucMec = "./temp_pars/nuclear_mechanics_32.txt"
		elseif any(i .== 129:132)
			nucMec = "./temp_pars/nuclear_mechanics_33.txt"
		elseif any(i .== 133:136)
			nucMec = "./temp_pars/nuclear_mechanics_34.txt"
		elseif any(i .== 137:140)
			nucMec = "./temp_pars/nuclear_mechanics_35.txt"
		elseif any(i .== 141:144)
			nucMec = "./temp_pars/nuclear_mechanics_36.txt"
		elseif any(i .== 145:148)
			nucMec = "./temp_pars/nuclear_mechanics_37.txt"
		elseif any(i .== 149:152)
			nucMec = "./temp_pars/nuclear_mechanics_38.txt"
		elseif any(i .== 153:156)
			nucMec = "./temp_pars/nuclear_mechanics_39.txt"
		elseif any(i .== 157:160)
			nucMec = "./temp_pars/nuclear_mechanics_40.txt"
		elseif any(i .== 161:164)
			nucMec = "./temp_pars/nuclear_mechanics_41.txt"
		elseif any(i .== 165:168)
			nucMec = "./temp_pars/nuclear_mechanics_42.txt"
		elseif any(i .== 169:172)
			nucMec = "./temp_pars/nuclear_mechanics_43.txt"
		elseif any(i .== 173:176)
			nucMec = "./temp_pars/nuclear_mechanics_44.txt"
		elseif any(i .== 177:180)
			nucMec = "./temp_pars/nuclear_mechanics_45.txt"
		elseif any(i .== 181:184)
			nucMec = "./temp_pars/nuclear_mechanics_46.txt"
		elseif any(i .== 185:188)
			nucMec = "./temp_pars/nuclear_mechanics_47.txt"
		elseif any(i .== 189:192)
			nucMec = "./temp_pars/nuclear_mechanics_48.txt"
		elseif any(i .== 193:196)
			nucMec = "./temp_pars/nuclear_mechanics_49.txt"
		elseif any(i .== 197:200)
			nucMec = "./temp_pars/nuclear_mechanics_50.txt"
		elseif any(i .== 201:204)
			nucMec = "./temp_pars/nuclear_mechanics_51.txt"
		elseif any(i .== 205:208)
			nucMec = "./temp_pars/nuclear_mechanics_52.txt"
		elseif any(i .== 209:212)
			nucMec = "./temp_pars/nuclear_mechanics_53.txt"
		elseif any(i .== 213:216)
			nucMec = "./temp_pars/nuclear_mechanics_54.txt"
		elseif any(i .== 217:220)
			nucMec = "./temp_pars/nuclear_mechanics_55.txt"
		elseif any(i .== 221:224)
			nucMec = "./temp_pars/nuclear_mechanics_56.txt"
		elseif any(i .== 225:228)
			nucMec = "./temp_pars/nuclear_mechanics_57.txt"
		elseif any(i .== 229:232)
			nucMec = "./temp_pars/nuclear_mechanics_58.txt"
		elseif any(i .== 233:236)
			nucMec = "./temp_pars/nuclear_mechanics_59.txt"
		elseif any(i .== 237:240)
			nucMec = "./temp_pars/nuclear_mechanics_60.txt"
		elseif any(i .== 241:244)
			nucMec = "./temp_pars/nuclear_mechanics_61.txt"
		elseif any(i .== 245:248)
			nucMec = "./temp_pars/nuclear_mechanics_62.txt"
		elseif any(i .== 249:252)
			nucMec = "./temp_pars/nuclear_mechanics_63.txt"
		elseif any(i .== 253:256)
			nucMec = "./temp_pars/nuclear_mechanics_64.txt"
		elseif any(i .== 257:260)
			nucMec = "./temp_pars/nuclear_mechanics_65.txt"
		elseif any(i .== 261:264)
			nucMec = "./temp_pars/nuclear_mechanics_66.txt"
		elseif any(i .== 265:268)
			nucMec = "./temp_pars/nuclear_mechanics_67.txt"
		elseif any(i .== 269:272)
			nucMec = "./temp_pars/nuclear_mechanics_68.txt"
		elseif any(i .== 273:276)
			nucMec = "./temp_pars/nuclear_mechanics_69.txt"
		elseif any(i .== 277:280)
			nucMec = "./temp_pars/nuclear_mechanics_70.txt"
		elseif any(i .== 281:284)
			nucMec = "./temp_pars/nuclear_mechanics_71.txt"
		elseif any(i .== 285:288)
			nucMec = "./temp_pars/nuclear_mechanics_72.txt"
		elseif any(i .== 289:292)
			nucMec = "./temp_pars/nuclear_mechanics_73.txt"
		elseif any(i .== 293:296)
			nucMec = "./temp_pars/nuclear_mechanics_74.txt"
		elseif any(i .== 297:300)
			nucMec = "./temp_pars/nuclear_mechanics_75.txt"
		elseif any(i .== 301:304)
			nucMec = "./temp_pars/nuclear_mechanics_76.txt"
		elseif any(i .== 305:308)
			nucMec = "./temp_pars/nuclear_mechanics_77.txt"
		elseif any(i .== 309:312)
			nucMec = "./temp_pars/nuclear_mechanics_78.txt"
		elseif any(i .== 313:316)
			nucMec = "./temp_pars/nuclear_mechanics_79.txt"
		elseif any(i .== 317:320)
			nucMec = "./temp_pars/nuclear_mechanics_80.txt"
		elseif any(i .== 321:324)
			nucMec = "./temp_pars/nuclear_mechanics_81.txt"
		elseif any(i .== 325:328)
			nucMec = "./temp_pars/nuclear_mechanics_82.txt"
		elseif any(i .== 329:332)
			nucMec = "./temp_pars/nuclear_mechanics_83.txt"
		elseif any(i .== 333:336)
			nucMec = "./temp_pars/nuclear_mechanics_84.txt"
		elseif any(i .== 337:340)
			nucMec = "./temp_pars/nuclear_mechanics_85.txt"
		elseif any(i .== 341:344)
			nucMec = "./temp_pars/nuclear_mechanics_86.txt"
		elseif any(i .== 345:348)
			nucMec = "./temp_pars/nuclear_mechanics_87.txt"
		elseif any(i .== 349:352)
			nucMec = "./temp_pars/nuclear_mechanics_88.txt"
		elseif any(i .== 353:356)
			nucMec = "./temp_pars/nuclear_mechanics_89.txt"
		elseif any(i .== 357:360)
			nucMec = "./temp_pars/nuclear_mechanics_90.txt"
		elseif any(i .== 361:364)
			nucMec = "./temp_pars/nuclear_mechanics_91.txt"
		elseif any(i .== 365:368)
			nucMec = "./temp_pars/nuclear_mechanics_92.txt"
		elseif any(i .== 369:372)
			nucMec = "./temp_pars/nuclear_mechanics_93.txt"
		elseif any(i .== 373:376)
			nucMec = "./temp_pars/nuclear_mechanics_94.txt"
		elseif any(i .== 377:380)
			nucMec = "./temp_pars/nuclear_mechanics_95.txt"
		elseif any(i .== 381:384)
			nucMec = "./temp_pars/nuclear_mechanics_96.txt"
		elseif any(i .== 385:388)
			nucMec = "./temp_pars/nuclear_mechanics_97.txt"
		elseif any(i .== 389:392)
			nucMec = "./temp_pars/nuclear_mechanics_98.txt"
		elseif any(i .== 393:396)
			nucMec = "./temp_pars/nuclear_mechanics_99.txt"
		elseif any(i .== 397:400)
			nucMec = "./temp_pars/nuclear_mechanics_100.txt"
		elseif any(i .== 401:404)
			nucMec = "./temp_pars/nuclear_mechanics_101.txt"
		elseif any(i .== 405:408)
			nucMec = "./temp_pars/nuclear_mechanics_102.txt"
		elseif any(i .== 409:412)
			nucMec = "./temp_pars/nuclear_mechanics_103.txt"
		elseif any(i .== 413:416)
			nucMec = "./temp_pars/nuclear_mechanics_104.txt"
		elseif any(i .== 417:420)
			nucMec = "./temp_pars/nuclear_mechanics_105.txt"
		elseif any(i .== 421:424)
			nucMec = "./temp_pars/nuclear_mechanics_106.txt"
		elseif any(i .== 425:428)
			nucMec = "./temp_pars/nuclear_mechanics_107.txt"
		elseif any(i .== 429:432)
			nucMec = "./temp_pars/nuclear_mechanics_108.txt"
		elseif any(i .== 433:436)
			nucMec = "./temp_pars/nuclear_mechanics_109.txt"
		elseif any(i .== 437:440)
			nucMec = "./temp_pars/nuclear_mechanics_110.txt"
		elseif any(i .== 441:444)
			nucMec = "./temp_pars/nuclear_mechanics_111.txt"
		elseif any(i .== 445:448)
			nucMec = "./temp_pars/nuclear_mechanics_112.txt"
		elseif any(i .== 449:452)
			nucMec = "./temp_pars/nuclear_mechanics_113.txt"
		elseif any(i .== 453:456)
			nucMec = "./temp_pars/nuclear_mechanics_114.txt"
		elseif any(i .== 457:460)
			nucMec = "./temp_pars/nuclear_mechanics_115.txt"
		elseif any(i .== 461:464)
			nucMec = "./temp_pars/nuclear_mechanics_116.txt"
		elseif any(i .== 465:468)
			nucMec = "./temp_pars/nuclear_mechanics_117.txt"
		elseif any(i .== 469:472)
			nucMec = "./temp_pars/nuclear_mechanics_118.txt"
		elseif any(i .== 473:476)
			nucMec = "./temp_pars/nuclear_mechanics_119.txt"
		elseif any(i .== 477:480)
			nucMec = "./temp_pars/nuclear_mechanics_120.txt"
		elseif any(i .== 481:484)
			nucMec = "./temp_pars/nuclear_mechanics_121.txt"
		elseif any(i .== 485:488)
			nucMec = "./temp_pars/nuclear_mechanics_122.txt"
		elseif any(i .== 489:492)
			nucMec = "./temp_pars/nuclear_mechanics_123.txt"
		elseif any(i .== 493:496)
			nucMec = "./temp_pars/nuclear_mechanics_124.txt"
		elseif any(i .== 497:500)
			nucMec = "./temp_pars/nuclear_mechanics_125.txt"
		elseif any(i .== 501:504)
			nucMec = "./temp_pars/nuclear_mechanics_126.txt"
		elseif any(i .== 505:508)
			nucMec = "./temp_pars/nuclear_mechanics_127.txt"
		elseif any(i .== 509:512)
			nucMec = "./temp_pars/nuclear_mechanics_128.txt"
		elseif any(i .== 513:516)
			nucMec = "./temp_pars/nuclear_mechanics_129.txt"
		elseif any(i .== 517:520)
			nucMec = "./temp_pars/nuclear_mechanics_130.txt"
		elseif any(i .== 521:524)
			nucMec = "./temp_pars/nuclear_mechanics_131.txt"
		elseif any(i .== 525:528)
			nucMec = "./temp_pars/nuclear_mechanics_132.txt"
		elseif any(i .== 529:532)
			nucMec = "./temp_pars/nuclear_mechanics_133.txt"
		elseif any(i .== 533:536)
			nucMec = "./temp_pars/nuclear_mechanics_134.txt"
		elseif any(i .== 537:540)
			nucMec = "./temp_pars/nuclear_mechanics_135.txt"
		elseif any(i .== 541:544)
			nucMec = "./temp_pars/nuclear_mechanics_136.txt"
		elseif any(i .== 545:548)
			nucMec = "./temp_pars/nuclear_mechanics_137.txt"
		elseif any(i .== 549:552)
			nucMec = "./temp_pars/nuclear_mechanics_138.txt"
		elseif any(i .== 553:556)
			nucMec = "./temp_pars/nuclear_mechanics_139.txt"
		elseif any(i .== 557:560)
			nucMec = "./temp_pars/nuclear_mechanics_140.txt"
		elseif any(i .== 561:564)
			nucMec = "./temp_pars/nuclear_mechanics_141.txt"
		elseif any(i .== 565:568)
			nucMec = "./temp_pars/nuclear_mechanics_142.txt"
		elseif any(i .== 569:572)
			nucMec = "./temp_pars/nuclear_mechanics_143.txt"
		elseif any(i .== 573:576)
			nucMec = "./temp_pars/nuclear_mechanics_144.txt"
		elseif any(i .== 577:580)
			nucMec = "./temp_pars/nuclear_mechanics_145.txt"
		elseif any(i .== 581:584)
			nucMec = "./temp_pars/nuclear_mechanics_146.txt"
		elseif any(i .== 585:588)
			nucMec = "./temp_pars/nuclear_mechanics_147.txt"
		elseif any(i .== 589:592)
			nucMec = "./temp_pars/nuclear_mechanics_148.txt"
		elseif any(i .== 593:596)
			nucMec = "./temp_pars/nuclear_mechanics_149.txt"
		elseif any(i .== 597:600)
			nucMec = "./temp_pars/nuclear_mechanics_150.txt"
		elseif any(i .== 601:604)
			nucMec = "./temp_pars/nuclear_mechanics_151.txt"
		elseif any(i .== 605:608)
			nucMec = "./temp_pars/nuclear_mechanics_152.txt"
		elseif any(i .== 609:612)
			nucMec = "./temp_pars/nuclear_mechanics_153.txt"
		elseif any(i .== 613:616)
			nucMec = "./temp_pars/nuclear_mechanics_154.txt"
		elseif any(i .== 617:620)
			nucMec = "./temp_pars/nuclear_mechanics_155.txt"
		elseif any(i .== 621:624)
			nucMec = "./temp_pars/nuclear_mechanics_156.txt"
		elseif any(i .== 625:628)
			nucMec = "./temp_pars/nuclear_mechanics_157.txt"
		elseif any(i .== 629:632)
			nucMec = "./temp_pars/nuclear_mechanics_158.txt"
		elseif any(i .== 633:636)
			nucMec = "./temp_pars/nuclear_mechanics_159.txt"
		elseif any(i .== 637:640)
			nucMec = "./temp_pars/nuclear_mechanics_160.txt"
		elseif any(i .== 641:644)
			nucMec = "./temp_pars/nuclear_mechanics_161.txt"
		elseif any(i .== 645:648)
			nucMec = "./temp_pars/nuclear_mechanics_162.txt"
		elseif any(i .== 649:652)
			nucMec = "./temp_pars/nuclear_mechanics_163.txt"
		elseif any(i .== 653:656)
			nucMec = "./temp_pars/nuclear_mechanics_164.txt"
		elseif any(i .== 657:660)
			nucMec = "./temp_pars/nuclear_mechanics_165.txt"
		elseif any(i .== 661:664)
			nucMec = "./temp_pars/nuclear_mechanics_166.txt"
		elseif any(i .== 665:668)
			nucMec = "./temp_pars/nuclear_mechanics_167.txt"
		elseif any(i .== 669:672)
			nucMec = "./temp_pars/nuclear_mechanics_168.txt"
		elseif any(i .== 673:676)
			nucMec = "./temp_pars/nuclear_mechanics_169.txt"
		elseif any(i .== 677:680)
			nucMec = "./temp_pars/nuclear_mechanics_170.txt"
		elseif any(i .== 681:684)
			nucMec = "./temp_pars/nuclear_mechanics_171.txt"
		elseif any(i .== 685:688)
			nucMec = "./temp_pars/nuclear_mechanics_172.txt"
		elseif any(i .== 689:692)
			nucMec = "./temp_pars/nuclear_mechanics_173.txt"
		elseif any(i .== 693:696)
			nucMec = "./temp_pars/nuclear_mechanics_174.txt"
		elseif any(i .== 697:700)
			nucMec = "./temp_pars/nuclear_mechanics_175.txt"
		elseif any(i .== 701:704)
			nucMec = "./temp_pars/nuclear_mechanics_176.txt"
		elseif any(i .== 705:708)
			nucMec = "./temp_pars/nuclear_mechanics_177.txt"
		elseif any(i .== 709:712)
			nucMec = "./temp_pars/nuclear_mechanics_178.txt"
		elseif any(i .== 713:716)
			nucMec = "./temp_pars/nuclear_mechanics_179.txt"
		elseif any(i .== 717:720)
			nucMec = "./temp_pars/nuclear_mechanics_180.txt"
		elseif any(i .== 721:724)
			nucMec = "./temp_pars/nuclear_mechanics_181.txt"
		elseif any(i .== 725:728)
			nucMec = "./temp_pars/nuclear_mechanics_182.txt"
		elseif any(i .== 729:732)
			nucMec = "./temp_pars/nuclear_mechanics_183.txt"
		elseif any(i .== 733:736)
			nucMec = "./temp_pars/nuclear_mechanics_184.txt"
		elseif any(i .== 737:740)
			nucMec = "./temp_pars/nuclear_mechanics_185.txt"
		elseif any(i .== 741:744)
			nucMec = "./temp_pars/nuclear_mechanics_186.txt"
		elseif any(i .== 745:748)
			nucMec = "./temp_pars/nuclear_mechanics_187.txt"
		elseif any(i .== 749:752)
			nucMec = "./temp_pars/nuclear_mechanics_188.txt"
		elseif any(i .== 753:756)
			nucMec = "./temp_pars/nuclear_mechanics_189.txt"
		elseif any(i .== 757:760)
			nucMec = "./temp_pars/nuclear_mechanics_190.txt"
		elseif any(i .== 761:764)
			nucMec = "./temp_pars/nuclear_mechanics_191.txt"
		elseif any(i .== 765:768)
			nucMec = "./temp_pars/nuclear_mechanics_192.txt"
		elseif any(i .== 769:772)
			nucMec = "./temp_pars/nuclear_mechanics_193.txt"
		elseif any(i .== 773:776)
			nucMec = "./temp_pars/nuclear_mechanics_194.txt"
		elseif any(i .== 777:780)
			nucMec = "./temp_pars/nuclear_mechanics_195.txt"
		elseif any(i .== 781:784)
			nucMec = "./temp_pars/nuclear_mechanics_196.txt"
		elseif any(i .== 785:788)
			nucMec = "./temp_pars/nuclear_mechanics_197.txt"
		elseif any(i .== 789:792)
			nucMec = "./temp_pars/nuclear_mechanics_198.txt"
		elseif any(i .== 793:796)
			nucMec = "./temp_pars/nuclear_mechanics_199.txt"
		elseif any(i .== 797:800)
			nucMec = "./temp_pars/nuclear_mechanics_200.txt"
		elseif any(i .== 801:804)
			nucMec = "./temp_pars/nuclear_mechanics_201.txt"
		elseif any(i .== 805:808)
			nucMec = "./temp_pars/nuclear_mechanics_202.txt"
		elseif any(i .== 809:812)
			nucMec = "./temp_pars/nuclear_mechanics_203.txt"
		elseif any(i .== 813:816)
			nucMec = "./temp_pars/nuclear_mechanics_204.txt"
		elseif any(i .== 817:820)
			nucMec = "./temp_pars/nuclear_mechanics_205.txt"
		elseif any(i .== 821:824)
			nucMec = "./temp_pars/nuclear_mechanics_206.txt"
		elseif any(i .== 825:828)
			nucMec = "./temp_pars/nuclear_mechanics_207.txt"
		elseif any(i .== 829:832)
			nucMec = "./temp_pars/nuclear_mechanics_208.txt"
		elseif any(i .== 833:836)
			nucMec = "./temp_pars/nuclear_mechanics_209.txt"
		elseif any(i .== 837:840)
			nucMec = "./temp_pars/nuclear_mechanics_210.txt"
		elseif any(i .== 841:844)
			nucMec = "./temp_pars/nuclear_mechanics_211.txt"
		elseif any(i .== 845:848)
			nucMec = "./temp_pars/nuclear_mechanics_212.txt"
		elseif any(i .== 849:852)
			nucMec = "./temp_pars/nuclear_mechanics_213.txt"
		elseif any(i .== 853:856)
			nucMec = "./temp_pars/nuclear_mechanics_214.txt"
		elseif any(i .== 857:860)
			nucMec = "./temp_pars/nuclear_mechanics_215.txt"
		elseif any(i .== 861:864)
			nucMec = "./temp_pars/nuclear_mechanics_216.txt"
		elseif any(i .== 865:868)
			nucMec = "./temp_pars/nuclear_mechanics_217.txt"
		elseif any(i .== 869:872)
			nucMec = "./temp_pars/nuclear_mechanics_218.txt"
		elseif any(i .== 873:876)
			nucMec = "./temp_pars/nuclear_mechanics_219.txt"
		elseif any(i .== 877:880)
			nucMec = "./temp_pars/nuclear_mechanics_220.txt"
		elseif any(i .== 881:884)
			nucMec = "./temp_pars/nuclear_mechanics_221.txt"
		elseif any(i .== 885:888)
			nucMec = "./temp_pars/nuclear_mechanics_222.txt"
		elseif any(i .== 889:892)
			nucMec = "./temp_pars/nuclear_mechanics_223.txt"
		elseif any(i .== 893:896)
			nucMec = "./temp_pars/nuclear_mechanics_224.txt"
		elseif any(i .== 897:900)
			nucMec = "./temp_pars/nuclear_mechanics_225.txt"
		elseif any(i .== 901:904)
			nucMec = "./temp_pars/nuclear_mechanics_226.txt"
		elseif any(i .== 905:908)
			nucMec = "./temp_pars/nuclear_mechanics_227.txt"
		elseif any(i .== 909:912)
			nucMec = "./temp_pars/nuclear_mechanics_228.txt"
		elseif any(i .== 913:916)
			nucMec = "./temp_pars/nuclear_mechanics_229.txt"
		elseif any(i .== 917:920)
			nucMec = "./temp_pars/nuclear_mechanics_230.txt"
		elseif any(i .== 921:924)
			nucMec = "./temp_pars/nuclear_mechanics_231.txt"
		elseif any(i .== 925:928)
			nucMec = "./temp_pars/nuclear_mechanics_232.txt"
		elseif any(i .== 929:932)
			nucMec = "./temp_pars/nuclear_mechanics_233.txt"
		elseif any(i .== 933:936)
			nucMec = "./temp_pars/nuclear_mechanics_234.txt"
		elseif any(i .== 937:940)
			nucMec = "./temp_pars/nuclear_mechanics_235.txt"
		elseif any(i .== 941:944)
			nucMec = "./temp_pars/nuclear_mechanics_236.txt"
		elseif any(i .== 945:948)
			nucMec = "./temp_pars/nuclear_mechanics_237.txt"
		elseif any(i .== 949:952)
			nucMec = "./temp_pars/nuclear_mechanics_238.txt"
		elseif any(i .== 953:956)
			nucMec = "./temp_pars/nuclear_mechanics_239.txt"
		elseif any(i .== 957:960)
			nucMec = "./temp_pars/nuclear_mechanics_240.txt"
		elseif any(i .== 961:964)
			nucMec = "./temp_pars/nuclear_mechanics_241.txt"
		elseif any(i .== 965:968)
			nucMec = "./temp_pars/nuclear_mechanics_242.txt"
		elseif any(i .== 969:972)
			nucMec = "./temp_pars/nuclear_mechanics_243.txt"
		end

		fileName1 = simulation("INIT",10000,"INIT_FITTING_"*lpad(i,3,"0"),"load"; returnFoldername = true, importFolder = "2024-04-08_125833_INIT_NI_"*lpad(cases[i],2,"0"), adherent = true, simulationDate = simulationDate, newTargetVolume = 720, nuclearMechPars = nucMec)

		fileName2 = simulation("AFM", 5  , "AFM_FITTING_"*lpad(i,3,"0"), "load"; returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName1, simulationDate = simulationDate, nuclearMechPars = nucMec)
	
		rm(".\\results\\"*fileName1; recursive = true)

		zip("./results/"*fileName2)

		rm(".\\results\\"*fileName2; recursive = true)
	end

elseif sim == 3 # initialize 8 hpi nucleus

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 24

	numInitSims = 64;

	local cases = []

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads 
	Threads.@threads for i = 1:numCases*simPerCase

		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName1 = simulation("INIT",10000,"INIT_NI_"*lpad(i,3,"0"),"load"; returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"), adherent = true, simulationDate = simulationDate, newTargetVolume = 720)

		rm(".\\results\\2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		fileName2 = simulation("AFM", 5, "AFM_NI_"*lpad(i,3,"0"), "load"; returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName1, simulationDate = simulationDate)
	
		rm(".\\results\\"*fileName1; recursive = true)

		zip("./results/"*fileName2)

		rm(".\\results\\"*fileName2; recursive = true)
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

	numCases = 9

	simPerCase = 12

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase

		if any(i .== 1:12)
			nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		elseif any(i .== 13:24)
			nucMec = "./temp_pars/nuclear_mechanics_002.txt"
		elseif any(i .== 25:36)
			nucMec = "./temp_pars/nuclear_mechanics_003.txt"
		elseif any(i .== 37:48)
			nucMec = "./temp_pars/nuclear_mechanics_004.txt"
		elseif any(i .== 49:60)
			nucMec = "./temp_pars/nuclear_mechanics_005.txt"
		elseif any(i .== 61:72)
			nucMec = "./temp_pars/nuclear_mechanics_006.txt"
		elseif any(i .== 73:84)
			nucMec = "./temp_pars/nuclear_mechanics_007.txt"
		elseif any(i .== 85:96)
			nucMec = "./temp_pars/nuclear_mechanics_008.txt"
		elseif any(i .== 97:108)
			nucMec = "./temp_pars/nuclear_mechanics_009.txt"
		end

		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_8hpi_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", newTargetVolume = 720, simulationDate = simulationDate, nuclearMechPars = nucMec)

		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		fileName6 = simulation("AFM", 5, "AFM_8hpi_2xlam_outward_changes_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)



	end
elseif sim == 10

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 4

	simPerCase = 12

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 110-72:110-72#1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		# if any(i .== 1:12)
		# 	lamDis = 0.6
		# elseif any(i .== 13:24)
		# 	lamDis = 0.7
		# elseif any(i .== 25:36)
		# 	lamDis = 0.8
		# elseif any(i .== 37:48)
		# 	lamDis = 0.9
		# # elseif any(i .== 49:60)
		# # 	lamDis = 0.40
		# # elseif any(i .== 61:72)
		# # 	lamDis = 0.50
		# end

		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_8hpi_"*lpad(i+72,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", newTargetVolume = 720, simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis)

		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		fileName6 = simulation("AFM", 5, "AFM_8hpi_2xlam_lamina_disintegration_"*lpad(i+72,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)



	end

elseif sim == 11

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 12

	local cases = []

	numInitSims = 12

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		
		unzip("./results/2024-05-07_131615_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_12hpi_"*lpad(i+72,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-05-07_131615_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", newTargetVolume = 890, simulationDate = simulationDate)

		rm("./results/2024-05-07_131615_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		fileName6 = simulation("AFM", 5, "AFM_8hpi_2xlam_lamina_disintegration_"*lpad(i+72,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", simulationDate = simulationDate)
		
		zip("./results/"*fileName5)
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)



	end

elseif sim == 12

    simulationDate = "2024-05-10_091503_" #Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	fileName1 = simulation("INIT",600,"ADHERENT_SHELL","new"; adherent = true, returnFoldername = true, simulationDate = simulationDate, noChromatin = true, newTargetVolume = 890, resetVertexDistancesTime = 400., nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt")

	# Threads.@threads 
	Threads.@threads for i = inds


		fileName2 = simulation("INIT" ,5, "init_P1_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true,  importFolder = fileName1, noEnveSolve = true, simPars = "./parameters/simulation_parameters_init_1.txt",
		returnFoldername = true, simulationDate = simulationDate, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt")

		# create the crosslinks
		fileName3 = simulation("INIT" ,1000, "init_P2_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true,  noEnveSolve = true, importFolder = fileName2, simPars = "./parameters/simulation_parameters_init_2.txt",
					returnFoldername = true, simulationDate = simulationDate, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt")

		# relax the whole system
		fileName4 = simulation("INIT", 20, "INIT_NI_12hpi_"*lpad(i,3,"0"), "load"; adherent = true, adherentStatic = true, importFolder = fileName3, simulationDate = simulationDate,returnFoldername = true,)

		rm(".\\results\\"*fileName2; recursive = true)
		rm(".\\results\\"*fileName3; recursive = true)

		zip("./results/"*fileName4)

		rm(".\\results\\"*fileName4; recursive = true)

    end

	rm(".\\results\\"*fileName1; recursive = true)


elseif sim == 13

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 9

	simPerCase = 12

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase

		if any(i .== 1:12)
			nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		elseif any(i .== 13:24)
			nucMec = "./temp_pars/nuclear_mechanics_002.txt"
		elseif any(i .== 25:36)
			nucMec = "./temp_pars/nuclear_mechanics_003.txt"
		elseif any(i .== 37:48)
			nucMec = "./temp_pars/nuclear_mechanics_004.txt"
		elseif any(i .== 49:60)
			nucMec = "./temp_pars/nuclear_mechanics_005.txt"
		elseif any(i .== 61:72)
			nucMec = "./temp_pars/nuclear_mechanics_006.txt"
		elseif any(i .== 73:84)
			nucMec = "./temp_pars/nuclear_mechanics_007.txt"
		elseif any(i .== 85:96)
			nucMec = "./temp_pars/nuclear_mechanics_008.txt"
		elseif any(i .== 97:108)
			nucMec = "./temp_pars/nuclear_mechanics_009.txt"
		end

		unzip("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_12hpi_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", newTargetVolume = 890, simulationDate = simulationDate, nuclearMechPars = nucMec)

		rm("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0"); recursive = true)

		fileName6 = simulation("AFM", 5, "AFM_12hpi_2xlam_outward_changes_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)



	end

elseif sim == 14

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 10

	simPerCase = 12

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		if any(i .== 1:12)
			lamDis = 0.0
		elseif any(i .== 13:24)
			lamDis = 0.1
		elseif any(i .== 25:36)
			lamDis = 0.2
		elseif any(i .== 37:48)
			lamDis = 0.3
		elseif any(i .== 49:60)
			lamDis = 0.4
		elseif any(i .== 61:72)
			lamDis = 0.5
		elseif any(i .== 73:84)
			lamDis = 0.6
		elseif any(i .== 85:96)
			lamDis = 0.7
		elseif any(i .== 97:108)
			lamDis = 0.8
		elseif any(i .== 109:120)
			lamDis = 0.9
		end

		unzip("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_12hpi_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", newTargetVolume = 890, simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis)

		fileName6 = simulation("AFM", 5, "AFM_12hpi_2xlam_lamina_disintegration_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)


	end

elseif sim == 15

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 10

	simPerCase = 3

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		if any(i .== 1:12)
			lamDis = 0.0
		elseif any(i .== 13:24)
			lamDis = 0.1
		elseif any(i .== 25:36)
			lamDis = 0.2
		elseif any(i .== 37:48)
			lamDis = 0.3
		elseif any(i .== 49:60)
			lamDis = 0.4
		elseif any(i .== 61:72)
			lamDis = 0.5
		elseif any(i .== 73:84)
			lamDis = 0.6
		elseif any(i .== 85:96)
			lamDis = 0.7
		elseif any(i .== 97:108)
			lamDis = 0.8
		elseif any(i .== 109:120)
			lamDis = 0.9
		end

		unzip("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"lamina_disintegration_nuclear_size_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis)

		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)

	end

elseif sim == 16

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 8

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_8hpi_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", newTargetVolume = 720, simulationDate = simulationDate, nuclearMechPars = nucMec)

		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		fileName6 = simulation("AFM", 5, "AFM_INF_enve_viscocity_2.0_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
	end

elseif sim == 17

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 8

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName6 = simulation("AFM", 5, "AFM_NI_chro_viscocity_10.0_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)
		
		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName6; recursive = true)
	end

elseif sim == 18

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 12

	simPerCase = 8

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		nucMec = "./temp_pars/nuclear_mechanics_001.txt"
		
		if any(i .== 1:24)
			nucMec = "./temp_pars/nuclear_mechanics_001.txt"
			if any(i .== 1:8)
				lamDis = 0.0
			elseif any(i .== 9:16)
				lamDis = 0.2
			else
				lamDis = 0.5
			end
		elseif any(i .== 25:48)
			nucMec = "./temp_pars/nuclear_mechanics_002.txt"
			if any(i .== 25:32)
				lamDis = 0.0
			elseif any(i .== 33:40)
				lamDis = 0.2
			else
				lamDis = 0.5
			end
		elseif any(i .== 49:72)
			nucMec = "./temp_pars/nuclear_mechanics_003.txt"
			if any(i .== 49:56)
				lamDis = 0.0
			elseif any(i .== 57:64)
				lamDis = 0.2
			else
				lamDis = 0.5
			end
		elseif any(i .== 73:96)
			nucMec = "./temp_pars/nuclear_mechanics_004.txt"
			if any(i .== 73:80)
				lamDis = 0.0
			elseif any(i .== 81:88)
				lamDis = 0.2
			else
				lamDis = 0.5
			end
		end

		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"lamina_disintegration_nuclear_size_"*lpad(i+192,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec, laminaDisintegration = lamDis, newTargetVolume = 720)

		fileName6 = simulation("AFM", 5, "AFM_INF_multi_sensitivity_"*lpad(i+192,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, nuclearMechPars = nucMec)

		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)

	end

elseif sim == 19 # NI AFM, size test

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 24

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		
		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName6 = simulation("AFM", 5, "AFM_NI_size_tests_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"), simulationDate = simulationDate)
		
		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName6; recursive = true)
	end

elseif sim == 20 # 8hpi AFM, size test

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 24

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		
		unzip("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_8hpi_no_LADS_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt", simulationDate = simulationDate, newTargetVolume = 720)

		fileName6 = simulation("AFM", 5, "AFM_8hpi_no_LADS_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, simulationDate = simulationDate, replPars = "./parameters/replication_compartment_mechanics_8hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_8hpi.txt")
		
		rm("./results/2024-04-19_090545_INIT_NI_"*lpad(cases[i],3,"0"); recursive = true)

		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)
	end

elseif sim == 21 # 12hpi AFM, size test

	simulationDate = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS_")

	numCases = 1

	simPerCase = 24

	local cases = []

	numInitSims = 64;

	for i = 1:numCases
		cases = vcat(cases, randperm(numInitSims)[1:simPerCase])
	end

	# Threads.@threads
	Threads.@threads for i = 1:numCases*simPerCase
		
		unzip("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0")*".zip")

		fileName5 = simulation("INIT",10000,"INIT_8hpi_size_tests_"*lpad(i,3,"0"),"load"; vrc = true, returnFoldername = true, importFolder = "2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,'0'), adherent = true, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt", simulationDate = simulationDate, newTargetVolume = 890)

		fileName6 = simulation("AFM", 5, "AFM_12hpi_size_tests_lam_x2_"*lpad(i,3,"0"), "load";  returnFoldername = true, adherentStatic = true, stickyBottom = true, importFolder = fileName5, simulationDate = simulationDate, replPars = "./parameters/replication_compartment_mechanics_12hpi.txt", nuclearPropPars = "./parameters/nuclear_properties_12hpi.txt")
		
		rm("./results/2024-05-10_091503_INIT_NI_12hpi_"*lpad(cases[i],3,"0"); recursive = true)
		
		zip("./results/"*fileName6)

		rm(".\\results\\"*fileName5; recursive = true)
		rm(".\\results\\"*fileName6; recursive = true)
	end

end

