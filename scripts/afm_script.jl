
include("..\\analysis_functions.jl")

for i = 1:48

    f,d = analyze_afm(;folder = ".\\results\\2024-02-26_075943_AFM_8hpi_lamina_disintegration_"*lpad(i,2,'0'))

    writedlm("AFM_8hpi_lamina_disintegration_"*lpad(i,2,'0')*".csv",[d[2:end] f[2:end]])

end