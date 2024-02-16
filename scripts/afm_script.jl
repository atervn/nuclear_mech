
include("..\\analysis_functions.jl")

for i = 1:12

    f,d = analyze_afm(;folder = ".\\results\\2024-02-13_093723_NI_AFM_"*lpad(i,2,'0'))

    writedlm("NI_AFM_"*lpad(i,2,'0')*".csv",[d[2:end] f[2:end]])

end