
for i = 1:12

    f,d = analyze_afm(;folder = ".\\results\\2024-01-08_234232_INF_12hpi_AFM_"*lpad(i,2,'0'))

    writedlm("inf_12hpi_afm_"*lpad(i,2,'0')*".csv",[d[2:end] f[2:end]])

end