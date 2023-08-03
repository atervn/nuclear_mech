using DelimitedFiles, Plots, Statistics
theme(:ggplot2)

include("setup_functions.jl")

if !(@isdefined envelopeType)
    include("NuclearMechTypes.jl")
    using .NuclearMechTypes
end

function analyze_strain_2()

    folder = pick_folder(pwd()*"\\results")

    if folder == ""
        return
    end

    ipar = inputParametersType();
    ipar = read_parameters(ipar,folder*"\\parameters.txt")


    F = readdlm(folder*"\\nuclearForces.csv");
    L = readdlm(folder*"\\nuclearLength.csv");

    F1 = mean(F[200:400])*10e-6*1e-6/60;
    F2 = mean(F[600:800])*10e-6*1e-6/60;
    F3 = mean(F[1000:1200])*10e-6*1e-6/60;
    F4 = mean(F[1400:1600])*10e-6*1e-6/60;
    F5 = mean(F[1800:2000])*10e-6*1e-6/60;

    p = plot((L .- L[1])./L[1],F.*ipar.viscosity.*ipar.scalingLength./ipar.scalingTime,legend=false,xaxis = "Strain", yaxis = "Force (N)", dpi=600, xlims = (0,0.8), ylims = (0,5.5e-9))

    annotate!(p,[(0.28,F1-3e-10, (string(round(F1*1e9;digits=3)), 8, :black, :center))])
    annotate!(p,[(0.45,F2-3e-10, (string(round(F2*1e9;digits=3)), 8, :black, :center))])
    annotate!(p,[(0.58,F3-3e-10, (string(round(F3*1e9;digits=3)), 8, :black, :center))])
    annotate!(p,[(0.68,F4-3e-10, (string(round(F4*1e9;digits=3)), 8, :black, :center))])
    annotate!(p,[(0.76,F5-3e-10, (string(round(F5*1e9;digits=3)), 8, :black, :center))])
    scatter!(p,[0.28],[F1]; color = :red)
    scatter!(p,[0.45],[F2]; color = :red)
    scatter!(p,[0.58],[F3]; color = :red)
    scatter!(p,[0.68],[F4]; color = :red)
    scatter!(p,[0.76],[F5]; color = :red)

    annotate!(p,[(0.01,5.3e-9, (basename(folder), 9, :black, :left))])
    annotate!(p,[(0.01,5.0e-9, ("laminaStiffness: "*string(round(ipar.laminaStiffness*1e4,digits = 1))*"e-4", 8, :black, :left))])
    annotate!(p,[(0.01,4.75e-9, ("laminaFriction: "*string(ipar.laminaFriction*1e5)*"e-5", 8, :black, :left))])
    annotate!(p,[(0.01,4.5e-9, ("areaCompressionModulus: "*string(ipar.areaCompressionModulus*1e-3)*"e3", 8, :black, :left))])
    annotate!(p,[(0.01,4.25e-9, ("bulkModulus: "*string(ipar.bulkModulus), 8, :black, :left))])

    return p

end