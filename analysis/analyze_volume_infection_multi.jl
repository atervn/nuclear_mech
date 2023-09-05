function analyze_volume_infection_multi(nSims::Int)

    v = Vector{Vector{Float64}}(undef,nSims)

    for i = 1:nSims
        v[i],p = analyze_volume_infection()
    end

    means = mean(v)
    stds = std(v)

    return means,stds

end