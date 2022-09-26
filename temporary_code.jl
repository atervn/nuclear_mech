# VARIABLE TIME STEPPING

enveVelocities = Vector{Vec{3,Float64}}(undef,length(nuc.vert))
for i = 1:length(nuc.vert)
    enveVelocities[i] = Vec(vX[i],vY[i],vZ[i])

end
chroVelocities = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)
for k = 1:spar.chromatinLength*spar.chromatinNumber
    chroVelocities[k] = Vec(vX[length(nuc.vert)+k],vY[length(nuc.vert)+k],vZ[length(nuc.vert)+k])
end

enveMovements = Vector{Vec{3,Float64}}(undef,length(nuc.vert))
chroMovements = Vector{Vec{3,Float64}}(undef,spar.chromatinLength*spar.chromatinNumber)

while true

    enveMovements = enveVelocities.*dt;
    if all(norm.(enveMovements) .< 0.05)
        chroMovements = chroVelocities.*dt;
        if all(norm.(chroMovements) .< 0.05)
            break
        else
            dt = dt/2;
        end
    else
        dt = dt/2;
    end
end

nuc.vert .+= enveMovements
chro.vert .+= chroMovements