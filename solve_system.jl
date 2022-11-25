function solve_system!(nuc,spar,simset,dt,ext)

    
    movements = Vector{Vec{3,Float64}}(undef,length(nuc.enve.vert)+length(nuc.chro.vert))

    if simset.simType == "VRC"
        replCompMovements = Vector{Vec{3,Float64}}(undef,length(nuc.repl.vert))
    end
    maxMovement::Float64 = 0;
 
    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = nuc.enve.forces.total .+ enveFlucs;
        totalChro = nuc.chro.forces.total .+ fluctuationForces;

        solX = cg(simset.frictionMatrix,[getindex.(totalEnve,1);getindex.(totalChro,1)],Pl = simset.iLU)
        solY = cg(simset.frictionMatrix,[getindex.(totalEnve,2);getindex.(totalChro,2)],Pl = simset.iLU)
        solZ = cg(simset.frictionMatrix,[getindex.(totalEnve,3);getindex.(totalChro,3)],Pl = simset.iLU)

        movements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

        maxMovement = 0;

        for i = eachindex(movements)
            
            movementNorm = norm(movements[i]);
            if movementNorm >= maxMovement
                maxMovement = movementNorm
            end

            if movementNorm >= 0.5
                everythingIsFine = false
                simset.timeStepMultiplier = simset.timeStepMultiplier/2
                break
            end
        end

        if simset.simType == "VRC"
            solX = cg(nuc.repl.frictionMatrix,getindex.(nuc.repl.forces.total,1),Pl = nuc.repl.iLU)
            solY = cg(nuc.repl.frictionMatrix,getindex.(nuc.repl.forces.total,2),Pl = nuc.repl.iLU)
            solZ = cg(nuc.repl.frictionMatrix,getindex.(nuc.repl.forces.total,3),Pl = nuc.repl.iLU)

            replCompMovements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

            for i = eachindex(replCompMovements)
            
                movementNorm = norm(replCompMovements[i]);
                if movementNorm >= maxMovement
                    maxMovement = movementNorm
                end
    
                if movementNorm >= 0.5
                    everythingIsFine = false
                    simset.timeStepMultiplier = simset.timeStepMultiplier/2
                    break
                end
            end
        end

        if everythingIsFine
            break
        end

    end

    nuc.enve.vert .+= movements[1:length(nuc.enve.vert)]
    nuc.chro.vert .+= movements[length(nuc.enve.vert)+1:end]

    if simset.simType == "VRC"
        nuc.repl.vert .+= replCompMovements
        nuc.repl.normalVolume += 5000*dt*simset.timeStepMultiplier
    elseif simset.simType == "PC"
        planeMovement = (9e5 - ext[4])*dt*simset.timeStepMultiplier
        if planeMovement > 0.001
            ext[1] -= 10*dt.*simset.timeStepMultiplier
        else
            println((9e5 - ext[4]))
            ext[1] -= planeMovement
        end
    end


    multiplier::Rational = 1;
    if simset.timeStepMultiplier != 1 && maxMovement <= 0.2
        multiplier = 2;
    end
    
    while true
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier
            simset.timeStepProgress += simset.timeStepMultiplier
            if simset.timeStepProgress == 1
                simset.timeStepProgress = 0
            end
            break
        else
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;
        end
    end
end

function solve_system_init!(nuc,spar,simset,dt,ext,noEnve)

    
    movements = Vector{Vec{3,Float64}}(undef,length(nuc.enve.vert)+length(nuc.chro.vert))

    maxMovement::Float64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = nuc.enve.forces.total .+ enveFlucs;
        totalChro = nuc.chro.forces.total .+ fluctuationForces;

        solX = cg(simset.frictionMatrix,[getindex.(totalEnve,1);getindex.(totalChro,1)],Pl = simset.iLU)
        solY = cg(simset.frictionMatrix,[getindex.(totalEnve,2);getindex.(totalChro,2)],Pl = simset.iLU)
        solZ = cg(simset.frictionMatrix,[getindex.(totalEnve,3);getindex.(totalChro,3)],Pl = simset.iLU)

        movements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

        maxMovement = 0;

        for i = eachindex(movements)
            
            movementNorm = norm(movements[i]);
            if movementNorm >= maxMovement
                maxMovement = movementNorm
            end

            if movementNorm >= 0.5
                everythingIsFine = false
                simset.timeStepMultiplier = simset.timeStepMultiplier/2
                break
            end
        end

        if everythingIsFine
            break
        end

    end
    if !noEnve
        nuc.enve.vert .+= movements[1:length(nuc.enve.vert)]
    end
    nuc.chro.vert .+= movements[length(nuc.enve.vert)+1:end]

    if simset.simType == "PC"
        planeMovement = (9e5 - ext[4])*dt.*simset.timeStepMultiplier
        if planeMovement > 0.001
            ext[1] -= 10*dt.*simset.timeStepMultiplier
        else
            println((9e5 - ext[4]))
            ext[1] -= planeMovement
        end
    end

    multiplier::Rational = 1;
    if simset.timeStepMultiplier != 1 && maxMovement <= 0.2
        multiplier = 2;
    end
    
    while true
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier
            simset.timeStepProgress += simset.timeStepMultiplier
            if simset.timeStepProgress == 1
                simset.timeStepProgress = 0
            end
            break
        else
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;
        end
    end
end

function solve_system_adh_init!(nuc,spar,simset,dt,ext,forceToPlate)

    movements = Vector{Vec{3,Float64}}(undef,length(nuc.enve.vert))

    maxMovement::Float64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        
        totalEnve = nuc.enve.forces.total .+ enveFlucs;
        
        solX = cg(simset.frictionMatrix,getindex.(totalEnve,1),Pl = simset.iLU)
        solY = cg(simset.frictionMatrix,getindex.(totalEnve,2),Pl = simset.iLU)
        solZ = cg(simset.frictionMatrix,getindex.(totalEnve,3),Pl = simset.iLU)

        movements = Vec.(solX,solY,solZ).*dt.*simset.timeStepMultiplier

        maxMovement = 0;

        for i = eachindex(movements)
            
            movementNorm = norm(movements[i]);
            if movementNorm >= maxMovement
                maxMovement = movementNorm
            end

            if movementNorm >= 0.5
                everythingIsFine = false
                simset.timeStepMultiplier = simset.timeStepMultiplier/2
                break
            end
        end

        if everythingIsFine
            break
        end

    end

    nuc.enve.vert .+= movements

        F = 50e-9/spar.viscosity*spar.scalingTime/spar.scalingLength;
        forceDifference = (F - forceToPlate)/F
        planeMovement = dt.*simset.timeStepMultiplier
        # if planeMovement > 0.001
        ext[1] -= 1*forceDifference*dt.*simset.timeStepMultiplier
        # else
        #     println((9e5 - ext[4]))
        #     ext[1] -= planeMovement
        # end


    multiplier::Rational = 1;
    if simset.timeStepMultiplier != 1 && maxMovement <= 0.2
        multiplier = 2;
    end
    
    while true
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier
            simset.timeStepProgress += simset.timeStepMultiplier
            if simset.timeStepProgress == 1
                simset.timeStepProgress = 0
            end
            break
        else
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;
        end
    end
end