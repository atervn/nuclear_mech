function solve_system!(enve, chro, spar,simset,dt,ext,noEnveSolve)

    
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert)+length(chro.vert))

    if simset.simType == "VRC"
        replCompMovements = Vector{Vec{3,Float64}}(undef,length(repl.vert))
    end
    maxMovement::Float64 = 0;
 
    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,enve,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = enve.forces.total .+ enveFlucs;
        totalChro = chro.forces.total .+ fluctuationForces;

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

    if !noEnveSolve
        enve.vert .+= movements[1:length(enve.vert)]
    end
    chro.vert .+= movements[length(enve.vert)+1:end]

    if simset.adh.adherent
        cellForcesOnPlane = sum(getindex.(enve.forces.volume[simset.adh.touchingTop],3))
        F = 50e-9/spar.viscosity*spar.scalingTime/spar.scalingLength;

        forceDifference = (F - cellForcesOnPlane)/F
        simset.adh.topPlane -= 1*forceDifference*dt.*simset.timeStepMultiplier
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

function solve_system!(enve::envelopeType, chro::chromatinType, repl::replicationCompartmentType, spar::scaledParametersType,simset::simulationSettingsType,dt::Number,ext)

    
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert)+length(chro.vert))

    if simset.simType == "VRC"
        replCompMovements = Vector{Vec{3,Float64}}(undef,length(repl.vert))
    end
    maxMovement::Float64 = 0;
 
    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,enve,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
                
        totalEnve = enve.forces.total .+ enveFlucs;
        totalChro = chro.forces.total .+ fluctuationForces;

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

        solX = cg(repl.frictionMatrix,getindex.(repl.forces.total,1),Pl = repl.iLU)
        solY = cg(repl.frictionMatrix,getindex.(repl.forces.total,2),Pl = repl.iLU)
        solZ = cg(repl.frictionMatrix,getindex.(repl.forces.total,3),Pl = repl.iLU)

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

        if everythingIsFine
            break
        end

    end

    enve.vert .+= movements[1:length(enve.vert)]
    chro.vert .+= movements[length(enve.vert)+1:end]

    repl.vert .+= replCompMovements
    repl.normalVolume += 5000*dt*simset.timeStepMultiplier
    
    if simset.adh.adherent
        cellForcesOnPlane = sum(getindex.(enve.forces.volume[simset.adh.touchingTop],3))
        F = 50e-9/spar.viscosity*spar.scalingTime/spar.scalingLength;
        forceDifference = (F - cellForcesOnPlane)/F
        simset.adh.topPlane -= 1*forceDifference*dt.*simset.timeStepMultiplier
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

function solve_system!(enve, spar,simset,dt,ext)

    
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert))

    maxMovement::Float64 = 0;
 
    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,enve,dt)
        
        totalEnve = enve.forces.total .+ enveFlucs;

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

    enve.vert .+= movements[1:length(enve.vert)]

    if simset.adh.adherent
        cellForcesOnPlane = sum(getindex.(enve.forces.volume[simset.adh.touchingTop],3))
        F = 50e-9/spar.viscosity*spar.scalingTime/spar.scalingLength;
        forceDifference = (F - cellForcesOnPlane)/F
        simset.adh.topPlane -= 1*forceDifference*dt.*simset.timeStepMultiplier
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

function solve_system_init!(enve,chro,spar,simset,dt,ext,noEnve)

    
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert)+length(chro.vert))

    maxMovement::Float64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,enve,dt)
        fluctuationForces = get_random_fluctuations(spar,dt)
        
        totalEnve = enve.forces.total .+ enveFlucs;
        totalChro = chro.forces.total .+ fluctuationForces;

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
        enve.vert .+= movements[1:length(enve.vert)]
    end
    chro.vert .+= movements[length(enve.vert)+1:end]

    if simset.simType == "PC"
        planeMovement = (9e5 - ext[4])*dt.*simset.timeStepMultiplier
        if planeMovement > 0.001
            ext[1] -= 10*dt.*simset.timeStepMultiplier
        else
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

    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert))

    maxMovement::Float64 = 0;

    while true

        everythingIsFine = true
        enveFlucs = get_random_enve_fluctuations(spar,nuc,dt)
        
        totalEnve = enve.forces.total .+ enveFlucs;
        
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

    enve.vert .+= movements

        F = 50e-9/spar.viscosity*spar.scalingTime/spar.scalingLength;
        forceDifference = (F - forceToPlate)/F
        planeMovement = dt.*simset.timeStepMultiplier
        # if planeMovement > 0.001
        ext[1] -= 1*forceDifference*dt.*simset.timeStepMultiplier
        # else
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