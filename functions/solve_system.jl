function solve_system!(enve::envelopeType, chro::chromatinType, spar::scaledParametersType ,simset::simulationSettingsType , ext, intTime)

    # initialize the movements array to store vertex movements
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert)+length(chro.vert))

    # initialize the maximum movement variable
    maxMovement::Float64 = 0;
 
    while true

        # generate random fluctuations for envelope and chromatin
        enveluctuation = get_random_fluctuations(spar,simset,length(enve.vert))
        chroFluctuation = get_random_fluctuations(spar,simset,spar.chromatinLength*spar.chromatinNumber)
        
        # compute total forces on envelope and chromatin
        totalEnve = enve.forces.total .+ enveluctuation;
        totalChro = 1/spar.chromatinViscosityMultiplier.*(chro.forces.total .+ chroFluctuation);

        # solve the linear system of equations using conjugate gradient method
        movements,maxMovement = run_cg!(spar,simset,simset.frictionMatrix,simset.iLU,[getindex.(totalEnve,1);getindex.(totalChro,1)],[getindex.(totalEnve,2);getindex.(totalChro,2)],[getindex.(totalEnve,3);getindex.(totalChro,3)])

        # check if the run_cg! function did not return false (indicating an issue)
        if !isa(movements,Bool)

            # exit the loop if the linear system is successfully solved
            break
        end

    end

    # move vertices of envelope and chromatin based on the computed movements
    move_vertices!(enve,chro,movements,simset)

    simset.totalEnvelopeSpeed = sum(norm.(movements[1:length(enve.vert)]))/(spar.dt.*simset.timeStepMultiplier)

    # update the position of the adherens plane
    move_adherens_plane!(enve,simset,spar)

    # perform time stepping based on the maximum movement
    time_stepping!(spar,simset,maxMovement)

    # move the micromanipulator
    move_micromanipulator!(ext, spar, simset, intTime)

    # move afm
    move_afm!(enve,ext,spar,simset)

end

# REPLICATION COMPARTMENT
function solve_system!(enve::envelopeType, chro::chromatinType, repl::replicationCompartmentType, spar::scaledParametersType,simset::simulationSettingsType,ext,intTime)

    # initialize the movement arrays to store vertex movements
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert)+length(chro.vert))
    replMovements = Vector{Vec{3,Float64}}(undef,length(repl.vert))

    # initialize the maximum movement variables
    maxMovement::Float64 = 0;
    replMaxMovement::Float64 = 0;
    
    while true

        # generate random fluctuations for envelope and chromatin
        enveluctuation = get_random_fluctuations(spar,simset,length(enve.vert))
        chroFluctuation = get_random_fluctuations(spar,simset,spar.chromatinLength*spar.chromatinNumber)
                
        # compute total forces on envelope and chromatin
        totalEnve = enve.forces.total .+ enveluctuation;
        totalChro = chro.forces.total .+ chroFluctuation;

        # solve the linear system of equations using conjugate gradient method
        movements,maxMovement = run_cg!(spar,simset,simset.frictionMatrix,simset.iLU,[getindex.(totalEnve,1);getindex.(totalChro,1)],[getindex.(totalEnve,2);getindex.(totalChro,2)],[getindex.(totalEnve,3);getindex.(totalChro,3)])

        # check if the run_cg! function did not return false (indicating an issue)
        if !isa(movements,Bool)

            # solve the linear system of equations using conjugate gradient method for the replication compartment
            replMovements,replMaxMovement = run_cg!(spar,simset,repl.frictionMatrix,repl.iLU,getindex.(repl.forces.total,1),getindex.(repl.forces.total,2),getindex.(repl.forces.total,3))

            # check if the run_cg! function did not return false (indicating an issue)
            if !isa(replMovements,Bool)

                # exit the loop if the linear system is successfully solved
                break

            end
        end
    end

    # check if replication compartment max movement is larger, and if yes, set it to the largest one
    if replMaxMovement > maxMovement
        maxMovement = replMaxMovement
    end

    # move vertices of envelope and chromatin based on the computed movements
    move_vertices!(enve,chro,movements,simset)

    # move vertices of replication compartment based on the computed movements
    move_repl_vertices!(repl,replMovements)

    # update the position of the adherens plane
    move_adherens_plane!(enve,simset,spar)
    
    # remodel lamina
    remodel_lamina!(enve,simset,spar)

    # perform time stepping based on the maximum movement
    time_stepping!(spar,simset,maxMovement)

    # move the micromanipulator
    move_micromanipulator!(ext, spar, simset, intTime)

    # move afm
    move_afm!(enve,ext,spar,simset)

end

function solve_system!(enve, spar,simset,ext,intTime)

    # initialize the movements array to store vertex movements
    movements = Vector{Vec{3,Float64}}(undef,length(enve.vert))

    # initialize the maximum movement variable
    maxMovement::Float64 = 0;
 
    while true

        # generate random fluctuations for envelope
        enveFluctuation = get_random_fluctuations(spar,simset,length(enve.vert))
        
        # compute total forces on envelope
        totalEnve = enve.forces.total .+ enveFluctuation;
        
        # solve the linear system of equations using conjugate gradient method
        movements,maxMovement = run_cg!(spar,simset,simset.frictionMatrix,simset.iLU,getindex.(totalEnve,1),getindex.(totalEnve,2),getindex.(totalEnve,3))


        # check if the run_cg! function did not return false (indicating an issue)
        if !isa(movements,Bool)

            # exit the loop if the linear system is successfully solved
            break
        end

    end

    # move vertices of envelope based on the computed movements
    move_vertices!(enve,movements)

    # update the position of the adherens plane
    move_adherens_plane!(enve,simset,spar)

    # perform time stepping based on the maximum movement
    time_stepping!(spar,simset,maxMovement)

    # move the micromanipulator
    move_micromanipulator!(ext, spar, simset, intTime)

    # move afm
    move_afm!(enve,ext,spar,simset)
end

function run_cg!(spar,simset,frictionMatrix,iLU,forcesX,forcesY,forcesZ)
    
    # flag for checking the status
    everythingIsFine = true

    # solve the linear equations using conjugate gradient method for the different components
    solX = cg(frictionMatrix,forcesX,Pl = iLU)
    solY = cg(frictionMatrix,forcesY,Pl = iLU)
    solZ = cg(frictionMatrix,forcesZ,Pl = iLU)

    
    # compute the movements based on the solutions
    movements = Vec.(solX,solY,solZ).*spar.dt.*simset.timeStepMultiplier

    # init the max movement
    maxMovement = 0;

    # iterate through the movements
    for i = eachindex(movements)
            
        # compute the norm of the current movement
        movementNorm = norm(movements[i]);

        # if too large movement
        if movementNorm >= spar.maxMovement

            # change flag status
            everythingIsFine = false

            # reduce the time step multiplier
            simset.timeStepMultiplier = simset.timeStepMultiplier/2
            
            break

        end
        
        # update the maximum movement if the current movement is larger
        if movementNorm >= maxMovement
            maxMovement = movementNorm
        end

    end

    if everythingIsFine

        # return the movements and maximum movement
        return movements,maxMovement

    else

        # return false to indicate an issue occurred
        return false,false

    end
end

function move_vertices!(enve,chro,movements,simset)

    # If the flag `noEnveSolve` is false, move the enve vertices
    if !simset.noEnveSolve
        enve.vert .+= movements[1:length(enve.vert)]
    end

    # move chro vertices 
    chro.vert .+= movements[length(enve.vert)+1:end]

end

function move_vertices!(enve,movements)

    # move enve vertices
    enve.vert .+= movements
    
end

function move_adherens_plane!(enve,simset,spar)

    # if adherent
    if simset.adh.adherent && !simset.adh.static

        nucleusHeight = simset.adh.topPlane + spar.cytoskeletonPlaneRadius - simset.adh.bottomPlane - 2*spar.repulsionDistance

        if abs(nucleusHeight - spar.nucleusHeight) > 0.1

            if nucleusHeight > spar.nucleusHeight
                simset.adh.topPlane -= 5*spar.dt.*simset.timeStepMultiplier
            elseif nucleusHeight < spar.nucleusHeight
                simset.adh.topPlane += 5*spar.dt.*simset.timeStepMultiplier
            end
        end
        
    end
end

function time_stepping!(spar,simset,maxMovement)

    # set multiplier to 1 by default
    multiplier::Rational = 1;

    # check if time step multiplier is not 1 and maxMovement is less than or equal to limit
    if simset.timeStepMultiplier != 1 && maxMovement <= spar.minMovement

        # double the multiplier
        multiplier = 2;
    end
    
    while true

        # check if increasing in time will not exceed 1
        if simset.timeStepProgress + simset.timeStepMultiplier*multiplier <= 1

            # update the time step size
            simset.timeStepMultiplier = simset.timeStepMultiplier*multiplier

            # progress the time
            simset.timeStepProgress += simset.timeStepMultiplier

            # if the current time step is 1 (full time step has been progressed)
            if simset.timeStepProgress == 1

                # set the time step to 0
                simset.timeStepProgress = 0
            end

            break

        # otherwise the current step size goes beyond the full time step
        else

            # half the time step
            simset.timeStepMultiplier = simset.timeStepMultiplier/2;

        end
    end

end

function move_repl_vertices!(repl,replMovements)

    # move the replication compartment vertices
    repl.vert .+= replMovements

end

function remodel_lamina!(enve,simset,spar)

    # if simset.laminaRemodel == "tension_remodeling"
    #     strains = (enve.edgeVectorNorms .- enve.normalLengths)./enve.normalLengths
    #     for i = eachindex(enve.edges)
    #         if strains[i] > 0.2
    #             enve.envelopeMultipliers[i] -= 0.5*spar.dt.*simset.timeStepMultiplier
    #         end
    #     end
    # elseif simset.laminaRemodel == "strain remodeling"

    # end

end

function move_afm!(enve,ext,spar,simset)

    if simset.simType == "AFM"

        forces = getindex.(enve.forces.total[ext.touchingIdx],3) - getindex.(enve.forces.afmRepulsion[ext.touchingIdx],3)

        forces = forces[findall(forces .> 0)]

        enveForcesOnBead = sum(forces)

        friction = 0.0001/spar.viscosity;

        distance = norm(ext.topPosition - ext.beadPosition);

        springConstant = 0.05/spar.viscosity*spar.scalingTime;

        cantileverForce = springConstant*(distance - ext.normDistance)

        if !ext.directionDown
            ext.topPosition = ext.topPosition + Vec(0.,0.,spar.cantileverSpeed*spar.dt*simset.timeStepMultiplier)
        elseif abs(cantileverForce) < 2e-9/spar.viscosity/spar.scalingLength*spar.scalingTime && ext.directionDown
            ext.topPosition = ext.topPosition - Vec(0.,0.,spar.cantileverSpeed*spar.dt*simset.timeStepMultiplier)
        else
            ext.directionDown = false
            simset.stopSimulation = true
        end

        ext.forceOnBead = enveForcesOnBead

        areaForceOnBead = getindex.(enve.forces.area[ext.touchingIdx],3)
        volumeForceOnBead = getindex.(enve.forces.volume[ext.touchingIdx],3)
        bendingForceOnBead = getindex.(enve.forces.bending[ext.touchingIdx],3)
        elasticForceOnBead = getindex.(enve.forces.elastic[ext.touchingIdx],3)
        chroRepulsionForceOnBead = getindex.(enve.forces.chromationRepulsion[ext.touchingIdx],3)
        ladForceOnBead = getindex.(enve.forces.ladEnveForces[ext.touchingIdx],3)
        cytoskeletonForceOnBead = getindex.(enve.forces.afmRepulsion[ext.touchingIdx],3)

        ext.areaForceOnBead = sum(areaForceOnBead[findall(areaForceOnBead  .> 0)])
        ext.volumeForceOnBead = sum(volumeForceOnBead[findall(volumeForceOnBead .> 0)])
        ext.bendingForceOnBead = sum(bendingForceOnBead[findall(bendingForceOnBead .> 0)])
        ext.elasticForceOnBead = sum(elasticForceOnBead[findall(elasticForceOnBead .> 0)])
        ext.chroRepulsionForceOnBead = sum(chroRepulsionForceOnBead[findall(chroRepulsionForceOnBead .> 0)])
        ext.ladForceOnBead = sum(ladForceOnBead[findall(ladForceOnBead .> 0)])
        ext.cytoskeletonForceOnBead = sum(cytoskeletonForceOnBead)

        ext.negareaForceOnBead = sum(areaForceOnBead[findall(areaForceOnBead  .< 0)])
        ext.negvolumeForceOnBead = sum(volumeForceOnBead[findall(volumeForceOnBead .< 0)])
        ext.negbendingForceOnBead = sum(bendingForceOnBead[findall(bendingForceOnBead .< 0)])
        ext.negelasticForceOnBead = sum(elasticForceOnBead[findall(elasticForceOnBead .< 0)])
        ext.negchroRepulsionForceOnBead = sum(chroRepulsionForceOnBead[findall(chroRepulsionForceOnBead .< 0)])
        ext.negladForceOnBead = sum(ladForceOnBead[findall(ladForceOnBead .< 0)])
        ext.negcytoskeletonForceOnBead = sum(cytoskeletonForceOnBead)

        forceOnBead = enveForcesOnBead + cantileverForce;
        # minus cantileverSpeed, since it is moving downward
        beadSpeed = (-friction*spar.cantileverSpeed + forceOnBead)./(1+friction)

        ext.beadPosition = ext.beadPosition + Vec(0.,0.,spar.dt*simset.timeStepMultiplier*beadSpeed);
        
    end

end