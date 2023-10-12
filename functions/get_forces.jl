function get_forces!(enve::envelopeType, chro::chromatinType ,spar::scaledParametersType, ext, simset::simulationSettingsType)

    # calculate the forces on the envelope
    get_envelope_forces!(enve,spar,ext,simset)    

    # calculate the forces on the chromatin
    get_chromatin_forces!(enve,chro,spar,simset)

    # calculate the total forces on the envelope
    get_total_envelope_forces!(enve,simset,"noRepl")

    # calculate the total forces on the chromatin
    get_total_chromatin_forces!(chro,"noRepl")
    
end

function get_forces!(enve, chro, repl, spar, ext, simset)

    # calculate the forces on the envelope
    get_envelope_forces!(enve,repl,spar,ext,simset)

    # calculate the forces on the chromatin
    get_chromatin_forces!(enve,chro,spar,simset)

    # calculate the forces on the replication compartment
    get_repl_forces!(enve,chro,repl,spar)

    # calculate the total forces on the envelope with replication compartment
    get_total_envelope_forces!(enve,simset,"repl")
    
    # calculate the total forces on the chromatin with replication compartment
    get_total_chromatin_forces!(chro,"repl")

    # calculate the total forces on the replication compartment
    get_total_repl_forces!(repl)
    
end

function get_forces!(enve::envelopeType, spar::scaledParametersType, ext, simset::simulationSettingsType)

    # calculate the forces on the envelope
    get_envelope_forces!(enve,spar,ext,simset)
 
    # calculate the total forces on the envelope       
    get_total_envelope_forces!(enve,simset,"noRepl")

end

function get_envelope_forces!(enve,spar,ext,simset)

    # calculate volume forces on the envelope
    get_volume_forces!(enve,spar);
     
    # calculate area forces on the envelope
    get_area_forces!(enve,spar);
    
    # calculate bending forces on the envelope
    get_bending_forces!(enve,spar);
      
    # calculate elastic forces on the envelope
    get_elastic_forces!(enve,spar);

    # Calculate repulsion forces within the envelope
    get_repulsion_forces!(enve,spar,simset);

    # check if the simulation involves adhesion
    if simset.adh.adherent

        # calculate repulsion forces between the envelope and an adherent plane
        get_plane_repulsion!(enve,simset,spar)

    end
    
    # if micropipette aspiration simulation
    if simset.simType =="MA" 

        # calculate repulsion forces due to aspiration in micropipette aspiration (MA)
        get_aspiration_repulsion_forces!(enve,ext[1],spar);

        # Calculate aspiration forces in micropipette aspiration (MA) simulation
        get_aspiration_forces!(enve,spar);

    # if micromanipulation simulation
    elseif simset.simType == "MM"

        # calculate micromanipulation forces in micromanipulation (MM) simulation
        get_micromanipulation_forces!(enve,ext[1],spar)

    elseif simset.simType == "AFM"
       
        get_afm_forces!(enve,ext,spar)

    end

end

function get_envelope_forces!(enve,repl,spar,ext,simset)

    # calculate volume forces on the envelope
    get_volume_forces!(enve,repl,spar);
     
    # calculate area forces on the envelope
    get_area_forces!(enve,spar);
    
    # calculate bending forces on the envelope
    get_bending_forces!(enve,spar);
      
    # calculate elastic forces on the envelope
    get_elastic_forces!(enve,spar);

    # Calculate repulsion forces within the envelope
    get_repulsion_forces!(enve,spar,simset);

    # check if the simulation involves adhesion
    if simset.adh.adherent

        # calculate repulsion forces between the envelope and an adherent plane
        get_plane_repulsion!(enve,simset,spar)

    end
    
    # if micropipette aspiration simulation
    if simset.simType =="MA" 

        # calculate repulsion forces due to aspiration in micropipette aspiration (MA)
        get_aspiration_repulsion_forces!(enve,ext[1],spar);

        # Calculate aspiration forces in micropipette aspiration (MA) simulation
        get_aspiration_forces!(enve,spar);

    # if micromanipulation simulation
    elseif simset.simType == "MM"

        # calculate micromanipulation forces in micromanipulation (MM) simulation
        get_micromanipulation_forces!(enve,ext[1],spar)

    elseif simset.simType == "AFM"
       
        get_afm_forces!(enve,ext,spar)

    end

end

function get_total_envelope_forces!(enve,simset,replStatus)

    # calculate the total forces on the envelope by summing individual force components
    enve.forces.total = enve.forces.elastic .+ enve.forces.bending .+ enve.forces.volume .+ enve.forces.area; # .+ enve.forces.envelopeRepulsion

    # check if chromatin forces are considered
    if !simset.noChromatin

        # add chromatin repulsion forces and lad-envelope forces to the total envelope forces
        enve.forces.total .+= enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces
    end
    
    # check if the simulation involves adhesion
    if simset.adh.adherent

        # add plane repulsion forces to the total envelope forces
        enve.forces.total .+= enve.forces.planeRepulsion;

    end

    # Check the simulation type
    if simset.simType == "MA" 

        # add pipette repulsion forces and aspiration forces to the total envelope forces in micropipette aspiration (MA) simulation
        enve.forces.total .+=  enve.forces.pipetteRepulsion .+ enve.forces.aspiration
 
    elseif simset.simType == "MM" 
        
        # add micromanipulation forces to the total envelope forces in micromanipulation (MM) simulation
        enve.forces.total .+= enve.forces.micromanipulation
   
    elseif simset.simType == "AFM" 


        enve.forces.total .+= enve.forces.afmRepulsion

    end

    # check the replication compartment status
    if replStatus == "repl"

        # Add replication compartment repulsion forces to the total envelope forces
        enve.forces.total .+= enve.forces.replRepulsion
    end

end

function get_chromatin_forces!(enve,chro,spar,simset)

    # calculate linear forces on the chromatin
    get_linear_chromatin_forces!(chro,spar)
    
    # calculate bending forces on the chromatin
    get_bending_chromatin_forces!(chro,spar)
    
    # calculate chromatin-chromatin repulsion forces
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    
    # calculate envelope-chromatin repulsion forces
    get_envelope_chromatin_repulsion_forces!(enve,chro,spar,simset.envelopeTree)
    
    # calculate crosslink forces within the chromatin
    get_crosslink_forces!(chro,spar)

    # Calculate forces between the envelope and chromatin at the lamin-associated domain (LAD) regions
    get_lad_forces!(enve,chro,spar)

end

function get_total_chromatin_forces!(chro,replStatus)

    # calculate the total forces on the chromatin by summing individual force components
    chro.forces.total = chro.forces.linear .+ chro.forces.chroRepulsion .+ chro.forces.crosslink .+ chro.forces.ladChroForces .+ chro.forces.enveRepulsion .+ chro.forces.bending

    # Check the replication compartment status
    if replStatus == "repl"

        # add replication compartment repulsion forces to the total chromatin forces
        chro.forces.total .+= chro.forces.replRepulsion
    end

end

function get_repl_forces!(enve,chro,repl,spar)
    
    # calculate replication compartment elastic forces
    get_repl_elastic_forces!(repl,spar)
    
    # calculate replication compartment volume forces
    get_repl_volume_forces!(repl,spar)
        
    # calculate replication compartment bending forces
    get_repl_bending_forces!(repl,spar)
    
    # calculate replication compartment envelope repulsion forces
    get_repl_comp_enve_repulsion_forces!(enve, repl, spar)
    
    # calculate replication compartment chromatin repulsion forces
    get_repl_chromatin_repulsion_forces!(chro,repl,spar)

end

function get_total_repl_forces!(repl)

    # calculate the total forces within the replication compartment
    repl.forces.total = repl.forces.elastic .+ repl.forces.volume .+ repl.forces.bending .+ repl.forces.chromationRepulsion .+ repl.forces.envelopeRepulsion;

end