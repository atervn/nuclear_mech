function get_forces!(enve, chro ,spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(enve,ext[1],spar);
        aspiration = get_aspiration_forces(enve,ext[1],spar);
    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)
    end

    get_linear_chromatin_forces!(chro,spar);
    get_bending_chromatin_forces!(chro,spar)
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(enve,chro,spar,simset.envelopeTree)
    get_crosslink_forces!(chro,spar)
    get_lad_forces!(enve,chro,spar)

    if simset.adh.adherent
        get_plane_repulsion(enve,simset,spar)
    end

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;# ;# .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion
    
    if simset.adh.adherent
        enve.forces.total .+= enve.forces.planeRepulsion;
    end

    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    end

    chro.forces.total = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ chro.forces.enveRepulsion .+ chro.forces.crosslink .+ chro.forces.ladChroForces
end

function get_forces!(enve, chro, repl, spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);
    get_repl_comp_enve_repulsion_forces!(enve, repl, spar)

    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)


    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);

    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)

    elseif cmp(simset.simType,"PC") == 0
        planeRepulsion = get_plane_repulsion(nuc,ext,spar)
    end

    get_linear_chromatin_forces!(chro,spar);
    get_bending_chromatin_forces!(chro,spar)
    get_chromation_chromation_repulsion_forces!(chro,spar,simset.chromatinTree)
    get_envelope_chromatin_repulsion_forces!(enve,chro,spar,simset.envelopeTree)
    get_crosslink_forces!(chro,spar)
    get_lad_forces!(enve,chro,spar)

    if simset.adh.adherent
        get_plane_repulsion(enve,simset,spar)
    end

    enve.forces.total = enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces .+ enve.forces.replCompRepulsion;  # .+ enve.forces.envelopeRepulsion
    
    if simset.adh.adherent
        enve.forces.total .+= enve.forces.planeRepulsion;
    end

    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    end

    get_repl_comp_elastic_forces!(repl,spar)
    get_repl_comp_volume_forces!(repl,spar)
    get_repl_comp_area_forces!(repl,spar)
    get_repl_comp_bending_forces!(repl,spar)
    get_repl_comp_chromatin_repulsion_forces!(chro,repl,spar)

    repl.forces.total = repl.forces.elastic .+ repl.forces.volume .+ repl.forces.area .+ repl.forces.bending .+ repl.forces.chromationRepulsion .+ repl.forces.envelopeRepulsion;

    chro.forces.total = chro.forces.linear .+ chro.forces.bending .+ chro.forces.chroRepulsion .+ chro.forces.enveRepulsion .+ chro.forces.crosslink .+ chro.forces.ladChroForces .+ chro.forces.replCompRepulsion
end

function get_forces!(enve, spar, ext, simset)

    get_volume_forces!(enve,spar);
    get_area_forces!(enve,spar);
    get_bending_forces!(enve,spar);
    get_elastic_forces!(enve,spar);

    
    # get_repulsion_forces!(nuc,spar,simset.envelopeTree);
    # enveFlucs = get_random_enve_fluctuations(spar,nuc)
    
    
    if cmp(simset.simType,"MA") == 0 
        repulsion = get_aspiration_repulsion_forces(nuc,ext[1],spar);
        aspiration = get_aspiration_forces(nuc,ext[1],spar);
    elseif cmp(simset.simType,"MM") == 0
        micromanipulation = get_micromanipulation_forces(enve,ext[1],spar)
    end
    
    if simset.adh.adherent
        get_plane_repulsion(enve,simset,spar)
    end
        
    enve.forces.total = enve.forces.elastic;#enve.forces.volume .+ enve.forces.area .+ enve.forces.elastic .+ enve.forces.bending; # .+ enve.forces.chromationRepulsion .+ enve.forces.ladEnveForces;  # .+ enve.forces.envelopeRepulsion

    if simset.adh.adherent
        enve.forces.total .+= enve.forces.planeRepulsion;
    end

    if cmp(simset.simType,"MA") == 0 
        enve.forces.total .+=  repulsion .+ aspiration
    elseif cmp(simset.simType,"MM") == 0
        enve.forces.total .+= micromanipulation
    end

end