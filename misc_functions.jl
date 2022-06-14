function get_friction_matrix(nucleus)

    frictionMatrix = sparse([],[],[],3*length(nucleus.x),3*length(nucleus.x));

    return frictionMatrix

end