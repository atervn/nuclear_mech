function plot_sphere!(nucleus)

    P = plot(legend = false; aspect_ratio=:equal);

    scatter!(P,nucleus.x,nucleus.y,nucleus.z)

    for i = 1:length(nucleus.edges[:,1])
        if nucleus.firstEdges[i] == 1
            plot!(P,[nucleus.x[nucleus.edges[i,1]] nucleus.x[nucleus.edges[i,2]]]',[nucleus.y[nucleus.edges[i,1]] nucleus.y[nucleus.edges[i,2]]]',[nucleus.z[nucleus.edges[i,1]] nucleus.z[nucleus.edges[i,2]]]',linecolor=:black);
        end
    end
    quiver!(P,nucleus.x,nucleus.y,nucleus.z,quiver=(nucleus.surfaceNormalUnitVectors[:,1],nucleus.surfaceNormalUnitVectors[:,2],nucleus.surfaceNormalUnitVectors[:,3]))
    P

end