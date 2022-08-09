function plot_sphere!(nuc,elasticForces,i)

    P = Plots.plot(legend = false, size = (600,600); aspect_ratio=1,xlim = (-1.5, 3),ylim = (-1.5, 1.5),zlim = (-1.5, 1.5), camera = (0, 0));

    for i = 1:length(nuc.edges[:,1])
        if nuc.firstEdges[i] == 1
            Plots.plot!(P,[nuc.x[nuc.edges[i,1]] nuc.x[nuc.edges[i,2]]]',[nuc.y[nuc.edges[i,1]] nuc.y[nuc.edges[i,2]]]',[nuc.z[nuc.edges[i,1]] nuc.z[nuc.edges[i,2]]]',linecolor=:black);
        end
    end

    a = [10];

    #Plots.scatter!(P,nuc.x,nuc.y,nuc.z,markercolor="black")
    Plots.scatter!(P,nuc.x[a],nuc.y[a],nuc.z[a],markercolor="red")

    # triCenters = zeros(Float64,size(nuc.tri,1),3)

    # for i = 1:size(nuc.tri,1)
        # triCenters[i,:] = [sum(nuc.x[nuc.tri[i,:]]) sum(nuc.y[nuc.tri[i,:]]) sum(nuc.z[nuc.tri[i,:]])]./3;
    # end
    Plots.plot(P, plot_title="$i");
    #quiver!(P,triCenters[:,1],triCenters[:,2],triCenters[:,3],quiver=(nuc.triangleNormalUnitVectors[:,1],nuc.triangleNormalUnitVectors[:,2],nuc.triangleNormalUnitVectors[:,3]))
    # quiver!(P,nuc.x,nuc.y,nuc.z,quiver=(elasticForces[:,1],elasticForces[:,2],elasticForces[:,3]))
    P

end