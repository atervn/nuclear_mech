function plot_sphere!(nuc,i)

    P = Plots.plot(legend = false, size = (600,600); aspect_ratio=1,xlim = (-5, 5),ylim = (-5, 5),zlim = (-5, 5), camera = (0, 0));

    for i = 1:length(nuc.enve.edges[:,1])
        if nuc.enve.firstEdges[i] == 1
            Plots.plot!(P,[nuc.enve.vert[nuc.enve.edges[i,1]][1] nuc.enve.vert[nuc.enve.edges[i,2]][1]]',[nuc.enve.vert[nuc.enve.edges[i,1]][2] nuc.enve.vert[nuc.enve.edges[i,2]][2]]',[nuc.enve.vert[nuc.enve.edges[i,1]][3] nuc.enve.vert[nuc.enve.edges[i,2]][3]]',linecolor=:black);
        end
    end

    #a = [10];

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    # Plots.plot!(rectangle(5,-0.2,1.0583,-0.3), color="black", opacity=.5)
    # Plots.plot!(rectangle(5,0.2,1.0583,0.3), color="black",opacity=.5)
    #Plots.scatter!(P,nuc.enve.x,nuc.enve.y,nuc.enve.z,markercolor="black")
    #Plots.scatter!(P,nuc.enve.x[a],nuc.enve.y[a],nuc.enve.z[a],markercolor="red")

    # triCenters = zeros(Float64,size(nuc.enve.tri,1),3)

    # for i = 1:size(nuc.enve.tri,1)
        # triCenters[i,:] = [sum(nuc.enve.x[nuc.enve.tri[i,:]]) sum(nuc.enve.y[nuc.enve.tri[i,:]]) sum(nuc.enve.z[nuc.enve.tri[i,:]])]./3;
    # end
    Plots.plot(P, plot_title="$i");
    #quiver!(P,triCenters[:,1],triCenters[:,2],triCenters[:,3],quiver=(nuc.enve.triangleNormalUnitVectors[:,1],nuc.enve.triangleNormalUnitVectors[:,2],nuc.enve.triangleNormalUnitVectors[:,3]))
    # quiver!(P,nuc.enve.x,nuc.enve.y,nuc.enve.z,quiver=(elasticForces[:,1],elasticForces[:,2],elasticForces[:,3]))
    P

end