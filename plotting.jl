function plot_sphere!(nuc,i)

    P = Plots.plot(legend = false, size = (600,600); aspect_ratio=1,xlim = (-1.5, 3),ylim = (-1.5, 1.5),zlim = (-1.5, 1.5), camera = (0, 0));

    for i = 1:length(nuc.edges[:,1])
        if nuc.firstEdges[i] == 1
            Plots.plot!(P,[nuc.vert[nuc.edges[i,1]][1] nuc.vert[nuc.edges[i,2]][1]]',[nuc.vert[nuc.edges[i,1]][2] nuc.vert[nuc.edges[i,2]][2]]',[nuc.vert[nuc.edges[i,1]][3] nuc.vert[nuc.edges[i,2]][3]]',linecolor=:black);
        end
    end

    #a = [10];

    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    Plots.plot!(rectangle(5,-0.2,1.0583,-0.3), color="black", opacity=.5)
    Plots.plot!(rectangle(5,0.2,1.0583,0.3), color="black",opacity=.5)
    #Plots.scatter!(P,nuc.x,nuc.y,nuc.z,markercolor="black")
    #Plots.scatter!(P,nuc.x[a],nuc.y[a],nuc.z[a],markercolor="red")

    # triCenters = zeros(Float64,size(nuc.tri,1),3)

    # for i = 1:size(nuc.tri,1)
        # triCenters[i,:] = [sum(nuc.x[nuc.tri[i,:]]) sum(nuc.y[nuc.tri[i,:]]) sum(nuc.z[nuc.tri[i,:]])]./3;
    # end
    Plots.plot(P, plot_title="$i");
    #quiver!(P,triCenters[:,1],triCenters[:,2],triCenters[:,3],quiver=(nuc.triangleNormalUnitVectors[:,1],nuc.triangleNormalUnitVectors[:,2],nuc.triangleNormalUnitVectors[:,3]))
    # quiver!(P,nuc.x,nuc.y,nuc.z,quiver=(elasticForces[:,1],elasticForces[:,2],elasticForces[:,3]))
    P

end