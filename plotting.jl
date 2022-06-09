function plot_sphere!(nucleus,elasticForces)

    P = plot(legend = false; aspect_ratio=:equal);

    scatter!(P,nucleus.x,nucleus.y,nucleus.z)

    for i = 1:length(nucleus.edges[:,1])
        if nucleus.firstEdges[i] == 1
            plot!(P,[nucleus.x[nucleus.edges[i,1]] nucleus.x[nucleus.edges[i,2]]]',[nucleus.y[nucleus.edges[i,1]] nucleus.y[nucleus.edges[i,2]]]',[nucleus.z[nucleus.edges[i,1]] nucleus.z[nucleus.edges[i,2]]]',linecolor=:black);
        end
    end

    triCenters = zeros(Float64,size(nucleus.tri,1),3)

    for i = 1:size(nucleus.tri,1)
        triCenters[i,:] = [sum(nucleus.x[nucleus.tri[i,:]]) sum(nucleus.y[nucleus.tri[i,:]]) sum(nucleus.z[nucleus.tri[i,:]])]./3;
    end

    #quiver!(P,triCenters[:,1],triCenters[:,2],triCenters[:,3],quiver=(nucleus.triangleNormalUnitVectors[:,1],nucleus.triangleNormalUnitVectors[:,2],nucleus.triangleNormalUnitVectors[:,3]))
    quiver!(P,nucleus.x,nucleus.y,nucleus.z,quiver=(elasticForces[:,1],elasticForces[:,2],elasticForces[:,3]))
    P

end