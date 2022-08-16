function main(maxT)

nuc = nucleusType();

ipar = inputParametersType();
ipar = read_parameters(ipar,"./parameters.txt");

radius = ipar.freeNucleusRadius/ipar.scalingLength;
nSubdivisions = 4;

nuc = create_icosahedron!(nuc,radius);
nuc = get_edges!(nuc);
nuc = subdivide_mesh!(nuc,radius,nSubdivisions)
get_triangle_normals!(nuc);

nuc.neighboringTriangles = zeros(Int64,size(nuc.edges,1),2)
nuc.edges3vertex = zeros(Int64,size(nuc.edges,1),2);
for i = 1:size(nuc.edges,1)
    temp = findall(sum(nuc.tri .== nuc.edges[i,1],dims=2) .> 0 .&& sum(nuc.tri .== nuc.edges[i,2],dims=2) .> 0);
    nuc.neighboringTriangles[i,:] = [j[1] for j in temp];
    thirdVertex1 = nuc.tri[nuc.neighboringTriangles[i,1],.!(nuc.tri[nuc.neighboringTriangles[i,1],:] .== nuc.edges[i,1] .|| nuc.tri[nuc.neighboringTriangles[i,1],:] .== nuc.edges[i,2])][1];
    thirdVertex2 = nuc.tri[nuc.neighboringTriangles[i,2],.!(nuc.tri[nuc.neighboringTriangles[i,2],:] .== nuc.edges[i,1] .|| nuc.tri[nuc.neighboringTriangles[i,2],:] .== nuc.edges[i,2])][1];

    nuc.edges3vertex[i,1] = thirdVertex1;
    nuc.edges3vertex[i,2] = thirdVertex2; 
end


get_triangle_normals!(nuc);
nuc.testview = Array{Any}(undef,length(nuc.edges))
nuc.p1 = Array{Any}(undef,length(nuc.tri))
nuc.p2 = Array{Any}(undef,length(nuc.tri))
nuc.p3 = Array{Any}(undef,length(nuc.tri))
nuc.trii = Array{Any}(undef,length(nuc.tri))
nuc.ep1 = Array{Any}(undef,length(nuc.edges))
nuc.ep2 = Array{Any}(undef,length(nuc.edges))
nuc.ep31 = Array{Any}(undef,length(nuc.edges))
nuc.ep32 = Array{Any}(undef,length(nuc.edges))
for i = 1:size(nuc.edges,1)
    nuc.testview[i] = @view nuc.triangleNormalUnitVectors[nuc.edgesTri[i,:]];
    nuc.ep1[i] = @view nuc.vert[nuc.edges[i,1]];
    nuc.ep2[i] = @view nuc.vert[nuc.edges[i,2]];
    nuc.ep31[i] = @view nuc.vert[nuc.edges3vertex[i,1]];
    nuc.ep32[i] = @view nuc.vert[nuc.edges3vertex[i,2]];
end
for i = 1:size(nuc.tri,1)
    nuc.p1[i] = @view nuc.vert[nuc.tri[i,1]];
    nuc.p2[i] = @view nuc.vert[nuc.tri[i,2]];
    nuc.p3[i] = @view nuc.vert[nuc.tri[i,3]];
    nuc.trii[i] = @view nuc.vert[nuc.tri[i,:]];
end

nuc.normalVolume = get_volume!(nuc);
nuc.normalTriangleAreas = get_area!(nuc);
nuc.normalArea = sum(nuc.normalTriangleAreas);
nuc.normalAngle = mean(get_triangle_angles(nuc));
lengths = zeros(Float64,Int64(size(nuc.edges,1)));

nuc.z .+= 0

for i = 1:size(nuc.edges,1)  
        lengths[i] = norm(nuc.vert[nuc.edges[i,2]] - nuc.vert[nuc.edges[i,1]]);
end
nuc.normalLengths = lengths;

spar = scaledParametersType();
spar = get_model_parameters(ipar,spar,nuc);

frictionMatrix = get_friction_matrix(nuc,spar);

# maxT = 1000;

p = Progress(maxT)
anim = Animation();





pip = generate_pipette_mesh();

for t = 0:maxT
    
    get_voronoi_areas!(nuc);
    get_local_curvatures!(nuc);
    get_area_unit_vectors!(nuc);
    get_triangle_normals!(nuc);

    get_volume_forces!(nuc,spar);
    get_area_forces!(nuc,spar);
    get_bending_forces!(nuc,spar);
    get_elastic_forces!(nuc,spar);

    # get_repulsion_forces!(nuc,spar);
    repulsion = get_aspiration_repulsion_forces(nuc,pip,spar);
    #=
    flatRepulsion = flat_repulsion_forces(nuc,spar);
    if t >= 0
        flatRepulsion2 = flat_repulsion_forces2(nuc,spar,t);
    else
        flatRepulsion2 = zeros(Float64,length(nuc.x),3);
    end
    =#

    local totalForces = nuc.forces.volume .+ nuc.forces.area .+ nuc.forces.elastic .+ repulsion;# .+ nuc.forces.repulsion# + repulsionX .+ nuc.forces.bending  ;# .+ nuc.forces.elastic .+ nuc.forces.repulsion# + repulsionX;
    # local totalForcesY = nuc.forces.volumeY .+ nuc.forces.areaY .+ nuc.forces.bendingY .+ nuc.forces.elasticY .+ nuc.forces.repulsionY + repulsionY;
    # local totalForcesZ = nuc.forces.volumeZ .+ nuc.forces.areaZ .+ nuc.forces.bendingZ .+ nuc.forces.elasticZ .+ nuc.forces.repulsionZ + repulsionZ;# .+ flatRepulsion2[:,3];
    
    # local totalForces[41] = totalForces[41] + Vec(1.,0.,0.);

    local vX = cg(frictionMatrix,getindex.(totalForces,1));
    local vY = cg(frictionMatrix,getindex.(totalForces,2));
    local vZ = cg(frictionMatrix,getindex.(totalForces,3));

    for i = 1:length(nuc.vert)
        nuc.vert[i] += Vec(vX[i],vY[i],vZ[i])*0.01
    end

    # nuc.x = nuc.x .+ vX.*0.01;
    # nuc.y = nuc.y .+ vY.*0.01;
    # nuc.z = nuc.z .+ vZ.*0.01;

    if mod(t,10) == 0
        plot_sphere!(nuc,t)
        Plots.frame(anim)
    end

    next!(p)

end

gif(anim, "test2.gif", fps = 30)

end