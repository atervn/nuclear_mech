using LinearAlgebra

n = 10000;

x = rand(n,1);
y = rand(n,1);
z = rand(n,1);


@time begin

    for i = 1:n

        norm.([x[i].-x, y[i].-y, z[i].-z],2);

    end

end

@time begin


    
    for i = 1:n

        sqrt.((x[i].-x).^2 .+ (y[i].-y).^2 .+ (z[i].-z).^2);

    end

end

xyz = [x y z];

@time begin


    
    @inbounds for i = 1:n
        sqrt.((xyz[i,1].-xyz[:,1]).^2 .+ (xyz[i,2].-xyz[:,2]).^2 .+ (xyz[i,3].-xyz[:,3]).^2);
    end

end

notify("Task completed")