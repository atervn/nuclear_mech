using LinearAlgebra

n = 10000;
a = rand(n,3);
b = rand(n,3);
d1 = zeros(Float64,n,1);
d2 = zeros(Float64,n,1);
@time begin


    for i = 1:n

        d1[i] = norm(a[i,:] .- b[i,:],2);

    end

end
x = a[:,1];
y = a[:,2];
z = a[:,3];
x2 = b[:,1];
y2 = b[:,2];
z2 = b[:,3];
@time begin


    
    for i = 1:n

        d2[i] = norm([x[i] - x2[i], y[i] - y2[i], z[i] - z2[i]],2);

    end

end

d1 == d2