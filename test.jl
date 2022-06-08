using LinearAlgebra

n = 1000000;

@time begin

    a = rand(n,3);
    b = rand(n,3);
    c = rand(n,3);

    d1 = zeros(Float64,n,3);
    d2 = zeros(Float64,n,3);

    for i = 1:n

        d1[i,:] = cross(b[i,:].-a[i,:],c[i,:].-a[i,:]);

    end

end

@time begin


    
    for i = 1:n

        d2[i,:] = [((b[i,2]-a[i,2])*(c[i,3]-a[i,3])) - ((b[i,3]-a[i,3])*(c[i,2]-a[i,2]))
        ((b[i,3]-a[i,3])*(c[i,1]-a[i,1])) - ((b[i,1]-a[i,1])*(c[i,3]-a[i,3]))
        ((b[i,1]-a[i,1])*(c[i,2]-a[i,2])) - ((b[i,2]-a[i,2])*(c[i,1]-a[i,1]))];

    end

end

d1 == d2