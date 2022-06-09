using LinearAlgebra

n = 10000;

@time begin

    a = rand(n,3);
    b = rand(n,3);
    c = rand(n,3);

    d1 = zeros(Float64,n,1);
    d2 = zeros(Float64,n,1);

    for i = 1:n

        norm(a[i,:],2);
        norm(b[i,:],2);
        norm(c[i,:],2);

    end

end

@time begin


    
    for i = 1:n

        d2[i] = sqrt(sum(a[i,:].^2));

    end

end




d1 == d2