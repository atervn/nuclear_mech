
n = 1000;
k = 1000;
x1 = rand(n);
y1 = rand(n);
z1 = rand(n);
x2 = rand(n);
y2 = rand(n);
z2 = rand(n);
A = Array{Point{3,Float64}}(undef,n)
B = Array{Point{3,Float64}}(undef,n)

for i = 1:n
    A[i] = Point(x1[i],y1[i],z1[i])
    B[i] = Point(x2[i],y2[i],z2[i])
end

@time begin
for j = 1:k
    for i = 1:n

        norm(A[i]-B[i]);

    end
end
end

@time begin
    for j = 1:k
    for i = 1:n

        norm([x1[i]-x2[i],y1[i]-y2[i],z1[i]-z2[i]]);

    end
end
end

@time begin
    for j = 1:k
    for i = 1:n

        sqrt((x1[i]-x2[i])^2+(y1[i]-y2[i])^2+(z1[i]-z2[i])^2);

    end
end
end