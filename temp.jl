function temp1(P, tri)

B = tri[1, :];
E0 = tri[2, :] - B;
   # E0 = E0/sqrt(sum(E0.^2)); %normalize vector
E1 = tri[3, :] - B;
   # E1 = E1/sqrt(sum(E1.^2)); %normalize vector
D = B - P;


a = E0[1]^2 + E0[2]^2 + E0[3]^2;
b = E0[1]*E1[1] + E0[2]*E1[2] + E0[3]*E1[3];
c = E1[1]^2 + E1[2]^2 + E1[3]^2;
d = E0[1]*D[1] + E0[2]*D[2] + E0[3]*D[3];
e = E1[1]*D[1] + E1[2]*D[2] + E1[3]*D[3];
f = D[1]^2 + D[2]^2 + D[3]^2;

#print "{0} {1} {2} ".format(B,E1,E0)
det = a * c - b * b
s = b * e - c * d
t = b * d - a * e

if (s + t) <= det
    if s < 0
        if t < 0
            # region4
            if d < 0
                t = 0
                if -d >= a
                    s = 1;
                    sqrdistance = a + 2*d + f;
                else
                    s = -d / a;
                    sqrdistance = d*s + f;
                end
            else
                s = 0;
                if e >= 0;
                    t = 0;
                    sqrdistance = f;
                else
                    if -e >= c
                        t = 1;
                        sqrdistance = c + 2*e + f;
                    else
                        t = -e/c;
                        sqrdistance = e*t + f;
                    end
                end
            end
        else
            # region 3
            s = 0;
            if e >= 0
                t = 0;
                sqrdistance = f;
            else
                if -e >= c
                    t = 1;
                    sqrdistance = c + 2*e + f;
                else
                    t = -e/c;
                    sqrdistance = e*t + f;
                end
            end
        end
    else
        if t < 0
            # region 5
            t = 0
            if d >= 0
                s = 0;
                sqrdistance = f;
            else
                if -d >= a
                    s = 1;
                    sqrdistance = a + 2*d + f;
                else
                    s = -d/a;
                    sqrdistance = d*s + f;
                end
            end
        else
            # region 0
            invDet = 1/det;
            s = s*invDet;
            t = t*invDet;
            sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
        end
    end
else
    if s < 0
        # region 2
        tmp0 = b + d;
        tmp1 = c + e;
        if tmp1 > tmp0  # minimum on edge s+t=1
            numer = tmp1 - tmp0;
            denom = a - 2*b + c;
            if numer >= denom
                s = 1;
                t = 0;
                sqrdistance = a + 2*d + f;
            else
                s = numer/denom;
                t = 1 - s;
                sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
            end

        else  # minimum on edge s=0
            s = 0.0
            if tmp1 <= 0
                t = 1;
                sqrdistance = c + 2*e + f;
            else
                if e >= 0
                    t = 0;
                    sqrdistance = f;
                else
                    t = -e/c;
                    sqrdistance = e*t + f;
                end
            end
        end
    else
        if t < 0
            # region6
            tmp0 = b + e;
            tmp1 = a + d;
            if tmp1 > tmp0
                numer = tmp1 - tmp0;
                denom = a - 2*b + c;
                if numer >= denom
                    t = 1;
                    s = 0;
                    sqrdistance = c + 2*e + f;
                else
                    t = numer / denom;
                    s = 1 - t;
                    sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                end
            else
                t = 0;
                if tmp1 <= 0
                    s = 1;
                    sqrdistance = a + 2*d + f;
                else
                    if d >= 0
                        s = 0;
                        sqrdistance = f;
                    else
                        s = -d/a;
                        sqrdistance = d*s + f;
                    end
                end
            end
        else
            # region 1
            numer = c + e - b - d;
            if numer <= 0
                s = 0
                t = 1
                sqrdistance = c + 2*e + f;
            else
                denom = a - 2*b + c;
                if numer >= denom
                    s = 1;
                    t = 0;
                    sqrdistance = a + 2*d + f;
                else
                    s = numer/denom;
                    t = 1 - s;
                    sqrdistance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                end
            end
        end
    end
end

# account for numerical round-off error
if sqrdistance < 0
    sqrdistance = 0
end

dist = sqrt(sqrdistance);

PP0 = Float64.(B + s*E0 + t*E1)
return dist, PP0

end

P = [1, 2, 1];

output = Matrix{Any}(undef,6,2);

@time begin

for i = 1:6
    
    local tri = rand(3,3);

    output[i,1],output[i,2] = temp1(P, tri)
    println(output[i,:])

end

end