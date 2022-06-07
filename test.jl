mutable struct nucleusType2
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

nucleus2 = nucleusType2([],[],[])

nucleus2.x

a = (1 + sqrt(5))/2;
x = [-1, 1, -1, 1, 0, 0, 0, 0, a, a, -a, -a]