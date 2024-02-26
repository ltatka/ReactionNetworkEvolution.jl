

struct A
    a::Int64
    b::Float64
end

st = A(5, 2.2)

test = Dict(fieldnames(A) .=> getfield.(Ref(test), fieldnames(S)))

println(test)