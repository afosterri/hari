using Distributed
@everywhere using Distributed, SharedArrays
y=ones(9)
x=SharedArray{Float64}(9)
x[:]=y
@sync @distributed for i=1:8
    x[i:(i+1)]=[i,-i]
end
println(x)