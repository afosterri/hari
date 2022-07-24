
using NLsolve, LinearAlgebra
function fg(x)
    if x>=0
    x^(1/2)
    else
    x=-x^2*1000
    end
end
function dfg(x)
  
    if x>=0
        1/2*x^(-1/2)
        else
        x=-x*2000
        end
    
end

function prepare()

    function f!(F,x)
        lambda=(x[nplots+1])
        s=0
        for i=1:nplots
        F[i]=dfg(x[i])-lambda*p[i]
        s=p[i]*x[i]+s
        end
        F[nplots+1]=(tB-s)
    end

nplots=100
p=[i for i=1:nplots]
tB=100
start=ones(nplots+1)
bottom=zeros(nplots+1)
top=[Inf for i=1:nplots+1]
res=nlsolve(f!,start)

#mcpsolve(f!,bottom, top,start)
end
#prepare()

#=
using NLsolve

function f!(F, x)
    F[1] = (x[1]+3)*(x[2]^3-7)+18
    F[2] = sin(x[2]*exp(x[1])-1)
end

function j!(J, x)
    J[1, 1] = x[2]^3-7
    J[1, 2] = 3*x[2]^2*(x[1]+3)
    u = exp(x[1])*cos(x[2]*exp(x[1])-1)
    J[2, 1] = x[2]*u
    J[2, 2] = u
end

nlsolve(f!, [ 0.1; 1.2])
=#

function g(x,A)
    y=x+A*x
end
function dg(x,A)
    (g.(x+I.*delta,A).-g(x,A))./delta
end
A=[1  2;  3 4]
x=[5 ;6]
delta=.001
dg.(x,A)


