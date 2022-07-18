using Distributed
@everywhere using Distributed,SharedArrays

function run()
min=1
max=20
step=.1
delta=.5
beta=1/2
xgrid=[i for i=min:step:max]
n=size(xgrid)[1]
n2=n*n
alpha=.5
k=3

function createdata(K,NN,BB,min,step,max) 
    #creates plots for K farmers, NN is distribution of number of plots, BB is range of input dist, min max is range of locations  
    N=rand(NN,K)
    #println("N ",N)
   
    idplot=zeros(Int16,K,20)

    i=0
    for k=1:K

        for n=1:N[k]
            i=i+1
            idplot[k,n]=i            
            #print("id" ,idplot)
        end
    end
    idplot
    nplots=sum(N)
    B=rand(BB,K)*K/nplots
    Ax=rand(nplots,2)*(max-min) .+ min
    input=rand(nplots)
    return K,nplots,Ax,N,idplot,B,input
end


function findneighbor(Ax)
    #finds neighboring plots
    mesh=delaunay(Ax)
    return mesh.vertex_neighbor_vertices
end



function clockwise(a1)
    #sorts neighbors so are clockwise
    nr=size(a1)[1]-1
    a=a1[1:nr,:]
    c=zeros(nr);
    b=sum(a1,dims=1)/(nr);
    for i=1:nr
        c[i]=atan(a[i,1]-b[1],a[i,2]-b[2])
    end
    d=sortperm(c);
    out=a[d,:];
    out=vcat(out,transpose(out[1,:]))
end

function spacecorr(Ax,sc,min,step,max)
    #perturbs locations so plots are more or less close depending on sc
    ngrid::Int16=(max-min)/step+1
    nplots=size(Ax)[1]
    Bx=zeros(nplots,2)
    ht=zeros(ngrid,ngrid,nplots)
    imin=Array{Int16}(undef,nplots)
    jmin=Array{Int16}(undef,nplots)
    for n=1:nplots
        #find closest grid point
        dmin=1000
        i=0
        for x=min:step:max
            i=i+1
            j=0
            for y=min:step:max
                j=j+1
                
                d=sqrt((Ax[n,1]-x)^2+(Ax[n,2]-y)^2)
                ht[i,j,n]=d
                if d<dmin 
                    dmin=d
                    imin[n]=i
                    jmin[n]=j
                end
            end
        end
    end
    #total distance of all plots to a given grid point
    ht2=sum(ht[:,:,i] for i=1:nplots)
    for n=1:nplots
        dx=0
        dy=0
        iv=imin[n]
        jv=jmin[n]
        try
        dx=ht2[iv+1,jv]-ht2[iv,jv]
        dy=ht2[iv,jv+1]-ht2[jv,jv]
        catch
            dx=0
            dy=0
        end    
        #println("dx ",dx," ",dy)
        #perturb location towards grid points that are farther from other points if sc>0

        Bx[n,1]=Ax[n,1]+dx*sc
        Bx[n,2]=Ax[n,2]+dy*sc
        #println("BB1 ",Bx[n,1])
        #println("BB2 ",Bx[n,2])
        #don't move outside boundaries
        if Bx[n,1]<min 
            Bx[n,1]=Ax[n,1]
        end
        if Bx[n,1]>max 
            Bx[n,1]=Ax[n,1]
        end
        if Bx[n,2]<min 
            Bx[n,2]=Ax[n,2]
        end
        if Bx[n,2]>max 
            Bx[n,2]=Ax[n,2]
        end
        #println("NBB1 ",Bx[n,1])
        #println("NBB2 ",Bx[n,2])
        
    end
    return(Bx)
end


function comparea(Ax,min,max)
    #computes the area of plots
    nr=size(Ax)[1]
    #points = [Point2(rand()*4+1, rand()*4+1) for _ in 1:nr]
    points=[Point2(Ax[i,1],Ax[i,2]) for i in 1:nr]
    #println(points')
    #println("min ",min,"max ",max)
    #println("min ",Point2(min,min),"max ",Point2(max,max))
    #outer rectangle
    xsquare=Rectangle(Point2(min,min), Point2(max,max))
    
    tess=voronoicells(points,xsquare)
    scatter(points, markersize = 6, label = "generators")
    annotate!([(points[n][1] + 0.02, points[n][2] + 0.03, Plots.text(n)) for n in 1:nr])
    zp=plot!(tess, legend = :topleft)
    display(zp)
    area=zeros(nr)
    Pmat=Vector{Vector{Point2{Float64}}}(undef,nr)
    for i=1:nr
        local P=tess.Cells[i]
        local P=cat(P,[P[1]],dims=1)
        Pmat[i]=P
        area[i]=0
        for j=1:(size(P)[1]-1) 
            area[i]=area[i]+P[j][1]*P[j+1][2]-P[j][2]*P[j+1][1]
        end
        area[i]=abs(area[i])/2
    end
   return area, Pmat
end


function xwalkfind(Ax,Pmat);
    #creates crosswalk of neighboring plots
    nplots=size(Pmat)[1]
    xwalk=zeros(nplots)
    Pmat2=Vector{Any}(undef,nplots)
    
    for j=1:nplots
        P=Vector{Any}(undef,nplots)
        for i=1:nplots
            k1=size(Pmat[j])[1]   
            P=[zeros(2) for _ in 1:k1]
            for k=1:k1
                P[k][1]=Pmat[j][k][1]
                P[k][2]=Pmat[j][k][2]
            end
                
            if PolygonOps.inpolygon(Ax[i,:],P)!=0 
                xwalk[i]=j 
            end
        end
        #println(Pmat[j])
        #println(P)
        Pmat2[j]=P
        #Pmat2[j]=[[P[k][1] P[k][2]] for k=1:size(P)[1] ]
      
    end
    return xwalk,Pmat2
end




function dist(x,y,delta)
    1/(1+exp(((x[1]-y[1])^2+(x[2]-y[2])^2)/delta))
end
function f(x)
    real((x.>0).*complex.(x).^beta).+(x.<0).*0
end

function inpoly(Pmat2)
    m=0
    for k=1:nplots
        for x=min:step:max
            for y=min:step:max
                m=m+1
                inpoly[m]=PolygonOps.inpolygon([x,y],Pmat2[k])
                if inpoly[m]
        if inpoly==1 | inpoly==-1 
            inpoly2[m]=k
        end
    end
    plotid=[i for i=1:nplots]
    inpoly3=inpoly2.==plotid
    return inpoly3
end

function weights()
    grid=[[xgrid[i] xgrid[j]] for i=1:n for  j=1:n]
    
    #@time distgrid=reshape([dist(grid[i],grid[j],delta) for i=1:n2 for j=1:n2],n2,n2)
    #@time distgrid=[dist(grid[i],grid[j],delta) for i=1:n2,  j=1:n2]
    distgrid=SharedArray{Float64}(n2,n2)
    @sync @distributed for i=1:n2
        for j=1:n2
            distgrid[i,j]=dist(grid[i],grid[j],delta)
        end
    end
    weights=distgrid./sum(distgrid,dims=2)
return weights
end

function agprod(weights,inputs,add)
    aginputs=(inputs.+alpha.*(weights*inputs))*add
    plotout=f.(aginputs)*add
    return aginputs,plotout
end

function self(inputi,kin,)
    for k=1:K
        agprod(weights,inputs)
