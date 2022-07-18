using Distributed
@everywhere using Delaunay, VoronoiCells, GeometryBasics, Random, PolygonOps, Plots, Optim
@everywhere using DataFrames, Econometrics, StatsBase, Distributed, SharedArrays

Random.seed!(1001)
#addprocs(3)
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





function aggproc(Pmat2,input,delta,nmat,min,step,max,maxin,inpoly)
    #computes aggregate production using a matrix of girid points accounting for inputs on other gridpoints and inverse weighted by distance
    
    nplots=size(input)[1]
    ngrid::Int64=(max-min)/step+1
    #println("ngrid ",ngrid)
    #println("nplots ", nplots)
    #=
    fz=zeros(ngrid,ngrid) #rupees per acre of inputs
    fz1=Array{Int16}(undef,ngrid,ngrid) #describes plot for each cell
    fz2=zeros(ngrid,ngrid) #inputs per acre among neighbors
    fz3=zeros(ngrid,ngrid) #fraction of spillover from self
    dsave=zeros(ngrid,ngrid)
    fz4=Vector{Any}(undef,nplots) #cell ids in plot
    fz5=zeros(nplots) #cells in plot
    =#
    fz=SharedArray{Float64}(ngrid,ngrid)
    fz1=SharedArray{Int16}(ngrid,ngrid)
    fz2=SharedArray{Float64}(ngrid,ngrid)
    fz3=SharedArray{Float64}(ngrid,ngrid)
    dsave=SharedArray{Float64}(ngrid,ngrid)
    fz4=SharedArray{Int64}(nplots,maxin,2)
    fz5=SharedArray{Float64}(nplots)
    xmat=SharedArray{Float64}(ngrid^2)
    ymat=SharedArray{Float64}(ngrid^2)
    imat=SharedArray{Int64}(ngrid^2)
    jmat=SharedArray{Int64}(ngrid^2)
    
   convert(SharedArray,inpoly)

   z=0
   for i=1:ngrid
        for j=1:ngrid
            z=z+1
            imat[z]=i
            jmat[z]=j
            xmat[z]=min+step*(i-1)
            ymat[z]=min+step*(j-1)
            
        end
    end
    @sync @distributed for k=1:nplots
        #println("Working on k=$k")
        flag=0
        tar=0
        i=0
        for x=min:step:max
            i=i+1
            j=0
            for y=min:step:max
                j=j+1                 
                #print(PointInPolygon([x,y],Pmat[k]));
                #inp=PolygonOps.inpolygon([x,y],Pmat2[k])
                inp=inpoly[i,j,k]
                #println("input ",inp)
                if (inp==1) | (inp==-1)
                    #println([i,j] ," in ",input[k] )
                    fz[i,j]=input[k]
                    fz1[i,j]=trunc(Int,k)
                    tar=1+tar
                    fz4[k,tar,1]=i
                    fz4[k,tar,2]=j
                    fz4[k,maxin,1]=tar
                end                         
            end
        end
        fz5[k]=tar
    end
    

    @sync @distributed for z=1:(ngrid^2)
        x=xmat[z]
        y=ymat[z]
        i=imat[z]
        j=jmat[z]
        
        
        i1=0
        sd=0
        sw=0
        fin=0
        
        for x1=min:step:max
            i1=i1+1
            j1=0
            for y1=min:step:max
                j1=j1+1
                d=1/(1+exp(sqrt((x1-x)^2+(y1-y)^2)/delta))
                #println(nmat)
                #if nmat[fz1[i,j],fz1[i1,j1]]==1 
                #    d=1
                #else
                #    d=0
                #end
                if fz1[i,j]==fz1[i1,j1]
                    dplot=0   #cell is in same plot
                #    d=1
                else
                    dplot=1  #cell not in same plot
                end
                #if (i==41) & (j==41) & (d==1) & (dplot==1)
                #    println("input ",fz[i1,j1])
                #end
                sd=(d*fz[i1,j1]*dplot+sd)
                sw=(d*dplot+sw)
                fin=(d*(1-dplot)+fin)
                
                if (i==20) & (j==20) 
                    #println(x," ",y," ",x1," ", y1)
                    dsave[i1,j1]=d 
                end
            end
        end
        fz2[i,j]=sd/sw   #inputs per acre
        fz3[i,j]=fin/(fin+sw)  #fraction spillover in 
    end
    return fz,fz1,fz2,fz3,fz4,fz5,dsave
end

function disttoplot(Ax)
    nplots=size(Ax)[1]
    distance=Ax[:,1].^2 .+ Ax[:,2].^2
    return distance
end
function f(x)
    #output of a plot
    #x is rupees per acre
    if x>0 
        y=x^(1/2)
    else
        y=-1000*x^2
    end
    return y
end

function inputmax(z,z2,z3,z4,z5,B,N,idplot,area,alpha,flag,starting,maxin,distance,bdist,cellsperacre)
    #determines inputs that maximize output for a farmer
    nplots=size(z4)[1]
    K=size(B)[1]
    pz=1  
    prplotv=zeros(nplots)
    prplotvso=zeros(nplots)
 
    function profit(x,xbar,wt)
        #output depends on own input and spillovers
        #x rupees per acre
        y=f(x+alpha*(wt*x+(1-wt)*xbar))
        #println(y," ",x," ",xbar)
        return y
    end
    function sumpro(x,k)
        #add up output 
        #x rupees per acre

        pr=0
        tinput=0
        
        #println("x ",x)
        for n=1:N[k]
            i=idplot[k,n]
            #println(" i ",i)
            #println(" x ",x)
            #println(area[i])
            #tinput=tinput+x[n]*area[i]
            tinput=tinput+x[n]*z5[i]/cellsperacre*exp(bdist*distance[i])
        end
        #i=idplot[k,N[k]]    
        #xlast=(B[k]-tinput)/(z5[i]/cellsperacre)
        #x1=vcat(x[:],[xlast])
        x1=x
        #println(x1)
        #println("x1 ",x1)
        for n=1:N[k]
            i=idplot[k,n]
            #println("x1 ",x1[n]," ",n)
            #println(i)
            prplot=0
            #println("z4 ",size(z4[i])[1])
            for m=1:z4[i,maxin,1]
                #println(" i ",i," m ", m)
                iwt=z4[i,m,:]
                #println(iwt)
                wt=z3[iwt[1],iwt[2]]
                #println(wt)
                z2wt=z2[iwt[1],iwt[2]]
                #println(z2wt)
                #println(x1[n]," ",z2wt," ",wt)
                
                prplot=prplot+profit(x1[n],z2wt,wt)/cellsperacre
            end
            prplotv[i]=prplot
            pr=pr+prplot
        end
        #println("tinput ",tinput," B ", B[k])
        pr=-(pr-1000000*((tinput-B[k])>0)*(tinput-B[k])^2)
        return pr
    end
    
    function totprof(x,kin,input)
        #calculates total profits across all farmers
        #println("N ",N[kin])
        count=count+1
        if count==1000
            println("1000")
            count=0
        end
        
        prplot=0
        #prplotv=Array{Float64}(undef,nplots)
        for k=1:K
            tinput=0
            prplot1=0
            for n=1:N[k]
                i=idplot[k,n]
                for m=1:z4[i,maxin,1]
                    iwt=z4[i,m,:]
                    #println(iwt)
                    wt=z3[iwt[1],iwt[2]]
                    #println(wt)
                    z2wt=z2[iwt[1],iwt[2]]                
                   
                    if k==kin
                     #println("x ",x,"i ",i,"kin ",kin)   
                    prplot1=prplot1+profit(x[n],z2wt,wt)/cellsperacre 
                    tinput=tinput+x[n]/cellsperacre*exp(bdist*distance[i])
                    else  
                    #println("input ",input[i])
                    prplot1=prplot1+profit(input[i],z2wt,wt)/cellsperacre
                    #tinput=tinput+input[i]/cellsperacre*exp(bdist*distance[i])
                    end  
                end
                #prplotv[i]=prplot1
                prplot=prplot+prplot1
            end
            prplot=prplot-1000000*((tinput-B[k])>0)*(tinput-B[k])^2
            #println("tinput ",tinput, " B ", B[k])
        end
        #println("pr profit ", prplot)
        prplot=-prplot
        return prplot
    end

    if flag==1
        input1=zeros(nplots)
        solvec=zeros(nplots)
        result=zeros(nplots)
        prplotv=zeros(nplots)
        for k=1:K

            #println("k ",k)
            #if N[k]==1
            #    id=idplot[k,1]
            #    result[id]=B[k]/z5[id]
            #else
                #println(starting)
                #println(idplot[k,:])
                startk=starting[[idplot[k,j] for j=1:N[k]]]
                #println(startk)
                res=optimize(b->sumpro(b,k),startk)
                #res=optimize(b->sumpro(b,k),ones(N[k]))
                res2=Optim.minimizer(res)
                
                #println("res ",res)
                #println("res2 ",res2)
                sx=0
                for j=1:N[k]
                    id=idplot[k,j]
                    result[id]=res2[j]
                    sx=sx+res2[j]*z5[id]
                end
                #id=idplot[k,N[k]]
                #result[id]=(B[k]-sx)/z5[id]
            #end

        
            
        
        end
    else
        count=0
        result=zeros(nplots)
        prplotv=zeros(nplots)
        for k=1:K
            res=optimize(b->totprof(b,k,z),ones(N[k]))
            res2=Optim.minimizer(res)
            res3=Optim.minimum(res)
            #println("plot v ",prplotv)
            #println(" optimal ", result)
            for j=1:N[k]
                id=idplot[k,j]
                result[id]=res2[j]
                prplotv[id]=-res3
                #println("profit ",res3)
            end
        end

    end
    
    #result=result.*tz/tarea
   
    return result,prplotv
end

function inputcomax(nplots,Btot,K,Pmat2,distmetric,nmat,min,step,max,maxin,inpoly,idplot,N,spill,cellsperacre,distance,bdist)
    function profit(x,xbar,wt)
        #output depends on own input and spillovers
        #x rupees per acre
        y=f(x+spill*(wt*x+(1-wt)*xbar))
        #println(y," ",x," ",xbar)
        return y
    end
    function totprof2(x,kin,input,lambda,Btot,z2,z3,z4)
        #calculates total profits across all farmers
        #println("N ",N[kin])
        count=count+1
        if count==1000
            println("1000")
            count=0
        end
        prplot1=0
        prplot=0
        tinput=0
        #prplotv=Array{Float64}(undef,nplots)
        for k=1:K            
            for n=1:N[k]
                i=idplot[k,n]
                for m=1:z4[i,maxin,1]
                    iwt=z4[i,m,:]
                    #println(iwt)
                    wt=z3[iwt[1],iwt[2]]
                    #println(wt)
                    z2wt=z2[iwt[1],iwt[2]]                
                   
                    if k==kin
                     #println("x ",x,"i ",i,"kin ",kin)   
                    prplot1=prplot1+profit(x[n],z2wt,wt)/cellsperacre  
                    tinput=tinput+x[n]/cellsperacre*exp(bdist*distance[i])
                    else  
                    #println("input ",input[i])
                    prplot1=prplot1+profit(input[i],z2wt,wt)/cellsperacre
                    tinput=tinput+input[i]/cellsperacre*exp(bdist*distance[i])
                    end  
                end
               
                
            end
            
            #println("tinput ",tinput, " B ", B[k])
        end
        prplot=prplot1+exp(lambda)*(Btot-tinput)
        #println("bc ",Btot," ",tinput)
        #println("pr profit ", prplot)
        prplot=-prplot
        #println("prplot ",kin," ",prplot)
        return prplot
    end

    function gn(lambda1,x,Btot)
   


        lambda=lambda1[1]
        del2=100
        xold=x[:]
       
        while del2>.0001
            z0,z1,z2,z3,z4,z5,dsave=aggproc(Pmat2,x,distmetric,nmat,min,step,max,maxin,inpoly)
            xold=x[:]
            for kin=1:K
                istart=idplot[kin,1]
                iend=istart+N[kin]-1
                xi=x[istart:iend]
                #println("istart",istart," ",iend)
                #println("kin ",kin)
                res=optimize(z->totprof2(z,kin,x,lambda,Btot,z2,z3,z4),xi)
                #println(res)
                x[istart:iend]=Optim.minimizer(res)
                #println(xold[istart:iend])
                #println(x[istart:iend])
            end
            del=x-xold
            del2=(del'del)/nplots
            println("del2 ", del2)
            
            #println("del2 ",del2)
            #println("optim ",Optim.minimum(res))
        end
        
        x1=x[:]
        #println("sumxi",sum(x1))
        z=Optim.minimum(res)
        #println("minimand ",z)
        return -z
    end
    count=0
    x1=zeros(nplots)
    lambda1=ones(1)*4
    res=optimize(z->gn(z,x1,Btot),lambda1,Optim.Options(iterations=30,g_abstol=.1))
    println(res)
    z=Optim.minimum(res)
    lambda1=Optim.minimizer(res)
    println("lambda1 ",lambda1)
    println("x1 ",x1)
    return lambda1,z,x1
end





function mnneighbors(input,area,nplots,K,N,idplot,nmat,B,prplotv,distance)
    function ln(x)
        log.(x)
    end
    #caluclates input of neighbors and runs regressions
    sparseinput=nmat.*ln.(input)'
    sparsedist=nmat.*ln.(distance)'
    mninput=[sum(sparseinput[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    mndist=[sum(sparsedist[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    tarea=zeros(K)
    herf=zeros(K)
    for k=1:K
        tarea[k]=(sum(area[idplot[k,n]] for n=1:N[k]))
        herf[k]=sum((area[idplot[k,n]]/tarea[k])^2 for n=1:N[k])
    end
    plottarea=zeros(nplots)
    plotherf=zeros(nplots)
    plotB=zeros(nplots)
    owner=zeros(nplots)
    for k=1:K
        for i=1:N[k]
            plottarea[idplot[k,i]]=tarea[k] 
            plotherf[idplot[k,i]]=herf[k] 
            plotB[idplot[k,i]]=B[k] 
            owner[idplot[k,i]]=k
        end    
    end
    println(plottarea)
    #println(plotB)
    sparseplottarea=nmat.*ln(plottarea)'
    mnareaf=[sum(sparseplottarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    sparseplotarea=nmat.*ln(area)'
    mnareap=[sum(sparseplotarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    sparseplotherf=nmat.*plotherf'
    mnherf=[sum(sparseplotherf[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    global data=DataFrame(input=input,area=area,mninput=mninput,mnareaf=mnareaf,mnareap=mnareap,plotB=plotB,plottarea=plottarea,owner=owner,prplotv=prplotv,distance=distance,mndist=mndist,plotherf=plotherf,mnherf=mnherf)
    mnarea=sum(Matrix(data[:,["area", "mnareap","input","prplotv"]]))/nplots
    areacov=cov(Matrix(data[:,["area", "mnareap","input","prplotv"]]))
    bhat1=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+distance+mndist+mnherf+absorb(owner)),data)
    bhat1a=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+plottarea +distance+mndist+plotherf+mnherf ),data)
    bhat2=fit(EconometricModel, @formula(input ~ area +distance+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
    bhat3=fit(EconometricModel, @formula(prplotv ~ area + distance+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
    bhat4=fit(EconometricModel, @formula(input ~ area +plottarea+ distance +plotherf+(mninput~mnareaf+mnareap+mndist+mnherf)),data)
    bhat5=fit(EconometricModel, @formula(prplotv ~ area +plottarea+distance+plotherf+ (mninput~mnareaf+mnareap+mndist+mnherf)),data)
    println(bhat1)
    println(bhat1a)
    println(bhat2)
    println(bhat3)
    println(bhat4)
    println(bhat5)
    return bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov
end

function run()
    #loops through procedures
    nfarmers=1
    min=1
    max::Int64=trunc(2*sqrt(nfarmers))
    #println(max)
    step=.1 #.1 #distance between grid points4
    plotdist=4:4 #1:3 #distriubtion of Plots
    Bdist=1:1 #distribution of endowments
    distmetric=1 #how fast to spillovers fall off
    spill=0#-.1 #spillover coefficient
    maxin=1200  #max number of neighbors 
    bdist=0#.0001 #coefficient on cost of distance from plot 
    scorr=0  #spatial correlation of plot sizes
    cellsperacre=((max-min)/step+1)^2/(max-min)^2
    #println("xxx" ,cellsperacre)
    #first create plots 
    K,nplots,Ax,N,idplot,B,input=createdata(nfarmers,plotdist,Bdist,min,step,max)

    Ax2=spacecorr(Ax,scorr,min,step,max)
    Ax=Ax2
    Btot=sum(B)
    println("Btot ",Btot)
    distance=disttoplot(Ax)
    #println(nplots)
    #println(Ax)
    #println(" N ",N)
    nmat=findneighbor(Ax)
    #println(nmat)
    #println("input ",input)
    area,Pmat=comparea(Ax,min,max)
    xwalk,Pmat2=xwalkfind(Ax,Pmat)
    svg5=0
    svg4=0
    svg3=0
    svg2=0
    svg=0
    dsave=0
    prplotv=zeros(nplots)
    prplotvso=zeros(nplots)        
    inpoly=[PolygonOps.inpolygon([x,y],Pmat2[k]) for x=min:step:max, y=min:step:max, k=1:nplots]
    #println("size ",size(inpoly))
    #println("input ", input)
    
    del=100
    while del>.0001 
        println("nash sqdif ", del)
        #println("input ",input)
        gz,gz1,gz2,gz3,gz4,gz5,dsave=aggproc(Pmat2,input,distmetric,nmat,min,step,max,maxin,inpoly)
        svg5=gz5
        svg4=gz4
        svg3=gz3
        svg2=gz2
        svg=gz
        ninput,prplotv=inputmax(gz,gz2,gz3,gz4,gz5,B,N,idplot,area,spill,1,input,maxin,distance,bdist,cellsperacre)
        del=(input.-ninput)'*(input.-ninput)/nplots
        input=.8*input+.2*ninput
    end

    #println(" iput ",input)  
    #println(input.*area)
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5=0,0,0,0,0,0
    #bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5=mnneighbors(input,area,nplots,K,N,idplot,nmat,B,prplotv,distance)
    del=0
    inputso=input[:]
    del=100
    while del>.000001
        println("coord sqdif ", del)
        #println("input ",input)
        gz,gz1,gz2,gz3,gz4,gz5,dsave=aggproc(Pmat2,inputso,distmetric,nmat,min,step,max,maxin,inpoly)
        svg5=gz5
        svg4=gz4
        svg3=gz3
        svg2=gz2
        svg=gz
        ninputso,prplotvso=inputmax(svg,svg2,svg3,svg4,svg5,B,N,idplot,area,spill,0,inputso,maxin,distance,bdist,cellsperacre)
      
        #println("in ",inputso, "nin ", ninputso)
        del=(inputso.-ninputso)'*(inputso.-ninputso)/nplots
        inputso=.8*inputso+.2*ninputso
    end
    
    
    
    lambda,prplotvco,inputco=inputcomax(nplots,Btot,K,Pmat2,distmetric,nmat,min,step,max,maxin,inpoly,idplot,N,spill,cellsperacre,distance,bdist)
    return bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,area,input,inputso,inputco,prplotv,prplotvso,prplotvco,N,B
end
bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,area,input,inputso,inputco,prplotv,prplotvso,prplotvco,N,B=run()


