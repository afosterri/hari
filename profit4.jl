using Distributed
@everywhere using Delaunay, VoronoiCells, GeometryBasics, Random, PolygonOps, Plots, Optim, Roots
@everywhere using DataFrames, Econometrics, StatsBase, Distributed, SharedArrays, ForwardDiff,NLsolve,LinearAlgebra
Random.seed!(1001)

function run()

    function disttoplot(Ax)
        nplots=size(Ax)[1]
        distance=sqrt.(Ax[:,1].^2 .+ Ax[:,2].^2)
        return distance
    end
    function curvy(beta,delta,y1,Ax)
        function f1(_Z)
           
            out=(-beta*delta*y1*z - beta*y1*z + _Z*delta + cos(beta*y1)*z - cos(_Z) + _Z - z + 1)
            #println(out)
            return out
        end
        z1=zeros(nplots,2)
        z=0
    for i=1:nplots
    z=Ax[i,1]
    
    z1[i,1]=find_zero(f1,y1/2)/beta
    z=Ax[i,2]
    z1[i,2]=find_zero(f1,y1/2)/beta
    #println(z1[i,:])
    end
    return z1 
    end

    function createdata(K,NN,BB,min,step,max) 
        #creates plots for K farmers, NN is distribution of number of plots, BB is range of input dist, min max is range of locations  
        N=rand(NN,K)
        #println("N ",N)
        #N=[2,2,2,2,1]    
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
        #B[5]=B[1]/2
        #println(" B ",B)
        #B[1]=2
        #B[2]=2
        #B[2]=4
        Ax1=rand(nplots,2)
        Ax=curvy(1.5,.1,max-min,Ax1).+min
        #=
        min=1
        max=5
        a=1.833333333333
        b=3.5
        c=5.16666666667
        Ax=[a a; a b; a c; b a; b b ; b c;c a ; c b ;c c]
        #Ax=[1 1 ; 1 3.5; 1 2 1.75 ; 2 2 ]
        =#
        
        inputs=ones(nplots)
        return K,nplots,Ax,N,idplot,B,inputs
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
    function spacecorr(A1x,sc,min,step,max)
        #perturbs locations so plots are more or less close depending on sc
        Ax=A1x[:,:]            
        Bx=A1x[:,:]
        
        for i=1:5
        
        #println("min ",min," ", max)
        
        wt=exp.(-((Ax[:,1].-Ax[:,1]').^2+(Ax[:,2].-Ax[:,2]').^2).^(1/2))
                dx=wt*Ax[:,1]./sum(wt,dims=2)-Ax[:,1]
        dy=wt*Ax[:,2]./sum(wt,dims=2)-Ax[:,2]
        
            for n=1:nplots
                if mod(n,2)==0
                Bx[n,1]=Ax[n,1]-dx[n]*sc
                Bx[n,2]=Ax[n,2]-dy[n]*sc
                else
                Bx[n,1]=Ax[n,1]+dx[n]*sc
                Bx[n,2]=Ax[n,2]+dy[n]*sc
                end
                #println(Bx[n,1]," ",Ax[n,1])
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
        Ax=Bx[:,:]
        #println("BX ",Bx)
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
        
        area=zeros(nr)
        center=zeros(nr,2)
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
            ta=zeros(2)
            twt=0
            for j=2:(size(P)[1]-1)
                ce=1/3*(P[j]+P[j+1]+P[1])
                a=[(P[j]-P[1])';(P[j+1]-P[1])']
                wt=a[1,1]*a[2,2]-a[1,2]*a[2,1]
                ta=ta+wt*ce
                twt=twt+wt
            end

            center[i,:]=ta/twt



            
        end

        points=[Point2(center[i,1],center[i,2]) for i in 1:nr]
        scatter(points, markersize = 6, label = "generators")
        annotate!([(points[n][1] + 0.02, points[n][2] + 0.03, Plots.text(n)) for n in 1:nr])
        zp=plot!(tess, legend = :topleft)
        display(zp)
        return area, Pmat,center
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
        
        
        z=1/(1+exp((((x[1]-y[1])^2+(x[2]-y[2])^2)^(1/2)/delta)-1))
        
        #z=0
        #if (sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)<1) 
        #    z=1
        #
        return z
    end
    function f(x)
        z=10*(x-beta*x^2)
        #=
        if x>0 
            z=10*x^beta
        else
            z=-1000*x^2
        end
        =#
        #println(x)
        #z=(real((x.>=0).*(complex.(x)).^beta)).+((x.<0).*(-10000))
        #println(z)
        #if sum(x.<0)>0 
        #println("prod",x," ",(z))
        #end
        return z
    end
    function g(x)
        z=(x+xi*x^2)
        #=
        if x>0 
            z=x^xi
            #z=x
        else
            z=0
            #z=-1000*x^2
            #z=x
        end
        =#
        #println(x)
        #z=(real((x.>=0).*(complex.(x)).^beta)).+((x.<0).*(-10000))
        #println(z)
        #if sum(x.<0)>0 
        #println("prod",x," ",(z))
        #end
        return z
    end
    function gprime(x)
        z=2*xi*x
        #=
        if x>0 
            z=xi*x^(xi-1)
            #z=x
        else
            z=0
            #z=-1000*x^2
            #z=x
        end
        =#
        #println(x)
        #z=(real((x.>=0).*(complex.(x)).^beta)).+((x.<0).*(-10000))
        #println(z)
        #if sum(x.<0)>0 
        #println("prod",x," ",(z))
        #end
        return z
    end
    function inpoly(Pmat2)
        
        inpoly2=zeros(n2)
        plotfarmer=[false for i=1:nplots , j=1:nfarmers]
      
        m=0
        for m=1:nplots
            for k=1:nplots
                inpoly=PolygonOps.inpolygon(Ax[m,1:2],Pmat2[k])                    
                if abs(inpoly)==1 
                    inpoly2[m]=k
                end
            end         
        end
        plotid=[i for i=1:nplots]
        for i=1:nplots
            for k=1:K
                for j=1:N[k]
                    if idplot[k,j]==plotid[i]
                        plotfarmer[i,k]=1
                    end
                end
            end
        end
        inplot=inpoly2.==plotid'
        infarmer=inplot*plotfarmer
        return inplot,infarmer,plotfarmer
    end
    function weights(centers)
        #grid=[[xgrid[i] xgrid[j]] for i=1:n for  j=1:n]
        grid=[centers[i,:] for i=1:nplots]
        n2=nplots
        #@time distgrid=reshape([dist(grid[i],grid[j],delta) for i=1:n2 for j=1:n2],n2,n2)
        #@time distgrid=[dist(grid[i],grid[j],delta) for i=1:n2,  j=1:n2]
        distgrid=SharedArray{Float64}(n2,n2)
        #println(centers)
        for i=1:n2
            for j=1:n2

                distgrid[i,j]=dist(grid[i],grid[j],delta)*area[j]
            end
        end
        #distgrid=(distgrid.>.5).*distgrid
        #den=maximum(sum(distgrid,dims=2))
        den=(sum(distgrid,dims=2))        
        #println(" dist grid ", distgrid)

        weights=distgrid./den
        #println(" weights ", weights)
        return weights
    end
   
    function agprod(inputs)
        inputsgrid=inputs[:,:]
        aginputs=(inputsgrid.+spill.*(weights*inputsgrid))          
        #println(" weights ", inputsgrid)
        plotout=(f.(aginputs).*area)'
        plotout=convert(Array{Float64},plotout)        
        
        price=exp.(gdist*distance).*area        
        bc=((price.*inputs)'*infarmer)'           
        bc=(bc.-B)
        #println("bc ", bc)
        cost=g.(bc)'
        return plotout,cost
    end
    

    function dagprod(inputs,add3)
        #println("hello2")
        dchange=.0000001
        
        eye=I(nplots)
        eye=eye[:,add3]
        #println("add3 ",add3)
        inputs1=inputs.+dchange.*eye                     
        p1,cost1=agprod(inputs)
        p2,cost2=agprod(inputs1)      
        #println("p1 ",p1)    
        #println("cost1 ",cost1)    
        dagprod2=(p2.-p1)./dchange        
        dcost=(cost2.-cost1)./dchange
        #println("transfer ",transfer)
        #println("inputs ",inputs)
        #println("dcost1 ", dcost)        
        return dagprod2,dcost
    end
    

    function indmaxso1(inputs,mu)
                
        function nlf!(Gfun,inplam)  
            #lambda=inplam[xnk+2]
            #c=inplam[xnk+1]                    
            xinputs3=xinputs2[:]
            xinputs3[xstart:xend]=inplam[1:xnk]        
 
            dplotout,dcost=dagprod(xinputs3,xstart:xend)   
            #println(" dplotout ",dplotout)
            #println(" dcost ",dcost)
            #println("xstart ",xstart)
            #println("dcost ",dcost)

            wtfarmer=(1-mu)*infarmer+mu*ones(nplots,nfarmers)            
            #println("wt farmer ",wtfarmer)
            dplotout=dplotout*wtfarmer
            #println(wtfarmercost)
            dp1=dplotout[:,xk]        
            pr1=dcost[:,xk]           
            deriv=dp1 .-pr1   #-1000*(inplam[1:xnk].<0).*inplam[1:xnk]
            #println("deriv ",deriv)
            Gfun[1:xnk]=deriv
        end
        
        xinputs1=SharedArray{Float64}(nplots)
        xinputs=inputs[:]
        xinputs2=xinputs[:]
        xinputs1[:]=xinputs
        xstart=0
        xend=0
        xk=0
        xnk=0
        del=100
        #lamvec1=ones(K)
        #lamvec2=lamvec1[:]
        while del>.00001
        @sync @distributed for k=1:K            
           # for k=1:K     
                xstart=idplot[k,1]
                xend=xstart+N[k]-1
                xk=k
                xnk=N[k]
                start3=xinputs[xstart:xend]                         
                z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)     
                xinputs1[xstart:xend]=z1.zero            
            end
            #lamvec1=lamvec2[:]        
            xinputs2=xinputs1[:]
            del1=xinputs2-xinputs
            xinputs=.8*xinputs2+.2*xinputs
            del=del1'del1/nplots  
            println("del ",del)
        end
        inpt=0
        #println(xinputs')   
        price=exp.(gdist*distance).*area
        
        
        prof,cost=agprod(xinputs)
        prof1=(prof*infarmer)
        util=prof1-cost        
        return xinputs,prof,util
    end


    function indmaxco1(inputs,mu)
                
        function nlf!(Gfun,inplam)  
            #lambda=inplam[xnk+2]
            #c=inplam[xnk+1]        
                   
            xinputs3=xinputs2[:]
            xinputs3[xstart:xend]=inplam[1:xnk]
            
            #println("inputs ",xinputs3)
            price=exp.(gdist*distance).*area
            tinputs=xinputs3'*price/nfarmers
            tstock=sum(B)/nfarmers
            lambda=-gprime(tinputs-tstock)
            
            dplotout,dcost=dagprod(xinputs3,xstart:xend)               
            wtfarmer=(1-mu)*infarmer+mu*ones(nplots,nfarmers)            
            dplotout=dplotout*wtfarmer
            #println(wtfarmercost)
            dp1=dplotout[:,xnk]
         
            #println("dplotout ",dp1)        
            #println("dlambda ",lambda)   
            #println("price ",price)
            #println(size(dp1),size(lambda),size(price[xstart:xend]))
            deriv=dp1 .+lambda.*price[xstart:xend]   #-1000*(inplam[1:xnk].<0).*inplam[1:xnk]
            #println("deriv ",deriv) 
            Gfun[1:xnk]=deriv
        end
        
 


        xinputs1=SharedArray{Float64}(nplots)
        xinputs=inputs[:]
        xinputs2=xinputs[:]
        xinputs1[:]=xinputs
        xtransferlam=zeros(nfarmers+1)
        #xtransferlam[nfarmers+1]=1
        xstart=0
        xend=0
        xk=0
        xnk=0
        del=100
        #lamvec1=ones(K)
        #lamvec2=lamvec1[:]
        while del>.00001
        # @sync @distributed for k=1:K            
            for k=1:K     
                xstart=idplot[k,1]
                xend=xstart+N[k]-1
                xk=k
                xnk=N[k]
                start3=xinputs[xstart:xend]            
                z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)     
                xinputs2[xstart:xend]=z1.zero 
                #println("xinputs1 ",xinputs1)           
            end
            #xinputs2=xinputs1[:]
            del1=xinputs2-xinputs
            xinputs=.8*xinputs2+.2*xinputs
            del=del1'del1/nplots  
            println("del ",del)
        end
        inpt=0                  
        prof,cost=agprod(xinputs)
        prof1=(prof*infarmer)        
        price=exp.(gdist*distance).*area
        tinputs=xinputs'*price/nfarmers
        tstock=sum(B)/nfarmers
        cost=g(tinputs-tstock)/nfarmers
        println(" cost ",cost)
        
        util=prof1.-cost        
        return xinputs,prof,util
    end
   
    
    function mnneighbors(input,prplotv)
        function ln(x)
            (x.>=0).*real.(log.(complex.(x)))+(x.<0).*0
        end
        #caluclates input of neighbors and runs regressions
        #println(size(input),size(nmat))
        #mninput=log.(weights*input-Diagonal(diag(weights))*input)
        #mndist=log.(weights*distance-Diagonal(diag(weights))*distance)
        sparseinput=nmat.*(input)'
        sparsedist=nmat.*(distance)'
        mninput=ln.([sum(sparseinput[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots])
        mndist=ln.([sum(sparsedist[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots])
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
        #println(plottarea)
        #println(plotB)
        sparseplottarea=nmat.*ln(plottarea)'
        mnareaf=[sum(sparseplottarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        sparseplotarea=nmat.*ln(area)'
        mnareap=[sum(sparseplotarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        sparseplotherf=nmat.*plotherf'
        mnherf=[sum(sparseplotherf[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
    
        areavec=vec(area)
        prplotvec=vec(prplotv)

        global data=DataFrame(input=input,area=areavec,mninput=mninput,mnareaf=mnareaf,mnareap=mnareap,plotB=plotB,plottarea=plottarea,owner=owner,prplotv=prplotvec,distance=distance,mndist=mndist,plotherf=plotherf,mnherf=mnherf)
        mnarea=sum(Matrix(data[:,["area", "mnareap","input","prplotv"]]))/nplots
        areacov=cov(Matrix(data[:,["area", "mnareap","input","prplotv"]]))
        bhat0a=fit(EconometricModel, @formula(input~area+mnareap+absorb(owner)),data)
        bhat0=fit(EconometricModel, @formula(input~area+distance+mnareap+mnareaf+mndist+mnherf+absorb(owner)),data)
        bhat1=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+distance+mndist+mnherf+absorb(owner)),data)
        bhat1a=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+plottarea +distance+mndist+plotherf+mnherf ),data)
        bhat2=fit(EconometricModel, @formula(input ~ area +distance+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
        bhat3=fit(EconometricModel, @formula(prplotv ~ area + distance+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
        bhat4=fit(EconometricModel, @formula(input ~ area +plottarea+ distance +plotherf+(mninput~mnareaf+mnareap+mndist+mnherf)),data)
        bhat5=fit(EconometricModel, @formula(prplotv ~ area +plottarea+distance+plotherf+ (mninput~mnareaf+mnareap+mndist+mnherf)),data)
        print(bhat0a)
        print(bhat0)
        println(bhat1)
        println(bhat1a)
        println(bhat2)
        println(bhat3)
        println(bhat4)
        println(bhat5)
        return bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov
    end



    nfarmers=50#number farmers
    
    min=1
    max=trunc(2*sqrt(nfarmers))
    max=convert(Int64,max)
    max=6
    #println(max)
    step=.1 #.1 #distance between grid points4
    plotdist=2:2 #1:3 #distriubtion of Plots
    Bdist=10:10 #distribution of endowments
    delta=.75 #how fast to spillovers fall off
    spill=-.5#-.1 #spillover coefficient
    maxin=1200  #max number of neighbors 
    gdist=0.001#.0001 #coefficient on cost of distance from plot 
    scorr=-.5  #spatial correlation of plot sizes
    cellsperacre=((max-min)/step+1)^2/(max-min)^2
      
    #NN=2:2
    xgrid=[i for i=min:step:max]
   
    #n=size(xgrid)[1]
    #n2=n*n
    #beta=1/2
    #xi=2
    beta=.01
    xi=.05
    K,nplots,Ax,N,idplot,B,inputs=createdata(nfarmers,plotdist,Bdist,min,step,max)
    n=nplots
    n2=nplots
    tB=sum(B) 
    Ax2=Ax[:,:]
    #Ax2=spacecorr(Ax,scorr,min,step,max)
    println("Ax ",Ax)
    println("Ax2 ",Ax2)
    Ax=Ax2[:,:]
    Btot=sum(B)
    println("Btot ",Btot)
    distance=disttoplot(Ax)
    global distsv=distance
    global Axsv=Ax
    nmat=findneighbor(Ax)
    global nmatsv=nmat
    area,Pmat,centers=comparea(Ax,min,max)
    Ax=centers[:,:]
    xwalk,Pmat2=xwalkfind(Ax,Pmat)
    #println(Pmat2)
    inplot,infarmer,plotfarmer=inpoly(Pmat2)
    global areasv=area
    global inplot1=inplot
    global infarm1=infarmer
    global polotfarm1=plotfarmer
    disttogrid=inplot*distance
    #println(centers)
    weights=weights(centers)
    global weightsv=weights
    #println(weights)
    #println(n)
    #println("sum grid ",sum(inplot,dims=2))
    #println("inp",inplot)
    #println("inf",infarmer)
    @time inputsmasv,profitsmasv,utilmasv=indmaxso1(inputs,0)
    #time inputsmasv,profitsmasv=indmaxma(inputs)
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=0,0,0,0,0,0,0,0
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=mnneighbors(inputsmasv[1:nplots],profitsmasv)
    
    #println(profitsv)
    println("now so")
    @time inputssosv,profitssosv,utilsosv=indmaxso1(inputs,1)
    #time inputssosv,profitssosv=indmaxso(inputs)
    #println(profitsosv)
    #println(infarmer)
    #println(B)
    #inputscosv,profitscosv=indmaxco(inputs)
    @time inputscosv,profitscosv,utilcosv=indmaxco1(inputs,1)
    #inputscosv=zeros(nplots)
    #profitscosv=zeros(nplots)
    global wtsv=weights*inplot
    return(inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv[1:nplots],inputssosv[1:nplots],inputscosv[1:nplots],profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv)
end   

inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv,inputssosv,inputscosv,profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv=run()
global output=DataFrame(ima=inputsmasv,iso=inputssosv,ico=inputscosv,pma=vec(profitsmasv),pso=vec(profitssosv),pco=vec(profitscosv))


println(sum(utilmasv))
println(sum(utilsosv))
println(sum(utilcosv))

scatter(inputsmasv,inputscosv)
scatter!(inputsmasv,inputssosv)