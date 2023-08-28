using Distributed, XLSX, JLD

@everywhere using Delaunay, VoronoiCells, GeometryBasics, Random, PolygonOps, Plots, Optim, Roots, JLD2
@everywhere using DataFrames, Econometrics, StatsBase, Distributed, SharedArrays, ForwardDiff,NLsolve,LinearAlgebra


function run(deltagrid,spillgrid,maonly)

    function disttoplot(Ax)
        nplots=size(Ax)[1]
        distance=sqrt.(Ax[:,1].^2 .+ Ax[:,2].^2)
        return distance
    end
    function curvy(beta1,delta1,y1,Ax)
        function f1(_Z)
           
            out=(-beta1*delta1*y1*z - beta1*y1*z + _Z*delta1 + cos(beta1*y1)*z - cos(_Z) + _Z - z + 1)
            #println(out)
            return out
        end
        z1=zeros(nplots,2)
        z=0
        for i=1:nplots
        z=Ax[i,1]
        
        z1[i,1]=find_zero(f1,y1/2)/beta1
        z=Ax[i,2]
        z1[i,2]=find_zero(f1,y1/2)/beta1
        #println(z1[i,:])
        end
        return z1 
    end

    function createdata(K,NN,BB,min,step,max) 
        Random.seed!(1001)
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
        return area, Pmat,center,tess
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
    
    function dist(x,y,delta2)
        
        
        z=1/(1+exp((((x[1]-y[1])^2+(x[2]-y[2])^2)^(1/2)/delta2)-1))
        
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
   
    function areaweights(centers,tess)
    
    
        ngrid=convert(Int64,ceil((max-min)/dgrid)+1)
        mat=SharedArray{Float64}(ngrid,ngrid)
        ivec=[i for i=min:dgrid:max]
        jvec=[j for j=min:dgrid:max]
        @sync @distributed for i=1:ngrid
            x=ivec[i]
            for j=1:ngrid
                y=jvec[j]    
                for n=1:nplots
                    cell=[tess.Cells[n];[tess.Cells[n][1]]]         
                    inthere=PolygonOps.inpolygon([x,y],cell)
                    #println("flag ",inthere," ",x," ",y)
                    if inthere!=0
                        #println("flag2 ",inthere," ",x," ",y)
                        mat[i,j]=n
                        break
                    end                        
                end
        
            end
        end
        
        weights=SharedArray{Float64}(nplots,nplots)
        onev=ones(ngrid^2)
        
        vmat=vec(mat)
        #println("ivec ",ivec)
        #println("j vec" ,jvec)
        @sync @distributed for k=1:nplots
            dist=vec(sqrt.((ivec.-centers[k,1]).^2 .+ (jvec'.-centers[k,2]).^2))
           
            dwt=(1 ./ (1 .+ exp.(dist./delta)))*dgrid^2
            
            for l=1:nplots    
                weights[k,l]=sum(dwt[vmat.==l]) 
            end
        
        end
        weights=weights./sum(weights,dims=2)
        #println("centers ",center')
        #z1=plot(tess)
        #display(z1)
        #println("mat ",weights)
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
        avdist=zeros(K)
        for k=1:K
            tarea[k]=(sum(area[idplot[k,n]] for n=1:N[k]))
            herf[k]=sum((area[idplot[k,n]]/tarea[k])^2 for n=1:N[k])
            avdist[k]=sum(distance[idplot[k,n]] for n=1:N[k])/N[k]
        end
        plottarea=zeros(nplots)
        plotherf=zeros(nplots)
        plotB=zeros(nplots)
        owner=zeros(nplots)
        avdistp=zeros(nplots)
        for k=1:K
            for i=1:N[k]
                plottarea[idplot[k,i]]=tarea[k] 
                plotherf[idplot[k,i]]=herf[k] 
                plotB[idplot[k,i]]=B[k] 
                owner[idplot[k,i]]=k
                avdistp[idplot[k,i]]=avdist[k]
            end    
        end
        mnn=sum(nmat,dims=2)
        #println(plottarea)
        #println(plotB)
        sparseplottarea=nmat.*ln(plottarea)'
        mnareaf=[sum(sparseplottarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        sparseplotarea=nmat.*ln(area)'
        mnareap=[sum(sparseplotarea[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        sparseplotherf=nmat.*plotherf'
        mnherf=[sum(sparseplotherf[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        println(size(avdistp));
        println(size(nmat));
        sparseplotdist=nmat.*avdistp'
        mnavdist=[sum(sparseplotdist[i,j] for j=1:nplots)/sum(nmat[i,j] for j=1:nplots) for i=1:nplots]
        areavec=vec(area)
        prplotvec=vec(prplotv)
        mnn=vec(mnn)
        xinput=ln(input)
        xareavec=ln(areavec)
        xareavec2=xareavec.^2
        xtarea=ln(plottarea)

        xprplotvec=ln(prplotvec)-xareavec
        global data=DataFrame(input=xinput,area=xareavec,tarea=xtarea,avdist=avdistp,mnavdist=mnavdist,area2=xareavec2,mninput=mninput,mnareaf=mnareaf,mnareap=mnareap,plotB=plotB,plottarea=plottarea,owner=owner,prplotv=xprplotvec,distance=distance,mndist=mndist,plotherf=plotherf,mnherf=mnherf,mnn=mnn)
        mnarea=sum(Matrix(data[:,["area", "mnareap","input","prplotv"]]),dims=1)/nplots
        areacov=cov(Matrix(data[:,["area", "mnareap","input","prplotv"]]))
        #global meansdata=[xinput,prplotv,area,area2,mnn,distance,mnareaf,lnmnsa,mnsher,mndist,mnavdist]
        bhat1=0
        bhat1a=0
        bhat4=0
        bhat5=0
        #bhat0a=fit(EconometricModel, @formula(input~area+area2+mnn+mnareap+absorb(owner)),data)
        #bhat0=fit(EconometricModel, @formula(input~area+area2+distance+mnn+mnareap+mnareaf+mndist+mnherf+absorb(owner)),data)
        #bhat1=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+area2+mnn+distance+mndist+mnherf+absorb(owner)),data)
        #bhat1a=fit(EconometricModel, @formula(mninput~mnareaf+mnareap+area+area2+mnn+plottarea +distance+mndist+plotherf+mnherf ),data)
        #bhat2=fit(EconometricModel, @formula(input ~ area+area2 +distance+mnn+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
        #bhat3=fit(EconometricModel, @formula(prplotv ~ area +area2+ distance+mnn+ absorb(owner)+ (mninput~mnareaf+mnareap+mndist + mnherf)),data)
        bhat2=fit(EconometricModel, @formula(input ~ (mninput~mnareaf+mnareap+mndist+mnavdist+mnherf)+area+ distance  +mnn+absorb(owner)),data)
        bhat3=fit(EconometricModel, @formula(prplotv ~ (mninput~mnareaf+mnareap+mndist+mnavdist+mnherf)+area+distance+ mnn+absorb(owner)),data)
        
        bhat4=fit(EconometricModel, @formula(input ~ (mninput~mnareaf+mnareap+mndist+mnavdist+mnherf)+area +tarea+plotherf+ distance +avdist +mnn),data)
        bhat5=fit(EconometricModel, @formula(prplotv ~ (mninput~mnareaf+mnareap+mndist+mnavdist+mnherf)+area+tarea +plotherf+distance+avdist+ mnn),data)
        #print(bhat0a)
        #print(bhat0)
        #println(bhat1)
        #println(bhat1a)
        println(bhat2)
        println(bhat3)
        #println(bhat4)
        #println(bhat5)
        
        return bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov
    end

  
    function heat(inputmat)
        function inputv(input,x,y)
            lev=0 
            for i=1:nplots
                cell=[tesssv.Cells[i];[tesssv.Cells[i][1]]] 
                incheck=PolygonOps.inpolygon([x,y],cell)
              
                if incheck!=0
                    lev=log(input[i])
                    break
                end
            end
            return lev
        end
        x=min:dgrid:max  
        z=inputmat
        zmax=ceil(log(maximum(z)))
        println("zmax ",zmax)
        heatmap(x,x,(x1,y1)->inputv(inputmat,x1,y1), c=:thermal, clims=(0,zmax))
        plot!(tesssv)
    end

    function tables()
    results=zeros(5,5,7)
        table1=zeros(15,5)
        table2=zeros(10,5)
        table3=zeros(10,5)
        global i1=0
        for delta=.25:.25:1.25
            global i1=i1+1
            global j1=0
            for spill=-.5:.25:.5 
                global j1=j1+1
                fit,inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv,inputssosv,inputscosv,profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv=run(delta,spill,false)
                global output=DataFrame(ima=inputsmasv,iso=inputssosv,ico=inputscosv,pma=vec(profitsmasv),pso=vec(profitssosv),pco=vec(profitscosv))
                c1=coef(bhat2)
                c2=coef(bhat3)
                results[i1,j1,:]=[sum(utilmasv),sum(utilsosv),sum(utilcosv),c1[2],c1[5],c2[2],c2[5]]
                table1[i1,j1]=sum(utilmasv)
                table1[i1+1,j1]=sum(utilsosv)
                table1[i1+2,j1]=sum(utilcosv)
                table2[i1,j1]=c1[2]
                table2[i1+1,j1]=c1[5]
                table3[i1,j1]=c2[2]
                table3[i1+1,j1]=c2[5]
                

                println(sum(utilmasv))
                println(sum(utilsosv))
                println(sum(utilcosv))
            end
        end
        JLD2.@save "results.jld2" table1 table2 table3 
    end



    nfarmers=50#number farmers
    min=1
    max=trunc(2*sqrt(nfarmers))
    max=convert(Int64,max)
    delta=deltagrid
    spill=spillgrid
    #max=6
    #println(max)
    dgrid=.01 #.1 #distance between grid points4
    plotdist=1:5 #1:3 #distriubtion of Plots
    Bdist=20:40 #distribution of endowments
    #delta=.75 #how fast to spillovers fall off
    #spill=-.5#-.1 #spillover coefficient
    maxin=1200  #max number of neighbors 
    gdist=0.00#.0001 #coefficient on cost of distance from plot 
    #cellsperacre=((max-min)/step+1)^2/(max-min)^2
      
    #NN=2:2
    #xgrid=[i for i=min:step:max]
   
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
    #println("Ax ",Ax)
    #println("Ax2 ",Ax2)
    Ax=Ax2[:,:]
    Btot=sum(B)
    #println("Btot ",Btot)
    distance=disttoplot(Ax)
    global distsv=distance
    global Axsv=Ax
    nmat=findneighbor(Ax)
    global nmatsv=nmat
    area,Pmat,centers,tess=comparea(Ax,min,max)
    global tesssv=tess
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
    weights1=weights(centers)
    weights2=areaweights(centers,tess)
    global weightsv1=weights1
    global weightsv2=weights2
    weights=weightsv2
    #println(weights)
    #println(n)
    #println("sum grid ",sum(inplot,dims=2))
    #println("inp",inplot)
    #println("inf",infarmer)
   
    @time inputsmasv,profitsmasv,utilmasv=indmaxso1(inputs,0)
    #time inputsmasv,profitsmasv=indmaxma(inputs)
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=0,0,0,0,0,0,0,0
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=mnneighbors(inputsmasv[1:nplots],profitsmasv)
    inputssosv,profitssosv,utilsosv=zeros(nplots),zeros(nplots),zeros(nplots)
    inputscosv,profitscosv,utilcosv=zeros(nplots),zeros(nplots),zeros(nplots)
    if maonly==false 
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

        global zma=heat(inputsmasv)
        global zso=heat(inputssosv)
        global zco=heat(inputscosv)
    end

    function fitdata(parm)

        x1=    XLSX.readxlsx("simom2.xlsx")
        println(x1)
        xs1=x1["Sheet1"]
        global inputcoefs=xs1["A1:A8"][[ 1 2 ],1]
        global yieldcoefs=xs1["A15:A22"][[1 2],1]
        global means=xs1["A30:A41"][[3 8 2 1],1]
        
        results=Optim.optimize(sse,parm,g_tol=.0001)
    end
    function sse(parm)
        delta=parm[1]
        spill=parm[2]    
        inputs=inputsmasv
        @time inputsmasv,profitsmasv,utilmasv=indmaxso1(inputs,0)
        bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=0,0,0,0,0,0,0,0
        bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=mnneighbors(inputsmasv[1:nplots],profitsmasv)
        vecsim1=bhat2.β[[5 2]] 
        vecsim2=bhat3.β[[5 2]]
        vecsim=vcat(vecsim1',vecsim2')
        vecdat=vcat(inputcoefs',yieldcoefs')
        dvec=vecsim-vecdat
        omega=[0 0 0 0 ;0 1 0 0;0 0 1 0;0 0 0 1]
        sse=(dvec'*omega*dvec)[1,1]
        println("sse momenges ",sse," delta ", delta,"spill",spill)
        println([vecsim vecdat])
        println([mnarea])
        println([means])
        
        return sse
    end
    fitresults=0
    if maonly==true
    parm=[deltagrid,spillgrid]
    fitresults=fitdata(parm)
    end
    tables()

    return(fitresults,inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv[1:nplots],inputssosv[1:nplots],inputscosv[1:nplots],profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv)
end   






parmout=[.75,-.5]
fit1,inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv,inputssosv,inputscosv,profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv=run(parmout[1],parmout[2],true)
println("fit ",fit1)
parmout=Optim.minimizer(fit1)
#parmout=[.75,-.5]
fit1,inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv,inputssosv,inputscosv,profitsmasv,profitssosv,profitscosv,utilmasv,utilsosv,utilcosv=run(parmout[1],parmout[2],false)
JLD2.@save "variables.jld2" fit1 inplot infarmer Pmat2 plotfarmer idplot bhat1 bhat1a bhat2 bhat3 bhat4 bhat5 inputsmasv inputssosv inputscosv profitsmasv profitssosv profitscosv utilmasv utilsosv utilcosv


scatter(areasv,inputssosv,xlabel="Plot Area",ylabel="Input/Area",label="Cooperative",legend_position=:topleft)
scatter!(areasv,inputscosv,label="Social Planner")
pinput=scatter!(areasv,inputsmasv,label="Nash")

scatter(areasv,vec(profitssosv)./areasv,xlabel="Plot Area",ylabel="Output/Area",label="Cooperative",legend_position=:topleft)
scatter!(areasv,vec(profitscosv)./areasv,label="Social Planner")
poutput=scatter!(areasv,vec(profitsmasv)./areasv,label="Nash")
areafarm=(areasv'*infarmer)'
scatter(areafarm,vec(utilsosv)./areafarm,xlabel="Farm Area",ylabel="Profit/Area",label="Cooperative",legend_position=:topleft)
scatter!(areafarm,vec(utilcosv)./areafarm,label="Social Planner")
pnet=scatter!(areafarm,vec(utilmasv)./areafarm,label="Nash")

display(zma)
display(zso)
display(zco)
savefig(zma,"Nashheat.png")
savefig(zso,"Coopheat.png")
savefig(zco,"Egalheat.png")
savefig(pinput,"Inputs.png")
savefig(poutput,"Outputs.png")
savefig(pnet,"Profits.png")
