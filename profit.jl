using Distributed
@everywhere using Delaunay, VoronoiCells, GeometryBasics, Random, PolygonOps, Plots, Optim
@everywhere using DataFrames, Econometrics, StatsBase, Distributed, SharedArrays, ForwardDiff,NLsolve,LinearAlgebra
Random.seed!(1001)

function run()

    function disttoplot(Ax)
        nplots=size(Ax)[1]
        distance=sqrt.(Ax[:,1].^2 .+ Ax[:,2].^2)
        return distance
    end

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
        B[1]=2
        B[2]=2
        #B[2]=4
        Ax=rand(nplots,2)*(max-min) .+ min
        Ax=[1 1 ; 1 1.25; 2 1.75 ; 2 2 ]
        inputs=rand(nplots)
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
        
        
        z=1/(1+exp(((x[1]-y[1])^2+(x[2]-y[2])^2)/delta))
        #=
        z=0
        if (sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)<1) 
            z=1
        end
        =#
        return z
    end
    function f(x)
        if x>=0 
            z=x^beta
        else
            z=-1000*x^2
        end
        #z=(real((x.>=0).*(complex.(x)).^beta)).+((x.<0).*(-10000))
        #if sum(x.<0)>0 
        #println("prod",x," ",(z))
        #end
        return z
    end
    function inpoly(Pmat2)
        
        inpoly2=zeros(n2)
        plotfarmer=[false for i=1:nplots , j=1:nfarmers]

        m=0
        for x=min:step:max
            for y=min:step:max
                m=m+1
                for k=1:nplots
                    inpoly=PolygonOps.inpolygon([x,y],Pmat2[k])                    
                    if abs(inpoly)==1 
                        inpoly2[m]=k
                    end
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
        #den=maximum(sum(distgrid,dims=2))
        den=1
        weights=distgrid./den
        return weights
    end
    function agprod(inputs,add,add2)
        #println(size(inputs))
        #println(size(inplot))
        #println(add)
        #println(add2)
        #println(weights)
        inputsgrid=inplot*inputs
        aginputs=(inputsgrid.+spill.*(weights*inputsgrid))
        #icost=exp.(disttogrid.*gdist)
        #println(icost," icos")
        #println(inputsgrid)
        #println(size(inputsgrid),size(add2),size(cellsperacre))
        tinputs=inputsgrid'*add2./cellsperacre
        plotout=f.(aginputs)'*add./cellsperacre
        return plotout,tinputs
    end

    function dagprod(inputs,add,add2)
        #println("hello2")
        
        dchange=.001
        inputs1=inputs.+dchange.*I(nplots)
        global dinpu1=(inputs)
        global dinpu2=(inputs1)
        p1,i1=agprod(inputs,add,add2)
        p2,i2=agprod(inputs1,add,add2)
        dagprod2=(p2.-p1)./dchange        
        return dagprod2'
    end
    function dagprod(inputs,add,add2,add3)
        #println("hello2")
        
        dchange=.0000001
        eye=I(nplots)
        eye=eye[:,add3]
        inputs1=inputs.+dchange.*eye
        
        global dinpu1=(inputs)
        global dinpu2=(inputs1)
        p1,i1=agprod(inputs,add,add2)
        p2,i2=agprod(inputs1,add,add2)
        #println(size(p1),size(p2),size(inputs1))
        dagprod2=(p2.-p1)./dchange        
        return dagprod2'
    end
#=
    function indmax(inputs)
        function profits(inputi,k,vec,inputs)
            #println("inputi ",inputi)
            x=inputs[:]
    
            x[vec]=inputi
            
            plotout,tinputs=agprod(x,infarmer[:,k],infarmer[:,k])
            #println("plotout ",plotout)
            #println("tinputs ",tinputs)
            tplotout=plotout
            profit=tplotout-10000*(tinputs>B[k])*(tinputs-B[k])^2-10000*sum((inputi.<0).*inputi.^2)
            #println("plot out ", k," ",tplotout," ",tinputs)
            return -profit
        end
        local xinputs1=SharedArray{Float64}(nplots)
        local profitsv=SharedArray{Float64}(nfarmers)
        xinputs=inputs[:]
        del=100
        while del>0.0001
            xinputs1[:]=xinputs.+0
            @sync @distributed  for k=1:nfarmers
                start=idplot[k,1]
                ends=start+N[k]-1
                vec=start:ends
                inputi=xinputs1[vec]    
                #println(vec," ",inputi)
                res=optimize(z->profits(z,k,vec,xinputs1),inputi)
                res2=Optim.minimizer(res)
                #println("res2 ",res2)
                #println("vec ",vec)
                xinputs1[vec]=res2
                #println("xinputs1 ",xinputs1)
                profitsv[k]=Optim.minimum(res)
                #println("profits 1",profitsv)
            end
            #println("xinputs1 ",xinputs1," ",xinputs)
            delv=xinputs1-xinputs
            xinputs=xinputs1[:]
            del=delv'delv/nplots
            println("del ",del)
        end
        println("profits2 ", profitsv)
        return xinputs,-profitsv
    end

    function indmaxso(inputs)

        function profitsso(inputi,k,vec,inputs)
            #println("inputi ",inputi)
            x=inputs[:]
    
            x[vec]=inputi
            
            plotout,tinputs=agprod(x,infarmer,infarmer[:,k])
            #println("plotout ",plotout)
            #println("tinputs ",tinputs)
            tinputs=tinputs
            tplotout=sum(plotout)
            profit=tplotout-10000*(tinputs>B[k])*(tinputs-B[k])^2-10000*sum((inputi.<0).*inputi.^2)
            #println("plot out ", k," ",tplotout," ",tinputs)
            return -profit
        end
        local xinputs1=SharedArray{Float64}(nplots)
        local profitsv=SharedArray{Float64}(nfarmers)
        
        xinputs=inputs[:]
        del=100
        while del>0.0001
            xinputs1[:]=xinputs.+0
            @sync @distributed  for k=1:nfarmers
                start=idplot[k,1]
                ends=start+N[k]-1
                vec=start:ends
                inputi=xinputs1[vec]    
                res=optimize(z->profitsso(z,k,vec,xinputs1),inputi)
                res2=Optim.minimizer(res)
                #println("res2 ",res2)
                #println("vec ",vec)
                xinputs1[vec]=res2
                #println("xinputs1 ",xinputs1)
                profitsv[k]=Optim.minimum(res)
                #println("profits 1",profitsv)
            end
            #println("xinputs1 ",xinputs1," ",xinputs)
            delv=xinputs1-xinputs
            xinputs=xinputs1[:]
            del=delv'delv/nplots
            println("del ",del)
        end
        #println("profits2 ", profitsv)
        return xinputs,-profitsv
    end

    function indmaxco(xinputs)
   

        function profitsco(inputi,k,vec,inputs,lambda)
            #println("inputi ",inputi)
            x=inputs[:]
    
            x[vec]=inputi
            
            plotout,tinputs=agprod(x,infarmer,infarmer)
            
            #println("plotout ",plotout)
            #println("tinputs ",tinputs)
            tinputs=sum(tinputs)
            tplotout=sum(plotout)
            profit=tplotout+exp(lambda)*(tB-tinputs)-10000*sum((inputi.<0).*inputi.^2)
            #println("plot out ", k," ",tplotout," ",tinputs)
            return -profit
        end

        function lambdaco(xinputs,lambda1)
          
            lambda=lambda1[1]
            println("lambda ",lambda)
            xinputs1[:]=xinputs.+0
            @sync @distributed  for k=1:nfarmers
                start=idplot[k,1]
                ends=start+N[k]-1
                vec=start:ends
                inputi=xinputs1[vec]    
                #println("check ",inputi)
                res=optimize(z->profitsco(z,k,vec,xinputs1,lambda),inputi)
                res2=Optim.minimizer(res)
                #println("res2 ",res2)
                #println("vec ",vec)
                xinputs1[vec]=res2
                #println("xinputs1 ",xinputs1)
                profitsv[k]=Optim.minimum(res)
                #println("profits 1",profitsv)
            end
            println("pr ",xinputs1)
            tinputs,tprofits=agprod(xinputs1,infarmer,infarmer)
            #println("tin ",tinputs,tprofits,lambda1)
            ttinputs=sum(tinputs)
            ttprofits=sum(tprofits)
      
            profits=ttprofits+exp(lambda)*(tB-ttinputs)
            return profits
        end
        local xinputs1=SharedArray{Float64}(nplots)
        local profitsv=SharedArray{Float64}(nfarmers)
        lambda1=zeros(1)
        res=optimize(lambda->lambdaco(xinputs,lambda),lambda1)
        res1=Optim.minimizer(res)
        pr1=lambdaco(xinputs,res1)
        tinputs=xinputs1
        tprofits=-profitsv
        return tinputs,tprofits
    end

=#

function indmaxma1(inputs)
            
    function nlf!(Gfun,inplam)              
        lambda=(inplam[xnk+1])
        xinputs3=xinputs2[:]
        xinputs3[xstart:xend]=inplam[1:xnk]
        #plotout,tinputs=agprod(xinputs,inplot,inplot)
        dplotout=dagprod(xinputs3,inplot,inplot,xstart:xend)
        #dplotout=diag(plotfarmer*plotfarmer'*dplotout) 
        #println(size(dplotout),size(plotfarmer))
        dplotout=plotfarmer'*dplotout
        #println(dplotout)
        #println(xk)
        dp1=dplotout[xk,:]
        #println("dp1 ",dp1)
        #dp1=diag(dp1)

        price=exp.(gdist*distance).*area'
        pr1=price[xstart:xend]
        tinputs=sum(pr1.*inplam[1:xnk])
        deriv=dp1-lambda.*pr1 #-1000*(inplam[1:xnk].<0).*inplam[1:xnk]
        Gfun[1:xnk]=deriv
        #println(size(tinputs),size(B),size(plotfarmer))
        bc=B[xk]-tinputs[1]
        #println(size(Gfun),size(tB),size(tinputs))
        Gfun[xnk+1]=bc
        #println(Gfun)
        #G1[:]=G
    # println("F" ,G)
    #println(Gfun)
    #println(typeof(Gfun))
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
    lamvec1=ones(K)
    lamvec2=lamvec1[:]
    while del>.0000001
        @sync @distributed for k=1:K            
        #for k=1:K            
            xstart=idplot[k,1]
            xend=xstart+N[k]-1
            xk=k
            xnk=N[k]
            start3=ones(xnk+1)
            start3[1:xnk]=xinputs[xstart:xend]
            start3[xnk+1]=lamvec1[xk]
            
            z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)
            println("z1 ",z1.zero)
            xinputs1[xstart:xend]=z1.zero[1:xnk]
            
            
            lamvec2[xk]=z1.zero[xnk+1]
        end
        lamvec1=lamvec2[:]
        xinputs2=xinputs1[:]
        del1=xinputs-xinputs1
        println("del1 ",del1)
        xinputs=.5*xinputs+.5*xinputs1
        xinputs1=xinputs[:]
        del=del1'del1/nplots  
        println("del ",del)
    end     
    println(xinputs')  
    prof,inpt=agprod(xinputs,inplot,inplot)
    return xinputs,prof
end

function indmaxso1(inputs)
            
    function nlf!(Gfun,inplam)              
        lambda=(inplam[xnk+1])
        
        xinputs3=xinputs2[:]
        xinputs3[xstart:xend]=inplam[1:xnk]
        #plotout,tinputs=agprod(xinputs,inplot,inplot)
        dplotout=dagprod(xinputs,inplot,inplot,xstart:xend)
        #dplotout=diag(plotfarmer*plotfarmer'*dplotout) 
        dplotout=sum(dplotout,dims=1)'
        #println(size(dplotout))
        dp1=dplotout
        price=exp.(gdist*distance).*area'
        pr1=price[xstart:xend]
        tinputs=sum(pr1.*inplam[1:xnk])
        #println(size(price),size(pr1),size(tinputs),size(lambda),size(dp1))
        deriv=dp1-lambda.*pr1 #-1000*(inplam[1:xnk].<0).*inplam[1:xnk]
        Gfun[1:xnk]=deriv
        #println(size(tinputs),size(B),size(plotfarmer))
        bc=B[xk]-tinputs[1]
        #println(size(Gfun),size(tB),size(tinputs))
        Gfun[xnk+1]=bc
        #println(Gfun)
        #G1[:]=G
    # println("F" ,G)
    #println(Gfun)
    #println(typeof(Gfun))
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
    lamvec1=ones(K)
    lamvec2=lamvec1[:]
    while del>.0000001
    @sync @distributed for k=1:K            
        #for k=1:K            
            xstart=idplot[k,1]
            xend=xstart+N[k]-1
            xk=k
            xnk=N[k]
            start3=ones(xnk+1)
            start3[1:xnk]=xinputs[xstart:xend]
            start3[xnk+1]=lamvec1[xk]
            z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)
            
            println("z1 ",z1.zero)
            xinputs1[xstart:xend]=z1.zero[1:xnk]
            lamvec2[xk]=z1.zero[xnk+1]
        end
        lamvec1=lamvec2[:]
        xinputs2=xinputs[:]
        del1=xinputs-xinputs1
        xinputs=.8*xinputs+.2*xinputs1
        del=del1'del1/nplots  
        println("del ",del)
    end
    println(xinputs')     
    prof,inpt=agprod(xinputs,inplot,inplot)
    return xinputs,prof
end
    function indmaxma(xinputs)
            
        function nlf!(Gfun,inplam)              
            lambda=(inplam[nplots+1:nplots+nfarmers])
            inputs=inplam[1:nplots]
            plotout,tinputs=agprod(inputs,inplot,inplot)
            #println(inputs," ",inplot)
        
            #println(inputs)
            #println(lambda)

            dplotout=dagprod(inputs,inplot,inplot)
            global dbase=dplotout
            #area=sum(inplot,dims=1)/cellsperacre
            #println(dplotout)
            #println(size(dplotout))
            #dplotout=sum(dplotout,dims=1)'  
            #println(size(plotfarmer),size(dplotout),size(plotfarmer))
            dplotout=diag(plotfarmer*plotfarmer'*dplotout) 
            global dplsv=dplotout     
            global areasv=area   
            global pricsv=exp.(gdist*distance)
            global tinpsv=tinputs
            global inpsv=inputs
            tinputs=tinputs'.*exp.(gdist*distance)
        
            #for i=1:nplots
            #Gfun[i]=dplotout[i]-lambda*exp.(gdist*distance[i])*area[i]
            #end
            #println(size(dplotout),size(lambda),size(gdist),size(distance),size(area))
            lamplot=plotfarmer*lambda
            price=exp.(gdist*distance).*area'
            Gfun[1:nplots]=dplotout-lamplot.*price
            #println(size(tinputs),size(B),size(plotfarmer))
            bc=B.-plotfarmer'*tinputs
            #println(size(Gfun),size(tB),size(tinputs))
            Gfun[nplots+1:nplots+nfarmers]=bc
            #G1[:]=G
        # println("F" ,G)
        #println(Gfun)
        #println(typeof(Gfun))
        end
        
        start3=ones(nplots+nfarmers)
        z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)
        inall=ones(n2,1)
        prof,inpt=agprod(z1.zero[1:nplots],inplot,inplot)
        println(z1.zero)
        return z1.zero,prof
    end

    
    function indmaxso(xinputs)
        
        function nlf!(Gfun,inplam)              
            lambda=(inplam[nplots+1:nplots+nfarmers])
            inputs=inplam[1:nplots]
            plotout,tinputs=agprod(inputs,inplot,inplot)
            
            #println(inputs," ",inplot)
          
            #println(inputs)
            #println(lambda)

            dplotout=dagprod(inputs,inplot,inplot)
            global dbase1=dplotout
            #area=sum(inplot,dims=1)/cellsperacre
            #println(dplotout)
            #println(size(dplotout))
            dplotout=sum(dplotout,dims=1)'           
            tinputs=tinputs'.*exp.(gdist*distance)
            #for i=1:nplots
            #Gfun[i]=dplotout[i]-lambda*exp.(gdist*distance[i])*area[i]
            #end
            #println(size(dplotout),size(lambda),size(gdist),size(distance),size(area))
            lamplot=plotfarmer*lambda
            price=exp.(gdist*distance).*area'
            global dplsv1=dplotout     
            global areasv1=area   
            global pricsv1=exp.(gdist*distance)
            global tinpsv1=tinputs
            global inpsv1=inputs
            Gfun[1:nplots]=dplotout-lamplot.*price
            #println(size(tinputs),size(B),size(plotfarmer))
            bc=B.-plotfarmer'*tinputs
            #println(size(Gfun),size(tB),size(tinputs))
            Gfun[nplots+1:nplots+nfarmers]=bc
            #G1[:]=G
           # println("F" ,G)
           # println(Gfun)
           #println(typeof(Gfun))
        end
 

        
        start3=ones(nplots+nfarmers)
        z1=nlsolve(nlf!,start3,show_trace=false,ftol=1e-8)
        inall=ones(n2,1)
        prof,inpt=agprod(z1.zero[1:nplots],inplot,inplot)
        println(z1.zero)
        return z1.zero,prof
    end


    function indmaxco(xinputs)
        
        function nlf!(Gfun,inplam)              
            lambda=(inplam[nplots+1])
            inputs=inplam[1:nplots]
            plotout,tinputs=agprod(inputs,inplot,inplot)
            #println(inputs," ",inplot)
          
            #println(inputs)
            #println(lambda)

            dplotout=dagprod(inputs,inplot,inplot)
            
            #println(dplotout)
            #println(size(dplotout))
            dplotout=sum(dplotout,dims=1)'           
            tinputs=tinputs*exp.(gdist*distance)
            #for i=1:nplots
            #Gfun[i]=dplotout[i]-lambda*exp.(gdist*distance[i])*area[i]
            #end
            #println(size(dplotout),size(lambda),size(gdist),size(distance),size(area))
            Gfun[1:nplots]=dplotout-lambda*exp.(gdist*distance).*area'
            
            bc=tB-tinputs[1]
            #println(size(Gfun),size(tB),size(tinputs))
            Gfun[nplots+1]=bc
            #G1[:]=G
           # println("F" ,G)
           # println(Gfun)
           #println(typeof(Gfun))
        end
        
        start3=ones(nplots+1)
        z1=nlsolve(nlf!,start3,show_trace=false,ftol=.001)
        inall=ones(n2,1)
        prof,inpt=agprod(z1.zero[1:nplots],inplot,inplot)
        println(z1.zero)
        return z1.zero,prof
    end

    
function mnneighbors(input,prplotv)
    function ln(x)
        (x.>=0).*real.(log.(complex.(x)))+(x.<0).*0
    end
    #caluclates input of neighbors and runs regressions
    #println(size(input),size(nmat))
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


 

    nfarmers=2#number farmers
    min=1
    max=trunc(2*sqrt(nfarmers))
    max=convert(Int64,max)
    
    #println(max)
    step=.1 #.1 #distance between grid points4
    plotdist=2:2 #1:3 #distriubtion of Plots
    Bdist=3:3 #distribution of endowments
    delta=.5 #how fast to spillovers fall off
    spill=-.01#-.1 #spillover coefficient
    maxin=1200  #max number of neighbors 
    gdist=0#.0001 #coefficient on cost of distance from plot 
    scorr=0  #spatial correlation of plot sizes
    cellsperacre=((max-min)/step+1)^2/(max-min)^2
      
    NN=2:2
    xgrid=[i for i=min:step:max]
    n=size(xgrid)[1]
    n2=n*n
    beta=1/2
    K,nplots,Ax,N,idplot,B,inputs=createdata(nfarmers,plotdist,Bdist,min,step,max)
    
    tB=sum(B) 
    Ax2=spacecorr(Ax,scorr,min,step,max)
    Ax=Ax2
    Btot=sum(B)
    println("Btot ",Btot)
    distance=disttoplot(Ax)
    
    nmat=findneighbor(Ax)
    area2,Pmat=comparea(Ax,min,max)
    xwalk,Pmat2=xwalkfind(Ax,Pmat)
    inplot,infarmer,plotfarmer=inpoly(Pmat2)
    global inplot1=inplot
    global infarm1=infarmer
    global polotfarm1=plotfarmer
    disttogrid=inplot*distance
    weights=weights()
    area=sum(inplot,dims=1)/cellsperacre
    println("area ",area2, " ",area)

    #println(weights)
    #println(n)
    #println("sum grid ",sum(inplot,dims=2))
    #println("inp",inplot)
    #println("inf",infarmer)
    @time inputsmasv,profitsmasv=indmaxma1(inputs)
    #time inputsmasv,profitsmasv=indmaxma(inputs)
    bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=0,0,0,0,0,0,0,0
    #bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,mnarea,areacov=mnneighbors(inputsmasv[1:nplots],profitsmasv)
    
    #println(profitsv)
    println("now so")
    @time inputssosv,profitssosv=indmaxso1(inputs)
    #time inputssosv,profitssosv=indmaxso(inputs)
    #println(profitsosv)
    #println(infarmer)
    #println(B)
    inputscosv,profitscosv=indmaxco(inputs)
    global wtsv=weights*inplot
    return(inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv[1:nplots],inputssosv[1:nplots],inputscosv[1:nplots],profitsmasv,profitssosv,profitscosv)
end   

inplot,infarmer,Pmat2,plotfarmer,idplot,bhat1,bhat1a,bhat2,bhat3,bhat4,bhat5,inputsmasv,inputssosv,inputscosv,profitsmasv,profitssosv,profitscosv=run()
println(sum(profitsmasv))
println(sum(profitssosv))
println(sum(profitscosv))