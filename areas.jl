using Distributed
@everywhere using Delaunay, VoronoiCells, GeometryBasics, Random, PolygonOps
function run()
max=10
min=1
nplots=10
del=1
Ax=rand(nplots,2).*(max-min).+min
points=Point2.(Ax[:,1],Ax[:,2])
xsquare=Rectangle(Point2(min,min), Point2(max,max))
tess=voronoicells(points,xsquare)
ngrid=convert(Int64,ceil((max-min)/del)+1)
mat=zeros(ngrid,ngrid)
i=0
 for x=min:del:max
    i=i+1
    j=0
    for y=min:del:max
        j=j+1
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
        #println(mat[i,j]," ",x," ",y)        
    end
end
        nr=nplots
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




weights=zeros(nplots,nplots)
onev=ones(ngrid^2)
ivec=[i for i=min:del:max]
jvec=[j for j=min:del:max]
vmat=vec(mat)
for k=1:nplots
    local dist=vec(sqrt.((ivec.-center[k,1]).^2 .+ (jvec'.-center[k,2]).^2))
    dwt=1 ./ (1 .+ exp.(dist./del))
    for l=1:nplots    
        weights[k,l]=sum(dwt[vmat.==l]) 
    end

end
weights=weights./sum(weights,dims=2)'
println("centers ",center')
plot(tess)
println("mat ",weights)
end
run()
