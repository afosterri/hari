using Plots,PolygonOps
function heat(inputmat)
function inputv(input,x,y)
    lev=0 
    for i=1:nplots
        cell=[tesssv.Cells[i];[tesssv.Cells[i][1]]] 
        incheck=PolygonOps.inpolygon([x,y],cell)
      
        if incheck!=0
            lev=input[i]
            break
        end
    end
    return lev
end
min=1
max=14
dgrid=.1
x=min:dgrid:max
nplots=149

heatmap(x,x,(x1,y1)->log(inputv(inputmat,x1,y1)), c=:thermal, clims=(0,4))
plot!(tesssv)
end
zma=heat(inputsmasv)

zso=heat(inputssosv)

zco=heat(inputscosv)



