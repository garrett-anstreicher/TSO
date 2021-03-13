Ndim=8
Nherm=5

using PyPlot




include("fherm.jl")
using Main.fherm

include("mkltab.jl")

using Optim


using DelimitedFiles
b=readdlm("bopt")

  bsat=b[1:2]
  bact=b[3:4]
  sigtest=exp.(b[5:7])
  k=7
  plow=0.05
  theta=b[k+1:k+Ndim]
  k+=Ndim
  pthind=exp.([b[k+1:k+Ndim-1];0.0])
  ptheta=(1-Ndim*plow)*pthind/sum(pthind).+plow
  k+=Ndim-1
  wtsat=zeros(3,Nherm)
  for j=1:3
    global k
    wtsat[j,:]=wtret(b[k+1:k+Nherm-1],Nherm)
    k+=Nherm-1
  end

btest=zeros(2,3)
btest[2,1]=1.0
btest[:,2]=bsat
btest[:,3]=bact
label=["Constant","\$\\theta\$"]
mkltab(label,btest)
itheta=sortperm(theta)
xtheta=[ptheta[itheta] theta[itheta]]
mknltab(xtheta)


label=["\$\\sigma\$"; "\$ x^0 \$"; "\$ x^1 \$"; "\$ x^2 \$"; "\$ x^3 \$"; "\$ x^4 \$"]
mkltab(label,[sigtest' ;wtsat'])
plot(theta,ptheta,linestyle="None",marker="o",color="blue")
vlines(theta,zeros(size(theta)),ptheta,color="blue")
hlines(0,minimum(theta),maximum(theta),color="red")
savefig("ddist.pdf")


NN=300
xfd=zeros(NN,5)
for i=1:NN
 xfd[i,1]=-2.5+5*i/NN
end
for j=1:3
   xfd[:,j+1]=fh(xfd[:,1],ptheta'*theta,sigtest[j],wtsat[j,:],Nherm)
 end

fig,figg=subplots()
figg[:plot](xfd[:,1],xfd[:,2],"b-",label="Math SAT")
figg[:plot](xfd[:,1],xfd[:,3],"r:",label="Verbal SAT") 
figg[:plot](xfd[:,1],xfd[:,4],"g--",label="ACT")
figg[:legend](loc="best")
title("Density of Test Measurement Error")
savefig("me.pdf")


bgps=readdlm("bgps")
pgps=zeros(Ndim,5)
for igp=1:4
  pthind2=exp.([bgps[1:Ndim-1,igp];0.0])
  ptheta2=pthind2/sum(pthind2)
  pgps[:,igp+1]=ptheta2[itheta]
end
pgps[:,1]=theta[itheta]
mknltab(pgps)
fig,figg=subplots()
figg[:plot](pgps[:,1],pgps[:,2],"b-",label="White/Asian Men")
figg[:plot](pgps[:,1],pgps[:,3],"r:",label="White/Asian Women")
figg[:plot](pgps[:,1],pgps[:,4],"g-",label="Black/Hispanic Men")
figg[:plot](pgps[:,1],pgps[:,5],"k-",label="Black/Hispanic Women")
figg[:legend](loc="best")
title("Distribution of Ability by race and sex")
savefig("all.pdf")

bst3=readdlm("bst3")

label=["Constant";"\$\\theta\$";;"\$\\sigma\$"; "\$ x^0 \$"; "\$ x^1 \$"; "\$ x^2 \$"; "\$ x^3 \$"; "\$ x^4 \$"]
Nherm2=5
bout=zeros(3+Nherm2,4)
for igp=1:4
 bout[4:3+Nherm2,igp]=wtret(bst3[4:3+Nherm2-1,igp],Nherm2)
end
bout[1:2,:]=bst3[1:2,:]
bout[3,:]=exp.(bst3[3,:])
mkltab(label,bout)

for igp=1:4
   xfd[:,igp+1]=fh(xfd[:,1],ptheta'*theta,bout[3,igp],bout[4:3+Nherm2,igp],Nherm2)
 end

fig,figg=subplots()
figg[:plot](xfd[:,1],xfd[:,2],"b-",label="White/Asian Men")
figg[:plot](xfd[:,1],xfd[:,3],"r:",label="White/Asian Women")
figg[:plot](xfd[:,1],xfd[:,4],"g--",label="Black/Hispanic Men")
figg[:plot](xfd[:,1],xfd[:,5],"k-",label="Black/Hispanic Women")
figg[:legend](loc="best")
title("Density of College Measurement Error")
savefig("cmeas.pdf")
