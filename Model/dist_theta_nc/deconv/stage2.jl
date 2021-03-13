
@everywhere include("getdat2.jl")
@everywhere using Main.dat2

@everywhere include("fherm.jl")
@everywhere using Main.fherm

@everywhere using Optim
@everywhere using DelimitedFiles

@everywhere function estmod(Ndim,Nherm,igp,nrun)

  b=readdlm("bopt")


   bsat=b[1:2]
   bact=b[3:4]
   sigtest=exp.(b[5:7])
   k=7
   theta=b[k+1:k+Ndim]
   k+=Ndim

   #next would come ptheta
   k+=Ndim-1

   wtsat=zeros(3,Nherm)
   for j=1:3
     wtsat[j,:]=wtret(b[k+1:k+Nherm-1],Nherm)
     k+=Nherm-1
   end

  N=nsize[igp]
  fliki=ones(N,Ndim)
  for i=1:N,ihet=1:Ndim
    if tstob[igp][i,1]==1
      fliki[i,ihet]*=fhi(sat[igp][i,1],theta[ihet],sigtest[1],wtsat[1,:],Nherm)*
         fhi(sat[igp][i,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],wtsat[2,:],Nherm)
    end
    if tstob[igp][i,2]==1
      fliki[i,ihet]*=fhi(act[igp][i],bact[1]+bact[2]*theta[ihet],sigtest[3],wtsat[3,:],Nherm)
    end
  end


  function ff(b)

  pthind=exp.([b[1:Ndim-1];0.0])
  ptheta=pthind/sum(pthind)

  fp=fliki*ptheta
  f=-sum(log.(fp))

  return f
  end



    k2=(Ndim-1)
    out=Array{Any}(undef,nrun)
    for irun=1:nrun
      b=randn(k2)
      bopt=optimize(ff,b,LBFGS(),Optim.Options(iterations=10000))
      out[irun]=bopt
    end
    return out
end
    

function runmod(Ndim,Nherm)
  nsim=4  
  out1=Array{Any}(undef,4,8)
  j=1
  for igp=1:4,inode=1:8
    j+=1
    out1[igp,inode]=remotecall(estmod,j,Ndim,Nherm,igp,nsim)
  end 
  out2=Array{Any}(undef,4,32)
  for igp=1:4,inode=1:8
    x=fetch(out1[igp,inode])
    for isim=1:nsim
      out2[igp,nsim*(inode-1)+isim]=x[isim]
    end     
println("end ",igp," ",inode)
  end 
  return out2
end



runs=runmod(8,5)
best=zeros(7,4)
for igp=1:4
  xout=zeros(32)
  for irun=1:32
    xout[irun]=runs[igp,irun].minimum
  end
  display(sort(xout))
  sx=sortperm(xout)
  best[:,igp]=runs[igp,sx[1]].minimizer
end

writedlm("bgps",best)
  
  
  
  
