
@everywhere include("getdat2.jl")
@everywhere using Main.dat2

@everywhere include("fherm.jl")
@everywhere using Main.fherm

@everywhere using Optim
@everywhere using DelimitedFiles

@everywhere function estmod(Ndim,Nherm1,Nherm2,igp,nrun)



  k=(7+Ndim+3*(Nherm1-1))
  b=readdlm("bopt")

  bsat=b[1:2]
  bact=b[3:4]
  sigtest=exp.(b[5:7])
  k=7
  theta=b[k+1:k+Ndim]
  k+=Ndim
  # next would come ptheta
  k+=Ndim-1
  wtsat=zeros(3,Nherm1)
  for j=1:3
    wtsat[j,:]=wtret(b[k+1:k+Nherm1-1],Nherm1)
    k+=Nherm1-1
  end

  bgps=readdlm("bgps")
  pgps=zeros(Ndim,4)
  for iigp=1:4 
    pthind2=exp.([bgps[1:Ndim-1,iigp];0.0])
    pgps[:,iigp]=pthind2/sum(pthind2)
  end


  N=nsize[igp]
  fliki=ones(N,Ndim)
  for i=1:N,ihet=1:Ndim
    if tstob[igp][i,1]==1
      fliki[i,ihet]*=fhi(sat[igp][i,1],theta[ihet],sigtest[1],wtsat[1,:],Nherm1)*
         fhi(sat[igp][i,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],wtsat[2,:],Nherm1)
    end
    if tstob[igp][i,2]==1
      fliki[i,ihet]*=fhi(act[igp][i],bact[1]+bact[2]*theta[ihet],sigtest[3],wtsat[3,:],Nherm1)
    end
  end


  function ff(b)
    bcq=b[1:2]
    sigcq=exp(b[3])
    wtcq=wtret(b[4:3+Nherm2-1],Nherm2)
    f=0.0
    for i=1:N
      if cqob[igp][i]==1
        Fi=zeros(Ndim)
        for ihet=1:Ndim
          Fi[ihet]=fliki[i,ihet].*fhi(cqual[igp][i],bcq[1]+bcq[2]*theta[ihet],sigcq,wtcq,Nherm2)
        end
        f-=sum(log.(max.(Fi'*pgps[:,igp],10^-300)))
      end
    end
    return f
  end 


    k2=(3+(Nherm2-1))
    out=Array{Any}(undef,nrun)
    for irun=1:nrun
      b=randn(k2)
      bopt=optimize(ff,b,LBFGS(),Optim.Options(iterations=10000))
      out[irun]=bopt
    end
    return out
end
    


  
  
  function runmod(Ndim,Nherm1,Nherm2)
  nsim=4
  out1=Array{Any}(undef,4,8)
  j=1
  for igp=1:4,inode=1:8
    j+=1
    out1[igp,inode]=remotecall(estmod,j,Ndim,Nherm1,Nherm2,igp,nsim)
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

Nherm2=5
runs=runmod(8,5,Nherm2)
k2=(3+(Nherm2-1))

best=zeros(k2,4)
for igp=1:4 
  xout=zeros(32)
  for irun=1:32
    xout[irun]=runs[igp,irun].minimum
  end
  display(sort(xout))
  sx=sortperm(xout)
  best[:,igp]=runs[igp,sx[1]].minimizer
end

writedlm("bst3",best)

