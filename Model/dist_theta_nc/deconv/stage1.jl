
@everywhere include("getdat.jl")
@everywhere using Main.dat




@everywhere include("fherm.jl")
@everywhere using Main.fherm

@everywhere using Optim


@everywhere function estmod(Ndim,Nherm,nrun)

  function ff(b)
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
      wtsat[j,:]=wtret(b[k+1:k+Nherm-1],Nherm)
      k+=Nherm-1
    end
    fiouts=zeros(Nsato,Ndim)
    fiouta=zeros(Nacto,Ndim)
    fioutsa=zeros(Nboth,Ndim)
      for ihet=1:Ndim
        fiouts[:,ihet]=fh(sato[:,1],theta[ihet],sigtest[1],wtsat[1,:],Nherm).*
         fh(sato[:,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],wtsat[2,:],Nherm)
        fiouta[:,ihet]=fh(acto,bact[1]+bact[2]*theta[ihet],sigtest[3],wtsat[3,:],Nherm)
        fioutsa[:,ihet]=fh(satact[:,1],theta[ihet],sigtest[1],wtsat[1,:],Nherm).*
         fh(satact[:,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],wtsat[2,:],Nherm).*
         fh(satact[:,3],bact[1]+bact[2]*theta[ihet],sigtest[3],wtsat[3,:],Nherm)
      end
    fi=[fiouts;fiouta;fioutsa]*ptheta
    f=-sum(log.(max.(fi,10^-300)))
    return f
  end

function dff!(storage,b)
    #translate b to parameters
    bsat=b[1:2]
    bact=b[3:4]
    sigtest=exp.(b[5:7])
    k=7
    plow=0.05
    theta=b[k+1:k+Ndim]
    k+=Ndim
    pthind=exp.([b[k+1:k+Ndim-1];0.0])
    pth1=pthind/sum(pthind)
    ptheta=(1-Ndim*plow)*pth1.+plow
    dpdb=zeros(Ndim,Ndim-1)
    for j1=1:Ndim
      if (j1<Ndim)
        dpdb[j1,j1]=pth1[j1]
      end
      for j2=1:Ndim-1
        dpdb[j1,j2]-=pth1[j1]*pth1[j2]
      end 
    end 
    dpdb=dpdb*(1-Ndim*plow)
    

    k+=Ndim-1
    wtsat=zeros(3,Nherm)
    dwtsat=zeros(3,Nherm,Nherm-1)
    for j=1:3
      (wtsat[j,:],dwtsat[j,:,:])=dwtret(b[k+1:k+Nherm-1],Nherm)
      k+=Nherm-1
    end 
    #declear key pieces
    fiouts=zeros(Nsato,Ndim)
    dfiouts=zeros(Nsato,Ndim,k)
    fiouta=zeros(Nacto,Ndim)
    dfiouta=zeros(Nacto,Ndim,k)
    fioutsa=zeros(Nboth,Ndim)
    dfioutsa=zeros(Nboth,Ndim,k)
    #calculate likelihood function for three different pieces
    for ihet=1:Ndim
      (fhm,dfdmum,dfdsigm,dfdwtm)=dfh(sato[:,1],theta[ihet],sigtest[1],
            wtsat[1,:],Nherm)
      (fhv,dfdmuv,dfdsigv,dfdwtv)=dfh(sato[:,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],
            wtsat[2,:],Nherm)
      fiouts[:,ihet]=fhm.*fhv
      dfiouts[:,ihet,1]=fhm.*dfdmuv
      dfiouts[:,ihet,2]=fhm.*dfdmuv*theta[ihet]
      dfiouts[:,ihet,7+ihet]=dfdmum.*fhv+fhm.*dfdmuv*bsat[2]
      dfiouts[:,ihet,5]=dfdsigm.*fhv*sigtest[1]
      dfiouts[:,ihet,6]=fhm.*dfdsigv*sigtest[2]
      dfiouts[:,ihet,6+2*Ndim+1:6+2*Ndim+(Nherm-1)]=(fhv*ones(1,Nherm-1)).*(dfdwtm*dwtsat[1,:,:])
      dfiouts[:,ihet,6+2*Ndim+(Nherm-1)+1:6+2*Ndim+2*(Nherm-1)]=
              (fhm*ones(1,Nherm-1)).*(dfdwtv*dwtsat[2,:,:])
      (fha,dfdmua,dfdsiga,dfdwta)=dfh(acto,bact[1]+bact[2]*theta[ihet],sigtest[3],
            wtsat[3,:],Nherm)
      fiouta[:,ihet]=fha
      dfiouta[:,ihet,3]=dfdmua
      dfiouta[:,ihet,4]=dfdmua*theta[ihet]
      dfiouta[:,ihet,7+ihet]=dfdmua*bact[2]
      dfiouta[:,ihet,7]=dfdsiga*sigtest[3]
      dfiouta[:,ihet,6+2*Ndim+2*(Nherm-1)+1:6+2*Ndim+3*(Nherm-1)]=dfdwta*dwtsat[3,:,:]
      (fhm,dfdmum,dfdsigm,dfdwtm)=dfh(satact[:,1],theta[ihet],sigtest[1],
            wtsat[1,:],Nherm)
      (fhv,dfdmuv,dfdsigv,dfdwtv)=dfh(satact[:,2],bsat[1]+bsat[2]*theta[ihet],sigtest[2],
            wtsat[2,:],Nherm)
      (fha,dfdmua,dfdsiga,dfdwta)=dfh(satact[:,3],bact[1]+bact[2]*theta[ihet],sigtest[3],
            wtsat[3,:],Nherm)
      fioutsa[:,ihet]=fhm.*fhv.*fha
      dfioutsa[:,ihet,1]=fhm.*dfdmuv.*fha
      dfioutsa[:,ihet,2]=fhm.*dfdmuv.*fha*theta[ihet]
      dfioutsa[:,ihet,3]=fhm.*fhv.*dfdmua
      dfioutsa[:,ihet,4]=fhm.*fhv.*dfdmua*theta[ihet]
      dfioutsa[:,ihet,7+ihet]=dfdmum.*fhv.*fha+fhm.*dfdmuv.*fha*bsat[2]+
                  fhm.*fhv.*dfdmua*bact[2]
      dfioutsa[:,ihet,5]=dfdsigm.*fhv.*fha*sigtest[1]
      dfioutsa[:,ihet,6]=fhm.*dfdsigv.*fha*sigtest[2]
      dfioutsa[:,ihet,7]=fhm.*fhv.*dfdsiga*sigtest[3]
      dfioutsa[:,ihet,6+2*Ndim+1:6+2*Ndim+(Nherm-1)]=((fhv.*fha)*ones(1,Nherm-1)).*(dfdwtm*dwtsat[1,:,:])
      dfioutsa[:,ihet,6+2*Ndim+(Nherm-1)+1:6+2*Ndim+2*(Nherm-1)]=
              ((fhm.*fha)*ones(1,Nherm-1)).*(dfdwtv*dwtsat[2,:,:])
      dfioutsa[:,ihet,6+2*Ndim+2*(Nherm-1)+1:6+2*Ndim+3*(Nherm-1)]=
              ((fhm.*fhv)*ones(1,Nherm-1)).*(dfdwta*dwtsat[3,:,:])
    end
    # integrate over heterogeneity and calculate
    fip=[fiouts;fiouta;fioutsa]
    fi=fip*ptheta
    dfi=[dfiouts;dfiouta;dfioutsa]
    f=-sum(log.(max.(fi,10^-300)))
    df=zeros(6+2*Ndim+3*(Nherm-1))
    dfdp=zeros(Ndim)
    for i=1:(Nsato+Nacto+Nboth),ihet=1:Ndim
      df-=ptheta[ihet]*dfi[i,ihet,:]/fi[i]
      dfdp[ihet]-=fip[i,ihet]/fi[i]
    end
    df[8+Ndim:6+2*Ndim]=dfdp'*dpdb
    k=(6+2*Ndim+3*(Nherm-1))
    for j=1:k
      storage[j]=df[j]
     end
  end
  
  
  
  
  
  k=(6+2*Ndim+3*(Nherm-1))
  out=Array{Any}(undef,nrun)
  for irun=1:nrun
    b=randn(k)
    println("Ndim: ",Ndim," Nherm: ",Nherm," run ",irun)
    println(b)
    bopt=optimize(ff,dff!,b,LBFGS(),Optim.Options(iterations=20000))
    out[irun]=bopt
  end
  return out
end

output=Array{Any}(undef,20,10)

function runmod(Ndim,Nherm)
  nsim=2
  out1=Array{Any}(undef,64)
  for j=2:33
    out1[j-1]=remotecall(estmod,j,Ndim,Nherm,nsim)
  end 
  out2=Array{Any}(undef,128)
  for j=1:32
    x=fetch(out1[j])
    for isim=1:nsim
      out2[nsim*(j-1)+isim]=x[isim]
    end 
println("j end ",j)
  end 
  return out2
end


@time a75=fetch(runmod(8,5))
for i=1:64
 ff[i]=a75[i].minimum
println(ff[i])
display(a75[i].minimizer)
end
display(sort(ff))
writedlm("a75",ff)


