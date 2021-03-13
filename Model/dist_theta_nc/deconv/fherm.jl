# Module to deal with hermite polynomials
module fherm

using LinearAlgebra
#Preliminares, Agh and Wgh are for guass hermite polynomials
using DelimitedFiles
Agh=readdlm("Agh")
Wgh=readdlm("Wgh")
Wgh=Wgh*sqrt(2*pi)

# K dimensional polynomial
function P(x,K)
  Nx=length(x)
  pp=ones(Nx,K)
  for j=2:K
    pp[:,j]=pp[:,j-1].*x
  end
  return pp
end

function dP(x,K)
  Nx=length(x)
  pp=ones(Nx,K)
  dpp=zeros(Nx,K)
  for j=2:K
    pp[:,j]=pp[:,j-1].*x
    dpp[:,j]=(j-1)*x.^(j-2)
  end
  return pp,dpp
end

# pdf of norm (mean zero)
function phi(x)
 return exp.(-x.*x/2)/sqrt(2*pi)
end

function dphi(x)
 pp=exp.(-x.*x/2)/sqrt(2*pi)
 dpp=-x.*pp
 return pp,dpp
end

# this is the extra noise parts
sige=0.05


# the weight vector
function wtret(theta,Nherm)
 ww=Wgh[1:2*Nherm,2*Nherm]
 wt=[1.0;theta[1:Nherm-1]]
 aa=Agh[1:2*Nherm,2*Nherm]
 ff=ww'*(P(aa,Nherm)*wt).^2
 wt=wt*sqrt((1-sige)/ff)
 return wt
end

function dwtret(theta,Nherm)
 ww=Wgh[1:2*Nherm,2*Nherm]
 wt=[1.0;theta[1:Nherm-1]]
 aa=Agh[1:2*Nherm,2*Nherm]
 Pol=P(aa,Nherm)
 ff=0.0
 dff=zeros(Nherm)
 for j1=1:2*Nherm
   fin=0.0
   dfin=zeros(Nherm)
   for j2=1:Nherm
     fin+=Pol[j1,j2]*wt[j2]
     dfin[j2]=Pol[j1,j2]
   end
   ff+=ww[j1]*fin^2
   dff+=2.0*ww[j1]*fin*dfin
 end
 wt=wt*sqrt((1-sige)/ff)
 dwt=Matrix(I,Nherm,Nherm)*sqrt((1-sige)/ff)-
      0.5*wt*dff'/ff
 return wt,dwt[:,2:Nherm]
end


# the  hermite polynomial likelihood (vector form)
function fh(x,mu,sig,wt,Nherm)
 xdiff=(x.-mu)/sig
 return (exp.(-xdiff.*xdiff/2.0).*((P(xdiff,Nherm)*wt).^2)+sige*phi(xdiff))/sig
end

function dfh(x,mu,sig,wt,Nherm)
 xdiff=(x.-mu)/sig
 dxdsig=-xdiff/sig
 (PP,dPP)=dP(xdiff,Nherm)
 (pphi,dpphi)=dphi(xdiff)
 f=(exp.(-xdiff.*xdiff/2.0).*((PP*wt).^2)+sige*pphi)/sig
 dfdxd=(-exp.(-xdiff.*xdiff/2.0).*((PP*wt).^2).*xdiff+
      2.0*exp.(-xdiff.*xdiff/2.0).*((PP*wt)).*(dPP*wt)+
      sige*dpphi)/sig
 ddd=0.0000001
 xdiffd=(x.+(ddd-mu))/sig
 fd=(exp.(-xdiffd.*xdiffd/2.0).*((P(xdiffd,Nherm)*wt).^2)+sige*phi(xdiffd))/sig
 dfdmu=-dfdxd/sig
 dfdsig=dfdxd.*dxdsig-f/sig
 dfdwt= (2.0*exp.(-xdiff.*xdiff/2.0).*((PP*wt))*ones(1,Nherm)).*PP/sig
 return f,dfdmu,dfdsig,dfdwt
end

# the  hermite polynomial likelihood (scalar form)
function fhi(x,mu,sig,wt,Nherm)
 xdiff=(x-mu)/sig
 return (exp(-xdiff*xdiff/2.0)*(((P(xdiff,Nherm)*wt).^2)[1])+sige*phi(xdiff))/sig
end

export  wtret,dwtret, fh,dfh,fhi

end

