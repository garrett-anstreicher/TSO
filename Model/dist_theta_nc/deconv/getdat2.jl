
module dat2
  using DelimitedFiles
  using Statistics
  xdat=readdlm("bbfact93.raw")
  const igp93=trunc.(Int64,xdat[:,2])
  id=xdat[:,1]
  const N93=length(id)
  satt93=xdat[:,3:4]
  actt93=xdat[:,5]
  cqualt93=xdat[:,6]
  const cert93=trunc.(Int64,xdat[:,7])
  const major93=trunc.(Int64,xdat[:,8])
  const majorname93=["STEM","Business","SS/Humanity","Education","Other"]

  xdat=readdlm("bbfact01.raw")
  const igp01=trunc.(Int64,xdat[:,2])
  id=xdat[:,1]
  const N01=length(id)
  actt01=xdat[:,3]
  satt01=xdat[:,4:5]
  cqualt01=xdat[:,6]

  xdat=readdlm("bbfact08.raw")
  const igp08=trunc.(Int64,xdat[:,2])
  id=xdat[:,1]
  const N08=length(id)
  satt08=xdat[:,3:4]
  actt08=xdat[:,5]
  cqualt08=xdat[:,6]

  satt=[satt93;satt01;satt08]
  actt=[actt93;actt01;actt08]
  igptt=[igp93;igp01;igp08]
  qualt=[cqualt93;cqualt01;cqualt08]
  
  norms=readdlm("svnorm")

  qualtx=qualt[qualt.>0]
  qm=mean(qualtx)
  qs=std(qualtx)
  
  testobs=zeros(Int64,length(actt),2)
  cqobs=zeros(Int64,length(actt))
  for i=1:length(actt)
    if satt[i,1]>0
      testobs[i,1]=1
      satt[i,1]=(satt[i,1]-norms[1])/norms[3]
      satt[i,2]=(satt[i,2]-norms[2])/norms[4]
    end
    if actt[i]>0
      testobs[i,2]=1
      actt[i]=(actt[i]-norms[5])/norms[6]
    end
    if qualt[i]>0
      cqobs[i]=1
      qualt[i]=(qualt[i]-qm)/qs
    end
  end

  sat=Array{Any,1}(undef,4)
  act=Array{Any,1}(undef,4)
  tstob=Array{Any,1}(undef,4)
  cqual=Array{Any,1}(undef,4)
  cqob=Array{Any,1}(undef,4)
  nsize=zeros(Int64,4)

  for i=1:4
    sat[i]=satt[igptt.==i,:]
    act[i]=actt[igptt.==i]
    tstob[i]=testobs[igptt.==i,:]
    nsize[i]=length(act[i])
    cqual[i]=qualt[igptt.==i]
    cqob[i]=cqobs[igptt.==i]
  end
 
  export sat,act,tstob, nsize,cqual,cqob
end




  


