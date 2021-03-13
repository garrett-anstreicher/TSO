
module dat
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

  xdat=readdlm("bbfact08.raw")
  const igp08=trunc.(Int64,xdat[:,2])
  id=xdat[:,1]
  const N08=length(id)
  satt08=xdat[:,3:4]
  actt08=xdat[:,5]

  satt=[satt93;satt01;satt08]
  actt=[actt93;actt01;actt08]
  

  satx=satt[satt[:,1].>0,:]
  satmmn=mean(satx[:,1])
  satvmn=mean(satx[:,2])
  satmstd=std(satx[:,1])
  satvstd=std(satx[:,2])



  actx=actt[actt.>=0]
  actm=mean(actx)
  acts=std(actx)

  svnorm=[satmmn;satvmn;satmstd;satvstd;actm;acts] 
  writedlm("svnorm",svnorm)

  const sato=satt[(satt[:,1].>0).&(actt.<0),:]
  sato[:,1]=(sato[:,1].-satmmn)/satmstd
  sato[:,2]=(sato[:,2].-satvmn)/satmstd
  const acto=(actt[(satt[:,1].<0).&(actt.>0),:].-actm)/acts
  const satact=[satt[(satt[:,1].>0).&(actt.>0),:] actt[(satt[:,1].>0).&(actt.>0),:]]
  satact[:,1]=(satact[:,1].-satmmn)/satmstd
  satact[:,2]=(satact[:,2].-satvmn)/satvstd
  satact[:,3]=(satact[:,3].-actm)/acts

  const Nsato=length(sato[:,1])
  const Nacto=length(acto)
  const Nboth=length(satact[:,1])
  export sato,acto,satact,Nsato,Nacto,Nboth
end
