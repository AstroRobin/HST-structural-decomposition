ellipsePlot = function(Data, modellist, model1='sersic', model2='sersic', loc1=1, loc2=2, pixscale=1, FWHM=1, SBlim=26, df=100, raw=FALSE){
  
  if(missing(Data)){stop('Data object of class profit.data must be provided!')}
  if(class(Data)!="profit.data"){stop("Data must be of class profit.data, as output by profitSetupData!")}
  if(missing(modellist)){modellist=Data$modellist}
  
  region = Data$region
  if(is.null(region)) {
    region = numeric(length(Data$image)) + 1
  }
  
  if (!is.null(model2)){
    comp1 = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list[[model1 = loc1), dim = Data$imagedim, psf=Data$psf)
    comp2 = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model2 = loc2), dim = Data$imagedim, psf=Data$psf)
    total = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model1 = 'all'), dim = Data$imagedim, psf=Data$psf)
    
    imageellipse=profitEllipse(Data$image*region, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    sigmaellipse=profitEllipse(Data$sigma*region, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    comp1ellipse=profitEllipse(comp1$z, xcen=modellist[[model1]]$xcen[loc1], ycen=modellist[[model1]]$ycen[loc1], ang=modellist[[model1]]$ang[loc1], axrat=modellist[[model1]]$axrat[loc1], box=modellist[[model1]]$box[loc1])
    comp2ellipse=profitEllipse(comp2$z, xcen=modellist[[model2]]$xcen[loc2], ycen=modellist[[model2]]$ycen[loc2], ang=modellist[[model2]]$ang[loc2], axrat=modellist[[model2]]$axrat[loc2], box=modellist[[model2]]$box[loc2])
    totalellipse=profitEllipse(total$z, xcen=modellist$sersic$xcen[loc2], ycen=modellist$sersic$ycen[loc2], ang=modellist$sersic$ang[loc2], axrat=modellist$sersic$axrat[loc2], box=modellist$sersic$box[loc2])
    psfellipse=cbind(seq(0,10*FWHM,len=1e3), dnorm(seq(0,10*FWHM,len=1e3),sd=FWHM/(2*sqrt(2*log(2)))))
    
    imageellipse[,1]=imageellipse[,1]*pixscale
    sigmaellipse[,1]=sigmaellipse[,1]*pixscale
    comp1ellipse[,1]=comp1ellipse[,1]*pixscale
    comp2ellipse[,1]=comp2ellipse[,1]*pixscale
    totalellipse[,1]=totalellipse[,1]*pixscale
    
    sigmaellipse=cbind(sigmaellipse,imageellipse[,2]-sigmaellipse[,2])
    sigmaellipse=cbind(sigmaellipse,imageellipse[,2]+sigmaellipse[,2])
    imageellipse[,2]=-2.5*suppressWarnings(log10(imageellipse[,2]))+5*log10(pixscale)+Data$magzero
    sigmaellipse[,2:4]=-2.5*suppressWarnings(log10(sigmaellipse[,2:4]))+5*log10(pixscale)+Data$magzero
    comp1ellipse[,2]=-2.5*suppressWarnings(log10(comp1ellipse[,2]))+5*log10(pixscale)+Data$magzero
    comp2ellipse[,2]=-2.5*suppressWarnings(log10(comp2ellipse[,2]))+5*log10(pixscale)+Data$magzero
    totalellipse[,2]=-2.5*suppressWarnings(log10(totalellipse[,2]))+5*log10(pixscale)+Data$magzero
    psfellipse[,2]=-2.5*suppressWarnings(log10(psfellipse[,2]))+5*log10(pixscale)+Data$magzero
    
    imageellipse[is.na(imageellipse)]=SBlim
    imageellipse[is.infinite(imageellipse)]=SBlim
    sigmaellipse[is.na(sigmaellipse)]=SBlim
    sigmaellipse[is.infinite(sigmaellipse)]=SBlim
    comp1ellipse[is.na(comp1ellipse)]=SBlim+5
    comp1ellipse[is.infinite(comp1ellipse)]=SBlim+5
    comp2ellipse[is.na(comp2ellipse)]=SBlim
    comp2ellipse[is.infinite(comp2ellipse)]=SBlim
    totalellipse[is.na(totalellipse)]=SBlim
    totalellipse[is.infinite(totalellipse)]=SBlim
    psfellipse[is.na(psfellipse)]=SBlim+5
    psfellipse[is.infinite(psfellipse)]=SBlim+5
    
    if(raw){
      predict.image=list(x=imageellipse[,1],y=imageellipse[,2])
      predict.sigma=list(x=sigmaellipse[,1],y=sigmaellipse[,2])
      predict.sigma.lo=list(x=sigmaellipse[,1],y=sigmaellipse[,3])
      predict.sigma.hi=list(x=sigmaellipse[,1],y=sigmaellipse[,4])
      predict.comp1=list(x=comp1ellipse[,1],y=comp1ellipse[,2])
      predict.comp2=list(x=comp2ellipse[,1],y=comp2ellipse[,2])
      predict.total=list(x=totalellipse[,1],y=totalellipse[,2])
      predict.psf=list(x=psfellipse[,1],y=psfellipse[,2])
    }else{
      
      smooth.image=smooth.spline(imageellipse,df=df)
      smooth.sigma.mid=smooth.spline(sigmaellipse[,c(1,2)],df=df)
      smooth.sigma.lo=smooth.spline(sigmaellipse[,c(1,3)],df=df)
      smooth.sigma.hi=smooth.spline(sigmaellipse[,c(1,4)],df=df)
      smooth.comp1=smooth.spline(comp1ellipse,df=df)
      smooth.comp2=smooth.spline(comp2ellipse,df=df)
      smooth.total=smooth.spline(totalellipse,df=df)
      
      predict.image=predict(smooth.image,imageellipse[,1])
      predict.sigma.mid=predict(smooth.sigma.mid,imageellipse[,1])
      predict.sigma.lo=predict(smooth.sigma.lo,imageellipse[,1])
      predict.sigma.hi=predict(smooth.sigma.hi,imageellipse[,1])
      predict.comp1=predict(smooth.comp1,imageellipse[,1])
      predict.comp2=predict(smooth.comp2,imageellipse[,1])
      predict.total=predict(smooth.total,imageellipse[,1])
    }
    
    psfellipse[,2]=psfellipse[,2]+min(predict.comp1$y)-min(psfellipse[1,2])
    
    smooth.total=smooth.spline(totalellipse,df=df)
    refpredict=predict(smooth.total,imageellipse[,1])
    
    xhi=refpredict$x[min(which(refpredict$y>SBlim))]
    yhi=min(predict.image$y,predict.total$y)
    
    sigma.polygon=rbind(cbind(predict.sigma.lo$x,predict.sigma.lo$y),
                        cbind(rev(predict.sigma.hi$x),rev(predict.sigma.hi$y))
    )
    
    layout(rbind(1,2), heights=c(0.7,0.4))
    par(oma=c(3.1,3.6,1.1,1.1))
    par(mar=c(0,0,0,0))
    
    if(pixscale==1){
      xlab='Project Major Axis / pix'
      ylab1=expression(mu*' (mag/pix'*''^2*')')
      ylab2=expression(Delta*mu*' (mag/pix'*''^2*')')
    }else{
      xlab='Project Major Axis / asec'
      ylab1=expression(mu*' (mag/asec'*''^2*')')
      ylab2=expression(Delta*mu*' (mag/asec'*''^2*')')
    }
    
    magplot(0,0,type='n',ylim=c(SBlim+1,yhi),xlim=c(0,xhi),col='black',xlab='', ylab=ylab1, grid=T,labels = c(F,T))
    polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
    lines(predict.image,col='black')
    lines(predict.comp1,col='red')
    lines(predict.comp2,col='blue')
    lines(predict.total,col='darkgreen')
    lines(psfellipse,col='purple')
    abline(h=SBlim,lty=3)
    abline(v=xhi,lty=3)
    abline(v=FWHM/2,lty=3)
    legend('topright',legend=c('Data','ProFit Total','ProFit Bulge','ProFit Disk','PSF'),lty=1,col=c('black','darkgreen','red','blue','purple'),bg='white')
  
  } else {
    total = profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(model1 = loc1), dim = Data$imagedim, psf=Data$psf)
  }
  
  magplot(predict.image$x, predict.image$y-predict.total$y, type='l', xlim=c(0,xhi), ylim=c(-0.5,0.5), col='darkgreen', xlab=xlab, ylab=ylab2, grid=T, labels=c(T,T))
  polygon(sigma.polygon[,1], sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)), col = hsv(v=0,alpha = 0.1), border=NA)
  abline(h=0)
  abline(v=FWHM/2,lty=3)
}